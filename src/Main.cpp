#include "GLApp/gl.h"
#include "GLApp/GLApp.h"
#include "GLApp/ViewBehavior.h"

#include "ImGuiCommon/ImGuiCommon.h"
#include "Parallel/Parallel.h"
#include "Tensor/Tensor.h"
#include "Common/File.h"

#include <memory>
#include <iostream>
#include <functional>
#include <fstream>
#include <algorithm>
#include <string>
#include <vector>
#include <list>
#include <regex>
#include <numeric>
#include <utility>

#include <cassert>

Parallel::Parallel parallel;

//https://stackoverflow.com/questions/16749069/c-split-string-by-regex
template<typename T>
T split(const std::string &string_to_split, const std::string& regexPattern) {
	std::regex rgx(regexPattern);
	std::sregex_token_iterator iter(string_to_split.begin(),
		string_to_split.end(),
		rgx,
		-1);
	std::sregex_token_iterator end;
	T result;
	for ( ; iter != end; ++iter) {
		result.push_back(*iter);
	}
	return result;
}

template<typename T>
std::string concat(const T& v, const std::string& sep) {
	bool first = true;
	std::string result = "";
	for (const auto& s : v) {
		if (!first) result += sep;
		result += s;
		first = false;
	}
	return result;
}

template<typename From, typename To>
To map(const From& from, std::function<typename To::value_type(typename From::value_type)> f) {
	To to;
	for (const auto& v : from) {
		to.push_back(f(v));
	}
	return to;
}

using real = double;

using int2 = Tensor::Vector<int, 2>;
using float2 = Tensor::Vector<float, 2>;
using real2 = Tensor::Vector<real, 2>;

//constexpr int vecdim = 2;
constexpr int vecdim = 3;
using vec = Tensor::Vector<real, vecdim>;	//n-dimensional vector of reals
using StateVec = Tensor::Vector<real, vecdim + 2>;

template<typename T>
typename T::value_type sum(const T& t) {
	return std::accumulate(t.begin(), t.end(), typename T::value_type());
}


struct Equations {
	virtual ~Equations() {}
	
	//TODO abstract away the arguments ...
	//or just turn all of Equation into a template parameter, and make everything compile-time
	//virtual std::pair<real, real> calcLambdaMinMax(vec normal, Prim W, real Cs) = 0;
};

struct EulerEquations : public Equations {
	real heatCapacityRatio = 1.4;

	struct Prim;

	struct Cons : public StateVec {
		real& rho() { return StateVec::v[0]; }
		vec& m() { return *(vec*)( StateVec::v + 1 ); }
		real& ETotal() { return StateVec::v[StateVec::size-1]; }
		
		const real& rho() const { return StateVec::v[0]; }
		const vec& m() const { return *(vec*)( StateVec::v + 1 ); }
		const real& ETotal() const { return StateVec::v[StateVec::size-1]; }

		Cons() {}
		
		Cons(const StateVec& v) : StateVec(v) {}

		Cons(real rho_, vec m_, real ETotal_) {
			rho() = rho_;
			m() = m_;
			ETotal() = ETotal_;
		}
	};

	struct Prim : public StateVec {
		real& rho() { return StateVec::v[0]; }
		vec& v() { return *(vec*)( StateVec::v + 1 ); }
		real& P() { return StateVec::v[StateVec::size-1]; }

		const real& rho() const { return StateVec::v[0]; }
		const vec& v() const { return *(vec*)( StateVec::v + 1 ); }
		const real& P() const { return StateVec::v[StateVec::size-1]; }

		Prim() {}
		
		Prim(const StateVec& v) : StateVec(v) {}
		
		Prim(real rho_, vec v_, real P_) {
			rho() = rho_;
			v() = v_;
			P() = P_;
		}
	};

	Cons consFromPrim(const Prim& W) {
		Cons U;
		U.rho() = W.rho();
		U.m() = W.v() * U.rho();
		U.ETotal() = W.P() / (heatCapacityRatio - 1.) + .5 * U.rho() * vec::lenSq(W.v());
		return U;
	}

	Prim primFromCons(const Cons& U) {
		Prim W;
		W.rho() = U.rho();
		W.v() = U.m() / W.rho();
		W.P() = (heatCapacityRatio - 1.) * (U.ETotal() - .5 * W.rho() * vec::lenSq(W.v()));
		return W;
	}

	real calc_hTotal(real rho, real P, real ETotal) {
		return (ETotal + P) / rho;
	}

	real calc_Cs_from_P_rho(real P, real rho) {
		return sqrt(heatCapacityRatio * P / rho);
	}

	real calc_Cs_from_v_hTotal(vec v, real hTotal) {
		return sqrt((heatCapacityRatio - 1) * (hTotal - .5 * vec::lenSq(v)));
	}
	
	virtual std::pair<real, real> calcLambdaMinMax(vec normal, Prim W, real Cs) {
		real v = vec::dot(normal, W.v());
		return std::make_pair<real, real>(v - Cs, v + Cs);
	}
};



using Cons = EulerEquations::Cons;
using Prim = EulerEquations::Prim;

EulerEquations eqn;

struct Vertex;
struct Edge;
struct Cell;

//zero-forms
struct Vertex {	//not required by finite volume algorithm
	vec pos;
	std::vector<int> edges;
	Vertex() {}
	Vertex(vec pos_) : pos(pos_) {}
};

//(n-1)-forms
struct Edge {
	vec pos;
	vec delta;
	vec normal;
	real length;
	real cellDist;
	StateVec flux;
	int vtxs[2];	//not required by finite volume algorithm
	int cells[2];	//there are always only 2 n-forms on either side of a (n-1)-form
	Edge(int va, int vb) : length(0), cellDist(0) {
		vtxs[0] = va;
		vtxs[1] = vb;
		cells[0] = -1;
		cells[1] = -1;
	}

	bool removeCell(int cellIndex) {
		if (cells[1] == cellIndex) cells[1] = -1;
		if (cells[0] == cellIndex) {
			cells[0] = cells[1];
			cells[1] = -1;
		}
		return cells[0] == -1 && cells[1] == -1;
	}
};

//n-forms
struct Cell {
	vec pos;
	real volume;
	Cons U;
	std::vector<int> edges;
	std::vector<int> vtxs;		//not required by finite volume algorithm
	Cell() : volume(0) {}
};

//2D polygon volume
real polyVol(const std::vector<vec>& vs) {
	size_t n = vs.size();
	real v = 0;
	for (size_t i = 0; i < n; ++i) {
		vec pi = vs[i];
		vec pj = vs[(i+1)%n];
		v += pi(0) * pj(1) - pi(1) * pj(0);
	}
	return .5 * v;
}

bool contains(const vec pos, std::vector<vec> vtxs) {
	size_t n = vtxs.size();
	for (size_t i = 0; i < n; ++i) {
		vec pi = vtxs[i] - pos;
		vec pj = vtxs[(i+1)%n] - pos;
		real vz = pi(0) * pj(1) - pi(1) * pj(0);
		if (vz > 0) return false;
	}
	return true;
}

struct MeshArgs {
	int2 size_ = int2(101, 101);
	MeshArgs& size(int2 x) { size_ = x; return *this; }
	
	real2 mins_ = real2(-1, -1);
	MeshArgs& mins(real2 x) { mins_ = x; return *this; }

	real2 maxs_ = real2(1, 1);
	MeshArgs& maxs(real2 x) { maxs_ = x; return *this; }

	using GridFunc = std::function<real2(real2)>;
	GridFunc grid_ = [](real2 x) -> real2 { return x; };
	MeshArgs& grid(GridFunc x) { grid_ = x; return *this; }
	
	int2 repeat_ = int2(0, 0);
	MeshArgs& repeat(int2 x) { repeat_ = x; return *this; }

	int2 capmin_ = int2(0, 0);
	MeshArgs& capmin(int2 x) { capmin_ = x; return *this; }
};

struct Mesh {
	std::vector<Vertex> vtxs;	//0-forms, which construct n and n-1 forms
	std::vector<Edge> edges;	//n-1-forms, hold flux information between cells
	std::vector<Cell> cells;	//n-forms, hold the info at each cell

//liu kang wins.  friendship.
private:
	Mesh() {}

public:

	static Mesh buildQuadChart(
		MeshArgs args = MeshArgs()
	) {
		int2 size = args.size_;
		real2 mins = args.mins_;
		real2 maxs = args.maxs_;
		auto grid = args.grid_;
		int2 repeat = args.repeat_;
		int2 capmin = args.capmin_;

		Mesh mesh;
		int m = size(0);
		int n = size(1);

		int vtxsize = m * n;
		if (capmin(0)) vtxsize++;
		mesh.vtxs.resize(vtxsize);
		for (int j = 0; j < n; ++j) {
			for (int i = 0; i < m; ++i) {
				real2 x = real2(
					((real)i + .5) / (real)size(0) * (maxs(0) - mins(0)) + mins(0),
					((real)j + .5) / (real)size(1) * (maxs(1) - mins(1)) + mins(1));
				
				real2 u = grid(x);
				std::function<real(int)> f = [&u](int i) -> real { return i < real2::size ? u(i) : 0.; };
				mesh.vtxs[i + m * j].pos = vec(f);
			}
		}
		
		int capindex = m * n;
		if (capmin(0)) {
			vec sum;
			for (int j = 0; j < n; ++j) {
				sum += mesh.vtxs[0 + m * j].pos;
			}
			mesh.vtxs[capindex].pos = sum / (real)n;
		}

		int imax = repeat(0) ? m : m-1;
		int jmax = repeat(1) ? n : n-1;
		for (int j = 0; j < jmax; ++j) {
			int jn = (j + 1) % n;
			for (int i = 0; i < imax; ++i) {
				int in = (i + 1) % m;
				mesh.addCell(std::vector<int>{i + m * j, in + m * j, in + m * jn, i + m * jn});
			}
		}

		if (capmin(0)) {
			for (int j = 0; j < jmax; ++j) {
				int jn = (j + 1) % n;
				mesh.addCell(std::vector<int>{ 0 + m * j, 0 + m * jn, capindex });
			}
		}


		mesh.calcAux();
		
		return mesh;
	}

	static Mesh buildTriChart(
		MeshArgs args = MeshArgs()
	) {
		int2 size = args.size_;
		real2 mins = args.mins_;
		real2 maxs = args.maxs_;
		auto grid = args.grid_;
		int2 repeat = args.repeat_;
		//int2 capmin = args.capmin_;
	
		Mesh mesh;
		int m = size(0);
		int n = size(1);
	
		mesh.vtxs.resize(m * n);
		for (int j = 0; j < n; ++j) {
			for (int i = 0; i < m; ++i) {
				real2 x = real2(
					((real)i + .5) / (real)size(0) * (maxs(0) - mins(0)) + mins(0),
					((real)j + .5) / (real)size(1) * (maxs(1) - mins(1)) + mins(1));
				
				real2 u = grid(x);
				std::function<real(int)> f = [&u](int i) -> real { return i < real2::size ? u(i) : 0.; };
				mesh.vtxs[i + m * j].pos = vec(f);
			}
		}
		
		int imax = repeat(0) ? m : m-1;
		int jmax = repeat(1) ? n : n-1;
		for (int j = 0; j < jmax; ++j) {
			int jn = (j + 1) % n;
			for (int i = 0; i < imax; ++i) {
				int in = (i + 1) % m;
				mesh.addCell(std::vector<int>{i + m * j, in + m * j, in + m * jn});
				mesh.addCell(std::vector<int>{in + m * jn, i + m * jn, i + m * j});
			}
		}
	
		mesh.calcAux();
		
		return mesh;
	}


	static Mesh buildFromFile(const std::string& fn) {
		Mesh mesh;	
		
		std::list<std::string> ls = split<std::list<std::string>>(Common::File::read(fn), "\n");
	
		std::string first = ls.front();
		ls.pop_front();
		std::vector<std::string> m_n = split<std::vector<std::string>>(ls.front(), "\\s+");
		ls.pop_front();
		int m = std::stoi(m_n[0]);
		int n = std::stoi(m_n[1]);
		std::list<std::string> _x = split<std::list<std::string>>(concat<std::list<std::string>>(ls, " "), "\\s+");
		if (_x.front() == "") _x.pop_front();
		if (_x.back() == "") _x.pop_back();
		std::function<real(const std::string&)> f = [](const std::string& s) -> real { return std::stod(s); };
		std::vector<real> x = map<std::list<std::string>, std::vector<real>>(_x, f);
		assert(x.size() == (size_t)(2 * m * n));
	
		auto us = std::vector(x.begin(), x.begin() + m*n);
		auto vs = std::vector(x.begin() + m*n, x.end());
		assert(us.size() == vs.size());

		mesh.vtxs.resize(m*n);
		for (int i = 0; i < (int)us.size(); ++i) {
			mesh.vtxs[i].pos = vec(us[i], vs[i]);
		}
	
		for (int j = 0; j < n-1; ++j) {
			for (int i = 0; i < m-1; ++i) {
				mesh.addCell(std::vector<int>{i + m * j, i + m * (j+1), i+1 + m * (j+1), i+1 + m * j});
			}
		}
	
		mesh.calcAux();
	
		return mesh;
	}

	int addEdge(int va, int vb) {
		int ei = 0;
		for (; ei < (int)edges.size(); ++ei) {
			Edge* e = &edges[ei];
			if ((e->vtxs[0] == va && e->vtxs[1] == vb) ||
				(e->vtxs[0] == vb && e->vtxs[1] == va)) 
			{
				return ei;
			}
		}
		assert(ei == (int)edges.size());
		edges.push_back(Edge(va, vb));
		vtxs[va].edges.push_back(ei);
		vtxs[vb].edges.push_back(ei);
		return ei;
	}

	void addCell(std::vector<int> vis) {
		Cell c;
		c.vtxs = vis;
		size_t n = vis.size();
		for (size_t i = 0; i < n; ++i) {
			c.edges.push_back(addEdge(vis[i], vis[(i+1)%n]));
		}

		c.volume = polyVol(map<std::vector<int>, std::vector<vec>>(c.vtxs, [this](int vi) -> vec { return vtxs[vi].pos; }));
		assert(c.volume > 0);

#if 1	//vertex average
		//TODO use COM
		c.pos = sum(map<
			std::vector<int>,
			std::vector<vec>
		>(c.vtxs, [this](int vi) {
			return vtxs[vi].pos;
		})) * (1. / (real)n);
#else	//COM
		for (int j = 0; j < 3; ++j) {
			for (int i = 0; i < n; ++i) {
				vec x1 = vtxs[vis[i]].pos;
				vec x2 = vtxs[vis[(i+1)%n]].pos;
				c.pos(j) += (x1(j) + x2(j)) / (x1(0) * x2(1) - x1(1) * x2(0));
			}
			c.pos(j) /= (6 * c.volume);
		}
#endif
		int ci = (int)cells.size();
		cells.push_back(c);
	
		for (int ei : c.edges) {
			Edge* e = &edges[ei];
			if (e->cells[0] == -1) {
				e->cells[0] = ci;
			} else if (e->cells[1] == -1) {
				e->cells[1] = ci;
			} else {
				throw Common::Exception() << "tried to add too many cells to an edge";
			}
		}
	}

	//calculate edge info
	//calculate cell volume info
	void calcAux() {		
		for (auto& e : edges) {
			{
				auto& a = vtxs[e.vtxs[0]];
				auto& b = vtxs[e.vtxs[1]];
				e.pos = (a.pos + b.pos) * .5;
				e.delta = a.pos - b.pos;
				e.length = vec::length(e.delta);
				e.normal = vec(e.delta(1), -e.delta(0));
				e.normal *= 1. / vec::length(e.normal);
			}

			{
				int a = e.cells[0];
				int b = e.cells[1];
				if (a != -1 && b != -1) {
					if (vec::dot(cells[a].pos - cells[b].pos, e.normal) < 0) {
						std::swap(a, b);
						e.cells[0] = a;
						e.cells[1] = b;
						e.normal *= -1;
					}
					//distance between cell centers
					//e.cellDist = vec::length(cells[b].pos - cells[a].pos);
					//distance projected to edge normal
					e.cellDist = fabs(vec::dot(e.normal, cells[b].pos - cells[a].pos));
				} else if (a != -1) {
					e.cellDist = vec::length(cells[a].pos - e.pos) * 2.;
				} else {
					throw Common::Exception() << "you are here";
				}
			}
		}
	}
};

// rotate vx,vy such that n now points along the x dir
vec rotateTo(vec v, vec n) {
	return vec(
		v(0) * n(0) + v(1) * n(1),
		v(1) * n(0) - v(0) * n(1)
	);
}

// rotate vx,vy such that the x dir now points along n 
vec rotateFrom(vec v, vec n) {
	return vec(
		v(0) * n(0) - v(1) * n(1),
		v(1) * n(0) + v(0) * n(1)
	);
}

template<typename T> void glVertex2v(const T* v);
template<> void glVertex2v<double>(const double* v) { glVertex2dv(v); }
template<> void glVertex2v<float>(const float* v) { glVertex2fv(v); }

StateVec matmul(const real* A, const StateVec& x) {
	StateVec y;
	for (int i = 0; i < StateVec::size; ++i) {
		real sum = 0;
		for (int j = 0; j < StateVec::size; ++j) {
			//C layout, so row-major
			sum += A[j + StateVec::size * i] * x(j);
		}
		y(i) = sum;
	}
	return y;
}

enum DisplayMethod {
	STATE,
	VOLUME,
	COUNT
};

static const char* displayMethodNames[DisplayMethod::COUNT] = {
	"state",
	"volume",
};

template<typename T> T cubed(const T& t) { return t * t * t; }

struct InitCond {
	virtual ~InitCond() {}
	virtual const char* name() const = 0;
	virtual Cons initCell(vec pos) const = 0;
	virtual void updateGUI() {}
};

struct InitCondConst : public InitCond {
	float rho = 1;
	float P = 1;
	float vx = 0;
	float vy = 0;

	virtual const char* name() const { return "constant"; }
	virtual Cons initCell(vec x) const {
		return eqn.consFromPrim(Prim(rho, vec(vx, vy), P));
	}
	
	virtual void updateGUI() {
		igInputFloat("rho", &rho, .1, 1, "%f", 0);
		igInputFloat("vx", &vx, .1, 1, "%f", 0);
		igInputFloat("vy", &vy, .1, 1, "%f", 0);
		igInputFloat("P", &P, .1, 1, "%f", 0);
	}
};

struct InitCondSod : public InitCond {
	virtual const char* name() const { return "Sod"; }
	virtual Cons initCell(vec x) const {
		bool lhs = x(0) < 0 && x(1) < 0;
		return eqn.consFromPrim(Prim(
			lhs ? 1. : .125,
			vec(),
			lhs ? 1. : .1
		));
	}
};

struct InitCondSpiral : public InitCond {
	virtual const char* name() const { return "Spiral"; }
	virtual Cons initCell(vec x) const {
		return eqn.consFromPrim(Prim(
			1,
			vec(-x(1), x(0)),
			1
		));
	}
};

std::vector<std::shared_ptr<InitCond>> initConds = {
	std::make_shared<InitCondConst>(),
	std::make_shared<InitCondSod>(),
	std::make_shared<InitCondSpiral>(),
};

std::vector<const char*> initCondNames = map<
	std::vector<std::shared_ptr<InitCond>>,
	std::vector<const char*>
>(
	initConds,
	[](std::shared_ptr<InitCond> ic) -> const char* { return ic->name(); }
);

std::vector<std::function<std::shared_ptr<Mesh>()>> meshes = {
	[]() -> std::shared_ptr<Mesh> {
		return std::make_shared<Mesh>(Mesh::buildFromFile("grids/n0012_113-33.p2dfmt"));
	},
	[]() -> std::shared_ptr<Mesh> {
		return std::make_shared<Mesh>(Mesh::buildQuadChart());
	},
	[]() -> std::shared_ptr<Mesh> {
		return std::make_shared<Mesh>(Mesh::buildTriChart());
	},
	[]() -> std::shared_ptr<Mesh> {
		return std::make_shared<Mesh>(Mesh::buildQuadChart(MeshArgs()
			.grid([](real2 v) -> real2 { return real2(cbrt(v(0)), cbrt(v(1))); })
		));
	},
	[]() -> std::shared_ptr<Mesh> {
		return std::make_shared<Mesh>(Mesh::buildQuadChart(MeshArgs()
			.grid([](real2 v) -> real2 { return real2(cubed(v(0)), cubed(v(1))); })
		));
	},
	[]() -> std::shared_ptr<Mesh> {
		return std::make_shared<Mesh>(Mesh::buildQuadChart(
			MeshArgs()
				.grid([](real2 v) -> real2 { 
					real r = real2::length(v);
					//real theta = std::max(0., 1. - r);
					real sigma = 3.;	//almost 0 at r=1
					const real rotationAmplitude = 3.;
					real theta = rotationAmplitude*sigma*r*exp(-sigma*sigma*r*r);
					real costh = cos(theta), sinth = sin(theta);
					return real2(
						costh * v(0) - sinth * v(1),
						sinth * v(0) + costh * v(1));
				})
		));
	},
	[]() -> std::shared_ptr<Mesh> {
		return std::make_shared<Mesh>(Mesh::buildQuadChart(
			MeshArgs()
				.size(int2(51, 201))
				.mins(real2(.5, 0))
				.maxs(real2(1, 2*M_PI))
				.grid([](real2 v) -> real2 { 
					return real2(cos(v(1)), sin(v(1))) * v(0);
				})
				.repeat(int2(0, 1))
				//.capmin(int2(1, 0))
		));
	},
};

std::vector<std::string> meshNameStrs = map<
	decltype(meshes),
	std::vector<std::string>
>(
	meshes,
	[](const std::function<std::shared_ptr<Mesh>()>& meshgen) -> std::string { return std::to_string((intptr_t)(void*)&meshgen); }
);

std::vector<const char*> meshNames = map<
	decltype(meshNameStrs),
	std::vector<const char*>
>(
	meshNameStrs,
	[](const std::string& s) -> const char* { return s.c_str(); }
);

struct CFDMeshApp : public ::GLApp::ViewBehavior<::GLApp::GLApp> {
	using Super = ::GLApp::ViewBehavior<::GLApp::GLApp>;
	
	std::shared_ptr<Mesh> m;
	std::function<Cons(vec)> initcond;

	double time = 0;
	bool running = false;
	bool singleStep = false;
	
	bool showVtxs = false;
	bool showEdges = false;
	bool showCellCenters = false;
	Cell* selectedCell = nullptr;

	int displayMethod = 0;
	float displayScalar = 1;

	int initCondIndex = 1;

	//1 = mirror boundary, -1 = freeflow boundary
	float restitution = 1;

	float2 mousepos;
		
	float cfl = .5;

	int meshIndex = 1;
	
	std::shared_ptr<ImGuiCommon::ImGuiCommon> gui = std::make_shared<ImGuiCommon::ImGuiCommon>(window, context);

	virtual const char* getTitle() {
		return "CFD Mesh";
	}
	
	CFDMeshApp(const Init& args) : Super(args) {
		
		view = viewOrtho;
		viewOrtho->zoom(0) = viewOrtho->zoom(1) = .5;

		glClearColor(.5, .75, .75, 1);
		
		resetMesh();	
		resetState();
	}

	void resetMesh() {
		m = meshes[meshIndex]();
	}

	void resetState() {
		running = false;
		singleStep = false;
		time = 0;

		InitCond* ic = initConds[initCondIndex].get();
		parallel.foreach(
			m->cells.begin(),
			m->cells.end(),
			[this, ic](auto& c) {
				c.U = ic->initCell(c.pos);
#ifdef DEBUG			
				Prim W = eqn.primFromCons(c.U);
				assert(W.rho() > 0);
				assert(vec::length(W.v()) >= 0);
				assert(W.P() > 0);
				real hTotal = calc_hTotal(W.rho(), W.P(), c.U.ETotal());
				assert(hTotal > 0);
#endif			
			}
		);
	}

	void draw() {
		glClear(GL_COLOR_BUFFER_BIT);
		for (const auto& c : m->cells) {
			Prim W = eqn.primFromCons(c.U);
			
			if (displayMethod == DisplayMethod::STATE) {
				glColor3f(displayScalar * W.rho(), displayScalar * vec::length(W.v()), displayScalar * W.P());
			} else if (displayMethod == DisplayMethod::VOLUME) {
				glColor3f(displayScalar * c.volume, .5, 1. - displayScalar * c.volume);
			}
			
			glBegin(GL_POLYGON);
			for (int vi : c.vtxs) {
				glVertex2v(m->vtxs[vi].pos.v);
			}
			glEnd();
		}
		if (selectedCell) {
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
			glLineWidth(3);
			glColor3f(1,1,0);
			glBegin(GL_POLYGON);
			for (int vi : selectedCell->vtxs) {
				glVertex2v(m->vtxs[vi].pos.v);
			}
			glEnd();
			glLineWidth(1);
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		}

		if (showVtxs || showCellCenters) {
			glColor3f(1,1,1);
			glPointSize(3);
			glBegin(GL_POINTS);
			if (showVtxs) {
				for (const auto& v : m->vtxs) {
					glVertex2v(v.pos.v);
				}
			}
			if (showCellCenters) {
				for (const auto& c : m->cells) {
					glVertex2v(c.pos.v);
				}
			}
			glEnd();
			glPointSize(1);
		}
		if (showEdges) {
			glColor3f(1,1,1);
			glBegin(GL_LINES);
			for (const auto& e : m->edges) {
				for (int vi : e.vtxs) {
					glVertex2v(m->vtxs[vi].pos.v);
				}
			}
			glEnd();
		}
	}

	std::pair<Cons, Cons> getEdgeStates(const Edge* e) {
		Cell* cL = e->cells[0] == -1 ? nullptr : &m->cells[e->cells[0]];
		Cell* cR = e->cells[1] == -1 ? nullptr : &m->cells[e->cells[1]];
		if (cL && cR) {
			return std::pair<Cons, Cons>{cL->U, cR->U};
		} else if (cL) {
			Cons UL, UR;
			UL = UR = cL->U;
			vec m = UR.m();
			UR.m() = m - e->normal * ((1 + restitution) * vec::dot(e->normal, m));
			return std::pair<Cons, Cons>{UL, UR};
		} else if (cR) {
			Cons UL, UR;
			UL = UR = cR->U;
			vec m = UL.m();
			UL.m() = m - e->normal * ((1 + restitution) * vec::dot(e->normal, m));
			return std::pair<Cons, Cons>{UL, UR};
		} 
		throw Common::Exception() << "here";
	}

	real calcDT() {
		//for (auto& c : m->cells) {
		real result = parallel.reduce(
			m->cells.begin(), 
			m->cells.end(), 
			[this](Cell& c) -> real {
				Prim W = eqn.primFromCons(c.U);
				real rho = W.rho();
				real P = W.P();
				real Cs = eqn.calc_Cs_from_P_rho(P, rho);
				real result = std::numeric_limits<real>::infinity();
#if 0 //check cartesian basis directions				
				for (int j = 0; j < vec::size; ++j) {
					vec normal;
					normal(j) = 1;
#endif
#if 1 //check mesh interfaces
				for (int ei : c.edges) {
					Edge* e = &m->edges[ei];
					vec normal = e->normal;
#endif
					real lambdaMin, lambdaMax;
					std::pair<real, real> lambdaMinMax = eqn.calcLambdaMinMax(normal, W, Cs);
					lambdaMin = lambdaMinMax.first;
					lambdaMax = lambdaMinMax.second;
					lambdaMin = std::min<real>(0, lambdaMin);
					lambdaMax = std::max<real>(0, lambdaMax);
					// TODO a better way to do this.  maybe use edges' lambdas?  maybe do this after calculating the eigenbasis?
					//real dx = sqrt(c.volume);
					real dx = e->length;
					real dum = dx / (abs(lambdaMax - lambdaMin) + 1e-9);
					result = std::min(result, dum);
				}
				return result;
			}, 
			std::numeric_limits<real>::infinity(), 
			[](real a, real b) -> real { return std::min(a,b); }
		);
		// calculate dt
		real dt = result * cfl;
		return dt;
	}

	Cons calcFluxRoe(Edge* e, Cons UL, Cons UR, real dt) {
		real ETotalL = UL.ETotal();
		Prim WL = eqn.primFromCons(UL);
		real rhoL = WL.rho();
assert(rhoL > 0);			
		vec vL = WL.v();
		real PL = WL.P();
assert(PL > 0);			
		real hTotalL = eqn.calc_hTotal(rhoL, PL, ETotalL);
assert(hTotalL > 0);			
		real sqrtRhoL = sqrt(rhoL);
		
		real ETotalR = UR.ETotal();
		Prim WR = eqn.primFromCons(UR);
		real rhoR = WR.rho();
assert(rhoR > 0);			
		vec vR = WR.v();
		real PR = WR.P();
assert(PR > 0);			
		real hTotalR = eqn.calc_hTotal(rhoR, PR, ETotalR);
assert(hTotalR > 0);			
		real sqrtRhoR = sqrt(rhoR);
			
		vec v = (vL * sqrtRhoL + vR * sqrtRhoR) / (sqrtRhoL + sqrtRhoR);
		real vx = v(0), vy = v(1), vz = v(2);
		real hTotal = (sqrtRhoL * hTotalL + sqrtRhoR * hTotalR) / (sqrtRhoL + sqrtRhoR);
		real Cs = eqn.calc_Cs_from_v_hTotal(v, hTotal);
		real CsSq = Cs * Cs;
		//e->roe = matrix{rho, vx, vy, vz, hTotal, Cs}

assert(isfinite(sqrtRhoL));
assert(isfinite(sqrtRhoR));
assert(isfinite(vx));
assert(isfinite(vy));
assert(isfinite(vz));
assert(isfinite(hTotal));
assert(isfinite(Cs));

		// 2) eigenbasis at interface
	
		real vSq = vec::lenSq(v);

		static_assert(StateVec::size == 5);

		//StateVec lambdas = {vx - Cs, vx, vx, vx, vx + Cs};
		StateVec lambdas;
		lambdas(0) = vx - Cs;
		lambdas(1) = vx;
		lambdas(2) = vx;
		lambdas(3) = vx;
		lambdas(4) = vx + Cs;
		
		real gamma_1 = eqn.heatCapacityRatio - 1;
		real evL[5][5] = {
			{(.5 * gamma_1 * vSq + Cs * vx) / (2 * CsSq),	(-Cs - gamma_1 * vx) / (2 * CsSq),	-gamma_1 * vy / (2 * CsSq),		-gamma_1 * vz / (2 * CsSq),	gamma_1 / (2 * CsSq),	},
			{1 - gamma_1 * vSq / (2 * CsSq),				gamma_1 * vx / CsSq,				gamma_1 * vy / CsSq,			gamma_1 * vz / CsSq,		-gamma_1 / CsSq,		},
			{-vy,											0,									1,								0,							0,						}, 
			{-vz,											0,									0,								1,							0,						},
			{(.5 * gamma_1 * vSq - Cs * vx) / (2 * CsSq),	(Cs - gamma_1 * vx) / (2 * CsSq),	-gamma_1 * vy / (2 * CsSq),		-gamma_1 * vz / (2 * CsSq),	gamma_1 / (2 * CsSq),	},
		};

		real evR[5][5] = {
			{1, 				1, 			0,		0,		1,				},
			{vx - Cs, 			vx, 		0,		0,		vx + Cs,		},
			{vy,				vy,			1,		0,		vy,				},
			{vz,				vz,			0,		1,		vz,				},
			{hTotal - Cs * vx, .5 * vSq, 	vy,		vz,		hTotal + Cs * vx},
		};
	
		Cons dU = UR - UL;
		Cons dUTilde = matmul(&evL[0][0], dU);
	
		Cons fluxTilde;
		for (int j = 0; j < 5; ++j) {
			real lambda = lambdas(j);
			real phi = 0;
			real sgnLambda = lambda >= 0 ? 1 : -1;
			real dx = e->cellDist;
			real epsilon = lambda * dt / dx;
			fluxTilde(j) = -.5 * lambda * dUTilde(j) * (sgnLambda + phi * (epsilon - sgnLambda));
		}
	
		Cons UAvg = (UR + UL) * .5;
		Cons UAvgTilde = matmul(&evL[0][0], UAvg);
		fluxTilde = fluxTilde + lambdas * UAvgTilde;
	
		Cons flux = matmul(&evR[0][0], fluxTilde);
		// here's the flux, aligned along the normal

		return flux;
	}

	Cons calcFluxHLL(Edge* e, Cons UL, Cons UR, real dt) {
		Prim WL = eqn.primFromCons(UL);
		Prim WR = eqn.primFromCons(UR);
		
		real gamma_1 = eqn.heatCapacityRatio - 1;
		
		real densityL = WL.rho();
		real invDensityL = 1. / densityL;
		vec velocityL = WL.v();
		real velocitySqL = vec::dot(velocityL, velocityL);
		real energyTotalL = UL.ETotal() * invDensityL;
		real energyKineticL = .5 * velocitySqL;
		real energyPotentialL = 0;//potentialBuffer[indexPrev];
		real energyInternalL = energyTotalL - energyKineticL - energyPotentialL;
		real pressureL = gamma_1 * densityL * energyInternalL;
		real enthalpyTotalL = energyTotalL + pressureL * invDensityL;
		real speedOfSoundL = sqrt(gamma_1 * (enthalpyTotalL - .5 * velocitySqL));
		real roeWeightL = sqrt(densityL);

		real densityR = WR.rho();
		real invDensityR = 1. / densityR;
		vec velocityR = WR.v();
		real velocitySqR = vec::dot(velocityR, velocityR);
		real energyTotalR = UR.ETotal() * invDensityR;
		real energyKineticR = .5 * velocitySqR;
		real energyPotentialR = 0;//potentialBuffer[index];
		real energyInternalR = energyTotalR - energyKineticR - energyPotentialR;
		real pressureR = gamma_1 * densityR * energyInternalR;
		real enthalpyTotalR = energyTotalR + pressureR * invDensityR;
		real speedOfSoundR = sqrt(gamma_1 * (enthalpyTotalR - .5 * velocitySqR));
		real roeWeightR = sqrt(densityR);
	
		real roeWeightNormalization = 1. / (roeWeightL + roeWeightR);
		vec velocity = (velocityL * roeWeightL + velocityR * roeWeightR) * roeWeightNormalization;
		real enthalpyTotal = (roeWeightL * enthalpyTotalL + roeWeightR * enthalpyTotalR) * roeWeightNormalization;
		real energyPotential = (roeWeightL * energyPotentialL + roeWeightR * energyPotentialR) * roeWeightNormalization; 
	
		real velocitySq = vec::dot(velocity, velocity);
		real speedOfSound = sqrt(gamma_1 * (enthalpyTotal - .5 * velocitySq - energyPotential));

		//eigenvalues

		real eigenvaluesMinL = velocityL(0) - speedOfSoundL;
		real eigenvaluesMaxR = velocityR(0) + speedOfSoundR;
		real eigenvaluesMin = velocity(0) - speedOfSound;
		real eigenvaluesMax = velocity(0) + speedOfSound;
	#if 0	//Davis direct
		real sl = eigenvaluesMinL;
		real sr = eigenvaluesMaxR;
	#endif
	#if 1	//Davis direct bounded
		real sl = std::min(eigenvaluesMinL, eigenvaluesMin);
		real sr = std::max(eigenvaluesMaxR, eigenvaluesMax);
	#endif

		Cons eigenvalues;
		eigenvalues(0) = sl;
		eigenvalues(1) = velocity(0);
		eigenvalues(2) = velocity(0);
		eigenvalues(3) = velocity(0);
		eigenvalues(4) = sr;

		//flux

		Cons fluxL;
		fluxL(0) = densityL * velocityL(0);
		fluxL(1) = densityL * velocityL(0) * velocityL(0) +  pressureL;
		fluxL(2) = densityL * velocityL(1) * velocityL(0);
		fluxL(3) = densityL * velocityL(2) * velocityL(0);
		fluxL(4) = densityL * enthalpyTotalL * velocityL(0);
		
		Cons fluxR;
		fluxR(0) = densityR * velocityR(0);
		fluxR(1) = densityR * velocityR(0) * velocityR(0) + pressureR;
		fluxR(2) = densityR * velocityR(1) * velocityR(0);
		fluxR(3) = densityR * velocityR(2) * velocityR(0);
		fluxR(4) = densityR * enthalpyTotalR * velocityR(0);	

		//HLL
		Cons flux;
		if (0. <= sl) {
			flux = fluxL;
		} else if (sl <= 0. && 0. <= sr) {
			//(sr * fluxL[j] - sl * fluxR[j] + sl * sr * (UR[j] - UL[j])) / (sr - sl)
			real invDenom = 1. / (sr - sl);
			for (int i = 0; i < StateVec::size; ++i) {
				flux(i) = (sr * fluxL(i) - sl * fluxR(i) + sl * sr * (UR(i) - UL(i))) * invDenom; 
			}
		} else if (sr <= 0.) {
			flux = fluxR;
		}
	
		return flux;
	}


	using CalcFluxFunc = Cons (CFDMeshApp::*)(Edge*, Cons, Cons, real);
	std::vector<CalcFluxFunc> calcFluxes = {
		&CFDMeshApp::calcFluxRoe,
		&CFDMeshApp::calcFluxHLL,
	};

	std::vector<const char*> calcFluxNames = {
		"Roe",
		"HLL",
	};

	int calcFluxIndex = 0;
	
	void step() {
		real dt = calcDT();

		//for (auto& e : m->edges) {
		parallel.foreach(
			m->edges.begin(),
			m->edges.end(),
			[this, dt](Edge& e) {
				std::pair<Cons, Cons> ULR = getEdgeStates(&e);
				Cons UL = ULR.first;
				Cons UR = ULR.second;
assert(e.cells[0] == -1 || UL == m->cells[e.cells[0]].U);
assert(e.cells[1] == -1 || UR == m->cells[e.cells[1]].U);
assert(UL(3) == 0);
assert(UR(3) == 0);

for (int i = 0; i < StateVec::size; ++i) {
	assert(isfinite(UL(0)));	
	assert(isfinite(UR(0)));	
}	

				// 1) roe values at edge 
				// rotate to align edge normal to x axis
				// so x-direction flux jacobian is good for calculating the flux 
				UL.m() = rotateTo(UL.m(), e.normal);
				UR.m() = rotateTo(UR.m(), e.normal);

				Cons flux = (this->*calcFluxes[calcFluxIndex])(&e, UL, UR, dt);
			
				flux.m() = rotateFrom(flux.m(), e.normal);
			
for (int i = 0; i < StateVec::size; ++i) {
	assert(isfinite(flux(0)));	
}	
				
				// here's the flux in underlying coordinates
				e.flux = flux;
assert(e.flux(3) == 0);
			}
		);

		//for (auto& c : m->cells) {
		parallel.foreach(
			m->cells.begin(),
			m->cells.end(),
			[this, dt](Cell& c) {
				Cons dU_dt;
				for (int ei : c.edges) {
					Edge* e = &m->edges[ei];
					//if (e.hasFlux) {
						if (&c == &m->cells[e->cells[0]]) {
							dU_dt -= e->flux * (e->length / c.volume);
						} else {
							dU_dt += e->flux * (e->length / c.volume);
						}
					//}
				}
				c.U += dU_dt * dt;
#ifdef DEBUG	
Prim W = eqn.primFromCons(c.U);
assert(W.rho() > 0);
assert(vec::length(W.v()) >= 0);
assert(W.P() > 0);
real hTotal = eqn.calc_hTotal(W.rho(), W.P(), c.U.ETotal());
assert(hTotal > 0);
Cons U2 = eqn.consFromPrim(W);
for (int i = 0; i < StateVec::size; ++i) {
	real delta = U2(i) - c.U(i);
	assert(fabs(delta) < 1e-9);
}
#endif			
			}
		);
	
		time += dt;
	}

	virtual void onUpdate() {
		Super::onUpdate();
		draw();
		
		gui->onUpdate([this](){
			igText("time: %f", time);
			igCheckbox("running", &running);
			
			if (igSmallButton("step")) singleStep = true;
		
			igInputFloat("cfl", &cfl, .1, 1, "%f", 0);

			igCheckbox("showVtxs", &showVtxs);
			igCheckbox("showCellCenters", &showCellCenters);
			igCheckbox("showEdges", &showEdges);
			//igCheckbox("ortho", &ortho);
		
			if (igSmallButton("reset state")) {
				resetState();
			}
	
			igCombo("display method", &displayMethod, displayMethodNames, DisplayMethod::COUNT, -1);
			
			igInputFloat("display scalar", &displayScalar, .1, 1., "%f", 0);
			igInputFloat("restitution", &restitution, .1, 1., "%f", 0);
		
			igCombo("flux", &calcFluxIndex, calcFluxNames.data(), calcFluxNames.size(), -1);
			
			igCombo("init cond", &initCondIndex, initCondNames.data(), initCondNames.size(), -1);

			igPushIDStr("init cond fields");
			initConds[initCondIndex]->updateGUI();
			igPopID();

			if (igCombo("mesh", &meshIndex, meshNames.data(), meshNames.size(), -1)) {
				resetMesh();
			}


			if (view == viewOrtho) {
				//find the view bounds
				//find the mouse position
				//find any cells at that position
				real2 pos2 = (real2)((mousepos - .5f) * float2(1, -1) * float2(aspectRatio, 1) / float2(viewOrtho->zoom(0), viewOrtho->zoom(1)) + viewOrtho->pos);
				std::function<real(int)> f = [&](int i) -> real { return i >= 2 ? 0 : pos2(i); };
				vec pos = vec(f);
					
				igBeginTooltip();
				igText("%f %f\n", pos(0), pos(1));
				igEndTooltip();

				selectedCell = nullptr;
				int selectedCellIndex = -1;
				for (int i = 0; i < (int)m->cells.size(); ++i) {
					Cell* c = &m->cells[i];
					if (contains(pos, map<std::vector<int>, std::vector<vec>>(c->vtxs, [this](int vi) -> vec { return m->vtxs[vi].pos; }))) {
						selectedCellIndex = i;
						break;
					}
				}

				if (selectedCellIndex != -1) {
					Cell* c = &m->cells[selectedCellIndex];
					selectedCell = c;
					
					igBeginTooltip();
					Prim W = eqn.primFromCons(c->U);
					igText("rho %f", W.rho());
					igText(" vx %f", W.v()(0));
					igText(" vy %f", W.v()(1));
					igText("  P %f", W.P());
					igEndTooltip();

#if 0 //erase
					if (leftButtonDown) {
						std::vector<int> toRemoveEdges;
						//remove c...
						for (int ei : c->edges) {
							Edge* e = &m->edges[ei];
							if (e->removeCell(i - m->cells.begin())) {
								//all cells on e are removed
								toRemoveEdges.push_back(ei);
							}
						}
						m->cells.erase(m->cells.begin() + selectedCellIndex);
					
						//sort largest -> smallest
						std::sort(toRemoveEdges.begin(), toRemoveEdges.end(), [](int a, int b) -> bool { return a > b; });
						for (int ei : toRemoveEdges) {
							m->edges.erase(m->edges.begin() + ei);
						}
					}
#endif			
				}
			}



		});
		
		if (running || singleStep) {
			step();
			if (singleStep) {
				running = false;
				singleStep = false;
			}
		}
		
	}

	virtual void onSDLEvent(SDL_Event& event) {
		if (gui) gui->onSDLEvent(event);
		bool canHandleMouse = !igGetIO()->WantCaptureMouse;
		bool canHandleKeyboard = !igGetIO()->WantCaptureKeyboard;
		
		if (canHandleMouse) {
			Super::onSDLEvent(event);
		}
		switch (event.type) {
		case SDL_KEYUP:
			if (canHandleKeyboard) {
				if (event.key.keysym.sym == SDLK_r) {
					resetState();
				} else if (event.key.keysym.sym == SDLK_u) {
					singleStep = true;
				} else if (event.key.keysym.sym == SDLK_SPACE) {
					running = !running;
				}
			}
			break;
		case SDL_MOUSEMOTION:
			mousepos = float2(
				(float)event.motion.x / (float)screenSize(0),
				(float)event.motion.y / (float)screenSize(1));
			break;
		}
	}
};

GLAPP_MAIN(CFDMeshApp)
