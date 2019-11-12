#include "GLApp/gl.h"
#include "ImGuiCommon/ImGuiCommon.h"
#include "GLApp/GLApp.h"

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
using real2 = Tensor::Vector<real, 2>;

//constexpr int vecdim = 2;
constexpr int vecdim = 3;
using vec = Tensor::Vector<real, vecdim>;	//n-dimensional vector of reals
using StateVec = Tensor::Vector<real, vecdim + 2>;

template<typename T>
typename T::value_type sum(const T& t) {
	return std::accumulate(t.begin(), t.end(), typename T::value_type());
}

real calc_hTotal(real rho, real P, real ETotal) {
	return (ETotal + P) / rho;
}

const real heatCapacityRatio = 1.4;
real calcSpeedOfSound(vec v, real hTotal) {
	return sqrt((heatCapacityRatio - 1) * (hTotal - .5 * vec::lenSq(v)));
}


struct Cons;
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
	
	Cons(const Prim& prim);
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
	
	Prim(const Cons& U);
};

Cons::Cons(const Prim& W) {
	rho() = W.rho();
	m() = W.v() * rho();
	ETotal() = W.P() / (heatCapacityRatio - 1.) + .5 * rho() * vec::lenSq(W.v());
}

Prim::Prim(const Cons& U) {
	rho() = U.rho();
	v() = U.m() / rho();
	P() = (heatCapacityRatio - 1.) * (U.ETotal() - .5 * rho() * vec::lenSq(v()));
}


struct Vertex;
struct Edge;
struct Cell;

struct Vertex {
	vec pos;
	std::vector<int> edges;
	Vertex() {}
	Vertex(vec pos_) : pos(pos_) {}
};

struct Edge {
	vec pos;
	vec delta;
	vec normal;
	real length;
	real cellDist;
	StateVec flux;
	int vtxs[2];	//this is always two
	int cells[2];	//this is always two
	Edge(int va, int vb) : length(0), cellDist(0) {
		vtxs[0] = va;
		vtxs[1] = vb;
		cells[0] = -1;
		cells[1] = -1;
	}
};

struct Cell {
	vec pos;
	real volume;
	Cons U;
	std::vector<int> edges;
	std::vector<int> vtxs;
	Cell() : volume(0) {}
};

//2D polygon volume
real polyVol(const std::vector<vec>& vs) {
	size_t n = vs.size();
	real v = 0;
	for (size_t i = 0; i < n; ++i) {
		vec pi = vs[i];
		vec pj = vs[(i+1)%n];
		v += .5 * (pi(0) * pj(1) - pi(1) * pj(0));
	}
	return fabs(v);
}

struct Mesh {
	std::vector<Vertex> vtxs;
	std::vector<Edge> edges;
	std::vector<Cell> cells;

	Mesh(
		int2 size,
		real2 mins,
		real2 maxs,
		std::function<real2(real2)> grid,
		Parallel::Parallel& parallel
	) {
		int m = size(0);
		int n = size(1);
	
		vtxs.resize(m * n);
		for (int j = 0; j < n; ++j) {
			for (int i = 0; i < m; ++i) {
				real2 x = real2(
					((real)i + .5) / (real)size(0),
					((real)j + .5) / (real)size(1)) * (maxs - mins) + mins;
				
				real2 u = grid(x);
				std::function<real(int)> f = [&u](int i) -> real { return i < real2::size ? u(i) : 0.; };
				vtxs[i + m * j].pos = vec(f);
			}
		}
		
		for (int i = 0; i < n-1; ++i) {
			for (int j = 0; j < m-1; ++j) {
				addCell(std::vector<int>{j + m * i, j + m * (i+1), j+1 + m * (i+1), j+1 + m * i});
			}
		}
	
		calcAux(parallel);
	}

	Mesh(const std::string& fn, Parallel::Parallel& parallel) {
		std::list<std::string> ls = split<std::list<std::string>>(Common::File::read("grids/n0012_113-33.p2dfmt"), "\n");
	
		std::string first = ls.front();
		ls.pop_front();
		std::vector<std::string> m_n = split<std::vector<std::string>>(ls.front(), "\\s+");
		ls.pop_front();
		int m = std::stoi(m_n[0]);
		int n = std::stoi(m_n[1]);
		std::list<std::string> _x = split<std::list<std::string>>(concat<std::list<std::string>>(ls, " "), "\\s+");
		if (_x.front() == "") _x.pop_front();
		if (_x.back() == "") _x.pop_back();
		std::vector<real> x = map<
			std::list<std::string>,
			std::vector<real>
		>(
			_x,
			[](const std::string& s) -> real {
				return std::stod(s);
			}
		);
		assert(x.size() == (size_t)(2 * m * n));
	
		auto us = std::vector(x.begin(), x.begin() + m*n);
		auto vs = std::vector(x.begin() + m*n, x.end());
		assert(us.size() == vs.size());

		vtxs.resize(m*n);
		for (int i = 0; i < (int)us.size(); ++i) {
			vtxs[i].pos = vec(us[i], vs[i]);
		}
	
		for (int i = 0; i < n-1; ++i) {
			for (int j = 0; j < m-1; ++j) {
				addCell(std::vector<int>{j + m * i, j + m * (i+1), j+1 + m * (i+1), j+1 + m * i});
			}
		}
	
		calcAux(parallel);
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
		
		//TODO com vs average pos
		c.pos = sum(map<
			std::vector<int>,
			std::vector<vec>
		>(c.vtxs, [this](int vi) {
			return vtxs[vi].pos;
		})) * (1. / (real)n);

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
	void calcAux(Parallel::Parallel& parallel) {		
		for (auto& e : edges) {
			{
				auto& a = vtxs[e.vtxs[0]];
				auto& b = vtxs[e.vtxs[1]];
				e.pos = (a.pos + b.pos) * .5;
				e.delta = a.pos - b.pos;
				e.length = vec::length(e.delta);
				e.normal = vec(-e.delta(1), e.delta(0));
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
					}
					e.cellDist = vec::length(cells[b].pos - cells[a].pos);
				} else if (a != -1) {
					e.cellDist = vec::length(cells[a].pos - e.pos) * 2.;
				} else {
					throw Common::Exception() << "you are here";
				}
			}
		}

		parallel.foreach(
			cells.begin(),
			cells.end(),
			[this](Cell& c) {
				c.volume = polyVol(map<
						std::vector<int>,
						std::vector<vec>
					>(c.vtxs, [this](int vi) -> vec {
						return vtxs[vi].pos;
					}));
			}
		);
		

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

struct CFDMeshApp : public GLApp::GLApp {
	using Super = ::GLApp::GLApp;
	
	std::shared_ptr<Mesh> m;
	std::shared_ptr<ImGuiCommon::ImGuiCommon> gui;
	std::function<Cons(vec)> initcond;

	bool running = false;
	bool singleStep = false;
	bool showVtxs = false;
	bool showEdges = false;
	bool showCellCenters = false;

	int displayMethod = 0;
	float displayScalar = 1;

	//1 = mirror boundary, -1 = freeflow boundary
	float restitution = 1;

	Parallel::Parallel parallel;

	virtual const char* getTitle() {
		return "CFD Mesh";
	}
	
	virtual void init() {
		Super::init();
		glClearColor(.5, .75, .75, 1);

		gui = std::make_shared<ImGuiCommon::ImGuiCommon>(window, context);
		
		//m = std::make_shared<Mesh>("grids/n0012_113-33.p2dfmt", parallel);
		m = std::make_shared<Mesh>(int2(101, 101), real2(-1), real2(1), [](real2 v) -> real2 { return v; }, parallel);
	
#if 0	//constant velocity
		initcond = [](vec) -> Cons {
			return Cons(Prim(1.,
				//behavior inconsistency:
				// vec(x) produces [x, x, ..., x]
				// vec(x, 0) produces [x, 0, ..., 0]
				vec(.1, 0.),
				1.));
		};
#endif
#if 1	//Sod
		initcond = [](vec x) -> Cons {
			bool lhs = x(0) < 0 && x(1) < 0;
			return Cons(Prim(
				lhs ? 1. : .125,
				vec(),
				lhs ? 1. : .1
			));
		};
#endif

		resetState();
	}

	void resetState() {
		running = false;
		singleStep = false;

		parallel.foreach(
			m->cells.begin(),
			m->cells.end(),
			[this](auto& c) {
				c.U = initcond(c.pos);
			
				Prim W(c.U);
				assert(W.rho() > 0);
				assert(vec::length(W.v()) >= 0);
				assert(W.P() > 0);
				real hTotal = calc_hTotal(W.rho(), W.P(), c.U.ETotal());
				assert(hTotal > 0);
			}
		);
	}

	virtual void shutdown() {
		gui = nullptr;
		Super::shutdown();
	}

	virtual void resize(int width, int height) {
		Super::resize(width, height);
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		
		float aspectRatio = (float)width / (float)height;
		
		//float zNear = .1f;
		//float zFar = 100;
		//glFrustum(-aspectRatio * zNear, aspectRatio * zNear, -zNear, zNear, zNear, zFar);
		
		float zNear = -1;
		float zFar = 1;
		float w = 1.;
		glOrtho(-aspectRatio * w, aspectRatio * w, -w, w, zNear, zFar);
	}
	
	void draw() {
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		//glTranslatef(0.f, 0.f, -2.f);
		
		glClear(GL_COLOR_BUFFER_BIT);
		glBegin(GL_QUADS);
		for (const auto& c : m->cells) {
			Prim W(c.U);
			
			if (displayMethod == DisplayMethod::STATE) {
				glColor3f(displayScalar * W.rho(), displayScalar * vec::length(W.v()), displayScalar * W.P());
			} else if (displayMethod == DisplayMethod::VOLUME) {
				glColor3f(displayScalar * c.volume, .5, 1. - displayScalar * c.volume);
			}
			
			for (int vi : c.vtxs) {
				glVertex2v(m->vtxs[vi].pos.v);
			}
		}
		glEnd();

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

	void step() {
		//for (auto& c : m->cells) {
		real result = parallel.reduce(
			m->cells.begin(), 
			m->cells.end(), 
			[this](Cell& c) -> real {
				Prim W(c.U);
				real rho = W.rho();
				real P = W.P();
				real Cs = sqrt(heatCapacityRatio * P / rho);
				real result = std::numeric_limits<real>::infinity();
				for (int j = 0; j < vec::size; ++j) {
					real v = W.v()(j);
					//local vx, vy = rotateTo(vx,vy, e->normal)
					//assert(vy == 0)
					real lambdaMin = v - Cs;
					real lambdaMax = v + Cs;
					lambdaMin = std::min<real>(0, lambdaMin);
					lambdaMax = std::max<real>(0, lambdaMax);
					// TODO a better way to do this.  maybe use edges' lambdas?  maybe do this after calculating the eigenbasis?
					real dx = sqrt(c.volume);
					real dum = dx / (abs(lambdaMax - lambdaMin) + 1e-9);
					result = std::min(result, dum);
				}
				return result;
			}, 
			std::numeric_limits<real>::infinity(), 
			[](real a, real b) -> real { return std::min(a,b); }
		);
		// calculate dt
		real cfl = .5;
		real dt = result * cfl;

		//for (auto& e : m->edges) {
		parallel.foreach(
			m->edges.begin(), 
			m->edges.end(), 
			[](auto& e) {
				e.flux = StateVec();
#if 0
for (int ci : e.cells) {
	Prim W(cells[ci].U);
	Cons U(W);
	std::cout << (U - cells[ci].U) << std::endl;
}
#endif
			}
		);

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

				real ETotalL = UL.ETotal();
				Prim WL(UL);
				real rhoL = WL.rho();
assert(rhoL > 0);			
				vec vL = WL.v();
				real PL = WL.P();
assert(PL > 0);			
				real hTotalL = calc_hTotal(rhoL, PL, ETotalL);
assert(hTotalL > 0);			
				real sqrtRhoL = sqrt(rhoL);
				
				real ETotalR = UR.ETotal();
				Prim WR(UR);
				real rhoR = WR.rho();
assert(rhoR > 0);			
				vec vR = WR.v();
				real PR = WR.P();
assert(PR > 0);			
				real hTotalR = calc_hTotal(rhoR, PR, ETotalR);
assert(hTotalR > 0);			
				real sqrtRhoR = sqrt(rhoR);
					
				vec v = (vL * sqrtRhoL + vR * sqrtRhoR) / (sqrtRhoL + sqrtRhoR);
				real vx = v(0), vy = v(1), vz = v(2);
				real hTotal = (sqrtRhoL * hTotalL + sqrtRhoR * hTotalR) / (sqrtRhoL + sqrtRhoR);
				real Cs = calcSpeedOfSound(v, hTotal);
				real CsSq = Cs * Cs;
				//e.roe = matrix{rho, vx, vy, vz, hTotal, Cs}

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
				
				real gamma_1 = heatCapacityRatio - 1;
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
					real dx = e.cellDist;
					real epsilon = lambda * dt / dx;
					fluxTilde(j) = -.5 * lambda * dUTilde(j) * (sgnLambda + phi * (epsilon - sgnLambda));
				}
			
				Cons UAvg = (UR + UL) * .5;
				Cons UAvgTilde = matmul(&evL[0][0], UAvg);
				fluxTilde = fluxTilde + lambdas * UAvgTilde;
			
				Cons flux = matmul(&evR[0][0], fluxTilde);
				// here's the flux, aligned along the normal
			
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
	
Prim W(c.U);
assert(W.rho() > 0);
assert(vec::length(W.v()) >= 0);
assert(W.P() > 0);
real hTotal = calc_hTotal(W.rho(), W.P(), c.U.ETotal());
assert(hTotal > 0);
Cons U2(W);
for (int i = 0; i < StateVec::size; ++i) {
	real delta = U2(i) - c.U(i);
	assert(fabs(delta) < 1e-9);
}
			}
		);
	}

	virtual void update() {
		Super::update();
		draw();
		
		gui->update([this](){
			igCheckbox("running", &running);
			
			if (igSmallButton("step")) singleStep = true;
			
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
		});
		
		if (running || singleStep) {
			step();
			if (singleStep) {
				running = false;
				singleStep = false;
			}
		}
	}

	virtual void sdlEvent(SDL_Event& event) {
		Super::sdlEvent(event);
		if (gui) gui->sdlEvent(event);
		//bool canHandleMouse = !igGetIO()->WantCaptureMouse;
		bool canHandleKeyboard = !igGetIO()->WantCaptureKeyboard;
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
		}
	}
};

GLAPP_MAIN(CFDMeshApp)
