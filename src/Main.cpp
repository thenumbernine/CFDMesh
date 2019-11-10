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
	std::vector<std::shared_ptr<Edge>> edges;
	Vertex(vec pos_) : pos(pos_) {}
};

struct Edge {
	vec pos;
	vec delta;
	vec normal;
	real length;
	real cellDist;
	StateVec flux;
	std::vector<std::shared_ptr<Vertex>> vtxs;
	std::vector<std::shared_ptr<Cell>> cells;
	Edge(std::vector<std::shared_ptr<Vertex>> vtxs_) : length(0), cellDist(0), vtxs(vtxs_) {}
};

struct Cell {
	vec pos;
	real volume;
	Cons U;
	std::vector<std::shared_ptr<Edge>> edges;
	std::vector<std::shared_ptr<Vertex>> vtxs;
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
	std::vector<std::shared_ptr<Vertex>> vtxs;
	std::vector<std::shared_ptr<Edge>> edges;
	std::vector<std::shared_ptr<Cell>> cells;

	Mesh(int2 size, real2 mins, real2 maxs, std::function<real2(real2)> grid) {
		int m = size(0);
		int n = size(1);
		for (int j = 0; j < n; ++j) {
			for (int i = 0; i < m; ++i) {
				real2 x = real2(
					((real)i + .5) / (real)size(0),
					((real)j + .5) / (real)size(1)) * (maxs - mins) + mins;
				
				//addVtx(vec(grid(u)));
				
				real2 u = grid(x);
				std::function<real(int)> f = [&u](int i) -> real { return i < vecdim ? u(i) : 0.; };
				addVtx(vec(f));
			}
		}
		
		for (int i = 0; i < n-1; ++i) {
			for (int j = 0; j < m-1; ++j) {
				addCell(std::vector<int>{j + m * i, j + m * (i+1), j+1 + m * (i+1), j+1 + m * i});
			}
		}
	
		calcAux();
		setInitState();
	}

	Mesh(const std::string& fn) {
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
	
		for (int i = 0; i < (int)us.size(); ++i) {
			addVtx(vec(us[i], vs[i]));
		}
		assert(vtxs.size() == (size_t)(m*n));
	
		for (int i = 0; i < n-1; ++i) {
			for (int j = 0; j < m-1; ++j) {
				addCell(std::vector<int>{j + m * i, j + m * (i+1), j+1 + m * (i+1), j+1 + m * i});
			}
		}

		calcAux();
		setInitState();
	}

	void addVtx(vec x) {
		vtxs.push_back(std::make_shared<Vertex>(x));
	}

	std::shared_ptr<Edge> addEdge(int i, int j) {
		std::shared_ptr<Vertex> a = vtxs[i];
		std::shared_ptr<Vertex> b = vtxs[j];
		for (const auto& e : edges) {
			if ((e->vtxs[0] == a && e->vtxs[1] == b) ||
				(e->vtxs[0] == b && e->vtxs[1] == a)) 
			{
				return e;
			}
		}
		std::shared_ptr<Edge> e = std::make_shared<Edge>(std::vector<std::shared_ptr<Vertex>>{a,b});
		edges.push_back(e);
		a->edges.push_back(e);
		b->edges.push_back(e);
		return e;
	}

	void addCell(std::vector<int> indexes) {
		std::shared_ptr<Cell> c = std::make_shared<Cell>();
		size_t n = indexes.size();
		for (size_t i = 0; i < n; ++i) {
			c->vtxs.push_back(vtxs[indexes[i]]);
		}
		for (size_t i = 0; i < n; ++i) {
			c->edges.push_back(addEdge(indexes[i], indexes[(i+1)%n]));
		}
		cells.push_back(c);
		c->pos = sum(map<
			std::vector<std::shared_ptr<Vertex>>,
			std::vector<vec>
		>(c->vtxs, [](std::shared_ptr<Vertex> v) { return v->pos; })) * (1. / (real)n);
		for (auto& e : c->edges) {
			e->cells.push_back(c);
		}
	}

	//calculate edge info
	//calculate cell volume info
	void calcAux() {
		for (auto& e : edges) {
			{
				auto a = e->vtxs[0];
				auto b = e->vtxs[1];
				e->pos = (a->pos + b->pos) * .5;
				e->delta = a->pos - b->pos;
				e->length = vec::length(e->delta);
				e->normal = vec(-e->delta(1), e->delta(0));
				e->normal *= 1. / vec::length(e->normal);
			}

			{
				std::shared_ptr<Cell> a, b;
				if (e->cells.size() > 0) a = e->cells[0];
				if (e->cells.size() > 1) b = e->cells[1];
				if (a && b) {
					if (vec::dot(a->pos - b->pos, e->normal) < 0) {
						std::swap(a, b);
						e->cells = std::vector<std::shared_ptr<Cell>>{a, b};
					}
					e->cellDist = vec::length(b->pos - a->pos);
				} else {
					e->cellDist = vec::length(a->pos - e->pos) * 2.;
				}
			}
		}
	
		for (auto& c : cells) {
			c->volume = polyVol(map<
					std::vector<std::shared_ptr<Vertex>>,
					std::vector<vec>
				>(c->vtxs, [](const std::shared_ptr<Vertex>& v) -> vec {
					return v->pos;
				}));
		}
	}

	void setInitState() {
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

template<typename T> void glVertex2v(T* v);
template<> void glVertex2v<double>(double* v) { glVertex2dv(v); }
template<> void glVertex2v<float>(float* v) { glVertex2fv(v); }

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
	bool showVtxs = true;
	bool showCellCenters = false;
	bool showEdges = true;

	int displayMethod = 0;
	float displayScalar = 1.;

	Parallel::Parallel parallel;

	virtual const char* getTitle() {
		return "CFD Mesh";
	}
	
	virtual void init() {
		Super::init();
		glClearColor(.5, .75, .75, 1.);

		gui = std::make_shared<ImGuiCommon::ImGuiCommon>(window, context);
		
		//m = std::make_shared<Mesh>("grids/n0012_113-33.p2dfmt");
		m = std::make_shared<Mesh>(
			int2(128, 128),
			real2(-1),
			real2(1),
			[](real2 v) { return v; }
		);
	
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
		//for (auto& c : m->cells) {
		parallel.foreach(
			m->cells.begin(),
			m->cells.end(),
			[this](std::shared_ptr<Cell> c) {
				c->U = initcond(c->pos);
//std::cout << c->U << std::endl;
			
				Prim W(c->U);
				assert(W.rho() > 0);
				assert(vec::length(W.v()) > 0);
				assert(W.P() > 0);
				real hTotal = calc_hTotal(W.rho(), W.P(), c->U.ETotal());
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
		float aspectRatio = (float)width / (float)height;
		float zNear = .1f;
		float zFar = 100.f;
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glFrustum(-aspectRatio * zNear, aspectRatio * zNear, -zNear, zNear, zNear, zFar);
	}
	
	void draw() {
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		glTranslatef(0.f, 0.f, -2.f);
		
		glClear(GL_COLOR_BUFFER_BIT);
		glBegin(GL_QUADS);
		for (const auto& c : m->cells) {
			Prim W(c->U);
			
			if (displayMethod == DisplayMethod::STATE) {
				glColor3f(displayScalar * W.rho(), displayScalar * vec::length(W.v()), displayScalar * W.P());
			} else if (displayMethod == DisplayMethod::VOLUME) {
				glColor3f(displayScalar * c->volume, .5, 1. - displayScalar * c->volume);
			}
			
			for (const auto& v : c->vtxs) {
				glVertex2v(v->pos.v);
			}
		}
		glEnd();

		if (showVtxs || showCellCenters) {
			glColor3f(1,1,1);
			glPointSize(3);
			glBegin(GL_POINTS);
			if (showVtxs) {
				for (const auto& v : m->vtxs) {
					glVertex2v(v->pos.v);
				}
			}
			if (showCellCenters) {
				for (const auto& c : m->cells) {
					glVertex2v(c->pos.v);
				}
			}
			glEnd();
			glPointSize(1);
		}
		if (showEdges) {
			glBegin(GL_LINES);
			for (const auto& e : m->edges) {
				for (const auto& v : e->vtxs) {
					glVertex2v(v->pos.v);
				}
			}
			glEnd();
		}
	}

	std::pair<Cons, Cons> getEdgeStates(std::shared_ptr<Edge> e) {
		std::shared_ptr<Cell> cL, cR;
		if (e->cells.size() > 0) cL = e->cells[0];
		if (e->cells.size() > 1) cR = e->cells[1];
		if (cL && cR) {
			Cons UL, UR;
			UL = cL->U;
			UR = cR->U;
			return std::pair<Cons, Cons>{UL, UR};
		} else if (cL) {
assert(false);			
			Cons UL, UR;
			UL = UR = cL->U;
			vec m = UR.m();
			UR.m() = m - e->normal * (2 * vec::dot(e->normal, m));
			return std::pair<Cons, Cons>{UL, UR};
		} else if (cR) {
assert(false);			
			Cons UL, UR;
			UL = UR = cR->U;
			vec m = UL.m();
			UL.m() = m - e->normal * (2 * vec::dot(e->normal, m));
			return std::pair<Cons, Cons>{UL, UR};
		} 
		throw Common::Exception() << "here";
	}

	void step() {
		//for (auto& c : m->cells) {
		real result = parallel.reduce(
			m->cells.begin(), 
			m->cells.end(), 
			[this](std::shared_ptr<Cell> c) -> real {
				Prim W(c->U);
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
					real dx = sqrt(c->volume);
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
			[](std::shared_ptr<Edge> e) {
				e->flux = StateVec();
#if 0
for (auto& c : e->cells) {
	Prim W(c->U);
	Cons U(W);
	std::cout << (U - c->U) << std::endl;
}
#endif
			}
		);

		//for (auto& e : m->edges) {
		parallel.foreach(
			m->edges.begin(),
			m->edges.end(),
			[this, dt](std::shared_ptr<Edge> e) {
//if (e->cells.size() != 2) continue;

				std::pair<Cons, Cons> ULR = getEdgeStates(e);
				Cons UL = ULR.first;
				Cons UR = ULR.second;
assert(UL == e->cells[0]->U);
assert(UR == e->cells[1]->U);
assert(UL(3) == 0);
assert(UR(3) == 0);

for (int i = 0; i < StateVec::size; ++i) {
	assert(isfinite(UL(0)));	
	assert(isfinite(UR(0)));	
}	

				// 1) roe values at edge 
				// rotate to align edge normal to x axis
				// so x-direction flux jacobian is good for calculating the flux 
				UL.m() = rotateTo(UL.m(), e->normal);
				UR.m() = rotateTo(UR.m(), e->normal);

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
					
//std::cout << (UL - Cons(WL)) << std::endl;
//std::cout << (UR - Cons(WR)) << std::endl;
				
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
					real dx = e->cellDist;
					real epsilon = lambda * dt / dx;
					fluxTilde(j) = -.5 * lambda * dUTilde(j) * (sgnLambda + phi * (epsilon - sgnLambda));
				}
			
				Cons UAvg = (UR + UL) * .5;
				Cons UAvgTilde = matmul(&evL[0][0], UAvg);
				fluxTilde = fluxTilde + lambdas * UAvgTilde;
			
				Cons flux = matmul(&evR[0][0], fluxTilde);
				// here's the flux, aligned along the normal
			
				flux.m() = rotateFrom(flux.m(), e->normal);
			
for (int i = 0; i < StateVec::size; ++i) {
	assert(isfinite(flux(0)));	
}	
				
				// here's the flux in underlying coordinates
				e->flux = flux;
assert(e->flux(3) == 0);
//std::cout << flux << std::endl;
			}
		);

		//for (auto& c : m->cells) {
		parallel.foreach(
			m->cells.begin(),
			m->cells.end(),
			[dt](std::shared_ptr<Cell> c) {
				Cons dU_dt;
				for (auto& e : c->edges) {
					//if (e.hasFlux) {
						if (c == e->cells.front()) {
							dU_dt -= e->flux * (e->length / c->volume);
						} else {
							dU_dt += e->flux * (e->length / c->volume);
						}
					//}
				}
				c->U += dU_dt * dt;
	
Prim W(c->U);
assert(W.rho() > 0);
assert(vec::length(W.v()) > 0);
assert(W.P() > 0);
real hTotal = calc_hTotal(W.rho(), W.P(), c->U.ETotal());
assert(hTotal > 0);
Cons U2(W);
for (int i = 0; i < StateVec::size; ++i) {
	real delta = U2(i) - c->U(i);
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
