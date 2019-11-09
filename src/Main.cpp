#include "GLApp/gl.h"
#include "ImGuiCommon/ImGuiCommon.h"
#include "GLApp/GLApp.h"

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
	for (const std::string& s : v) {
		if (!first) result += sep;
		result += s;
		first = false;
	}
	return result;
}

template<typename From, typename To>
To map(const From& from, std::function<typename To::value_type(typename From::value_type)> f) {
	To to;
	for (const typename From::value_type& v : from) {
		to.push_back(f(v));
	}
	return to;
}

using real = double;
//constexpr int vecdim = 2;
constexpr int vecdim = 3;
using vecnr = Tensor::Vector<real, vecdim>;
constexpr int statedim = vecdim + 2;
using StateVec = Tensor::Vector<real, statedim>;

template<typename T>
typename T::value_type sum(const T& t) {
	return std::accumulate(t.begin(), t.end(), typename T::value_type());
}


real calc_hTotal(real rho, real P, real ETotal) {
	return (ETotal + P) / rho;
}

const real heatCapacityRatio = 1.4;
real calcSpeedOfSound(vecnr v, real hTotal) {
	return sqrt((heatCapacityRatio - 1) * (hTotal - .5 * vecnr::lenSq(v)));
}


struct Cons;
struct Prim;

struct Cons : public StateVec {
	real& rho() { return StateVec::v[0]; }
	vecnr& m() { return *(vecnr*)( StateVec::v + 1 ); }
	real& ETotal() { return StateVec::v[statedim-1]; }
	
	const real& rho() const { return StateVec::v[0]; }
	const vecnr& m() const { return *(vecnr*)( StateVec::v + 1 ); }
	const real& ETotal() const { return StateVec::v[statedim-1]; }

	Cons() {}
	
	Cons(const StateVec& v) : StateVec(v) {}

	Cons(real rho_, vecnr m_, real ETotal_) {
		rho() = rho_;
		m() = m_;
		ETotal() = ETotal_;
	}
	
	Cons(const Prim& prim);
};

struct Prim : public StateVec {
	real& rho() { return StateVec::v[0]; }
	vecnr& v() { return *(vecnr*)( StateVec::v + 1 ); }
	real& P() { return StateVec::v[statedim-1]; }

	const real& rho() const { return StateVec::v[0]; }
	const vecnr& v() const { return *(vecnr*)( StateVec::v + 1 ); }
	const real& P() const { return StateVec::v[statedim-1]; }

	Prim() {}
	
	Prim(const StateVec& v) : StateVec(v) {}
	
	Prim(real rho_, vecnr v_, real P_) {
		rho() = rho_;
		v() = v_;
		P() = P_;
	}
	
	Prim(const Cons& U);
};

Cons::Cons(const Prim& W) {
	rho() = W.rho();
	m() = W.v() * rho();
	ETotal() = W.P() / (heatCapacityRatio - 1.) + .5 * rho() * vecnr::lenSq(W.v());
}

Prim::Prim(const Cons& U) {
	rho() = U.rho();
	v() = U.m() / rho();
	P() = (heatCapacityRatio - 1.) * (U.ETotal() - .5 * rho() * vecnr::lenSq(v()));
}


struct Vertex;
struct Edge;
struct Cell;

struct Vertex {
	vecnr pos;
	std::vector<std::shared_ptr<Edge>> edges;
	Vertex(vecnr pos_) : pos(pos_) {}
};

struct Edge {
	vecnr pos;
	vecnr delta;
	vecnr normal;
	real length;
	real cellDist;
	StateVec flux;
	std::vector<std::shared_ptr<Vertex>> vtxs;
	std::vector<std::shared_ptr<Cell>> cells;
	Edge(std::vector<std::shared_ptr<Vertex>> vtxs_) : length(0), cellDist(0), vtxs(vtxs_) {}
};

struct Cell {
	vecnr pos;
	real volume;
	Cons U;
	std::vector<std::shared_ptr<Edge>> edges;
	std::vector<std::shared_ptr<Vertex>> vtxs;
	Cell() : volume(0) {}
};

//2D polygon volume
real polyVol(const std::vector<vecnr>& vs) {
	size_t n = vs.size();
	real v = 0;
	for (size_t i = 0; i < n; ++i) {
		vecnr pi = vs[i];
		vecnr pj = vs[(i+1)%n];
		v += .5 * (pi(0) * pj(1) - pi(1) * pj(0));
	}
	return fabs(v);
}

struct Mesh {
	std::vector<std::shared_ptr<Vertex>> vtxs;
	std::vector<std::shared_ptr<Edge>> edges;
	std::vector<std::shared_ptr<Cell>> cells;

	void addVtx(vecnr x) {
		vtxs.push_back(std::make_shared<Vertex>(x));
	}

	std::shared_ptr<Edge> addEdge(int i, int j) {
		std::shared_ptr<Vertex> a = vtxs[i];
		std::shared_ptr<Vertex> b = vtxs[j];
		for (std::shared_ptr<Edge> e : edges) {
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

	std::shared_ptr<Cell> addCell(std::vector<int> indexes) {
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
			std::vector<vecnr>
		>(c->vtxs, [](std::shared_ptr<Vertex> v) { return v->pos; })) * (1. / (real)n);
		for (auto& e : c->edges) {
			e->cells.push_back(c);
		}
		return c;
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
			addVtx(vecnr(us[i], vs[i]));
		}
		assert(vtxs.size() == (size_t)(m*n));
	
		for (int i = 0; i < n-1; ++i) {
			for (int j = 0; j < m-1; ++j) {
				std::shared_ptr<Cell> c = addCell(std::vector<int>{
					j + m * i,
					j + m * (i+1),
					j+1 + m * (i+1),
					j+1 + m * i
				});
			}
		}
	
		for (auto& e : edges) {
			{
				auto a = e->vtxs[0];
				auto b = e->vtxs[1];
				e->pos = (a->pos + b->pos) * .5;
				e->delta = a->pos - b->pos;
				e->length = vecnr::length(e->delta);
				e->normal = vecnr(-e->delta(1), e->delta(0));
			}

			{
				std::shared_ptr<Cell> a, b;
				if (e->cells.size() > 0) a = e->cells[0];
				if (e->cells.size() > 1) a = e->cells[1];
				if (a && b) {
					if (vecnr::dot(a->pos - b->pos, e->normal) < 0) {
						std::swap(a, b);
						e->cells = std::vector<std::shared_ptr<Cell>>{a, b};
					}
					e->cellDist = vecnr::length(b->pos - a->pos);
				} else {
					e->cellDist = vecnr::length(a->pos - e->pos) * 2.;
				}
			}
		}
	
		for (auto& c : cells) {
			c->volume = polyVol(map<
					std::vector<std::shared_ptr<Vertex>>,
					std::vector<vecnr>
				>(c->vtxs, [](const std::shared_ptr<Vertex>& v) -> vecnr {
					return v->pos;
				}));
		}
	
		resetState();
	}

	void resetState() {
		for (auto& c : cells) {	
			c->U = Cons(Prim(1., vecnr(.1), 1.));
		}
	}
};

// rotate vx,vy such that n now points along the x dir
vecnr rotateTo(vecnr v, vecnr n) {
	return vecnr(
		v(0) * n(0) + v(1) * n(1),
		v(1) * n(0) - v(0) * n(1)
	);
}

// rotate vx,vy such that the x dir now points along n 
vecnr rotateFrom(vecnr v, vecnr n) {
	return vecnr(
		v(0) * n(0) - v(1) * n(1),
		v(1) * n(0) + v(0) * n(1)
	);
}

template<typename T> void glVertex2v(T* v);
template<> void glVertex2v<double>(double* v) { glVertex2dv(v); }
template<> void glVertex2v<float>(float* v) { glVertex2fv(v); }

StateVec matmul(const real* A, const StateVec& x) {
	StateVec y;
	for (int i = 0; i < statedim; ++i) {
		real sum = 0;
		for (int j = 0; j < statedim; ++j) {
			//C layout, so row-major
			sum += A[j + statedim * i] * x(j);
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

	virtual const char* getTitle() {
		return "CFD Mesh";
	}
	
	virtual void init() {
		Super::init();
		glClearColor(.5, .75, .75, 1.);

		gui = std::make_shared<ImGuiCommon::ImGuiCommon>(window, context);
		m = std::make_shared<Mesh>("grids/n0012_113-33.p2dfmt");
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
	
	bool running = false;
	bool singleStep = false;
	bool showVtxs = true;
	bool showCellCenters = false;
	bool showEdges = true;

	int displayMethod = 0;
	float displayScalar = 1.;

	void draw() {
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		glTranslatef(0.f, 0.f, -2.f);
		
		glClear(GL_COLOR_BUFFER_BIT);
		glBegin(GL_QUADS);
		for (const auto& c : m->cells) {
			Prim W(c->U);
			
			if (displayMethod == DisplayMethod::STATE) {
				glColor3f(displayScalar * W.rho(), displayScalar * vecnr::length(W.v()), displayScalar * W.P());
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
		if (e->cells.size() > 1) cL = e->cells[1];
		Cons UL, UR;
		if (cL && cR) {
			UL = cL->U;
			UR = cR->U;
		} else if (cL) {
			UL = UR = cL->U;
			vecnr m = UR.m();
			UR.m() = m - e->normal * (2 * vecnr::dot(e->normal, m));
		} else if (cR) {
			UL = UR = cR->U;
			vecnr m = UL.m();
			UL.m() = m - e->normal * (2 * vecnr::dot(e->normal, m));
		} else {
			throw Common::Exception() << "here";
		}
		return std::pair<Cons, Cons>{UL, UR};
	}

	void step() {
		// calculate dt
		real result = std::numeric_limits<real>::infinity();
		real cfl = .5;
		for (auto& c : m->cells) {
			Prim W(c->U);
			real rho = W.rho();
			real P = W.P();
			real Cs = sqrt(heatCapacityRatio * P / rho);
			for (int j = 0; j < vecnr::size; ++j) {
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
		}
		real dt = result * cfl;
		
		for (auto& e : m->edges) {
			std::pair<Cons, Cons> ULR = getEdgeStates(e);
			Cons UL = ULR.first;
			Cons UR = ULR.second;

			// 1) roe values at edge 
			// rotate to align edge normal to x axis
			// so x-direction flux jacobian is good for calculating the flux 
			UL.m() = rotateTo(UL.m(), e->normal);
			UR.m() = rotateTo(UR.m(), e->normal);

			real ETotalL = UL.ETotal();
			Prim WL(UL);
			real rhoL = WL.rho();
			vecnr vL = WL.v();
			real PL = WL.P();
			real hTotalL = calc_hTotal(rhoL, PL, ETotalL);
			real sqrtRhoL = sqrt(rhoL);
			
			real ETotalR = UR.ETotal();
			Prim WR(UR);
			real rhoR = WR.rho();
			vecnr vR = WR.v();
			real PR = WR.P();
			real hTotalR = calc_hTotal(rhoR, PR, ETotalR);
			real sqrtRhoR = sqrt(rhoR);
			
			vecnr v = (vL * sqrtRhoL + vR * sqrtRhoR) / (sqrtRhoL + sqrtRhoR);
			real vx = v(0), vy = v(1), vz = v(2);
			real hTotal = (sqrtRhoL * hTotalL + sqrtRhoR * hTotalR) / (sqrtRhoL + sqrtRhoR);
			real Cs = calcSpeedOfSound(v, hTotal);
			real CsSq = Cs * Cs;
			//e.roe = matrix{rho, vx, vy, vz, hTotal, Cs}
		
			// 2) eigenbasis at interface
		
			real vSq = vecnr::lenSq(v);

			static_assert(statedim == 5);

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
			// here's the flux in underlying coordinates
			e->flux = flux;
		}

		for (auto& c : m->cells) {
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
		}
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
				m->resetState();
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
		gui->sdlEvent(event);
	}
};

GLAPP_MAIN(CFDMeshApp)
