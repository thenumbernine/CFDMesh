#include "CFDMesh/Equation/EulerEquation.h"
#include "CFDMesh/Mesh/Mesh.h"
#include "CFDMesh/Util.h"
#include "CFDMesh/Vector.h"

#include "GLApp/gl.h"
#include "GLApp/GLApp.h"
#include "GLApp/ViewBehavior.h"

#include "ImGuiCommon/ImGuiCommon.h"
#include "Parallel/Parallel.h"
#include "Tensor/Tensor.h"
#include "Common/File.h"

#include <algorithm>
#include <memory>
#include <iostream>
#include <fstream>
#include <vector>
#include <list>

#include <cassert>

using namespace CFDMesh;
using namespace CFDMesh::Equation;


//config for everything, to hold everything in one place
template<typename real_, int vecdim_>
struct ConfigTemplate {
	using real = real_;
	enum { vecdim = vecdim_ };
	
	using real2 = Tensor::Vector<real, 2>;
	
	using vec = Tensor::Vector<real, vecdim>;	//n-dimensional vector of reals
	using StateVec = Tensor::Vector<real, vecdim + 2>;
};


//using Config = ConfigTemplate<double, 2>;
using Config = ConfigTemplate<double, 3>;


using real = Config::real;
constexpr int vecdim = Config::vecdim;
using real2 = Config::real2;
using vec = Config::vec;
using StateVec = Config::StateVec;



using ThisEquation = EulerEquation<Config>;
using Cons = ThisEquation::Cons;
using Prim = ThisEquation::Prim;

struct MeshConfig : public Config {
	using Cons = ThisEquation::Cons;
};

using ThisMeshNamespace = CFDMesh::MeshNamespace<MeshConfig>;
using Mesh = ThisMeshNamespace::Mesh;
using MeshFactory = ThisMeshNamespace::MeshFactory;
using Cell = ThisMeshNamespace::Cell;
using Face = ThisMeshNamespace::Face;



static Parallel::Parallel parallel;
static ThisEquation eqn;

// rotate vx,vy such that n now points along the x dir
static vec rotateTo(vec v, vec n) {
	return vec(
		v(0) * n(0) + v(1) * n(1),
		v(1) * n(0) - v(0) * n(1)
	);
}

// rotate vx,vy such that the x dir now points along n 
static vec rotateFrom(vec v, vec n) {
	return vec(
		v(0) * n(0) - v(1) * n(1),
		v(1) * n(0) + v(0) * n(1)
	);
}

template<typename T> T cubed(const T& t) { return t * t * t; }

template<typename T> void glVertex2v(const T* v);
template<> void glVertex2v<double>(const double* v) { glVertex2dv(v); }
template<> void glVertex2v<float>(const float* v) { glVertex2fv(v); }

enum DisplayMethod {
	STATE,
	VOLUME,
	COUNT
};

static const char* displayMethodNames[DisplayMethod::COUNT] = {
	"state",
	"volume",
};

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
	decltype(initConds),
	std::vector<const char*>
>(
	initConds,
	[](std::shared_ptr<InitCond> ic) -> const char* { return ic->name(); }
);


struct FileMeshFactory : public MeshFactory {
	std::string filename = {"grids/n0012_113-33.p2dfmt"};
	
	FileMeshFactory() : MeshFactory("p2dfmt mesh") {}

	virtual std::shared_ptr<Mesh> createMesh() const {
		std::shared_ptr<Mesh> mesh = createMeshSuper();
		
		std::list<std::string> ls = split<std::list<std::string>>(Common::File::read(filename), "\n");
	
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

		mesh->vtxs.resize(m*n);
		for (int i = 0; i < (int)us.size(); ++i) {
			mesh->vtxs[i].pos = vec(us[i], vs[i]);
		}
	
		for (int j = 0; j < n-1; ++j) {
			for (int i = 0; i < m-1; ++i) {
				mesh->addCell(std::vector<int>{i + m * j, i + m * (j+1), i+1 + m * (j+1), i+1 + m * j});
			}
		}
	
		mesh->calcAux();
		return mesh;
	}

	virtual void updateGUI() {
		//TODO filename popup
	}
};

struct ChartMeshFactory : public MeshFactory {
	int2 size = int2(101, 101);
	float2 mins = real2(-1, -1);
	float2 maxs = real2(1, 1);
	int2 repeat = int2(0, 0);
	int2 capmin = int2(0, 0);

	ChartMeshFactory(const char* name_) : MeshFactory(name_) {}

	virtual real2 grid(real2 x) const { return x; }
	
	virtual void updateGUI() {
		igInputInt2("size", size.v, 0);
		igInputFloat2("mins", mins.v, "%f", 0);
		igInputFloat2("maxs", maxs.v, "%f", 0);
		igInputInt2("repeat", repeat.v, 0);
		igInputInt2("capmin", capmin.v, 0);
	}
};

struct TriUnitMeshFactory : public ChartMeshFactory {
	TriUnitMeshFactory() : ChartMeshFactory("unit square of triangles") {}
	
	virtual std::shared_ptr<Mesh> createMesh() const {
		std::shared_ptr<Mesh> mesh = createMeshSuper();
		
		int m = size(0);
		int n = size(1);
	
		mesh->vtxs.resize(m * n);
		for (int j = 0; j < n; ++j) {
			for (int i = 0; i < m; ++i) {
				real2 x = real2(
					((real)i + .5) / (real)size(0) * (maxs(0) - mins(0)) + mins(0),
					((real)j + .5) / (real)size(1) * (maxs(1) - mins(1)) + mins(1));
				
				real2 u = grid(x);
				std::function<real(int)> f = [&u](int i) -> real { return i < real2::size ? u(i) : 0.; };
				mesh->vtxs[i + m * j].pos = vec(f);
			}
		}
		
		int imax = repeat(0) ? m : m-1;
		int jmax = repeat(1) ? n : n-1;
		for (int j = 0; j < jmax; ++j) {
			int jn = (j + 1) % n;
			for (int i = 0; i < imax; ++i) {
				int in = (i + 1) % m;
				mesh->addCell(std::vector<int>{i + m * j, in + m * j, in + m * jn});
				mesh->addCell(std::vector<int>{in + m * jn, i + m * jn, i + m * j});
			}
		}
	
		mesh->calcAux();
		return mesh;
	}
};

struct QuadUnitMeshFactory : public ChartMeshFactory {
	QuadUnitMeshFactory(const char* name_ = "unit square of quads") : ChartMeshFactory(name_) {}
	
	virtual std::shared_ptr<Mesh> createMesh() const {
		std::shared_ptr<Mesh> mesh = createMeshSuper();

		int m = size(0);
		int n = size(1);

		int vtxsize = m * n;
		if (capmin(0)) vtxsize++;
		mesh->vtxs.resize(vtxsize);
		for (int j = 0; j < n; ++j) {
			for (int i = 0; i < m; ++i) {
				real2 x = real2(
					((real)i + .5) / (real)size(0) * (maxs(0) - mins(0)) + mins(0),
					((real)j + .5) / (real)size(1) * (maxs(1) - mins(1)) + mins(1));
				
				real2 u = grid(x);
				std::function<real(int)> f = [&u](int i) -> real { return i < real2::size ? u(i) : 0.; };
				mesh->vtxs[i + m * j].pos = vec(f);
			}
		}
		
		int capindex = m * n;
		if (capmin(0)) {
			vec sum;
			for (int j = 0; j < n; ++j) {
				sum += mesh->vtxs[0 + m * j].pos;
			}
			mesh->vtxs[capindex].pos = sum / (real)n;
		}

		int imax = repeat(0) ? m : m-1;
		int jmax = repeat(1) ? n : n-1;
		for (int j = 0; j < jmax; ++j) {
			int jn = (j + 1) % n;
			for (int i = 0; i < imax; ++i) {
				int in = (i + 1) % m;
				mesh->addCell(std::vector<int>{i + m * j, in + m * j, in + m * jn, i + m * jn});
			}
		}

		if (capmin(0)) {
			for (int j = 0; j < jmax; ++j) {
				int jn = (j + 1) % n;
				mesh->addCell(std::vector<int>{ 0 + m * j, 0 + m * jn, capindex });
			}
		}

		mesh->calcAux();
		return mesh;
	}
};

struct QuadUnitCbrtMeshFactory : public QuadUnitMeshFactory {
	QuadUnitCbrtMeshFactory() : QuadUnitMeshFactory("unit square of quads, cbrt mapping") {}
	virtual real2 grid(real2 v) const {
		return real2(cbrt(v(0)), cbrt(v(1)));
	}
};

struct QuadUnitCubedMeshFactory : public QuadUnitMeshFactory {
	QuadUnitCubedMeshFactory() : QuadUnitMeshFactory("unit square of quads, cubed mapping") {}
	virtual real2 grid(real2 v) const {
		return real2(cubed(v(0)), cubed(v(1)));
	}
};

struct TwistQuadUnitMeshFactory : public QuadUnitMeshFactory {
	TwistQuadUnitMeshFactory() : QuadUnitMeshFactory("unit square of quads, twist in the middle") {}
	virtual real2 grid(real2 v) const {
		real r = real2::length(v);
		//real theta = std::max(0., 1. - r);
		real sigma = 3.;	//almost 0 at r=1
		const real rotationAmplitude = 3.;
		real theta = rotationAmplitude*sigma*r*exp(-sigma*sigma*r*r);
		real costh = cos(theta), sinth = sin(theta);
		return real2(
			costh * v(0) - sinth * v(1),
			sinth * v(0) + costh * v(1));	
	}
};

struct DonutQuadUnitMeshFactory : public QuadUnitMeshFactory {
	DonutQuadUnitMeshFactory() : QuadUnitMeshFactory("polar") {
		size = int2(51, 201);
		mins = real2(.5, 0);
		maxs = real2(1, 2*M_PI);
		repeat = int2(0, 1);
		//capmin = int2(1, 0);
	}
	virtual real2 grid(real2 v) const {
		return real2(cos(v(1)), sin(v(1))) * v(0);
	}
};

std::vector<std::shared_ptr<MeshFactory>> meshGenerators = {
	std::make_shared<FileMeshFactory>(),
	std::make_shared<QuadUnitMeshFactory>(),
	std::make_shared<TriUnitMeshFactory>(),
	std::make_shared<QuadUnitCbrtMeshFactory>(),
	std::make_shared<QuadUnitCubedMeshFactory>(),
	std::make_shared<TwistQuadUnitMeshFactory>(),
	std::make_shared<DonutQuadUnitMeshFactory>(),
};

std::vector<const char*> meshGenerationNames = map<
	decltype(meshGenerators),
	std::vector<const char*>
>(
	meshGenerators,
	[](const std::shared_ptr<MeshFactory>& m) -> const char* { return m->name; }
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

	int meshGenerationIndex = 1;
	
	std::shared_ptr<ImGuiCommon::ImGuiCommon> gui = std::make_shared<ImGuiCommon::ImGuiCommon>(window, context);

	virtual const char* getTitle() {
		return "CFD Mesh";
	}
	
	CFDMeshApp(const Init& args) : Super(args) {
		
		view = viewOrtho;
		viewOrtho->zoom(0) = viewOrtho->zoom(1) = .5;

		glClearColor(.5, .75, .75, 1);
		
		resetMesh();	//which calls resetState()
	}

	void resetMesh() {
		m = meshGenerators[meshGenerationIndex]->createMesh();
		resetState();
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
				real hTotal = eqn.calc_hTotal(W.rho(), W.P(), c.U.ETotal());
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
			for (const auto& e : m->faces) {
				for (int vi : e.vtxs) {
					glVertex2v(m->vtxs[vi].pos.v);
				}
			}
			glEnd();
		}
	}

	std::pair<Cons, Cons> getEdgeStates(const Face* e) {
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
				for (int ei : c.faces) {
					Face* e = &m->faces[ei];
					vec normal = e->normal;
#endif
					real lambdaMin, lambdaMax;
					std::pair<real, real> lambdaMinMax = eqn.calcLambdaMinMax(normal, W, Cs);
					lambdaMin = lambdaMinMax.first;
					lambdaMax = lambdaMinMax.second;
					lambdaMin = std::min<real>(0, lambdaMin);
					lambdaMax = std::max<real>(0, lambdaMax);
					// TODO a better way to do this.  maybe use faces' lambdas?  maybe do this after calculating the eigenbasis?
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

	Cons calcFluxRoe(Cons UL, Cons UR, real dx, real dt) {
		static_assert(StateVec::size == 5);
		
		ThisEquation::RoeAvg vars = eqn.calcRoeAvg(UL, UR);

		StateVec lambdas = eqn.getEigenvalues(vars);

		Cons dU = UR - UL;
		Cons dUTilde = eqn.apply_evL(dU, vars);
	
		Cons fluxTilde;
		for (int j = 0; j < 5; ++j) {
			real lambda = lambdas(j);
			real phi = 0;
			real sgnLambda = lambda >= 0 ? 1 : -1;
			real epsilon = lambda * dt / dx;
			fluxTilde(j) = -.5 * lambda * dUTilde(j) * (sgnLambda + phi * (epsilon - sgnLambda));
		}
	
		Cons UAvg = (UR + UL) * .5;
		Cons UAvgTilde = eqn.apply_evL(UAvg, vars);
		fluxTilde = fluxTilde + lambdas * UAvgTilde;
	
		Cons flux = eqn.apply_evR(fluxTilde, vars);
		// here's the flux, aligned along the normal

		return flux;
	}

	Cons calcFluxHLL(Cons UL, Cons UR, real dx, real dt) {
		ThisEquation::RoeAvg vars = eqn.calcRoeAvg(UL, UR);

		real CsL = eqn.calc_Cs_from_P_rho(vars.WL.P(), vars.WL.rho());
		real CsR = eqn.calc_Cs_from_P_rho(vars.WR.P(), vars.WR.rho());

		real lambdaMinL = eqn.calcLambdaMin(vars.WL.v()(0), CsL);
		real lambdaMaxR = eqn.calcLambdaMax(vars.WR.v()(0), CsR);
		real lambdaMin = eqn.calcLambdaMin(vars.v(0), vars.Cs);
		real lambdaMax = eqn.calcLambdaMax(vars.v(0), vars.Cs);
#if 0	//Davis direct
		real sl = lambdaMinL;
		real sr = lambdaMaxR;
#endif
#if 1	//Davis direct bounded
		real sl = std::min(lambdaMinL, lambdaMin);
		real sr = std::max(lambdaMaxR, lambdaMax);
#endif
		
		Cons flux;
		if (0. <= sl) {
			Cons fluxL = eqn.calcFluxFromCons(UL);
			flux = fluxL;
		} else if (sl <= 0. && 0. <= sr) {
			Cons fluxL = eqn.calcFluxFromCons(UL);
			Cons fluxR = eqn.calcFluxFromCons(UR);
			//(sr * fluxL[j] - sl * fluxR[j] + sl * sr * (UR[j] - UL[j])) / (sr - sl)
			real invDenom = 1. / (sr - sl);
			for (int i = 0; i < StateVec::size; ++i) {
				flux(i) = (sr * fluxL(i) - sl * fluxR(i) + sl * sr * (UR(i) - UL(i))) * invDenom; 
			}
		} else if (sr <= 0.) {
			Cons fluxR = eqn.calcFluxFromCons(UR);
			flux = fluxR;
		}
	
		return flux;
	}


	using CalcFluxFunc = Cons (CFDMeshApp::*)(Cons, Cons, real, real);
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

		//for (auto& e : m->faces) {
		parallel.foreach(
			m->faces.begin(),
			m->faces.end(),
			[this, dt](Face& e) {
				std::pair<Cons, Cons> ULR = getEdgeStates(&e);
				Cons UL = ULR.first;
				Cons UR = ULR.second;
assert(e.cells[0] == -1 || UL == m->cells[e.cells[0]].U);
assert(e.cells[1] == -1 || UR == m->cells[e.cells[1]].U);
assert(UL(3) == 0);
assert(UR(3) == 0);

for (int i = 0; i < StateVec::size; ++i) {
	assert(std::isfinite(UL(0)));	
	assert(std::isfinite(UR(0)));	
}	

				// roe values at edge 
				// rotate to align edge normal to x axis
				// so x-direction flux jacobian is good for calculating the flux 
				UL.m() = rotateTo(UL.m(), e.normal);
				UR.m() = rotateTo(UR.m(), e.normal);
				
				real dx = e.cellDist;
				Cons flux = (this->*calcFluxes[calcFluxIndex])(UL, UR, dx, dt);
		
				// rotate back to normal
				flux.m() = rotateFrom(flux.m(), e.normal);
			
for (int i = 0; i < StateVec::size; ++i) {
	assert(std::isfinite(flux(0)));	
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
				for (int ei : c.faces) {
					Face* e = &m->faces[ei];
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

			if (igCombo("mesh", &meshGenerationIndex, meshGenerationNames.data(), meshGenerationNames.size(), -1)) {
				resetMesh();
			}
			
			igPushIDStr("mesh generation fields");
			meshGenerators[meshGenerationIndex]->updateGUI();
			igPopID();
			
			if (igSmallButton("reset mesh")) {
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
					if (ThisMeshNamespace::contains(pos, map<std::vector<int>, std::vector<vec>>(c->vtxs, [this](int vi) -> vec { return m->vtxs[vi].pos; }))) {
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
						for (int ei : c->faces) {
							Face* e = &m->faces[ei];
							if (e->removeCell(i - m->cells.begin())) {
								//all cells on e are removed
								toRemoveEdges.push_back(ei);
							}
						}
						m->cells.erase(m->cells.begin() + selectedCellIndex);
					
						//sort largest -> smallest
						std::sort(toRemoveEdges.begin(), toRemoveEdges.end(), [](int a, int b) -> bool { return a > b; });
						for (int ei : toRemoveEdges) {
							m->faces.erase(m->faces.begin() + ei);
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
