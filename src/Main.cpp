#include "CFDMesh/Equation/Euler.h"
#include "CFDMesh/Equation/GLMMaxwell.h"

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


//config for everything, to hold everything in one place
using real = double;
//using real = float;


using real2 = Tensor::Vector<real, 2>;
using real3 = Tensor::Vector<real, 3>;

//using ThisEquation = Equation::EulerNamespace<real>::Euler;
using ThisEquation = Equation::GLMMaxwellNamespace<real>::GLMMaxwell;

using WaveVec = ThisEquation::WaveVec;
using Cons = ThisEquation::Cons;

struct MeshConfig {
	using real = ::real;
	using real2 = ::real2;
	using real3 = ::real3;
	using Cons = ThisEquation::Cons;
};

using ThisMeshNamespace = CFDMesh::MeshNamespace<MeshConfig>;
using Mesh = ThisMeshNamespace::Mesh;
using MeshFactory = ThisMeshNamespace::MeshFactory;
using Cell = ThisMeshNamespace::Cell;
using Face = ThisMeshNamespace::Face;



static Parallel::Parallel parallel;
static ThisEquation eqn;

template<typename T> void glVertex2v(const T* v);
template<> void glVertex2v<double>(const double* v) { glVertex2dv(v); }
template<> void glVertex2v<float>(const float* v) { glVertex2fv(v); }


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
			mesh->vtxs[i].pos = real3(us[i], vs[i]);
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
				mesh->vtxs[i + m * j].pos = real3(f);
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
				mesh->vtxs[i + m * j].pos = real3(f);
			}
		}
		
		int capindex = m * n;
		if (capmin(0)) {
			real3 sum;
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

template<typename T> T cubed(const T& t) { return t * t * t; }

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
	using Parent = ::GLApp::ViewBehavior<::GLApp::GLApp>;
	
	std::shared_ptr<Mesh> m;
	std::function<Cons(real3)> initcond;

	double time = 0;
	bool running = false;
	bool singleStep = false;
	
	bool showVtxs = false;
	bool showEdges = false;
	bool showCellCenters = false;
	Cell* selectedCell = nullptr;

	bool displayAutomaticRange = true;
	int displayMethodIndex = 0;
	
	using ValueRange = std::pair<float, float>; //min, max
	ValueRange displayValueRange = ValueRange(0, 1);

	int initCondIndex = 0;

	//1 = mirror boundary, -1 = freeflow boundary
	float restitution = 1;

	float2 mousepos;
		
	float cfl = .5;

	int meshGenerationIndex = 1;
	
	std::shared_ptr<ImGuiCommon::ImGuiCommon> gui = std::make_shared<ImGuiCommon::ImGuiCommon>(window, context);

	virtual const char* getTitle() {
		return "CFD Mesh";
	}

	GLuint gradientTex = {};

	CFDMeshApp(const Init& args) : Parent(args) {
		view = viewOrtho;
		viewOrtho->zoom(0) = viewOrtho->zoom(1) = .5;

		glClearColor(.5, .75, .75, 1);

		std::vector<float4> gradientColors = {
			{1,0,0,1},
			{1,1,0,1},
			{0,1,0,1},
			{0,1,1,1},
			{0,0,1,1},
			{1,0,1,1},
		};

		int gradientTexWidth = 256;
		std::vector<uchar4> gradientTexData(gradientTexWidth );
		for (int i = 0; i < gradientTexWidth ; ++i) {
			float f = (float)(i+.5)/(float)gradientTexWidth * (1 - 1e-7);
			f *= (float)gradientColors.size();
			int ip = (int)floor(f);
			float fn = f - (float)ip;
			float fp = 1 - fn;
			int n1 = ip;
			int n2 = (n1 + 1) % gradientColors.size();
			const float4& c1 = gradientColors[n1];
			const float4& c2 = gradientColors[n2];
			float4 c = (c1 * fp + c2 * fn) * 255;
			for (int j = 0; j < 4; ++j) {
				c(j) = std::clamp<float>(c(j), 0, 255);
			}
			gradientTexData[i] = (uchar4)c;
		}
		glGenTextures(1, &gradientTex);
		glBindTexture(GL_TEXTURE_1D, gradientTex);
		glTexImage1D(GL_TEXTURE_1D, 0, GL_RGBA, gradientTexData.size(), 0, GL_RGBA, GL_UNSIGNED_BYTE, gradientTexData.data()->v);
		glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glBindTexture(GL_TEXTURE_1D, 0);

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

		auto ic = eqn.initConds[initCondIndex].get();
		parallel.foreach(
			m->cells.begin(),
			m->cells.end(),
			[this, ic](auto& c) {
				c.U = ic->initCell(&eqn, c.pos);
			}
		);
		
		refreshDisplayValues();
	}

	void draw() {
		glClear(GL_COLOR_BUFFER_BIT);

		glEnable(GL_TEXTURE_1D);
		glBindTexture(GL_TEXTURE_1D, gradientTex);

		for (const auto& c : m->cells) {
			real f = (c.displayValue - displayValueRange.first) / (displayValueRange.second - displayValueRange.first);
			glTexCoord1f(f);
			glBegin(GL_POLYGON);
			for (int vi = 0; vi < c.vtxCount; ++vi) {
				glVertex2v(m->vtxs[m->cellVtxIndexes[vi + c.vtxOffset]].pos.v);
			}
			glEnd();
		}
		
		glBindTexture(GL_TEXTURE_1D, 0);
		glDisable(GL_TEXTURE_1D);
		
		if (selectedCell) {
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
			glLineWidth(3);
			glColor3f(1,1,0);
			glBegin(GL_POLYGON);
			for (int vi = 0; vi < selectedCell->vtxCount; ++vi) {
				glVertex2v(m->vtxs[m->cellVtxIndexes[vi + selectedCell->vtxOffset]].pos.v);
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
			return std::pair<Cons, Cons>{cL->U, eqn.reflect(cL->U, e->normal, restitution)};
		} else if (cR) {
			return std::pair<Cons, Cons>{eqn.reflect(cR->U, e->normal, restitution), cR->U};
		} 
		throw Common::Exception() << "here";
	}

	real calcDT() {
		//for (auto& c : m->cells) {
		real result = parallel.reduce(
			m->cells.begin(), 
			m->cells.end(), 
			[this](Cell& c) -> real {
				
				real result = std::numeric_limits<real>::infinity();
				for (int ei = 0; ei < c.faceCount; ++ei) {
					Face* e = &m->faces[m->cellFaceIndexes[ei+c.faceOffset]];
					
					Cons U = eqn.rotateTo(c.U, e->normal);
					ThisEquation::CalcLambdaVars vars(eqn, U);
					
					real lambdaMin, lambdaMax;
					std::pair<real, real> lambdaMinMax = eqn.calcLambdaMinMax(vars);
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
		ThisEquation::Eigen vars = eqn.calcRoeAvg(UL, UR);

		WaveVec lambdas = eqn.getEigenvalues(vars);

		Cons dU = UR - UL;
		WaveVec dUTilde = eqn.apply_evL(dU, vars);
	
		WaveVec fluxTilde;
		for (int j = 0; j < ThisEquation::numWaves; ++j) {
			real lambda = lambdas(j);
			real phi = 0;
			real sgnLambda = lambda >= 0 ? 1 : -1;
			real epsilon = lambda * dt / dx;
			fluxTilde(j) = -.5 * lambda * dUTilde(j) * (sgnLambda + phi * (epsilon - sgnLambda));
		}
	
		Cons UAvg = (UR + UL) * .5;
		WaveVec UAvgTilde = eqn.apply_evL(UAvg, vars);
		fluxTilde = fluxTilde + lambdas * UAvgTilde;
	
		Cons flux = eqn.apply_evR(fluxTilde, vars);
		// here's the flux, aligned along the normal

		return flux;
	}

	Cons calcFluxHLL(Cons UL, Cons UR, real dx, real dt) {
		real lambdaMinL = eqn.calcLambdaMin(ThisEquation::CalcLambdaVars(eqn, UL));
		
		real lambdaMaxR = eqn.calcLambdaMax(ThisEquation::CalcLambdaVars(eqn, UR));
		
		ThisEquation::Eigen vars = eqn.calcRoeAvg(UL, UR);
		std::pair<real, real> lambdaMinMax = eqn.calcLambdaMinMax(ThisEquation::CalcLambdaVars(vars));
		real lambdaMin = lambdaMinMax.first;
		real lambdaMax = lambdaMinMax.second;

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
			for (int i = 0; i < Cons::size; ++i) {
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
				//roe averaged values at edge 
				std::pair<Cons, Cons> ULR = getEdgeStates(&e);
				
				//rotate to align edge normal to x axis
				//so x-direction flux jacobian is good for calculating the flux 
				Cons UL = eqn.rotateTo(ULR.first, e.normal);
				Cons UR = eqn.rotateTo(ULR.second, e.normal);

for (int i = 0; i < Cons::size; ++i) {
	assert(std::isfinite(UL(i)));	
	assert(std::isfinite(UR(i)));	
}	
			
				Cons F = (this->*calcFluxes[calcFluxIndex])(UL, UR, e.cellDist, dt);
		
				// rotate back to normal
				e.flux = eqn.rotateFrom(F, e.normal);

for (int i = 0; i < Cons::size; ++i) {
	assert(std::isfinite(e.flux(i)));
}		
			}
		);

		//for (auto& c : m->cells) {
		parallel.foreach(
			m->cells.begin(),
			m->cells.end(),
			[this, dt](Cell& c) {
				Cons dU_dt;
				//for (int ei : c.faces) {
				//	Face* e = &m->faces[ei];
				for (int ei = 0; ei < c.faceCount; ++ei) {
					Face* e = &m->faces[m->cellFaceIndexes[ei+c.faceOffset]];
					//if (e.hasFlux) {
						if (&c == &m->cells[e->cells[0]]) {
							dU_dt -= e->flux * (e->length / c.volume);
						} else {
							dU_dt += e->flux * (e->length / c.volume);
						}
					//}
				}
				c.U += dU_dt * dt;
			}
		);

		refreshDisplayValues();

		time += dt;
	}
	
	void refreshDisplayValues() {
		parallel.foreach(
			m->cells.begin(),
			m->cells.end(),
			[this](Cell& c) {
				c.displayValue = eqn.displayMethods[displayMethodIndex]->f(&eqn, c.U);
			}
		);

		if (displayAutomaticRange) {
			displayValueRange = parallel.reduce(
				m->cells.begin(),
				m->cells.end(),
				[](const Cell& c) -> ValueRange { 
					return ValueRange(c.displayValue, c.displayValue);
				},
				ValueRange(INFINITY, -INFINITY),
				[](ValueRange a, ValueRange b) -> ValueRange {
					return ValueRange(
						std::min<float>(a.first, b.first),
						std::max<float>(a.second, b.second)
					);
				}
			);
		}
	}
	
	virtual void onUpdate() {
		Parent::onUpdate();
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
		
			igSeparator();
	
			if (igCombo("display method", &displayMethodIndex, eqn.displayMethodNames.data(), eqn.displayMethodNames.size(), -1)) {
				refreshDisplayValues();
			}
		
			igCheckbox("auto display range", &displayAutomaticRange);
			igInputFloat("display range min", &displayValueRange.first, .1, 1, "%f", 0);
			igInputFloat("display range max", &displayValueRange.second, .1, 1, "%f", 0);
			
			igSeparator();
			
			igInputFloat("restitution", &restitution, .1, 1., "%f", 0);
			
			igSeparator();
		
			igCombo("flux", &calcFluxIndex, calcFluxNames.data(), calcFluxNames.size(), -1);
			
			igSeparator();
			
			eqn.updateGUI();
			
			igSeparator();
			
			igCombo("init cond", &initCondIndex, eqn.initCondNames.data(), eqn.initCondNames.size(), -1);
			igPushIDStr("init cond fields");
			eqn.initConds[initCondIndex]->updateGUI();
			igPopID();
			
			if (igSmallButton("reset state")) {
				resetState();
			}
			
			igSeparator();
			
			if (igCombo("mesh", &meshGenerationIndex, meshGenerationNames.data(), meshGenerationNames.size(), -1)) {
				resetMesh();
			}
			igPushIDStr("mesh generation fields");
			meshGenerators[meshGenerationIndex]->updateGUI();
			igPopID();
			if (igSmallButton("reset mesh")) {
				resetMesh();
			}
			
			igSeparator();

			if (view == viewOrtho) {
				//find the view bounds
				//find the mouse position
				//find any cells at that position
				real2 pos2 = (real2)((mousepos - .5f) * float2(1, -1) * float2(aspectRatio, 1) / float2(viewOrtho->zoom(0), viewOrtho->zoom(1)) + viewOrtho->pos);
				std::function<real(int)> f = [&](int i) -> real { return i >= 2 ? 0 : pos2(i); };
				real3 pos = real3(f);
					
				bool canHandleMouse = !igGetIO()->WantCaptureMouse;
				
				if (canHandleMouse) {
					igBeginTooltip();
					igText("%f %f\n", pos(0), pos(1));
					igEndTooltip();
				}

				selectedCell = nullptr;
				int selectedCellIndex = -1;
				for (int i = 0; i < (int)m->cells.size(); ++i) {
					Cell* c = &m->cells[i];
					if (ThisMeshNamespace::contains(
						pos, 
						map<
							std::vector<int>, 
							std::vector<real3>
						>(
							//c->vtxs,
							//TODO change function to use an iterator
							std::vector<int>(m->cellVtxIndexes.begin() + c->vtxOffset, m->cellVtxIndexes.begin() + c->vtxOffset + c->vtxCount), 
							
							[this](int vi) -> real3 { 
								return m->vtxs[vi].pos; 
							}
						))) {
						selectedCellIndex = i;
						break;
					}
				}

				if (selectedCellIndex != -1) {
					Cell* c = &m->cells[selectedCellIndex];
					selectedCell = c;
					
					if (canHandleMouse) {
						igBeginTooltip();
						igText("%f", c->displayValue);
						igEndTooltip();
					}

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
			Parent::onSDLEvent(event);
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
