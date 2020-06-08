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

#include "Common/Exception.h"
#include "Common/File.h"

#include <algorithm>
#include <memory>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <list>



namespace CFDMesh {

//config for everything, to hold everything in one place
using real = double;
//using real = float;

using real2 = Tensor::Vector<real, 2>;
using real3 = Tensor::Vector<real, 3>;

static Parallel::Parallel parallel;//(1);

struct CFDMeshApp;

struct ISimulation {
	//written by app SDL events:
	bool running = true;
	bool singleStep = false;
	float2 mousepos;

	CFDMeshApp* app = {};

	ISimulation(CFDMeshApp* app_) : app(app_) {}

	virtual void draw() = 0;
	virtual void updateGUI() = 0;
	virtual void onUpdate() = 0;
	virtual void resetState() = 0;
};

template<typename real, int dim, typename ThisEquation>
struct Simulation : public ISimulation {
	using This = Simulation;
	using Super = ISimulation;

	using WaveVec = typename ThisEquation::WaveVec;
	using Cons = typename ThisEquation::Cons;

	using ThisMeshNamespace = CFDMesh::MeshNamespace<real, dim, Cons>;
	using Mesh = typename ThisMeshNamespace::Mesh;
	using MeshFactory = typename CFDMesh::MeshFactory<real, dim, Cons>;
	using Cell = typename ThisMeshNamespace::Cell;
	using Face = typename ThisMeshNamespace::Face;
	
	using DisplayMethod = typename ThisEquation::DisplayMethod;


	std::vector<std::shared_ptr<MeshFactory>> meshGenerators;
	std::vector<const char*> meshGenerationNames;

	std::shared_ptr<Mesh> m;
	std::function<Cons(real3)> initcond;

	ThisEquation eqn;

	double time = 0;

	//passed to mesh.draw
	typename Mesh::DrawArgs drawArgs;

	bool displayAutomaticRange = true;
	int displayMethodIndex = 0;
	
	int initCondIndex = 0;

	//1 = mirror boundary, -1 = freeflow boundary
	float restitution = -1;
	
	real cfl = .5;

	real dt = {};
	bool useFixedDT = false;

	int meshGenerationIndex = 0;

	//do we rotate to align fluxes with x-axis and only use the x-axis flux,
	// or do we use the general normal-based flux computation
	bool rotateToAlign = false;

	
	std::vector<std::shared_ptr<DisplayMethod>> displayMethods;
	std::vector<const char*> displayMethodNames;



	//TODO
	//this is not just reflection, but also has a lot of gui-specific data tied into it
	//maybe I should call it a separate variable, like 'gui =' instead of 'fields =' ... ?
	static constexpr auto fields = std::make_tuple(
		std::make_tuple("time", &This::time, GUIReadOnly()),
		
		//these are in Super ... inherit them automatically plz
		std::make_pair("running", &Super::running),
		std::make_pair("step", &Super::singleStep),	//TODO as a button? with lambdas maybe?
		//std::make_pair("mousepos", &Super::mousepos),		//TODO as a hover text / readonly

		std::make_pair("cfl", &This::cfl),
		std::make_pair("dt", &This::dt),
		std::make_pair("useFixedDT", &This::useFixedDT),
		
		GUISeparator(),
		
		std::make_pair("display method", &This::displayMethodIndex),	//TODO combo from displayMethodNames
		std::make_pair("auto display range", &This::displayAutomaticRange),

		std::make_pair("", &This::drawArgs),

		std::make_pair("restitution", &This::restitution),
		std::make_pair("rotateToAlign", &This::rotateToAlign),
		GUISeparator(),
		
		std::make_pair("flux", &This::calcFluxIndex),	//TODO combo from calcFluxNames
		GUISeparator(),

//		GUICall<This>([](This* sim) { sim->eqn.updateGUI(); }),

		GUISeparator(),
		
		std::make_pair("init cond", &This::initCondIndex),	//TODO combo from initCondNames
		//TODO eqn.initConds[initCondIndex]->updateGUI() here
		
		std::make_pair("reset state", &This::resetState),	//button, not a field ... hmm
		GUISeparator(),
		
		std::make_pair("mesh generation", &This::meshGenerationIndex),	//TODO combo from meshGenerationNames
		//TODO meshGenerators[meshGenerationIndex]->updateGUI() here
		
		std::make_pair("reset mesh", &This::resetMesh),	//button
		GUISeparator()
	);


	Simulation(CFDMeshApp* app_) : Super(app_) {


/*
test case:
between UL and UL, examine the flux for varying angles
results: looks like the flux is aligned with the normal, and the flux magnitude is constant for all angles
*/
#if 0
	using Prim = typename ThisEquation::Prim;
	constexpr int n = 500;
	Cons UL = eqn.consFromPrim(Prim( 1, real3(), 1 ));
	Cons UR = eqn.consFromPrim(Prim( .125, real3(), .1 ));
	for (int i = 0; i < n; ++i) {
		real theta = ((real)i + .5) / (real)n * 2. * M_PI;
		real3 n(cos(theta), sin(theta), 0);	
		Cons flux = (this->*calcFluxes[0])(UL, UR, .1, .1, n);
		std::cout << theta << "\t" << flux << "\t" << normal << std::endl;
	}
exit(0);
#endif


		//first copy display methods from equation
		displayMethods = eqn.displayMethods;

		//next incorporate any simulation-based display methods
		//this is a face-based display.  for now I'll just use the first face.
		//TODO gui fields for DisplayMethod ?
		displayMethods.push_back(std::make_shared<DisplayMethod>(
			"eigenbasis orthogonality error",
			[this](const ThisEquation* eqn, const Cell* c) -> float {
				const Face& face = m->faces[m->cellFaceIndexes[c->faceOffset]];
				auto [UL, UR] = getEdgeStates(&face);
				
				real3 fluxNormal;
				if (rotateToAlign) {	//rotate normal to x-axis
					UL = eqn->rotateFrom(UL, face.normal);
					UR = eqn->rotateFrom(UR, face.normal);
					fluxNormal = real3(1, 0, 0);
				} else {
					fluxNormal = face.normal;
				}
			

				Eigen vars = eqn->calcRoeAvg(UL, UR);
				
				real value = 0;
				for (int k = 0; k < eqn->numWaves; ++k) {
					Cons basis;
					for (int j = 0; j < eqn->numStates; ++j) {
						basis.ptr[j] = k == j ? 1 : 0;
					}
					
					WaveVec charVars = eqn->applyEigL(basis, vars, face.normal);
					Cons newbasis = eqn->applyEigR(charVars, vars, face.normal);
				
					for (int j = 0; j < eqn->numStates; ++j) {
						value += fabs(newbasis.ptr[j] - basis.ptr[j]);
					}
				}

				return (float)value;
			}
		));

		//last update the names
		displayMethodNames = map<
			decltype(displayMethods),
			std::vector<const char*>
		>(
			displayMethods,
			[](const std::shared_ptr<DisplayMethod>& m) -> const char* {
				return m->name.c_str();
			}
		);



		meshGenerators = ThisMeshNamespace::getGens();
		
		meshGenerationNames = map<
			decltype(meshGenerators),
			std::vector<const char*>
		>(
			meshGenerators,
			[](const std::shared_ptr<MeshFactory>& m) -> const char* { return m->name; }
		);
	

		resetMesh();	//which calls resetState()
	}

	virtual void draw();

	void resetMesh() {
		timeFunc("building mesh", [this]() {
			m = meshGenerators[meshGenerationIndex]->createMesh();
		});
		resetState();
	}

	virtual void resetState() {
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

	std::pair<Cons, Cons> getEdgeStates(const Face* e) {
		Cell* cL = e->cells(0) == -1 ? nullptr : &m->cells[e->cells(0)];
		Cell* cR = e->cells(1) == -1 ? nullptr : &m->cells[e->cells(1)];
		if (cL && cR) {
			return std::pair<Cons, Cons>{cL->U, cR->U};
		} else if (cL) {
			return std::pair<Cons, Cons>{cL->U, eqn.reflect(cL->U, e->normal, restitution)};
		} else if (cR) {
			return std::pair<Cons, Cons>{eqn.reflect(cR->U, e->normal, restitution), cR->U};
		}
		throw Common::Exception() << "here";
	}

	using CalcLambdaVars = typename ThisEquation::CalcLambdaVars;
	using Eigen = typename ThisEquation::Eigen;

	real calcDT() {
		//for (auto& c : m->cells) {
		real result = parallel.reduce(
			m->cells.begin(),
			m->cells.end(),
			[this](Cell& c) -> real {
				const Cons& U = c.U;
				real result = std::numeric_limits<real>::infinity();
				for (int i = 0; i < c.faceCount; ++i) {
					Face* f = &m->faces[m->cellFaceIndexes[i+c.faceOffset]];
					real dx = f->area;
					if (dx > 1e-7 && f->cells(0) != -1 && f->cells(1) != -1) {
						CalcLambdaVars vars(eqn, U, f->normal);
						auto [lambdaMin, lambdaMax] = eqn.calcLambdaMinMax(vars);
						result = std::min(result, dx / (std::max(fabs(lambdaMin), fabs(lambdaMax)) + 1e-9));
					}
				}
				return result;
			},
			std::numeric_limits<real>::infinity(),
			[](real a, real b) -> real { return std::min(a,b); }
		);
		// calculate dt
		return result * cfl;
	}

	//NOTICE: dt used for dt/dx with the flux limiter
	//which needs slope information
	//which I don't yet have
	Cons calcFluxRoe(Cons UL, Cons UR, real dx, real dt, real3 n) {
		Eigen vars = eqn.calcRoeAvg(UL, UR);

		WaveVec lambdas = eqn.getEigenvalues(vars, n);

		Cons dU = UR - UL;
		WaveVec dUTilde = eqn.applyEigL(dU, vars, n);
	
		WaveVec fluxTilde;
		for (int j = 0; j < ThisEquation::numWaves; ++j) {
			real lambda = lambdas(j);
			real phi = 0;
			real sgnLambda = lambda >= 0 ? 1 : -1;
			real epsilon = lambda * dt / dx;
			fluxTilde(j) = -.5 * lambda * dUTilde(j) * (sgnLambda + phi * (epsilon - sgnLambda));
		}
	
		Cons UAvg = (UR + UL) * .5;
		WaveVec UAvgTilde = eqn.applyEigL(UAvg, vars, n);
		fluxTilde = fluxTilde + lambdas * UAvgTilde;
	
		Cons flux = eqn.applyEigR(fluxTilde, vars, n);
		// here's the flux, aligned along the normal

#if 0	//debug print eigenbasis error
		real value = 0;
		for (int k = 0; k < eqn.numWaves; ++k) {
			Cons basis;
			for (int j = 0; j < eqn.numStates; ++j) {
				basis.ptr[j] = k == j ? 1 : 0;
			}
			
			WaveVec charVars = eqn.applyEigL(basis, vars, n);
			Cons newbasis = eqn.applyEigR(charVars, vars, n);
		
			for (int j = 0; j < eqn.numStates; ++j) {
				value += fabs(newbasis.ptr[j] - basis.ptr[j]);
			}
		}

		std::cout << "vars=" << vars << std::endl;
		std::cout << "n=" << n << std::endl;
		std::cout << "evL rows" << std::endl;
		for (int k = 0; k < eqn.numWaves; ++k) {
			Cons basis;
			for (int j = 0; j < eqn.numStates; ++j) {
				basis.ptr[j] = k == j ? 1 : 0;
			}
			std::cout << eqn.applyEigL(basis, vars, n) << std::endl;
		}

		std::cout << "evR rows" << std::endl;
		for (int k = 0; k < eqn.numWaves; ++k) {
			WaveVec basis;
			for (int j = 0; j < eqn.numWaves; ++j) {
				basis.ptr[j] = k == j ? 1 : 0;
			}
			std::cout << eqn.applyEigR(basis, vars, n) << std::endl;
		}
		std::cout << "ortho error = " << value << std::endl;
#endif

		return flux;
	}

	Cons calcFluxHLL(Cons UL, Cons UR, real dx, real dt, real3 n) {
		real lambdaMinL = eqn.calcLambdaMin(CalcLambdaVars(eqn, UL, n));
		real lambdaMaxR = eqn.calcLambdaMax(CalcLambdaVars(eqn, UR, n));
		
		Eigen vars = eqn.calcRoeAvg(UL, UR);
		std::pair<real, real> lambdaMinMax = eqn.calcLambdaMinMax(CalcLambdaVars(vars, n));
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
			Cons fluxL = eqn.calcFluxFromCons(UL, n);
			flux = fluxL;
		} else if (sl <= 0. && 0. <= sr) {
			Cons fluxL = eqn.calcFluxFromCons(UL, n);
			Cons fluxR = eqn.calcFluxFromCons(UR, n);
			//(sr * fluxL[j] - sl * fluxR[j] + sl * sr * (UR[j] - UL[j])) / (sr - sl)
			real invDenom = 1. / (sr - sl);
			for (int i = 0; i < Cons::size; ++i) {
				flux(i) = (sr * fluxL(i) - sl * fluxR(i) + sl * sr * (UR(i) - UL(i))) * invDenom;
			}
		} else if (sr <= 0.) {
			Cons fluxR = eqn.calcFluxFromCons(UR, n);
			flux = fluxR;
		}
	
		return flux;
	}

	using CalcFluxFunc = Cons (Simulation::*)(Cons, Cons, real, real, real3);

	int calcFluxIndex = 0;
	
	std::vector<CalcFluxFunc> calcFluxes = {
		&Simulation::calcFluxRoe,
		&Simulation::calcFluxHLL,
	};

	std::vector<const char*> calcFluxNames = {
		"Roe",
		"HLL",
	};
	
	void step() {
		if (!useFixedDT) {
			dt = calcDT();
		}
//std::cout << "dt " << dt << std::endl;

		parallel.foreach(
			m->faces.begin(),
			m->faces.end(),
			[this](Face& face) {
				if (face.area <= 1e-7) return;

				auto [UL, UR] = getEdgeStates(&face);

#if DEBUG
for (int i = 0; i < Cons::size; ++i) {
	if (!std::isfinite(UL(i))) {
		std::cerr
			<< "got non-finite left edge state " << UL
			<< " for face " << face
			<< std::endl;
		running = false;
	}
	if (!std::isfinite(UR(i))) {
		std::cerr
			<< "got non-finite right edge state " << UR
			<< " for face " << face;
		running = false;
	}
}
#endif

				real3 fluxNormal;
				if (rotateToAlign) {	//rotate normal to x-axis
					UL = eqn.rotateFrom(UL, face.normal);
					UR = eqn.rotateFrom(UR, face.normal);
					fluxNormal = real3(1, 0, 0);
				} else {
					fluxNormal = face.normal;
				}

#if DEBUG
for (int i = 0; i < Cons::size; ++i) {
	if (!std::isfinite(UL(i))) { 
		std::cerr << "got non-finite post-rotate left state " << UL << std::endl;
		running = false;
	}
	if (!std::isfinite(UR(i))) { 
		std::cerr << "got non-finite post-rotate right state " << UR << std::endl;
		running = false;
	}
}
#endif
				face.flux = (this->*calcFluxes[calcFluxIndex])(UL, UR, face.cellDist, dt, fluxNormal);
#if DEBUG
for (int i = 0; i < Cons::size; ++i) {
	if (!std::isfinite(face.flux(i))) { 
		std::cerr
			<< "got non-finite flux " << face.flux << "\n"
			<< " face=" << face << "\n"
			<< " cellL=" << m->cells[face.cells(0)] << "\n"
			<< " cellR=" << m->cells[face.cells(1)] << "\n"
			<< " UL=" << UL << "\n"
			<< " UR=" << UR << "\n"
			<< " WL=" << eqn.primFromCons(UL) << "\n"
			<< " WR=" << eqn.primFromCons(UR) << "\n"
			<< " vars=" << eqn.calcRoeAvg(UL, UR) << "\n"
			<< " lambdas=" << eqn.getEigenvalues(eqn.calcRoeAvg(UL, UR), fluxNormal) << "\n"
			<< " dUTilde=" << eqn.applyEigL(UR - UL, eqn.calcRoeAvg(UL, UR), fluxNormal) << "\n"
			<< " dt=" << dt << "\n"
			<< " fluxNormal=" << fluxNormal << "\n"
			<< std::endl;
		running = false;
		break;
	}
}
#endif

				if (rotateToAlign) {
					face.flux = eqn.rotateTo(face.flux, face.normal);
				}
			}
		);

		//TODO instead of f.e. directly into the cell, write to a derivBuf, and make the integrator modular to work with Runge-Kutta etc 
		parallel.foreach(
			m->cells.begin(),
			m->cells.end(),
			[this](Cell& c) {
//Cons U = c.U;
//std::cout << "integrating cell " << (&c - &*m->cells.begin()) << " = " << c << std::endl;				
				Cons dU_dt;
				for (int i = 0; i < c.faceCount; ++i) {
					int ei = m->cellFaceIndexes[i+c.faceOffset];
					Face* e = &m->faces[ei];
					if (&c == &m->cells[e->cells(0)]) {
//std::cout << " ... - " << *e << std::endl;
						dU_dt -= e->flux * (e->area / c.volume);
					} else {
//std::cout << " ... + " << *e << std::endl;
						dU_dt += e->flux * (e->area / c.volume);
					}
				}
				c.U += dU_dt * dt;
			
//std::cout << "cell before " << U << " after " << c << std::endl;

#if 0	//Euler-only
//TODO is template of ThisEquation, Equation::Euler::Euler
typename ThisEquation::Prim W = eqn.primFromCons(c.U);
if (W.P <= 0) {
	std::cerr << "pressure negative. "
		<< " W=" << W
		//<< " oldU=" << U
		<< " c=" << c
		<< std::endl;
	running = false;
}
#endif
			}
		);
		
		refreshDisplayValues();

		time += dt;
	
#if 0	//display v len min max over time
		auto iter = std::find(displayMethodNames.begin(), displayMethodNames.end(), std::string("v len"));
assert(iter != displayMethodNames.end());	
		int pushDisplayIndex = displayMethodIndex;
		displayMethodIndex = iter - displayMethodNames.begin();
		refreshDisplayValues();
		displayMethodIndex = pushDisplayIndex;
std::cout << "time=" << time 
	<< "\tvlenMin=" << drawArgs.displayValueRange.first
	<< "\tvlenMax=" << drawArgs.displayValueRange.second
	<< std::endl;
#endif
	}
	
	void refreshDisplayValues() {
		parallel.foreach(
			m->cells.begin(),
			m->cells.end(),
			[this](Cell& c) {
				c.displayValue = displayMethods[displayMethodIndex]->f(&eqn, &c);
			}
		);

		if (displayAutomaticRange) {
			using ValueRange = typename Mesh::DrawArgs::ValueRange;
			drawArgs.displayValueRange = parallel.reduce(
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

	virtual void updateGUI();
	
	virtual void onUpdate() {
		if (running || singleStep) {
			step();
			if (singleStep) {
				running = false;
				singleStep = false;
			}
		}
	}
};

std::vector<std::pair<const char*, std::function<std::shared_ptr<ISimulation>(CFDMeshApp*)>>> simGens = {
	{"2D Euler", [](CFDMeshApp* app) -> std::shared_ptr<ISimulation> { return std::make_shared<Simulation<real, 2, Equation::Euler::Euler<real, 2>>>(app); }},
	{"2D GLM-Maxwell", [](CFDMeshApp* app) -> std::shared_ptr<ISimulation> { return std::make_shared<Simulation<real, 2, Equation::GLMMaxwell::GLMMaxwell<real, 2>>>(app); }},
	
	{"3D Euler", [](CFDMeshApp* app) -> std::shared_ptr<ISimulation> { return std::make_shared<Simulation<real, 3, Equation::Euler::Euler<real, 3>>>(app); }},
};

std::vector<const char*> simGenNames = map<
	decltype(simGens),
	std::vector<const char*>
>(
	simGens,
	[](decltype(simGens)::value_type p) -> const char* { return p.first; }
);

struct CFDMeshApp : public ::GLApp::ViewBehavior<::GLApp::GLApp> {
	using Super = ::GLApp::ViewBehavior<::GLApp::GLApp>;

	std::shared_ptr<ISimulation> sim;

	std::shared_ptr<ImGuiCommon::ImGuiCommon> gui;

	GLuint gradientTex = {};

	int simGenIndex = 0;

	int viewIndex = 1;
	std::vector<std::shared_ptr<::GLApp::View>> views;
	std::vector<const char*> viewNames = {"ortho", "frustum"};


	virtual const char* getTitle() { return "CFD Mesh"; }

	virtual void init(const Init& args) {
		Super::init(args);

		gui = std::make_shared<ImGuiCommon::ImGuiCommon>(window, context);
		
		viewOrtho->zoom(0) = viewOrtho->zoom(1) = .5;
		viewFrustum->dist = 3;
		views = decltype(views){viewOrtho, viewFrustum};
		view = views[viewIndex];

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
			f *= (float)(gradientColors.size()-1);
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
		glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP);
		glBindTexture(GL_TEXTURE_1D, 0);

		glEnable(GL_DEPTH_TEST);

		resetSimulation();
	}
	
	void resetSimulation() {
		sim = simGens[simGenIndex].second(this);
	}

	virtual void onUpdate() {
		Super::onUpdate();
	
		sim->draw();
		
		gui->onUpdate([this](){

			if (igCombo("view method", &viewIndex, viewNames.data(), viewNames.size(), -1)) {
				view = views[viewIndex];
			}

			if (igCombo("simulation", &simGenIndex, simGenNames.data(), simGenNames.size(), -1)) {
				resetSimulation();
				return;
			}
			
			sim->updateGUI();
		});
	
		sim->onUpdate();
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
					sim->resetState();
				} else if (event.key.keysym.sym == SDLK_u) {
					sim->singleStep = true;
				} else if (event.key.keysym.sym == SDLK_SPACE) {
					sim->running = !sim->running;
				}
			}
			break;
		case SDL_MOUSEMOTION:
			sim->mousepos = float2(
				(float)event.motion.x / (float)screenSize(0),
				(float)event.motion.y / (float)screenSize(1));
			break;
		}
	}
};

template<typename real, int dim, typename Equation>
void Simulation<real, dim, Equation>::updateGUI() {

#if 1
	igPushIDStr("mesh generation");
	if (igCombo("", &meshGenerationIndex, meshGenerationNames.data(), meshGenerationNames.size(), -1)) {
		resetMesh();
	}
	igPopID();
	igSameLine(0, 0);
	igPushIDStr("mesh generation fields");
	if (igCollapsingHeader("", 0)) {
		meshGenerators[meshGenerationIndex]->updateGUI();
	}
	igPopID();
	if (igSmallButton("reset mesh")) {
		resetMesh();
	}
	
	igSeparator();
	
	igText("time: %e", time);

	//TODO put the fps calculation somewhere else, maybe so it is not based on gui updates but on simulation updates
	{
		static std::chrono::high_resolution_clock::time_point lasttime = std::chrono::high_resolution_clock::now();
		std::chrono::high_resolution_clock::time_point thistime = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> timespan = std::chrono::duration_cast<std::chrono::duration<double>>(thistime - lasttime);
		double dt = timespan.count();
		double fps = 1. / dt;
		lasttime = thistime;
		
		igSameLine(0, 0);
		igText(" fps: %.2f", fps);
	}

	igCheckbox("running", &running);
	
	if (igSmallButton("step")) singleStep = true;

	CFDMesh::updateGUI(&cfl, "cfl");
	CFDMesh::updateGUI(&dt, "dt");
	CFDMesh::updateGUI(&useFixedDT, "use fixed dt");

	CFDMesh::updateGUI(&drawArgs);

	igSeparator();

	if (igCombo("display method", &displayMethodIndex, displayMethodNames.data(), displayMethodNames.size(), -1)) {
		refreshDisplayValues();
	}

	igCheckbox("auto display range", &displayAutomaticRange);
	
	igSeparator();

	igCheckbox("rotate to align", &rotateToAlign);

	igInputFloat("restitution", &restitution, .1, 1., "%e", 0);
	
	igSeparator();

	igCombo("flux", &calcFluxIndex, calcFluxNames.data(), calcFluxNames.size(), -1);
	
	igSeparator();
	
	eqn.updateGUI();
	
	igSeparator();
	
	igPushIDStr("init cond");
	igCombo("", &initCondIndex, eqn.initCondNames.data(), eqn.initCondNames.size(), -1);
	igPopID();
	igPushIDStr("init cond fields");
	igSameLine(0, 0);
	if (igCollapsingHeader("", 0)) {
		eqn.initConds[initCondIndex]->updateGUI();
	}
	igPopID();
	
	if (igSmallButton("reset state")) {
		resetState();
	}
	
	igSeparator();
	
#else
	CFDMesh::updateGUI(this);
#endif

	if (app->view == app->viewOrtho) {
		//find the view bounds
		//find the mouse position
		//find any cells at that position
		real2 pos = (real2)((mousepos - .5f) * float2(1, -1) * float2(app->getAspectRatio(), 1) / float2(app->viewOrtho->zoom(0), app->viewOrtho->zoom(1)) + app->viewOrtho->pos);
			
		bool canHandleMouse = !igGetIO()->WantCaptureMouse;
		
		if (canHandleMouse) {
			igBeginTooltip();
			igText("%e %e\n", pos(0), pos(1));
			igEndTooltip();
		}

		drawArgs.selectedCellIndex = -1;
		for (int i = 0; i < (int)m->cells.size(); ++i) {
			Cell* c = &m->cells[i];
			if (ThisMeshNamespace::polygonContains(
				pos,
				m->cellVtxIndexes.begin() + c->vtxOffset,
				m->cellVtxIndexes.begin() + c->vtxOffset + c->vtxCount,
				[this](int i) -> real2 {
					typename ThisMeshNamespace::Vertex& vi = m->vtxs[i];
					return real2([&vi](int j) -> real { return vi.pos(j); });
				}
			)) {
				drawArgs.selectedCellIndex = i;
				break;
			}
		}

		if (drawArgs.selectedCellIndex != -1) {
			Cell* c = &m->cells[drawArgs.selectedCellIndex];
			
			if (canHandleMouse) {
				igBeginTooltip();
				igText("%e", c->displayValue);
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
				m->cells.erase(m->cells.begin() + drawArgs.selectedCellIndex);
			
				//sort largest -> smallest
				std::sort(toRemoveEdges.begin(), toRemoveEdges.end(), [](int a, int b) -> bool { return a > b; });
				for (int ei : toRemoveEdges) {
					m->faces.erase(m->faces.begin() + ei);
				}
			}
#endif
		}
	}
}


template<typename real, int dim, typename ThisEquation>
void Simulation<real, dim, ThisEquation>::draw() {
	drawArgs.gradientTex = app->gradientTex;
	m->draw(drawArgs);
}

}

GLAPP_MAIN(CFDMesh::CFDMeshApp)
