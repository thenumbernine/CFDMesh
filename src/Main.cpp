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

namespace CFDMesh {

//config for everything, to hold everything in one place
using real = double;
//using real = float;

using real2 = Tensor::Vector<real, 2>;
using real3 = Tensor::Vector<real, 3>;

static Parallel::Parallel parallel;

struct CFDMeshApp;

struct ISimulation {
	//written by app SDL events:
	bool running = {};
	bool singleStep = {};
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
	using MeshFactory = typename ThisMeshNamespace::MeshFactory;
	using Cell = typename ThisMeshNamespace::Cell;
	using Face = typename ThisMeshNamespace::Face;

	std::vector<std::shared_ptr<MeshFactory>> meshGenerators;
	std::vector<const char*> meshGenerationNames;

	std::shared_ptr<Mesh> m;
	std::function<Cons(real3)> initcond;

	ThisEquation eqn;

	double time = 0;
	
	using ValueRange = std::pair<float, float>; //min, max
	ValueRange displayValueRange = ValueRange(0, 1);
	
	int selectedCellIndex = -1;
	
	bool showCells = true;
	bool showVtxs = false;
	bool showFaces = false;
	bool showFaceCenters = false;
	bool showCellCenters = false;

	float cellScale = 1;

	bool displayAutomaticRange = true;
	int displayMethodIndex = 0;
	
	int initCondIndex = 0;

	//1 = mirror boundary, -1 = freeflow boundary
	float restitution = 1;
		
	float cfl = .5;

	int meshGenerationIndex = 0;

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
		
		std::make_pair("show cells", &This::showCells),
		std::make_pair("show vtxs", &This::showVtxs),
		std::make_pair("show faces", &This::showFaces),
		std::make_pair("show face centers", &This::showFaceCenters),
		std::make_pair("show cell centers", &This::showCellCenters),
		std::make_pair("cell scale", &This::cellScale),
		GUISeparator(),
		
		std::make_pair("display method", &This::displayMethodIndex),	//TODO combo from displayMethodNames
		std::make_pair("auto display range", &This::displayAutomaticRange),
		std::make_pair("display value range", &This::displayValueRange),
		GUISeparator(),
		
		std::make_pair("restitution", &This::restitution),
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
		m = meshGenerators[meshGenerationIndex]->createMesh();
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
			
				if constexpr (std::is_same_v<ThisEquation, CFDMesh::Equation::Euler::Euler<real>>) {
					assert(std::isfinite(c.U.rho) && c.U.rho > 0);
					assert(std::isfinite(c.U.m(0)));
					assert(std::isfinite(c.U.m(1)));
					assert(std::isfinite(c.U.m(2)));
					assert(std::isfinite(c.U.ETotal) && c.U.ETotal > 0);
				}
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
				
				real result = std::numeric_limits<real>::infinity();
				for (int ei = 0; ei < c.faceCount; ++ei) {
					Face* e = &m->faces[m->cellFaceIndexes[ei+c.faceOffset]];
					
					Cons U = eqn.rotateTo(c.U, e->normal);
					CalcLambdaVars vars(eqn, U);
					
					real lambdaMin, lambdaMax;
					std::pair<real, real> lambdaMinMax = eqn.calcLambdaMinMax(vars);
					lambdaMin = lambdaMinMax.first;
					lambdaMax = lambdaMinMax.second;
					lambdaMin = std::min<real>(0, lambdaMin);
					lambdaMax = std::max<real>(0, lambdaMax);
					// TODO a better way to do this.  maybe use faces' lambdas?  maybe do this after calculating the eigenbasis?
					//real dx = sqrt(c.volume);
					real dx = e->area;
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
		Eigen vars = eqn.calcRoeAvg(UL, UR);

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
		real lambdaMinL = eqn.calcLambdaMin(CalcLambdaVars(eqn, UL));
		
		real lambdaMaxR = eqn.calcLambdaMax(CalcLambdaVars(eqn, UR));
		
		Eigen vars = eqn.calcRoeAvg(UL, UR);
		std::pair<real, real> lambdaMinMax = eqn.calcLambdaMinMax(CalcLambdaVars(vars));
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

	using CalcFluxFunc = Cons (Simulation::*)(Cons, Cons, real, real);
	
	std::vector<CalcFluxFunc> calcFluxes = {
		&Simulation::calcFluxRoe,
		&Simulation::calcFluxHLL,
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
				auto [UL, UR] = getEdgeStates(&e);
				
				//rotate to align edge normal to x axis
				//so x-direction flux jacobian is good for calculating the flux 
				UL = eqn.rotateTo(UL, e.normal);
				UR = eqn.rotateTo(UR, e.normal);

#ifdef DEBUG
for (int i = 0; i < Cons::size; ++i) {
	assert(std::isfinite(UL(i)));	
	assert(std::isfinite(UR(i)));	
}	
#endif
				Cons F = (this->*calcFluxes[calcFluxIndex])(UL, UR, e.cellDist, dt);
		
				// rotate back to normal
				e.flux = eqn.rotateFrom(F, e.normal);
#ifdef DEBUG
for (int i = 0; i < Cons::size; ++i) {
	assert(std::isfinite(e.flux(i)));
}
#endif
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
						if (&c == &m->cells[e->cells(0)]) {
							dU_dt -= e->flux * (e->area / c.volume);
						} else {
							dU_dt += e->flux * (e->area / c.volume);
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
	{"3D Euler", [](CFDMeshApp* app) -> std::shared_ptr<ISimulation> { return std::make_shared<Simulation<real, 3, Equation::Euler::Euler<real>>>(app); }},
	
	{"2D Euler", [](CFDMeshApp* app) -> std::shared_ptr<ISimulation> { return std::make_shared<Simulation<real, 2, Equation::Euler::Euler<real>>>(app); }},
	{"2D GLM-Maxwell", [](CFDMeshApp* app) -> std::shared_ptr<ISimulation> { return std::make_shared<Simulation<real, 2, Equation::GLMMaxwell::GLMMaxwell<real>>>(app); }},
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

	int viewIndex = 0;
	std::vector<std::shared_ptr<::GLApp::View>> views;
	std::vector<const char*> viewNames = {"ortho", "frustum"};


	virtual const char* getTitle() { return "CFD Mesh"; }

	virtual void init(const Init& args) {
		Super::init(args);

		gui = std::make_shared<ImGuiCommon::ImGuiCommon>(window, context);
		
		view = viewOrtho;
		viewOrtho->zoom(0) = viewOrtho->zoom(1) = .5;

		views = decltype(views){viewOrtho, viewFrustum};

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
	igText("time: %f", time);
	igCheckbox("running", &running);
	
	if (igSmallButton("step")) singleStep = true;

	igInputFloat("cfl", &cfl, .1, 1, "%f", 0);

	igCheckbox("showCells", &showCells);
	igCheckbox("showVtxs", &showVtxs);
	igCheckbox("showFaceCenters", &showFaceCenters);
	igCheckbox("showCellCenters", &showCellCenters);
	igCheckbox("showFaces", &showFaces);
	igInputFloat("cell scale", &cellScale, .1, 1, "%f", 0);

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
			igText("%f %f\n", pos(0), pos(1));
			igEndTooltip();
		}

		selectedCellIndex = -1;
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
				selectedCellIndex = i;
				break;
			}
		}

		if (selectedCellIndex != -1) {
			Cell* c = &m->cells[selectedCellIndex];
			
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
}


template<typename real, int dim, typename ThisEquation>
void Simulation<real, dim, ThisEquation>::draw() {
	typename Mesh::DrawArgs args;
	args.gradientTex = app->gradientTex;
	args.displayValueRange = displayValueRange;
	args.selectedCellIndex = selectedCellIndex;
	args.showCells = showCells;
	args.showVtxs = showVtxs;
	args.showFaces = showFaces;
	args.showFaceCenters = showFaceCenters;
	args.showCellCenters = showCellCenters;
	args.cellScale = cellScale;
	m->draw(args);
}

}

GLAPP_MAIN(CFDMesh::CFDMeshApp)
