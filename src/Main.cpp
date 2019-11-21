#include "CFDMesh/Equation/Euler.h"
#include "CFDMesh/Equation/GLMMaxwell.h"

#include "CFDMesh/Mesh/Mesh.h"
#include "CFDMesh/Util.h"
#include "CFDMesh/Vector.h"

#include "GLApp/gl.h"
#include "GLApp/GLApp.h"
#include "GLApp/ViewBehavior.h"

#include "ImGuiCommon/ImGuiCommon.h"

#include "Image/Image.h"

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


template<typename T> T cubed(const T& t) { return t * t * t; }


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
	using Super = ISimulation;

	using WaveVec = typename ThisEquation::WaveVec;
	using Cons = typename ThisEquation::Cons;

	using ThisMeshNamespace = CFDMesh::MeshNamespace<real, dim, Cons>;
	using Mesh = typename ThisMeshNamespace::Mesh;
	using MeshFactory = typename ThisMeshNamespace::MeshFactory;
	using Cell = typename ThisMeshNamespace::Cell;
	using Face = typename ThisMeshNamespace::Face;


	struct FileMeshFactory : public MeshFactory {
		std::string filename = {"grids/n0012_113-33.p2dfmt"};
		
		FileMeshFactory() : MeshFactory("p2dfmt mesh") {}

		virtual std::shared_ptr<Mesh> createMesh() {
			std::shared_ptr<Mesh> mesh = MeshFactory::createMeshSuper();
			
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
			std::vector<real> x = map<std::list<std::string>, std::vector<real>>(_x, [](const std::string& s) -> real { return std::stod(s); });
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

	struct Chart2DMeshFactory : public MeshFactory {
		using This = Chart2DMeshFactory; 
		
		int2 size = int2(101, 101);
		float2 mins = real2(-1, -1);
		float2 maxs = real2(1, 1);
		int2 repeat = int2(0, 0);
		int2 capmin = int2(0, 0);
		
		static constexpr auto fields = std::make_tuple(
			std::make_pair("size", &This::size),
			std::make_pair("mins", &This::mins),
			std::make_pair("maxs", &This::maxs),
			std::make_pair("repeat", &This::repeat),
			std::make_pair("capmin", &This::capmin)
		);

		Chart2DMeshFactory(const char* name_) : MeshFactory(name_) {}

		virtual real2 coordChart(real2 x) const { return x; }
		
		virtual void updateGUI() {
			updateGUIForFields(this);
		}
	};

	struct TriUnitMeshFactory : public Chart2DMeshFactory {
		using Super = Chart2DMeshFactory;
		TriUnitMeshFactory() : Super("unit square of triangles") {}
		
		virtual std::shared_ptr<Mesh> createMesh() {
			std::shared_ptr<Mesh> mesh = MeshFactory::createMeshSuper();
			
			int m = Super::size(0);
			int n = Super::size(1);
		
			mesh->vtxs.resize(m * n);
			for (int j = 0; j < n; ++j) {
				for (int i = 0; i < m; ++i) {
					real2 x = real2(
						((real)i + .5) / (real)Super::size(0) * (Super::maxs(0) - Super::mins(0)) + Super::mins(0),
						((real)j + .5) / (real)Super::size(1) * (Super::maxs(1) - Super::mins(1)) + Super::mins(1));
					
					real2 u = Super::coordChart(x);
					mesh->vtxs[i + m * j].pos = real3([&u](int i) -> real { return i < real2::size ? u(i) : 0.; });
				}
			}
			
			int imax = Super::repeat(0) ? m : m-1;
			int jmax = Super::repeat(1) ? n : n-1;
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

	struct QuadUnitMeshFactory : public Chart2DMeshFactory {
		using Super = Chart2DMeshFactory;
		QuadUnitMeshFactory(const char* name_ = "unit square of quads") : Super(name_) {}
		
		virtual std::shared_ptr<Mesh> createMesh() {
			std::shared_ptr<Mesh> mesh = MeshFactory::createMeshSuper();

			int2 n = Super::size;
			int2 step(1, n(0));
			int vtxsize = n.volume();
			if (Super::capmin(0)) vtxsize++;
			mesh->vtxs.resize(vtxsize);
			int2 i;
			for (i(1) = 0; i(1) < n(1); ++i(1)) {
				for (i(0) = 0; i(0) < n(0); ++i(0)) {
					real2 x = ((real2)i + .5) / (real2)Super::size * (Super::maxs - Super::mins) + Super::mins;
					real2 u = Super::coordChart(x);
					mesh->vtxs[int2::dot(i, step)].pos = real3([&u](int i) -> real { return i < real2::size ? u(i) : 0.; });
				}
			}
			
			int capindex = n.volume();
			if (Super::capmin(0)) {
				real3 sum;
				for (int j = 0; j < n(1); ++j) {
					sum += mesh->vtxs[0 + n(0) * j].pos;
				}
				mesh->vtxs[capindex].pos = sum / (real)n(1);
			}

			int2 imax;
			for (int j = 0; j < 2; ++j) {
				imax(j) = Super::repeat(j) ? n(j) : n(j)-1;
			}
			int2 in;
			for (i(1) = 0; i(1) < imax(1); ++i(1)) {
				in(1) = (i(1) + 1) % n(1);
				for (i(0) = 0; i(0) < imax(0); ++i(0)) {
					in(0) = (i(0) + 1) % n(0);
					mesh->addCell(std::vector<int>{
						i(0) + n(0) * i(1),
						in(0) + n(0) * i(1),
						in(0) + n(0) * in(1),
						i(0) + n(0) * in(1)
					});
				}
			}

			if (Super::capmin(0)) {
				for (int j = 0; j < imax(1); ++j) {
					int jn = (j + 1) % n(1);
					mesh->addCell(std::vector<int>{ 0 + n(0) * j, 0 + n(0) * jn, capindex });
				}
			}

			mesh->calcAux();
			return mesh;
		}
	};

	struct QuadUnitCbrtMeshFactory : public QuadUnitMeshFactory {
		QuadUnitCbrtMeshFactory() : QuadUnitMeshFactory("unit square of quads, cbrt mapping") {}
		virtual real2 coordChart(real2 v) const {
			return real2(cbrt(v(0)), cbrt(v(1)));
		}
	};

	struct QuadUnitCubedMeshFactory : public QuadUnitMeshFactory {
		QuadUnitCubedMeshFactory() : QuadUnitMeshFactory("unit square of quads, cubed mapping") {}
		virtual real2 coordChart(real2 v) const {
			return real2(cubed(v(0)), cubed(v(1)));
		}
	};

	struct TwistQuadUnitMeshFactory : public QuadUnitMeshFactory {
		TwistQuadUnitMeshFactory() : QuadUnitMeshFactory("unit square of quads, twist in the middle") {}
		virtual real2 coordChart(real2 v) const {
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
		using Super = QuadUnitMeshFactory;
		DonutQuadUnitMeshFactory() : Super("polar") {
			Super::size = int2(51, 201);
			Super::mins = real2(.5, 0);
			Super::maxs = real2(1, 2*M_PI);
			Super::repeat = int2(0, 1);
			//Super::capmin = int2(1, 0);
		}
		virtual real2 coordChart(real2 v) const {
			return real2(cos(v(1)), sin(v(1))) * v(0);
		}
	};

	//not inheriting from QuadUnitMeshFactory because it has variable size and we want fixed size (based on image size)
	struct QuadUnitBasedOnImageMeshFactory : public Chart2DMeshFactory {
		using Super = Chart2DMeshFactory;
		QuadUnitBasedOnImageMeshFactory() : Super("unit based on image") {}

		std::string imageFilename = "layout.png";
		//std::string imageFilename = "layout.bmp";
		//std::string imageFilename = "layout.tiff";

		virtual std::shared_ptr<Mesh> createMesh() {
			std::shared_ptr<Mesh> mesh = MeshFactory::createMeshSuper();
			
			//auto iimg = ::Image::system->read(imageFilename);	//TODO fixme, it's not working
			auto iimg = ::Image::pngIO->read(imageFilename);
			
			auto img = std::dynamic_pointer_cast<Image::Image>(iimg);

			Super::size = img->getSize();
			
			int m = Super::size(0) + 1;
			int n = Super::size(1) + 1;
		
			mesh->vtxs.resize(m * n);
			for (int j = 0; j < n; ++j) {
				for (int i = 0; i < m; ++i) {
					real2 x = real2(
						((real)i + .5) / (real)Super::size(0) * (Super::maxs(0) - Super::mins(0)) + Super::mins(0),
						((real)j + .5) / (real)Super::size(1) * (Super::maxs(1) - Super::mins(1)) + Super::mins(1));
					
					real2 u = Super::coordChart(x);
					mesh->vtxs[i + m * j].pos = real3([&u](int i) -> real { return i < real2::size ? u(i) : 0.; });
				}
			}

			int imax = m-1;
			int jmax = n-1;
			for (int j = 0; j < jmax; ++j) {
				int jn = (j + 1) % n;
				for (int i = 0; i < imax; ++i) {
					int in = (i + 1) % m;
					if ((*img)(i, jmax-1-j)) {
						mesh->addCell(std::vector<int>{i + m * j, in + m * j, in + m * jn, i + m * jn});
					}
				}
			}

			mesh->calcAux();
			return mesh;
		}

		//override Chart2DMeshFactory and get rid of size
		//TODO add in imageFilename
		virtual void updateGUI() {
			//igInputInt2("size", size.v, 0);
			igInputFloat2("mins", Super::mins.v, "%f", 0);
			igInputFloat2("maxs", Super::maxs.v, "%f", 0);
			igInputInt2("repeat", Super::repeat.v, 0);
			igInputInt2("capmin", Super::capmin.v, 0);
		}
	};

	struct Chart3DMeshFactory : public MeshFactory {
		using This = Chart3DMeshFactory;
		
		int3 size = int3(31,31,31);
		float3 mins = float3(-1, -1, -1);
		float3 maxs = float3(1, 1, 1);

		//TODO support for inheritence and reflection
		static constexpr auto fields = std::make_tuple(
			std::make_pair("size", &This::size),
			std::make_pair("mins", &This::mins),
			std::make_pair("maxs", &This::maxs)
		);

		Chart3DMeshFactory() : MeshFactory("3D chart mesh") {}

		virtual real3 coordChart(real3 x) const { return x; }

		virtual void updateGUI() {
			updateGUIForFields(this);
		}
	};

	struct CubeUnitMeshFactory : public Chart3DMeshFactory {
		using Super = Chart3DMeshFactory;
		CubeUnitMeshFactory(const char* name_ = "unit cube of cubes") : Super(name_) {}
		
		virtual std::shared_ptr<Mesh> createMesh() {
			std::shared_ptr<Mesh> mesh = MeshFactory::createMeshSuper();
		
			int3 n = Super::size;
			int vtxsize = n.volume();
			mesh->vtxs.resize(vtxsize);
			int3 i;
			for (i(2) = 0; i(2) < n(2); ++i(2)) {
				for (i(1) = 0; i(1) < n(1); ++i(1)) {
					for (i(0) = 0; i(0) < n(0); ++i(0)) {
						real3 x = ((real3)i + .5) / (real3)Super::size * (Super::maxs - Super::mins) + Super::mins;
						mesh->vtxs[i(0) + n(0) * (i(1) + n(1) * i(2))].pos = Super::coordChart(x);
					}
				}
			}

			int3 imax;
			for (int j = 0; j < 3; ++j) {
				imax(j) = Super::repeat(j) ? n(j) : n(j)-1;
			}
			int3 in;
			for (i(2) = 0; i(2) < imax(2); ++i(2)) {
				in(2) = (i(2) + 1) % n(1);
				for (i(1) = 0; i(1) < imax(1); ++i(1)) {
					in(1) = (i(1) + 1) % n(1);
					for (i(0) = 0; i(0) < imax(0); ++i(0)) {
						in(0) = (i(0) + 1) % n(0);
						mesh->addCube(std::vector<int>{
							//using z-order
							i(0) + n(0) * (i(1) + n(1) * i(2)),
							in(0) + n(0) * (i(1) + n(1) * i(2)),
							i(0) + n(0) * (in(1) + n(1) * i(2)),
							in(0) + n(0) * (in(1) + n(1) * i(2)),
							
							i(0) + n(0) * (i(1) + n(1) * in(2)),
							in(0) + n(0) * (i(1) + n(1) * in(2)),
							i(0) + n(0) * (in(1) + n(1) * in(2)),
							in(0) + n(0) * (in(1) + n(1) * in(2)),
						});
					}
				}
			}
		
			mesh->calcAux();
		}
	};

	std::vector<std::shared_ptr<MeshFactory>> meshGenerators;
	std::vector<const char*> meshGenerationNames;

	std::shared_ptr<Mesh> m;
	std::function<Cons(real3)> initcond;

	ThisEquation eqn;

	double time = 0;
	
	bool showVtxs = false;
	bool showEdges = false;
	bool showCellCenters = false;
	int selectedCellIndex = -1;

	bool displayAutomaticRange = true;
	int displayMethodIndex = 0;
	
	using ValueRange = std::pair<float, float>; //min, max
	ValueRange displayValueRange = ValueRange(0, 1);

	int initCondIndex = 0;

	//1 = mirror boundary, -1 = freeflow boundary
	float restitution = 1;
		
	float cfl = .5;

	int meshGenerationIndex = 0;

	Simulation(CFDMeshApp* app_) : Super(app_) {
		
		meshGenerators = {
			std::make_shared<QuadUnitMeshFactory>(),
			std::make_shared<TriUnitMeshFactory>(),
			std::make_shared<QuadUnitCbrtMeshFactory>(),
			std::make_shared<QuadUnitCubedMeshFactory>(),
			std::make_shared<TwistQuadUnitMeshFactory>(),
			std::make_shared<DonutQuadUnitMeshFactory>(),
			std::make_shared<QuadUnitBasedOnImageMeshFactory>(),
			std::make_shared<FileMeshFactory>(),
		};

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
				std::pair<Cons, Cons> ULR = getEdgeStates(&e);
				
				//rotate to align edge normal to x axis
				//so x-direction flux jacobian is good for calculating the flux 
				Cons UL = eqn.rotateTo(ULR.first, e.normal);
				Cons UR = eqn.rotateTo(ULR.second, e.normal);

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
	{"Euler 2D", [](CFDMeshApp* app) -> std::shared_ptr<ISimulation> { return std::make_shared<Simulation<real, 2, Equation::Euler::Euler<real>>>(app); }},
	{"Euler 3D", [](CFDMeshApp* app) -> std::shared_ptr<ISimulation> { return std::make_shared<Simulation<real, 3, Equation::Euler::Euler<real>>>(app); }},
	{"GLM-Maxwell 2D", [](CFDMeshApp* app) -> std::shared_ptr<ISimulation> { return std::make_shared<Simulation<real, 2, Equation::GLMMaxwell::GLMMaxwell<real>>>(app); }},
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

	virtual const char* getTitle() { return "CFD Mesh"; }

	virtual void init(const Init& args) {
		Super::init(args);

		gui = std::make_shared<ImGuiCommon::ImGuiCommon>(window, context);
		
		view = viewOrtho;
		viewOrtho->zoom(0) = viewOrtho->zoom(1) = .5;

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

		resetSimulation();
	}
	
	int simGenIndex = 0;
	
	void resetSimulation() {
		sim = simGens[simGenIndex].second(this);
	}

	virtual void onUpdate() {
		Super::onUpdate();
	
		sim->draw();
		
		gui->onUpdate([this](){
			
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
			if (ThisMeshNamespace::contains(
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
	m->draw(
		app->gradientTex,
		displayValueRange.first,
		displayValueRange.second,
		selectedCellIndex,
		showVtxs,
		showEdges,
		showCellCenters
	);
}

}

GLAPP_MAIN(CFDMesh::CFDMeshApp)
