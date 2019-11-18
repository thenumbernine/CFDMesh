#pragma once

#include "CFDMesh/Equation/Equation.h"
#include "CFDMesh/Util.h"
#include "Tensor/Vector.h"
#include "cimgui.h"
#include <utility>
#include <cmath>
#include <cassert>

namespace CFDMesh {
namespace Equation {

template<typename real>
struct GLMMaxwellNamespace {

using real3 = Tensor::Vector<real, 3>;

enum { numCons = 12 };

struct Cons : public Tensor::GenericVector<real, numCons, real, Cons> {
	using Parent = Tensor::GenericVector<real, numCons, real, Cons>;
	
	real3& D() { return *(real3*)(Parent::v); }
	real3& B() { return *(real3*)(Parent::v + 3); }
	real& phi() { return Parent::v[6]; }
	real& psi() { return Parent::v[7]; }
	real& rhoCharge() { return Parent::v[8]; }
	real& sigma() { return Parent::v[9]; }
	real& sqrt_1_eps() { return Parent::v[10]; }
	real& sqrt_1_mu() { return Parent::v[11]; }

	const real3& D() const { return *(real3*)(Parent::v); }
	const real3& B() const { return *(real3*)(Parent::v + 3); }
	const real& phi() const { return Parent::v[6]; }
	const real& psi() const { return Parent::v[7]; }
	const real& rhoCharge() const { return Parent::v[8]; }
	const real& sigma() const { return Parent::v[9]; }
	const real& sqrt_1_eps() const { return Parent::v[10]; }
	const real& sqrt_1_mu() const { return Parent::v[11]; }

	using Parent::Parent;

	Cons(real3 D_, real3 B_, real rhoCharge_, real sigma_, real sqrt_1_eps_, real sqrt_1_mu_) {
		D() = D_;
		B() = B_;
		rhoCharge() = rhoCharge_;
		sigma() = sigma_;
		sqrt_1_eps() = sqrt_1_eps_;
		sqrt_1_mu() = sqrt_1_mu_;
	}
};

struct GLMMaxwell : public Equation<real, Cons, GLMMaxwell> {
	using Parent = Equation<real, Cons, GLMMaxwell>;

	enum { numWaves = 8 };
	using WaveVec = Tensor::Vector<real, 8>;
	
	using InitCond = typename Parent::InitCond;
	using DisplayMethod = typename Parent::DisplayMethod;

	struct InitCondDefault : public InitCond {
		float3 DL = float3(1, 1, 0);
		float3 DR = float3(-1, 1, 0);
		float3 BL = float3(0, -1, 1);
		float3 BR = float3(0, -1, -1);
	
		virtual const char* name() const { return "default"; }
		virtual Cons initCell(const GLMMaxwell* eqn, real3 x) const {
			bool lhs = x(0) < 0 && x(1) < 0;
			return Cons(
				lhs ? DL : DR,
				lhs ? BL : BR,
				0,	//rho
				0,	//sigma
				1,	//sqrt(1/eps)
				1	//sqrt(1/mu)
			);
		}
		virtual void updateGUI() {
			igInputFloat("DLx", &DL(0), .1, 1, "%f", 0);
			igInputFloat("DLy", &DL(1), .1, 1, "%f", 0);
			igInputFloat("DLz", &DL(2), .1, 1, "%f", 0);
		
			igInputFloat("DRx", &DR(0), .1, 1, "%f", 0);
			igInputFloat("DRy", &DR(1), .1, 1, "%f", 0);
			igInputFloat("DRz", &DR(2), .1, 1, "%f", 0);
		
			igInputFloat("BLx", &BL(0), .1, 1, "%f", 0);
			igInputFloat("BLy", &BL(1), .1, 1, "%f", 0);
			igInputFloat("BLz", &BL(2), .1, 1, "%f", 0);
		
			igInputFloat("BRx", &BR(0), .1, 1, "%f", 0);
			igInputFloat("BRy", &BR(1), .1, 1, "%f", 0);
			igInputFloat("BRz", &BR(2), .1, 1, "%f", 0);
		}
	};
	
	using Parent::Parent;

	void buildInitCondsAndDisplayVars() {
		Parent::initConds = {
			std::make_shared<InitCondDefault>(),
		};

		Parent::addDisplayVector("D", [](const GLMMaxwell* eqn, const Cons& U) -> real3 { return U.D(); });
		Parent::addDisplayVector("B", [](const GLMMaxwell* eqn, const Cons& U) -> real3 { return U.B(); });
		Parent::addDisplayScalar("phi", [](const GLMMaxwell* eqn, const Cons& U) -> real { return U.phi(); });
		Parent::addDisplayScalar("psi", [](const GLMMaxwell* eqn, const Cons& U) -> real { return U.psi(); });
		Parent::addDisplayScalar("rhoCharge", [](const GLMMaxwell* eqn, const Cons& U) -> real { return U.rhoCharge(); });
		Parent::addDisplayScalar("sigma", [](const GLMMaxwell* eqn, const Cons& U) -> real { return U.sigma(); });
		Parent::addDisplayScalar("sqrt_1_eps", [](const GLMMaxwell* eqn, const Cons& U) -> real { return U.sqrt_1_eps(); });
		Parent::addDisplayScalar("sqrt_1_mu", [](const GLMMaxwell* eqn, const Cons& U) -> real { return U.sqrt_1_mu(); });
	}

	real3 calc_E(const Cons& U) const {
		return U.D() * (U.sqrt_1_eps() * U.sqrt_1_eps());
	}

	real3 calc_H(const Cons& U) const {
		return U.B() * (U.sqrt_1_mu() * U.sqrt_1_mu());
	}

	using Prim = Cons;

	const Cons& consFromPrim(const Prim& W) { return W; }
	const Prim& primFromCons(const Cons& U) { return U; }

	//variables used by the roe avg, eigen basis, etc
	struct Eigen {
		real sqrt_1_eps;
		real sqrt_1_mu;
	};

	Eigen calcRoeAvg(Cons UL, Cons UR) {
		Eigen vars;
		vars.sqrt_1_eps = .5 * (UL.sqrt_1_eps() + UR.sqrt_1_eps());
		vars.sqrt_1_mu = .5 * (UL.sqrt_1_mu() + UR.sqrt_1_mu());
		return vars;
	}

	static constexpr float sqrt_2 = sqrt(2);
	static constexpr float sqrt_1_2 = 1. / sqrt_2;
	static constexpr float speedOfLight = 1;
	float divPhiWavespeed = speedOfLight;
	float divPsiWavespeed = speedOfLight;

	WaveVec getEigenvalues(const Eigen& vars) {
		real v_p = vars.sqrt_1_eps * vars.sqrt_1_mu;
		WaveVec lambdas;
		lambdas(0) = divPhiWavespeed;
		lambdas(1) = divPsiWavespeed;
		lambdas(2) = -v_p;
		lambdas(3) = -v_p;
		lambdas(4) = v_p;
		lambdas(5) = v_p;
		lambdas(6) = divPhiWavespeed;
		lambdas(7) = divPsiWavespeed;
		return lambdas;
	}

	struct CalcLambdaVars : public Eigen {
		CalcLambdaVars(const GLMMaxwell& eqn, const Cons& U) {
			Eigen::sqrt_1_eps = U.sqrt_1_eps();
			Eigen::sqrt_1_mu = U.sqrt_1_mu();
		}

		CalcLambdaVars(const Eigen& vars) : Eigen(vars) {}
	};

	std::pair<real, real> calcLambdaMinMax(const CalcLambdaVars& vars) {
		real v_p = vars.sqrt_1_eps * vars.sqrt_1_mu;
		real lambda = std::max((real)std::max<float>(divPsiWavespeed, divPhiWavespeed), v_p);
		return std::pair<real, real>(-lambda, lambda);
	}

	real calcLambdaMin(const CalcLambdaVars& vars) {
		return -calcLambdaMax(vars);
	}

	real calcLambdaMax(const CalcLambdaVars& vars) {
		real v_p = vars.sqrt_1_eps * vars.sqrt_1_mu;
		return std::max((real)std::max<float>(divPsiWavespeed, divPhiWavespeed), v_p);
	}

	WaveVec apply_evL(const Cons& x, const Eigen& vars) {
		real sqrt_eps = 1. / vars.sqrt_1_eps;
		real sqrt_mu = 1. / vars.sqrt_1_mu;
		real v_p = vars.sqrt_1_eps * vars.sqrt_1_mu;

		WaveVec y;
		y(0) = (x(6) - x(0)) * sqrt_eps * sqrt_1_2;
		y(1) = (x(7) - x(3)) * sqrt_eps * sqrt_1_2;
		y(2) = (x(4) * sqrt_eps + x(2) * sqrt_mu) * v_p * sqrt_1_2;
		y(3) = (x(5) * sqrt_eps - x(1) * sqrt_mu) * v_p * sqrt_1_2;
		y(4) = (x(4) * sqrt_eps - x(2) * sqrt_mu) * v_p * sqrt_1_2;
		y(5) = (x(5) * sqrt_eps + x(1) * sqrt_mu) * v_p * sqrt_1_2;
		y(6) = (x(0) + x(6)) * sqrt_eps * sqrt_1_2;
		y(7) = (x(3) + x(7)) * sqrt_eps * sqrt_1_2;
		return y;
	}

	Cons apply_evR(const WaveVec& x, const Eigen& vars) {
		real sqrt_eps = 1. / vars.sqrt_1_eps;
		real sqrt_mu = 1. / vars.sqrt_1_mu;

		Cons y;
		y(0) = ((-(x(0) - x(6))) / (sqrt_2 * sqrt_eps));
		y(1) = ((-(sqrt_eps * (x(3) - x(5)))) / sqrt_2);
		y(2) = ((sqrt_eps * (x(2) - x(4))) / sqrt_2);
		y(3) = ((-(x(1) - x(7))) / (sqrt_2 * sqrt_eps));
		y(4) = ((sqrt_mu * (x(2) + x(4))) / sqrt_2);
		y(5) = ((sqrt_mu * (x(3) + x(5))) / sqrt_2);
		y(6) = ((x(0) + x(6)) / (sqrt_2 * sqrt_eps));
		y(7) = ((x(1) + x(7)) / (sqrt_2 * sqrt_eps));
		return y;
	}

	Cons calcFluxFromCons(const Cons& U) {
		real3 E = calc_E(U);
		real3 H = calc_H(U);
		Cons F;
		F.D() = real3(U.phi() * divPhiWavespeed,  H(2), -H(1));
		F.B() = real3(U.psi() * divPsiWavespeed, -E(2),  E(1));
		F.phi() = U.D()(0) * divPhiWavespeed;
		F.psi() = U.B()(0) * divPsiWavespeed;
		return F;
	}
	
	void updateGUI() {
		igInputFloat("divPhiWavespeed", &divPhiWavespeed, .1, 1, "%f", 0);
		igInputFloat("divPsiWavespeed", &divPsiWavespeed, .1, 1, "%f", 0);
	}

	Cons rotateTo(Cons U, real3 normal) {
		U.D() = CFDMesh::rotateTo<real3>(U.D(), normal);
		U.B() = CFDMesh::rotateTo<real3>(U.B(), normal);
		return U;
	}
	
	Cons rotateFrom(Cons U, real3 normal) {
		U.D() = CFDMesh::rotateFrom<real3>(U.D(), normal);
		U.B() = CFDMesh::rotateFrom<real3>(U.B(), normal);
		return U;
	}

	Cons reflect(Cons U, real3 normal, real restitution) const {
		U.D() = U.D() - normal * ((1 + restitution) * real3::dot(normal, U.D()));
		U.B() = U.B() - normal * ((1 + restitution) * real3::dot(normal, U.B()));
		return U;
	}

};

};

}
}
