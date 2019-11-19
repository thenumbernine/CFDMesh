#pragma once

#include "CFDMesh/Equation/Equation.h"
#include "CFDMesh/Util.h"
#include "Tensor/Vector.h"
#include "cimgui.h"
#include <utility>
#include <tuple>
#include <cmath>
#include <cassert>

namespace CFDMesh {
namespace Equation {
namespace GLMMaxwell {

enum { numCons = 12 };

template<typename real>
union Cons_ {
	using real3 = Tensor::Vector<real, 3>;
	enum { size = numCons };
	real ptr[size];
	struct {
		real3 D = {};
		real3 B = {};
		real phi = {};
		real psi = {};
		real rhoCharge = {};
		real sigma = {};
		real sqrt_1_eps = {};
		real sqrt_1_mu = {};
	};

	Cons_() {}

	Cons_(real3 D_, real3 B_, real rhoCharge_, real sigma_, real sqrt_1_eps_, real sqrt_1_mu_) {
		D = D_;
		B = B_;
		rhoCharge = rhoCharge_;
		sigma = sigma_;
		sqrt_1_eps = sqrt_1_eps_;
		sqrt_1_mu = sqrt_1_mu_;
	}
	
	ADD_OPS(Cons_)

	static constexpr auto fields = std::make_tuple(
		std::make_pair("D", &Cons_::D),
		std::make_pair("B", &Cons_::B),
		std::make_pair("phi", &Cons_::phi),
		std::make_pair("psi", &Cons_::psi),
		std::make_pair("rhoCharge", &Cons_::rhoCharge),
		std::make_pair("sigma", &Cons_::sigma),
		std::make_pair("sqrt_1_eps", &Cons_::sqrt_1_eps),
		std::make_pair("sqrt_1_mu", &Cons_::sqrt_1_mu)
	);
};

template<typename T>
ADD_OSTREAM(Cons_<T>)


template<typename T>
using Prim_ = Cons_<T>;


template<typename real>
struct GLMMaxwell : public Equation<GLMMaxwell<real>, real, Cons_<real>, Prim_<real>> {
	using Parent = Equation<GLMMaxwell<real>, real, Cons_<real>, Prim_<real>>;
	using Cons = typename Parent::Cons;
	using Prim = typename Parent::Prim;

	enum { numWaves = 8 };
	using WaveVec = Tensor::Vector<real, 8>;
	
	using InitCond = typename Parent::InitCond;
	using DisplayMethod = typename Parent::DisplayMethod;
	
	using real3 = Tensor::Vector<real, 3>;

	struct InitCondDefault : public InitCond {
		Cons_<float> UL = Cons_<float>(float3(1, 1, 0), float3(0, -1, 1), 0, 0, 1, 1);
		Cons_<float> UR = Cons_<float>(float3(-1, 1, 0), float3(0, -1, -1), 0, 0, 1, 1);
		using InitCond::InitCond;
		virtual const char* name() const { return "default"; }
		virtual Cons initCell(const GLMMaxwell* eqn, real3 x) const {
			bool lhs = x(0) < 0 && x(1) < 0;
			return lhs ? UL : UR;
		}
		virtual void updateGUI() {
			updateGUIForFields(&UL, "UL");
			updateGUIForFields(&UR, "UR");
		}
	};
	
	using Parent::Parent;

	void buildInitCondsAndDisplayVars() {
		Parent::initConds = {
			std::make_shared<InitCondDefault>(),
		};
	}

	real3 calc_E(const Cons& U) const {
		return U.D * (U.sqrt_1_eps * U.sqrt_1_eps);
	}

	real3 calc_H(const Cons& U) const {
		return U.B * (U.sqrt_1_mu * U.sqrt_1_mu);
	}

	const Cons& consFromPrim(const Prim& W) { return W; }
	const Prim& primFromCons(const Cons& U) { return U; }

	//variables used by the roe avg, eigen basis, etc
	struct Eigen {
		real sqrt_1_eps;
		real sqrt_1_mu;
	};

	Eigen calcRoeAvg(Cons UL, Cons UR) {
		Eigen vars;
		vars.sqrt_1_eps = .5 * (UL.sqrt_1_eps + UR.sqrt_1_eps);
		vars.sqrt_1_mu = .5 * (UL.sqrt_1_mu + UR.sqrt_1_mu);
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
			Eigen::sqrt_1_eps = U.sqrt_1_eps;
			Eigen::sqrt_1_mu = U.sqrt_1_mu;
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
		F.D = real3(U.phi * divPhiWavespeed,  H(2), -H(1));
		F.B = real3(U.psi * divPsiWavespeed, -E(2),  E(1));
		F.phi = U.D(0) * divPhiWavespeed;
		F.psi = U.B(0) * divPsiWavespeed;
		return F;
	}
	
	void updateGUI() {
		igInputFloat("divPhiWavespeed", &divPhiWavespeed, .1, 1, "%f", 0);
		igInputFloat("divPsiWavespeed", &divPsiWavespeed, .1, 1, "%f", 0);
	}
};

}
}
}
