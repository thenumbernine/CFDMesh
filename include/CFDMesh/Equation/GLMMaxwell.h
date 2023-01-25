#pragma once

#include "CFDMesh/Equation/Equation.h"
#include "ImGuiCommon/Reflect.h"
#include "Tensor/Tensor.h"
#include "Common/String.h"
#include "Common/Macros.h"
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
	using This = Cons_;
	using real3 = Tensor::vec<real, 3>;
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
		std::make_pair("D", &This::D),
		std::make_pair("B", &This::B),
		std::make_pair("phi", &This::phi),
		std::make_pair("psi", &This::psi),
		std::make_pair("rhoCharge", &This::rhoCharge),
		std::make_pair("sigma", &This::sigma),
		std::make_pair("sqrt_1_eps", &This::sqrt_1_eps),
		std::make_pair("sqrt_1_mu", &This::sqrt_1_mu)
	);
};


template<typename T>
using Prim_ = Cons_<T>;


template<typename real, int dim_>
struct GLMMaxwell : public Equation<GLMMaxwell<real, dim_>, real, Cons_<real>, Prim_<real>> {
	using This = GLMMaxwell;
	using Super = Equation<GLMMaxwell<real, dim_>, real, Cons_<real>, Prim_<real>>;
	using Cons = typename Super::Cons;
	using Prim = typename Super::Prim;

	enum { numWaves = 8 };
	using WaveVec = Tensor::vec<real, numWaves>;
	
	using InitCond = typename Super::InitCond;
	using DisplayMethod = typename Super::DisplayMethod;
	
	using real3 = Tensor::vec<real, 3>;
	using real3x3 = Tensor::mat<real, 3, 3>;

	struct InitCondDefault : public InitCond {
		using InitCond::InitCond;
		
		Cons UL = Cons(real3(1, 1, 0), real3(0, -1, 1), 0, 1, 1, 1);
		Cons UR = Cons(real3(-1, 1, 0), real3(0, -1, -1), 0, 1, 1, 1);
	
		static constexpr auto fields = std::make_tuple(
			std::make_pair("UL", &InitCondDefault::UL),
			std::make_pair("UR", &InitCondDefault::UR)
		);

		virtual char const * name() const { return "default"; }
		
		virtual Cons initCell(GLMMaxwell const * eqn, real3 x) const {
			bool lhs = x(0) < 0 && x(1) < 0;
			return lhs ? UL : UR;
		}
		virtual void updateGUI() {
			ImGuiCommon::updateGUI(this);
		}
	};
	
	using Super::Super;

	void buildInitCondsAndDisplayVars() {
		Super::initConds = {
			std::make_shared<InitCondDefault>(),
		};
	}

	real3 calc_E(Cons const & U) const {
		return U.D * (U.sqrt_1_eps * U.sqrt_1_eps);
	}

	real3 calc_H(Cons const & U) const {
		return U.B * (U.sqrt_1_mu * U.sqrt_1_mu);
	}

	Cons const & consFromPrim(Prim const & W) { return W; }
	Prim const & primFromCons(Cons const & U) { return U; }

	//variables used by the roe avg, eigen basis, etc
	struct Eigen {
		real sqrt_1_eps;
		real sqrt_1_mu;
		
		static constexpr auto fields = std::make_tuple(
			std::make_pair("sqrt_1_eps", &Eigen::sqrt_1_eps),
			std::make_pair("sqrt_1_mu", &Eigen::sqrt_1_mu)
		);
	};

	Eigen calcRoeAvg(Cons UL, Cons UR) const {
		Eigen vars;
		vars.sqrt_1_eps = .5 * (UL.sqrt_1_eps + UR.sqrt_1_eps);
		vars.sqrt_1_mu = .5 * (UL.sqrt_1_mu + UR.sqrt_1_mu);
#if 0	//TODO map from eigen field to Cons field
		Common::TupleForEach(Eigen::fields, [&UL, &UR](auto x, size_t i) constexpr -> bool {
			auto field = std::get<1>(x);
			vars.*field = (UL.*field + UR.*field) * .5;
			return false;
		});
#endif
		return vars;
	}

	//static constexpr real sqrt_2 = sqrt(2);
	//clang++: error: constexpr variable 'sqrt_1_2' must be initialized by a constant expression
	//https://stackoverflow.com/a/8622277 : std::sqrt is not defined as constexpr, according to section 26.8 of N3291: the C++11 FDIS
	//doh.
	static constexpr real sqrt_2 = (real)1.41421356237309514547462185873882845044136047363281;
	static constexpr real sqrt_1_2 = 1. / sqrt_2;
	static constexpr real speedOfLight = 1;
	real divPhiWavespeed = speedOfLight;
	real divPsiWavespeed = speedOfLight;

	WaveVec getEigenvalues(Eigen const & vars, real3x3 const & n) {
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
		CalcLambdaVars(GLMMaxwell const & eqn, Cons const & U, real3x3 const & n_) {
			Eigen::sqrt_1_eps = U.sqrt_1_eps;
			Eigen::sqrt_1_mu = U.sqrt_1_mu;
		}

		CalcLambdaVars(Eigen const & vars, real3x3 const & n_) : Eigen(vars) {}
	};

	std::pair<real, real> calcLambdaMinMax(CalcLambdaVars const & vars) {
		real v_p = vars.sqrt_1_eps * vars.sqrt_1_mu;
		real lambda = std::max(std::max(divPsiWavespeed, divPhiWavespeed), v_p);
		return std::make_pair(-lambda, lambda);
	}

	real calcLambdaMin(CalcLambdaVars const & vars) const {
		return -calcLambdaMax(vars);
	}

	real calcLambdaMax(CalcLambdaVars const & vars) const {
		real v_p = vars.sqrt_1_eps * vars.sqrt_1_mu;
		return std::max(std::max(divPsiWavespeed, divPhiWavespeed), v_p);
	}

	WaveVec applyEigL(Cons const & x, Eigen const & vars, real3x3 const & n) const {
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

	Cons applyEigR(WaveVec const & x, Eigen const & vars, real3x3 const & n) const {
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

	Cons calcFluxFromCons(Cons const & U, real3x3 const & n) const {
		real3 E = calc_E(U);
		real3 H = calc_H(U);
		Cons F;
		F.D = real3(U.phi * divPhiWavespeed,  H(2), -H(1));
		F.B = real3(U.psi * divPsiWavespeed, -E(2),  E(1));
		F.phi = U.D(0) * divPhiWavespeed;
		F.psi = U.B(0) * divPsiWavespeed;
		return F;
	}

	static constexpr auto fields = std::make_tuple(
		std::make_pair("divPhiWavespeed", &This::divPhiWavespeed),
		std::make_pair("divPsiWavespeed", &This::divPsiWavespeed)
	);

	void updateGUI() {
		ImGuiCommon::updateGUI(this);
	}
};

}
}
}
