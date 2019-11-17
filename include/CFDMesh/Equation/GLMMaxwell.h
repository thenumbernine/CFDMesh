#pragma once

#include "CFDMesh/Equation/Equation.h"
#include "CFDMesh/Util.h"
#include "Tensor/Vector.h"
#include <utility>
#include <cmath>
#include <cassert>

namespace CFDMesh {
namespace Equation {

template<typename real>
struct GLMMaxwell : public Equation<GLMMaxwell<real>> {
	using real = real;
	using real2 = Tensor::Vector<real, 2>;
	using real3 = Tensor::Vector<real, 3>;

	//notice we only integrate 8 vars
	using StateVec = Tensor::Vector<real, 12>;

	enum { numWaves = 8 };
	using WaveVec = Tensor::Vector<real, 8>;

	struct Cons : public StateVec {
		real3& D() { return *(real3*)(StateVec::v); }
		real3& B() { return *(real3*)(StateVec::v + 3); }
		real& phi() { return StateVec::v[6]; }
		real& psi() { return StateVec::v[7]; }
		real& rhoCharge() { return StateVec::v[8]; }
		real& sigma() { return StateVec::v[9]; }
		real& sqrt_1_eps() { return StateVec::v[10]; }
		real& sqrt_1_mu() { return StateVec::v[11]; }

		using StateVec::StateVec;
	};

	struct DisplayMethod {
		using DisplayFunc = std::function<float(const GLMMaxwell*, const Cons&)>;
		std::string name;
		DisplayFunc f;
		DisplayMethod(const std::string& name_, DisplayFunc f_) : name(name_), f(f_) {}
	};

	std::vector<std::shared_ptr<DisplayMethod>> displayMethods;

	GLMMaxwell() {

		Super::displayMethods = {
			std::make_shared<DisplayMethod>("rho", [](const Euler* eqn, const Cons& U) -> float { return U.rho(); }),
			
			std::make_shared<DisplayMethod>("m", [](const Euler* eqn, const Cons& U) -> float { return real3::length(U.m()); }),
			std::make_shared<DisplayMethod>("mx", [](const Euler* eqn, const Cons& U) -> float { return U.m()(0); }),
			std::make_shared<DisplayMethod>("my", [](const Euler* eqn, const Cons& U) -> float { return U.m()(1); }),
			std::make_shared<DisplayMethod>("mz", [](const Euler* eqn, const Cons& U) -> float { return U.m()(2); }),
			
			std::make_shared<DisplayMethod>("ETotal", [](const Euler* eqn, const Cons& U) -> float { return U.ETotal(); }),
			
			std::make_shared<DisplayMethod>("v", [](const Euler* eqn, const Cons& U) -> float { return real3::length(U.m()) / U.rho(); }),
			std::make_shared<DisplayMethod>("vx", [](const Euler* eqn, const Cons& U) -> float { return U.m()(0) / U.rho(); }),
			std::make_shared<DisplayMethod>("vy", [](const Euler* eqn, const Cons& U) -> float { return U.m()(1) / U.rho(); }),
			std::make_shared<DisplayMethod>("vz", [](const Euler* eqn, const Cons& U) -> float { return U.m()(2) / U.rho(); }),
			
			std::make_shared<DisplayMethod>("P", [](const Euler* eqn, const Cons& U) -> float { return eqn->primFromCons(U).P(); }),
		};

		displayMethodNames = map<
			decltype(displayMethods),
			std::vector<const char*>
		>(
			displayMethods,
			[](const std::shared_ptr<DisplayMethod>& m) -> const char* { return m->name.c_str(); }
		);
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

	static constexpr real sqrt_1_2 = sqrt(.5);
	static constexpr real speedOfLight = 1;
	static constexpr real divPhiWavespeed = speedOfLight;
	static constexpr real divPsiWavespeed = speedOfLight;

	WaveVec getEigenvalues(const Eigen& vars) {
		real v_p = vars.sqrt_1_eps * vars.sqrt_1_mu;
		StateVec lambdas;
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

	std::pair<real, real> calcLambdaMinMax(const Cons& U) {
		real v_p = U.sqrt_1_eps * U.sqrt_1_mu;
		real lambda = std::max(std::max(divPsiWavespeed, divPhiWavespeed), v_p);
		return std::pair<real, real>(-lambda, lambda);
	}

	real calcLambdaMin(const Cons& U) {
		return -calcLambdaMax(U);
	}

	real calcLambdaMax(const Cons& U) {
		real v_p = U.sqrt_1_eps * U.sqrt_1_mu;
		return std::max(std::max(divPsiWavespeed, divPhiWavespeed), v_p);
	}

	static StateVec matmul(const real* A, const StateVec& x) {
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

	WaveVec apply_evL(Cons x, const Eigen& vars) {
		real sqrt_1_eps = vars.sqrt_1_eps;
		real sqrt_eps = 1. / sqrt_1_eps;
		real sqrt_1_mu = vars.sqrt_1_mu;
		real sqrt_mu = 1. / sqrt_1_mu;
		real v_p = sqrt_1_eps * sqrt_1_mu;

		Cons y;
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

	Cons apply_evR(WaveVec x, const Eigen& vars) {
		real sqrt_1_eps = vars.sqrt_1_eps;
		real sqrt_1_mu = vars.sqrt_1_mu;
		real sqrt_eps = 1. / vars.sqrt_1_eps;
		real sqrt_mu = 1. / vars.sqrt_1_mu;
		real sqrt_2 = 1. / sqrt_1_2;

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

	Cons calcFluxFromCons(Cons U) {
		real3 E = calc_E(U);
		real3 H = calc_H(U);
		Cons F;
		F.D = real3(U.phi * divPhiWavespeed,  H.z, -H.y);
		F.B = real3(U.psi * divPsiWavespeed, -E.z,  E.y);
		F.phi = U.D(0), divPhiWavespeed;
		F.psi = U.B(0), divPsiWavespeed;
		return F;
	}
};

}
}
