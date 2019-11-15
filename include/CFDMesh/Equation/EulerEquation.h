#pragma once

#include "CFDMesh/Equation/Equation.h"

#include <utility>
#include <cmath>
#include <cassert>

namespace CFDMesh {
namespace Equation {

template<typename Config>
struct EulerEquation : public Equation<Config> {
	using real = typename Config::real;
	using vec = typename Config::vec;
	using StateVec = typename Config::StateVec;

	real heatCapacityRatio = 1.4;

	struct Prim;

	struct Cons : public StateVec {
		real& rho() { return StateVec::v[0]; }
		vec& m() { return *(vec*)( StateVec::v + 1 ); }
		real& ETotal() { return StateVec::v[StateVec::size-1]; }
		
		const real& rho() const { return StateVec::v[0]; }
		const vec& m() const { return *(vec*)( StateVec::v + 1 ); }
		const real& ETotal() const { return StateVec::v[StateVec::size-1]; }

		Cons() {}
		
		Cons(const StateVec& v) : StateVec(v) {}

		Cons(real rho_, vec m_, real ETotal_) {
			rho() = rho_;
			m() = m_;
			ETotal() = ETotal_;
		}
	};

	struct Prim : public StateVec {
		real& rho() { return StateVec::v[0]; }
		vec& v() { return *(vec*)( StateVec::v + 1 ); }
		real& P() { return StateVec::v[StateVec::size-1]; }

		const real& rho() const { return StateVec::v[0]; }
		const vec& v() const { return *(vec*)( StateVec::v + 1 ); }
		const real& P() const { return StateVec::v[StateVec::size-1]; }

		Prim() {}
		
		Prim(const StateVec& v) : StateVec(v) {}
		
		Prim(real rho_, vec v_, real P_) {
			rho() = rho_;
			v() = v_;
			P() = P_;
		}
	};

	Cons consFromPrim(const Prim& W) {
		Cons U;
		U.rho() = W.rho();
		U.m() = W.v() * U.rho();
		U.ETotal() = W.P() / (heatCapacityRatio - 1.) + .5 * U.rho() * vec::lenSq(W.v());
		return U;
	}

	Prim primFromCons(const Cons& U) {
		Prim W;
		W.rho() = U.rho();
		W.v() = U.m() / W.rho();
		W.P() = (heatCapacityRatio - 1.) * (U.ETotal() - .5 * W.rho() * vec::lenSq(W.v()));
		return W;
	}

	real calc_hTotal(real rho, real P, real ETotal) {
		return (ETotal + P) / rho;
	}

	real calc_Cs_from_P_rho(real P, real rho) {
		return sqrt(heatCapacityRatio * P / rho);
	}

	real calc_Cs_from_vSq_hTotal(real vSq, real hTotal) {
		return sqrt((heatCapacityRatio - 1) * (hTotal - .5 * vSq));
	}

	//variables used by the roe avg, eigen basis, etc
	struct RoeAvg {
		Prim WL, WR;
		vec v;
		real vSq;
		real hTotal;
		real Cs;
		real CsSq;
	};

	RoeAvg calcRoeAvg(Cons UL, Cons UR) {
		RoeAvg vars;
		Prim& WL = vars.WL;
		Prim& WR = vars.WR;
		
		real ETotalL = UL.ETotal();
		WL = primFromCons(UL);
		real rhoL = WL.rho();
assert(rhoL > 0);			
		vec vL = WL.v();
		real PL = WL.P();
assert(PL > 0);			
		real hTotalL = calc_hTotal(rhoL, PL, ETotalL);
assert(hTotalL > 0);			
		real sqrtRhoL = sqrt(rhoL);
assert(std::isfinite(sqrtRhoL));
		
		real ETotalR = UR.ETotal();
		WR = primFromCons(UR);
		real rhoR = WR.rho();
assert(rhoR > 0);			
		vec vR = WR.v();
		real PR = WR.P();
assert(PR > 0);			
		real hTotalR = calc_hTotal(rhoR, PR, ETotalR);
assert(hTotalR > 0);			
		real sqrtRhoR = sqrt(rhoR);
assert(std::isfinite(sqrtRhoR));
		
		real invDenom = 1. / (sqrtRhoL + sqrtRhoR);

		vars.v = (vL * sqrtRhoL + vR * sqrtRhoR) * invDenom;
assert(std::isfinite(vars.v(0)));
assert(std::isfinite(vars.v(1)));
assert(std::isfinite(vars.v(2)));
		vars.vSq = vec::lenSq(vars.v);
		vars.hTotal = (sqrtRhoL * hTotalL + sqrtRhoR * hTotalR) * invDenom;
assert(std::isfinite(vars.hTotal));
		vars.Cs = calc_Cs_from_vSq_hTotal(vars.vSq, vars.hTotal);
assert(std::isfinite(vars.Cs));
		vars.CsSq = vars.Cs * vars.Cs;
		return vars;
	}

	StateVec getEigenvalues(const RoeAvg& vars) {
		const real& vx = vars.v(0);
		const real& Cs = vars.Cs;
		StateVec lambdas;
		lambdas(0) = vx - Cs;
		lambdas(1) = vx;
		lambdas(2) = vx;
		lambdas(3) = vx;
		lambdas(4) = vx + Cs;
		return lambdas;
	}

	std::pair<real, real> calcLambdaMinMax(vec normal, Prim W, real Cs) {
		real v = vec::dot(normal, W.v());
		return std::make_pair<real, real>(v - Cs, v + Cs);
	}

	real calcLambdaMin(real v, real Cs) {
		return v - Cs;
	}

	real calcLambdaMax(real v, real Cs) {
		return v + Cs;
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

	Cons apply_evL(Cons x, const RoeAvg& vars) {
		const real& Cs = vars.Cs;
		const vec& v = vars.v;
		const real& vSq = vars.vSq;
		const real& vx = v(0);
		const real& vy = v(1);
		const real& vz = v(2);
		real CsSq = Cs * Cs;
		real gamma_1 = heatCapacityRatio - 1;
		real evL[5][5] = {
			{(.5 * gamma_1 * vSq + Cs * vx) / (2 * CsSq),	(-Cs - gamma_1 * vx) / (2 * CsSq),	-gamma_1 * vy / (2 * CsSq),		-gamma_1 * vz / (2 * CsSq),	gamma_1 / (2 * CsSq),	},
			{1 - gamma_1 * vSq / (2 * CsSq),				gamma_1 * vx / CsSq,				gamma_1 * vy / CsSq,			gamma_1 * vz / CsSq,		-gamma_1 / CsSq,		},
			{-vy,											0,									1,								0,							0,						}, 
			{-vz,											0,									0,								1,							0,						},
			{(.5 * gamma_1 * vSq - Cs * vx) / (2 * CsSq),	(Cs - gamma_1 * vx) / (2 * CsSq),	-gamma_1 * vy / (2 * CsSq),		-gamma_1 * vz / (2 * CsSq),	gamma_1 / (2 * CsSq),	},
		};
		return matmul(&evL[0][0], x);
	}

	Cons apply_evR(Cons x, const RoeAvg& vars) {
		const real& Cs = vars.Cs;
		const vec& v = vars.v;
		const real& vSq = vars.vSq;
		const real& hTotal = vars.hTotal;
		const real& vx = v(0);
		const real& vy = v(1);
		const real& vz = v(2);
		real evR[5][5] = {
			{1, 				1, 			0,		0,		1,				},
			{vx - Cs, 			vx, 		0,		0,		vx + Cs,		},
			{vy,				vy,			1,		0,		vy,				},
			{vz,				vz,			0,		1,		vz,				},
			{hTotal - Cs * vx, .5 * vSq, 	vy,		vz,		hTotal + Cs * vx},
		};
		return matmul(&evR[0][0], x);
	}

	Cons calcFluxFromCons(Cons U) {
		Prim W = primFromCons(U);
		real hTotal = calc_hTotal(W.rho(), W.P(), U.ETotal());
		Cons flux;
		flux(0) = U.m()(0);
		flux(1) = U.m()(0) * W.v()(0) + W.P();
		flux(2) = U.m()(0) * W.v()(1);
		flux(3) = U.m()(0) * W.v()(2);
		flux(4) = U.m()(0) * hTotal;
		return flux;
	}
};

}
}
