#pragma once

#include "CFDMesh/Equation/Equation.h"
#include "CFDMesh/Util.h"
#include "Tensor/Vector.h"
#include <cimgui.h>
#include <utility>
#include <cmath>
#include <cassert>

namespace CFDMesh {
namespace Equation {

template<typename real>
struct EulerNamespace {

using real3 = Tensor::Vector<real, 3>;

enum { numCons = 5 };

using StateVec_ = Tensor::Vector<real, numCons>;

struct Cons : public StateVec_ {
	real& rho() { return StateVec_::v[0]; }
	real3& m() { return *(real3*)( StateVec_::v + 1 ); }
	real& ETotal() { return StateVec_::v[StateVec_::size-1]; }
	
	const real& rho() const { return StateVec_::v[0]; }
	const real3& m() const { return *(real3*)( StateVec_::v + 1 ); }
	const real& ETotal() const { return StateVec_::v[StateVec_::size-1]; }

	Cons() {}
	
	Cons(const StateVec_& v) : StateVec_(v) {}

	Cons(real rho_, real3 m_, real ETotal_) {
		rho() = rho_;
		m() = m_;
		ETotal() = ETotal_;
	}
};

struct Prim_ : public StateVec_ {
	real& rho() { return StateVec_::v[0]; }
	real3& v() { return *(real3*)( StateVec_::v + 1 ); }
	real& P() { return StateVec_::v[StateVec_::size-1]; }

	const real& rho() const { return StateVec_::v[0]; }
	const real3& v() const { return *(real3*)( StateVec_::v + 1 ); }
	const real& P() const { return StateVec_::v[StateVec_::size-1]; }

	Prim_() {}
	
	Prim_(const StateVec_& v) : StateVec_(v) {}
	
	Prim_(real rho_, real3 v_, real P_) {
		rho() = rho_;
		v() = v_;
		P() = P_;
	}
};

struct Euler : public Equation<real, Cons, Euler> {
	using Super = Equation<real, Cons, Euler>;

	using Prim = Prim_;

	using StateVec = StateVec_;

	enum { numWaves = numCons };
	using WaveVec = StateVec;

	using InitCond = typename Super::InitCond;
	using DisplayMethod = typename Super::DisplayMethod;

	real heatCapacityRatio = 1.4;


	struct InitCondConst : public InitCond {
		using InitCond::InitCond;
		float rho = 1;
		float P = 1;
		float vx = 0;
		float vy = 0;

		virtual const char* name() const { return "constant"; }
		virtual Cons initCell(const Euler* eqn, real3 x) const {
			return eqn->consFromPrim(Prim_(rho, real3(vx, vy), P));
		}
		
		virtual void updateGUI() {
			igInputFloat("rho", &rho, .1, 1, "%f", 0);
			igInputFloat("vx", &vx, .1, 1, "%f", 0);
			igInputFloat("vy", &vy, .1, 1, "%f", 0);
			igInputFloat("P", &P, .1, 1, "%f", 0);
		}
	};

	struct InitCondSod : public InitCond {
		using InitCond::InitCond;
		float rhoL = 1;
		float vxL = 0;
		float vyL = 0;
		float vzL = 0;
		float PL = 1;
		float rhoR = .125;
		float vxR = 0;
		float vyR = 0;
		float vzR = 0;
		float PR = .1;
		virtual const char* name() const { return "Sod"; }
		virtual Cons initCell(const Euler* eqn, real3 x) const {
			bool lhs = x(0) < 0 && x(1) < 0;
			return eqn->consFromPrim(Prim_(
				lhs ? rhoL : rhoR,
				lhs ? real3(vxL, vyL, vzL) : real3(vxR, vyR, vzR),
				lhs ? PL : PR
			));
		}

		virtual void updateGUI() {
			igInputFloat("rhoL", &rhoL, .1, 1, "%f", 0);
			igInputFloat("vxL", &vxL, .1, 1, "%f", 0);
			igInputFloat("vyL", &vyL, .1, 1, "%f", 0);
			igInputFloat("vzL", &vzL, .1, 1, "%f", 0);
			igInputFloat("PL", &PL, .1, 1, "%f", 0);
			igInputFloat("rhoR", &rhoR, .1, 1, "%f", 0);
			igInputFloat("vxR", &vxR, .1, 1, "%f", 0);
			igInputFloat("vyR", &vyR, .1, 1, "%f", 0);
			igInputFloat("vzR", &vzR, .1, 1, "%f", 0);
			igInputFloat("PR", &PR, .1, 1, "%f", 0);
		}
	};

	struct InitCondSpiral : public InitCond {
		using InitCond::InitCond;
		virtual const char* name() const { return "Spiral"; }
		virtual Cons initCell(const Euler* eqn, real3 x) const {
			return eqn->consFromPrim(Prim_(
				1,
				real3(-x(1), x(0)),
				1
			));
		}
	};


	Euler() {
		Super::initConds = {
			std::make_shared<InitCondConst>(),
			std::make_shared<InitCondSod>(),
			std::make_shared<InitCondSpiral>(),
		};

		Super::addDisplayScalar("rho", [](const Euler* eqn, const Cons& U) -> float { return U.rho(); });
		Super::addDisplayVector("m", [](const Euler* eqn, const Cons& U) -> float3 { return (float3)U.m(); });
		Super::addDisplayScalar("ETotal", [](const Euler* eqn, const Cons& U) -> float { return U.ETotal(); });
		Super::addDisplayVector("v", [](const Euler* eqn, const Cons& U) -> float3 { return (float3)(U.m() / U.rho()); });
		Super::addDisplayScalar("P", [](const Euler* eqn, const Cons& U) -> float { return eqn->primFromCons(U).P(); });


		Super::getNames();
	}

	Cons consFromPrim(const Prim& W) const {
		Cons U;
		U.rho() = W.rho();
		U.m() = W.v() * U.rho();
		U.ETotal() = W.P() / (heatCapacityRatio - 1.) + .5 * U.rho() * real3::lenSq(W.v());
		return U;
	}

	Prim primFromCons(const Cons& U) const {
		Prim W;
		W.rho() = U.rho();
		W.v() = U.m() / W.rho();
		W.P() = (heatCapacityRatio - 1.) * (U.ETotal() - .5 * W.rho() * real3::lenSq(W.v()));
		return W;
	}

	real calc_hTotal(real rho, real P, real ETotal) const {
		return (ETotal + P) / rho;
	}

	real calc_Cs_from_P_rho(real P, real rho) const {
		return sqrt(heatCapacityRatio * P / rho);
	}

	real calc_Cs_from_vSq_hTotal(real vSq, real hTotal) const {
		return sqrt((heatCapacityRatio - 1) * (hTotal - .5 * vSq));
	}

	//variables used by the roe avg, eigen basis, etc
	struct Eigen {
		Prim WL, WR;
		real3 v;
		real vSq;
		real hTotal;
		real Cs;
		real CsSq;
	};

	Eigen calcRoeAvg(Cons UL, Cons UR) {
		Eigen vars;
		Prim& WL = vars.WL;
		Prim& WR = vars.WR;
		
		real ETotalL = UL.ETotal();
		WL = primFromCons(UL);
		real rhoL = WL.rho();
assert(rhoL > 0);			
		real3 vL = WL.v();
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
		real3 vR = WR.v();
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
		vars.vSq = real3::lenSq(vars.v);
		vars.hTotal = (sqrtRhoL * hTotalL + sqrtRhoR * hTotalR) * invDenom;
assert(std::isfinite(vars.hTotal));
		vars.Cs = calc_Cs_from_vSq_hTotal(vars.vSq, vars.hTotal);
assert(std::isfinite(vars.Cs));
		vars.CsSq = vars.Cs * vars.Cs;
		return vars;
	}

	StateVec getEigenvalues(const Eigen& vars) {
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

	struct CalcLambdaVars {
		real v;
		real Cs;
		
		CalcLambdaVars(const Euler& eqn, const Cons& U) {
			Prim W = eqn.primFromCons(U);
			v = real3::length(W.v());
			Cs = eqn.calc_Cs_from_P_rho(W.P(), W.rho());
		}
	};

	std::pair<real, real> calcLambdaMinMax(const CalcLambdaVars& vars) {
		return std::make_pair<real, real>(vars.v - vars.Cs, vars.v + vars.Cs);
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

	Cons apply_evL(Cons x, const Eigen& vars) {
		const real& Cs = vars.Cs;
		const real3& v = vars.v;
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

	Cons apply_evR(Cons x, const Eigen& vars) {
		const real& Cs = vars.Cs;
		const real3& v = vars.v;
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

	Cons rotateTo(Cons U, real3 normal) {
		U.m() = CFDMesh::rotateTo<real3>(U.m(), normal);
		return U;
	}
	
	Cons rotateFrom(Cons U, real3 normal) {
		U.m() = CFDMesh::rotateFrom<real3>(U.m(), normal);
		return U;
	}
};

};

}
}
