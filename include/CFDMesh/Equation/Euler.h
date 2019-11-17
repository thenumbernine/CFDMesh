#pragma once

#include "CFDMesh/Equation/Equation.h"
#include "CFDMesh/Util.h"
#include "Tensor/Vector.h"
#include <cimgui.h>
#include <utility>
#include <memory>
#include <cmath>
#include <cassert>

namespace CFDMesh {
namespace Equation {

template<typename Config>
struct Euler : public Equation<Config> {
	using real = typename Config::real;
	using real3 = typename Config::real3;
	using StateVec = Tensor::Vector<real, 5>;

	enum { numWaves = 5 };
	using WaveVec = StateVec;

	real heatCapacityRatio = 1.4;

	struct Cons : public StateVec {
		real& rho() { return StateVec::v[0]; }
		real3& m() { return *(real3*)( StateVec::v + 1 ); }
		real& ETotal() { return StateVec::v[StateVec::size-1]; }
		
		const real& rho() const { return StateVec::v[0]; }
		const real3& m() const { return *(real3*)( StateVec::v + 1 ); }
		const real& ETotal() const { return StateVec::v[StateVec::size-1]; }

		Cons() {}
		
		Cons(const StateVec& v) : StateVec(v) {}

		Cons(real rho_, real3 m_, real ETotal_) {
			rho() = rho_;
			m() = m_;
			ETotal() = ETotal_;
		}
	};

	struct Prim : public StateVec {
		real& rho() { return StateVec::v[0]; }
		real3& v() { return *(real3*)( StateVec::v + 1 ); }
		real& P() { return StateVec::v[StateVec::size-1]; }

		const real& rho() const { return StateVec::v[0]; }
		const real3& v() const { return *(real3*)( StateVec::v + 1 ); }
		const real& P() const { return StateVec::v[StateVec::size-1]; }

		Prim() {}
		
		Prim(const StateVec& v) : StateVec(v) {}
		
		Prim(real rho_, real3 v_, real P_) {
			rho() = rho_;
			v() = v_;
			P() = P_;
		}
	};

	struct InitCond {
		virtual ~InitCond() {}
		virtual const char* name() const = 0;
		virtual Cons initCell(const Euler* eqn, real3 pos) const = 0;
		virtual void updateGUI() {}
	};

	struct InitCondConst : public InitCond {
		using InitCond::InitCond;
		float rho = 1;
		float P = 1;
		float vx = 0;
		float vy = 0;

		virtual const char* name() const { return "constant"; }
		virtual Cons initCell(const Euler* eqn, real3 x) const {
			return eqn->consFromPrim(Prim(rho, real3(vx, vy), P));
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
			return eqn->consFromPrim(Prim(
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
			return eqn->consFromPrim(Prim(
				1,
				real3(-x(1), x(0)),
				1
			));
		}
	};

	std::vector<std::shared_ptr<InitCond>> initConds;
	std::vector<const char*> initCondNames;


	struct DisplayMethod {
		using DisplayFunc = std::function<float(const Euler*, const Cons&)>;
		std::string name;
		DisplayFunc f;
		DisplayMethod(const std::string& name_, DisplayFunc f_) : name(name_), f(f_) {}
	};

	std::vector<std::shared_ptr<DisplayMethod>> displayMethods;
	std::vector<const char*> displayMethodNames;


	Euler() {
		initConds = {
			std::make_shared<InitCondConst>(),
			std::make_shared<InitCondSod>(),
			std::make_shared<InitCondSpiral>(),
		};

		initCondNames = map<
			decltype(initConds),
			std::vector<const char*>
		>(
			initConds,
			[](std::shared_ptr<InitCond> ic) -> const char* { return ic->name(); }
		);
	
		displayMethods = {
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

	std::pair<real, real> calcLambdaMinMax(real3 normal, Prim W, real Cs) {
		real v = real3::dot(normal, W.v());
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
};

}
}
