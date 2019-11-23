#pragma once

#include "CFDMesh/Equation/Equation.h"
#include "CFDMesh/GUI.h"
#include "CFDMesh/Util.h"
#include "Tensor/Vector.h"
#include "Common/Macros.h"
#include <utility>
#include <tuple>
#include <cmath>
#include <cassert>

namespace CFDMesh {
namespace Equation {
namespace Euler {

enum { numCons = 5 };

template<typename real>
union Cons_ {
	using This = Cons_;
	using real3 = Tensor::Vector<real, 3>;
	enum { size = numCons };
	real ptr[size];
	struct {
		real rho = {};
		real3 m = {};
		real ETotal = {};
	};

	Cons_() {}

	Cons_(real rho_, real3 m_, real ETotal_) {
		rho = rho_;
		m = m_;
		ETotal = ETotal_;
	}

	ADD_OPS(Cons_)

	static constexpr auto fields = std::make_tuple(
		std::make_pair("rho", &This::rho),
		std::make_pair("m", &This::m),
		std::make_pair("ETotal", &This::ETotal)
	);
};

template<typename T>
ADD_OSTREAM(Cons_<T>)


template<typename real>
union Prim_ {
	using This = Prim_;
	using real3 = Tensor::Vector<real, 3>;
	enum { size = numCons };
	real ptr[size];
	struct {
		real rho = {};
		real3 v = {};
		real P = {};
	};

	Prim_() {}

	Prim_(real rho_, real3 v_, real P_) {
		rho = rho_;
		v = v_;
		P = P_;
	}

	ADD_OPS(Prim_)

	static constexpr auto fields = std::make_tuple(
		std::make_pair("rho", &This::rho),
		std::make_pair("v", &This::v),
		std::make_pair("P", &This::P)
	);
};

template<typename T>
ADD_OSTREAM(Prim_<T>)


template<typename real, int dim_>
struct Euler : public Equation<Euler<real, dim_>, real, Cons_<real>, Prim_<real>> {
	using This = Euler;
	using Super = Equation<Euler<real, dim_>, real, Cons_<real>, Prim_<real>>;
	using Cons = typename Super::Cons;
	using Prim = typename Super::Prim;

	enum { dim = dim_ };

	enum { numWaves = numCons };
	using WaveVec = Cons;

	using InitCond = typename Super::InitCond;
	using DisplayMethod = typename Super::DisplayMethod;
	
	using real3 = Tensor::Vector<real, 3>;

	float heatCapacityRatio = 1.4;

	struct InitCondConst : public InitCond {
		using InitCond::InitCond;
		
		Prim W = Prim(1, float3(), 1);
		
		static constexpr auto fields = std::make_tuple(
			std::make_pair("W", &InitCondConst::W)
		);

		virtual const char* name() const { return "constant"; }
		virtual Cons initCell(const This* eqn, real3 x) const {	
			
			assert(std::isfinite(W.rho) && W.rho > 0);
			assert(std::isfinite(W.v(0)));
			assert(std::isfinite(W.v(1)));
			assert(std::isfinite(W.v(2)));
			assert(std::isfinite(W.P) && W.P > 0);
			
			Cons U = eqn->consFromPrim((Prim)W);
			
			assert(std::isfinite(U.rho) && U.rho > 0);
			assert(std::isfinite(U.m(0)));
			assert(std::isfinite(U.m(1)));
			assert(std::isfinite(U.m(2)));
			assert(std::isfinite(U.ETotal) && U.ETotal > 0);
			
			return U;
		}
		
		virtual void updateGUI() {
			CFDMesh::updateGUI(this);
		}
	};

	struct InitCondSod : public InitCond {
		using InitCond::InitCond;
		
		Prim WL = Prim(1, float3(), 1);
		Prim WR = Prim(.125, float3(), .1);

		static constexpr auto fields = std::make_tuple(
			std::make_pair("WL", &InitCondSod::WL),
			std::make_pair("WR", &InitCondSod::WR)
		);

		virtual const char* name() const { return "Sod"; }
		virtual Cons initCell(const This* eqn, real3 x) const {
			bool lhs = x(0) < 0 && x(1) < 0 && (
				This::dim == 3 ? x(2) < 0 : true
			);
			return eqn->consFromPrim(lhs ? WL : WR);
		}

		virtual void updateGUI() {
			CFDMesh::updateGUI(this);
		}
	};

	struct InitCondKelvinHelmholtz : public InitCond {
		using InitCond::InitCond;

		float rhoIn = 2;
		float rhoOut = 1;
		float velYAmp = 1e-2;
		float frequency = 2;
		float randomNoiseAmp = 1e-4;
		float backgroundPressure = 2.5;
		float velocity = .5;
		
		Prim Win = Prim(2, float3(-.5, 0), 2.5);
		Prim Wout = Prim(1, float3(.5, 0), 2.5);
	
		static constexpr auto fields = std::make_tuple(
			std::make_pair("rhoIn", &InitCondKelvinHelmholtz::rhoIn),
			std::make_pair("rhoOut", &InitCondKelvinHelmholtz::rhoOut),
			std::make_pair("velYAmp", &InitCondKelvinHelmholtz::velYAmp),
			std::make_pair("frequency", &InitCondKelvinHelmholtz::frequency),
			std::make_pair("randomNoiseAmp", &InitCondKelvinHelmholtz::randomNoiseAmp),
			std::make_pair("backgroundPressure", &InitCondKelvinHelmholtz::backgroundPressure),
			std::make_pair("velocity", &InitCondKelvinHelmholtz::velocity)
		);

		virtual const char* name() const { return "Kelvin-Helmholtz"; }

		virtual Cons initCell(const This* eqn, real3 x) const {
			//TODO get mins & maxs from mesh?
			real3 mins = real3(-1);
			real3 maxs = real3(1);
			bool inside = x(1) > -.5 && x(1) < .5;
			Prim W;
			W.rho = inside ? rhoIn : rhoOut;
			float theta = frequency * 2 * M_PI;
			theta *= (x(0) - mins(0)) / (maxs(0) - mins(0));
			W.v(0) = (inside ? -1 : 1) * velocity;
			float noise = (maxs(0) - mins(1)) * velYAmp;
			W.v(1) = noise * sin(theta);
			W.P = backgroundPressure;
			for (int j = 0; j < Prim::size; ++j) {
				W.ptr[j] += randomNoiseAmp * crand();
			}
			return eqn->consFromPrim(W);
		}
		
		virtual void updateGUI() {
			CFDMesh::updateGUI(this);
		}
	};

	struct InitCondSpiral : public InitCond {
		using InitCond::InitCond;
		
		virtual const char* name() const { return "Spiral"; }
		virtual Cons initCell(const This* eqn, real3 x) const {
			return eqn->consFromPrim(Prim(1, real3(-x(1), x(0)), 1));
		}
	};

	using Super::Super;

	void buildInitCondsAndDisplayVars() {
		Super::initConds = {
			std::make_shared<InitCondSod>(),
			std::make_shared<InitCondConst>(),
			std::make_shared<InitCondKelvinHelmholtz>(),
			std::make_shared<InitCondSpiral>(),
		};
	}

	Cons consFromPrim(const Prim& W) const {
		Cons U;
		U.rho = W.rho;
		U.m = W.v * U.rho;
		U.ETotal = W.P / (heatCapacityRatio - 1.) + .5 * U.rho * real3::lenSq(W.v);
		return U;
	}

	Prim primFromCons(const Cons& U) const {
		Prim W;
		W.rho = U.rho;
		W.v = U.m / W.rho;
		W.P = (heatCapacityRatio - 1.) * (U.ETotal - .5 * W.rho * real3::lenSq(W.v));
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
		
		real ETotalL = UL.ETotal;
		WL = primFromCons(UL);
		real rhoL = WL.rho;
		real3 vL = WL.v;
		real PL = WL.P;
		real hTotalL = calc_hTotal(rhoL, PL, ETotalL);
		real sqrtRhoL = sqrt(rhoL);
		
		real ETotalR = UR.ETotal;
		WR = primFromCons(UR);
		real rhoR = WR.rho;
		real3 vR = WR.v;
		real PR = WR.P;
		real hTotalR = calc_hTotal(rhoR, PR, ETotalR);
		real sqrtRhoR = sqrt(rhoR);
		
		real invDenom = 1. / (sqrtRhoL + sqrtRhoR);

		vars.v = (vL * sqrtRhoL + vR * sqrtRhoR) * invDenom;
		vars.vSq = real3::lenSq(vars.v);
		vars.hTotal = (sqrtRhoL * hTotalL + sqrtRhoR * hTotalR) * invDenom;
		vars.Cs = calc_Cs_from_vSq_hTotal(vars.vSq, vars.hTotal);
		vars.CsSq = vars.Cs * vars.Cs;
		return vars;
	}

	WaveVec getEigenvalues(const Eigen& vars) {
		const real& vx = vars.v(0);
		const real& Cs = vars.Cs;
		WaveVec lambdas;
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
		
		CalcLambdaVars(const This& eqn, const Cons& U) {
			Prim W = eqn.primFromCons(U);
			v = real3::length(W.v);
			Cs = eqn.calc_Cs_from_P_rho(W.P, W.rho);
		}
	
		CalcLambdaVars(const Eigen& vars) : v(vars.v(0)), Cs(vars.Cs) {}
	};

	std::pair<real, real> calcLambdaMinMax(const CalcLambdaVars& vars) const {
		return std::make_pair<real, real>(vars.v - vars.Cs, vars.v + vars.Cs);
	}

	real calcLambdaMin(const CalcLambdaVars& vars) const {
		return vars.v - vars.Cs;
	}

	real calcLambdaMax(const CalcLambdaVars& vars) const {
		return vars.v + vars.Cs;
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
		Cons y;
		y.ptr[0] = x.ptr[0]*((.5 * gamma_1 * vSq + Cs * vx) / (2 * CsSq)) + x.ptr[1]*((-Cs - gamma_1 * vx) / (2 * CsSq)) + x.ptr[2]*(-gamma_1 * vy / (2 * CsSq)) + x.ptr[3]*(-gamma_1 * vz / (2 * CsSq)) + x.ptr[4]*(gamma_1 / (2 * CsSq));
		y.ptr[1] = x.ptr[0]*(1 - gamma_1 * vSq / (2 * CsSq)) + x.ptr[1]*(gamma_1 * vx / CsSq) + x.ptr[2]*(gamma_1 * vy / CsSq) + x.ptr[3]*(gamma_1 * vz / CsSq) + x.ptr[4]*(-gamma_1 / CsSq);
		y.ptr[2] = x.ptr[0]*(-vy) + x.ptr[1]*(0) + x.ptr[2]*(1) + x.ptr[3]*(0) + x.ptr[4]*(0); 
		y.ptr[3] = x.ptr[0]*(-vz) + x.ptr[1]*(0) + x.ptr[2]*(0) + x.ptr[3]*(1) + x.ptr[4]*(0);
		y.ptr[4] = x.ptr[0]*((.5 * gamma_1 * vSq - Cs * vx) / (2 * CsSq)) + x.ptr[1]*((Cs - gamma_1 * vx) / (2 * CsSq)) + x.ptr[2]*(-gamma_1 * vy / (2 * CsSq)) + x.ptr[3]*(-gamma_1 * vz / (2 * CsSq)) + x.ptr[4]*(gamma_1 / (2 * CsSq));
		return y;
	}

	Cons apply_evR(Cons x, const Eigen& vars) {
		const real& Cs = vars.Cs;
		const real3& v = vars.v;
		const real& vSq = vars.vSq;
		const real& hTotal = vars.hTotal;
		const real& vx = v(0);
		const real& vy = v(1);
		const real& vz = v(2);
		Cons y;
		y.ptr[0] = x.ptr[0]*(1) + x.ptr[1]*(1) + x.ptr[2]*(0) + x.ptr[3]*(0) + x.ptr[4]*(1);
		y.ptr[1] = x.ptr[0]*(vx - Cs) + x.ptr[1]*(vx) + x.ptr[2]*(0) + x.ptr[3]*(0) + x.ptr[4]*(vx + Cs);
		y.ptr[2] = x.ptr[0]*(vy) + x.ptr[1]*(vy) + x.ptr[2]*(1) + x.ptr[3]*(0) + x.ptr[4]*(vy);
		y.ptr[3] = x.ptr[0]*(vz) + x.ptr[1]*(vz) + x.ptr[2]*(0) + x.ptr[3]*(1) + x.ptr[4]*(vz);
		y.ptr[4] = x.ptr[0]*(hTotal - Cs * vx) + x.ptr[1]*(.5 * vSq) + x.ptr[2]*(vy) + x.ptr[3]*(vz) + x.ptr[4]*(hTotal + Cs * vx);
		return y;
	}

	Cons calcFluxFromCons(Cons U) {
		Prim W = primFromCons(U);
		real hTotal = calc_hTotal(W.rho, W.P, U.ETotal);
		Cons flux;
		flux(0) = U.m(0);
		flux(1) = U.m(0) * W.v(0) + W.P;
		flux(2) = U.m(0) * W.v(1);
		flux(3) = U.m(0) * W.v(2);
		flux(4) = U.m(0) * hTotal;
		return flux;
	}

	static constexpr auto fields = std::make_tuple(
		std::make_pair("heatCapacityRatio", &This::heatCapacityRatio)
	);

	void updateGUI() {
		CFDMesh::updateGUI(this);
	}
};

}
}
}
