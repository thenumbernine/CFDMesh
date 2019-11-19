#pragma once

#include "CFDMesh/Equation/Equation.h"
#include "CFDMesh/Util.h"
#include "Tensor/Vector.h"
#include <cimgui.h>
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
		std::make_pair("rho", &Cons_::rho),
		std::make_pair("m", &Cons_::m),
		std::make_pair("ETotal", &Cons_::ETotal)
	);
};

template<typename T>
ADD_OSTREAM(Cons_<T>)


template<typename real>
union Prim_ {
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

	//TODO autogen from fields
	void updateGUI(std::string suffix = {}) {
		igInputFloat(("rho" + suffix).c_str(), &rho, .1, 1, "%f", 0);
		igInputFloat(("vx" + suffix).c_str(), &v(0), .1, 1, "%f", 0);
		igInputFloat(("vy" + suffix).c_str(), &v(1), .1, 1, "%f", 0);
		igInputFloat(("vz" + suffix).c_str(), &v(2), .1, 1, "%f", 0);
		igInputFloat(("P" + suffix).c_str(), &P, .1, 1, "%f", 0);
	}

	static constexpr auto fields = std::make_tuple(
		std::make_pair("rho", &Prim_::rho),
		std::make_pair("v", &Prim_::v),
		std::make_pair("P", &Prim_::P)
	);
};

template<typename T>
ADD_OSTREAM(Prim_<T>)


template<typename real>
struct Euler : public Equation<real, Cons_<real>, Euler<real>> {
	using Parent = Equation<real, Cons_<real>, Euler<real>>;
	using Cons = typename Parent::Cons;
	using Prim = Prim_<real>;

	enum { numWaves = numCons };
	using WaveVec = Cons;

	using InitCond = typename Parent::InitCond;
	using DisplayMethod = typename Parent::DisplayMethod;
	
	using real3 = Tensor::Vector<real, 3>;

	float heatCapacityRatio = 1.4;

	struct InitCondConst : public InitCond {
		Prim_<float> W = Prim_<float>(1, float3(), 1);
		using InitCond::InitCond;
		virtual const char* name() const { return "constant"; }
		virtual Cons initCell(const Euler* eqn, real3 x) const {	
			
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
			W.updateGUI();
		}
	};

	struct InitCondSod : public InitCond {
		using InitCond::InitCond;
		Prim_<float> WL = Prim_<float>(1, float3(), 1);
		Prim_<float> WR = Prim_<float>(.125, float3(), .1);
		
		virtual const char* name() const { return "Sod"; }
		virtual Cons initCell(const Euler* eqn, real3 x) const {
			bool lhs = x(0) < 0 && x(1) < 0;
			return eqn->consFromPrim(lhs ? WL : WR);
		}

		virtual void updateGUI() {
			WL.updateGUI("L");
			WR.updateGUI("R");
		}
	};

	struct InitCondSpiral : public InitCond {
		using InitCond::InitCond;
		virtual const char* name() const { return "Spiral"; }
		virtual Cons initCell(const Euler* eqn, real3 x) const {
			return eqn->consFromPrim(Prim(1, real3(-x(1), x(0)), 1));
		}
	};

	using Parent::Parent;

	void buildInitCondsAndDisplayVars() {
		Parent::initConds = {
			std::make_shared<InitCondConst>(),
			std::make_shared<InitCondSod>(),
			std::make_shared<InitCondSpiral>(),
		};

		
		//cycle through the Cons::fields tuple and add each of these
		tuple_for_each(Cons::fields, [this](auto x, size_t i) constexpr {
			auto field = x.second;
			using FieldType = typename MemberPointerInfo<decltype(field)>::FieldType;
			Parent::template addDisplayForType<FieldType>(
				x.first, 
				[field](const Euler* eqn, const Cons& U) -> typename FloatTypeForType<FieldType>::Type { 
					return U.*field; 
				}
			);
		});

		if constexpr (!std::is_same_v<Cons, Prim>) {
			tuple_for_each(Prim::fields, [this](auto x, size_t i) constexpr {
				auto field = x.second;
				using FieldType = typename MemberPointerInfo<decltype(field)>::FieldType;
				Parent::template addDisplayForType<FieldType>(
					x.first, 
					[field](const Euler* eqn, const Cons& U) -> typename FloatTypeForType<FieldType>::Type { 
						Prim W = eqn->primFromCons(U);
						return W.*field; 
					}
				);
			});	
		}
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
		
		CalcLambdaVars(const Euler& eqn, const Cons& U) {
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

	static Cons matmul(const real* A, const Cons& x) {
		Cons y;
		for (int i = 0; i < Cons::size; ++i) {
			real sum = 0;
			for (int j = 0; j < Cons::size; ++j) {
				//C layout, so row-major
				sum += A[j + Cons::size * i] * x(j);
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
		real hTotal = calc_hTotal(W.rho, W.P, U.ETotal);
		Cons flux;
		flux(0) = U.m(0);
		flux(1) = U.m(0) * W.v(0) + W.P;
		flux(2) = U.m(0) * W.v(1);
		flux(3) = U.m(0) * W.v(2);
		flux(4) = U.m(0) * hTotal;
		return flux;
	}

	void updateGUI() {
		igInputFloat("gamma", &heatCapacityRatio, .1, 1, "%f", 0);
	}

	//TODO enumerate vectors in Cons struct
	//automatically rotate and reflect them here
	//also automatically add them as vectors to display vars

	Cons rotateTo(Cons U, real3 normal) {
		U.m = CFDMesh::rotateTo<real3>(U.m, normal);
		return U;
	}
	
	Cons rotateFrom(Cons U, real3 normal) {
		U.m = CFDMesh::rotateFrom<real3>(U.m, normal);
		return U;
	}

	Cons reflect(Cons U, real3 normal, real restitution) const {
		U.m = U.m - normal * ((1 + restitution) * real3::dot(normal, U.m));
		return U;
	}
};

}
}
}
