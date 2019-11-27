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
			return eqn->consFromPrim((Prim)W);
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
		U.ETotal = W.P / (heatCapacityRatio - 1.) + .5 * U.rho * W.v.lenSq();
		return U;
	}

	Prim primFromCons(const Cons& U) const {
		Prim W;
		W.rho = U.rho;
		W.v = U.m / W.rho;
		W.P = (heatCapacityRatio - 1.) * (U.ETotal - .5 * W.rho * W.v.lenSq());
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
		vars.vSq = vars.v.lenSq();
		vars.hTotal = (sqrtRhoL * hTotalL + sqrtRhoR * hTotalR) * invDenom;
		vars.Cs = calc_Cs_from_vSq_hTotal(vars.vSq, vars.hTotal);
		vars.CsSq = vars.Cs * vars.Cs;
		return vars;
	}

	WaveVec getEigenvalues(const Eigen& vars, real3 n) {
		real v = real3::dot(vars.v, n);
		const real& Cs = vars.Cs;
		WaveVec lambdas;
		lambdas(0) = v - Cs;
		lambdas(1) = v;
		lambdas(2) = v;
		lambdas(3) = v;
		lambdas(4) = v + Cs;
		return lambdas;
	}

	struct CalcLambdaVars {
		real3 n;
		real v;
		real Cs;
		
		CalcLambdaVars(const This& eqn, const Cons& U, const real3& n_) {
			Prim W = eqn.primFromCons(U);
			n = n_;
			v = real3::dot(W.v, n);
			Cs = eqn.calc_Cs_from_P_rho(W.P, W.rho);
		}
	
		CalcLambdaVars(const Eigen& vars, const real3& n_) : n(n_), v(real3::dot(vars.v, n)), Cs(vars.Cs) {}
	};

	std::pair<real, real> calcLambdaMinMax(const CalcLambdaVars& vars) const {
		return std::pair<real, real>(vars.v - vars.Cs, vars.v + vars.Cs);
	}

	real calcLambdaMin(const CalcLambdaVars& vars) const {
		return vars.v - vars.Cs;
	}

	real calcLambdaMax(const CalcLambdaVars& vars) const {
		return vars.v + vars.Cs;
	}


	struct BuildPerpendicularBasis {
		static void go(
			real3 normal, 
			Tensor::Vector<real3, 2> &tangents) 
		{
			//1) pick normal's max abs component
			//2) fill in all axii but that component
			//3) apply Graham-Schmidt
			// what about coordinate system handedness?

			int maxAxis = -1;
			real maxValue = -HUGE_VAL;
			for (int k = 0; k < dim; ++k) {
				real absNormal = fabs(normal(k));
				if (absNormal > maxValue) {
					maxValue = absNormal;
					maxAxis = k;
				}
			}

			for (int j = 0; j < dim-1; ++j) {
				if (j < maxAxis) {
					tangents(j)(j) = 1;
				} else {
					tangents(j)(j+1) = 1;
				}
			
				for (int k = j-1; k >= 0; --k) {
					real num = real(0), denom = real(0);
					for (int i = 0; i < dim; ++i) {
						num += tangents(j)(i) * tangents(k)(i);
						denom += tangents(j)(i) * tangents(j)(i);
					}
					tangents(j) -= tangents(j) * (num / denom);
				}
				{
					real num = real(0), denom = real(0);
					for (int i = 0; i < dim; ++i) {
						num += tangents(j)(i) * normal(i);
						denom += tangents(j)(i) * tangents(j)(i);
					}
					tangents(j) -= tangents(j) * (num / denom);
				}
				{
					real len = real(0);
					for (int k = 0; k < dim; ++k) {
						len += sqrt(tangents(j)(k) * tangents(j)(k));
					}
					tangents(j) *= real(1) / len;
				}
			}
		}
	};

	Cons apply_evL(Cons x, const Eigen& vars, real3 n) {
		const real& Cs = vars.Cs;
		const real3& v = vars.v;
		const real& vSq = v.lenSq();
	
		real vn = real3::dot(v, n);

		real CsSq = Cs * Cs;
		real gamma_1 = heatCapacityRatio - 1;

		real denom = 2. * CsSq;
		real invDenom = 1. / denom;

		Tensor::Vector<real3, 2> t;
		BuildPerpendicularBasis::go(n, t);

		Cons y;
		//min row	
		y.ptr[0] = (
				x.ptr[0] * (.5 * gamma_1 * vSq + Cs * vn)
				+ x.ptr[1] * -(n(0) * Cs + gamma_1 * v(0))
				+ x.ptr[2] * -(n(1) * Cs + gamma_1 * v(1))
				+ x.ptr[3] * -(n(2) * Cs + gamma_1 * v(2))
				+ x.ptr[4] * gamma_1
			) * invDenom;
		//mid normal row
		y.ptr[1] =
			x.ptr[0] * (1. - gamma_1 * vSq * invDenom)
			+ x.ptr[1] * (gamma_1 * v(0) * 2. * invDenom)
			+ x.ptr[2] * (gamma_1 * v(1) * 2. * invDenom)
			+ x.ptr[3] * (gamma_1 * v(2) * 2. * invDenom)
			+ x.ptr[4] * (-gamma_1 * 2. * invDenom);
		//mid tangent row
		y.ptr[2] =
			x.ptr[0] * -real3::dot(v, t(0))
			+ x.ptr[1] * t(0)(0)
			+ x.ptr[2] * t(0)(1)
			+ x.ptr[3] * t(0)(2);
		y.ptr[3] =
			x.ptr[0] * -real3::dot(v, t(1))
			+ x.ptr[1] * t(1)(0)
			+ x.ptr[2] * t(1)(1)
			+ x.ptr[3] * t(1)(2);
		//max row
		y.ptr[4] = (
				x.ptr[0] * (.5 * gamma_1 * vSq - Cs * vn)
				+ x.ptr[1] * (n(0) * Cs - gamma_1 * v(0))
				+ x.ptr[2] * (n(1) * Cs - gamma_1 * v(1))
				+ x.ptr[3] * (n(2) * Cs - gamma_1 * v(2))
				+ x.ptr[4] * gamma_1
			) * invDenom;

		return y;
	}

	Cons apply_evR(Cons x, const Eigen& vars, real3 n) {
		const real& Cs = vars.Cs;
		const real3& v = vars.v;
		const real& vSq = v.lenSq();
		const real& hTotal = vars.hTotal;

		real vn = real3::dot(v, n);
		
		Tensor::Vector<real3, 2> t;
		BuildPerpendicularBasis::go(n, t);

		Cons y;
		
		//min eigenvector
		y.ptr[0] =
			x.ptr[0]
			+ x.ptr[0] * (v(0) - Cs * n(0))
			+ x.ptr[0] * (v(1) - Cs * n(1))
			+ x.ptr[0] * (v(2) - Cs * n(2))
			+ x.ptr[0] * (hTotal - Cs * vn);
		//mid eigenvectors (n)
		y.ptr[1] =
			x.ptr[1]
			+ x.ptr[1] * (v(0))
			+ x.ptr[1] * (v(1))
			+ x.ptr[1] * (v(2))
			+ x.ptr[1] * (.5 * vSq);
		//mid eigenvectors (tangents)
		y.ptr[2] =
			x.ptr[2] * t(0)(0)
			+ x.ptr[2] * t(0)(1)
			+ x.ptr[2] * t(0)(2)
			+ x.ptr[2] * real3::dot(v, t(0));
		y.ptr[3] =
			x.ptr[3] * t(1)(0)
			+ x.ptr[3] * t(1)(1)
			+ x.ptr[3] * t(1)(2)
			+ x.ptr[3] * real3::dot(v, t(1));
		//max eigenvector
		y.ptr[4] =
			x.ptr[4]
			+ x.ptr[4] * (v(0) + Cs * n(0))
			+ x.ptr[4] * (v(1) + Cs * n(1))
			+ x.ptr[4] * (v(2) + Cs * n(2))
			+ x.ptr[4] * (hTotal + Cs * vn);
		
		return y;
	}

	Cons calcFluxFromCons(Cons U, real3 n) {
		Prim W = primFromCons(U);
		real hTotal = calc_hTotal(W.rho, W.P, U.ETotal);
		Cons flux;
		real m = real3::dot(U.m, n);
		flux(0) = m;
		flux(1) = m * W.v(0) + W.P * n(0);
		flux(2) = m * W.v(1) + W.P * n(1);
		flux(3) = m * W.v(2) + W.P * n(2);
		flux(4) = m * hTotal;
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
