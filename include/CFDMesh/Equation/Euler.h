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

//hmm, fields enable_if detect_v works for Eigen struct below, but not for this ... why so?
template<typename real>
std::ostream& operator<<(std::ostream& o, const Cons_<real>& b) {
	o << "[";
	Common::TupleForEach(Cons_<real>::fields, [&o, &b](auto x, size_t i) constexpr {
		if (i > 0) o << ", ";
		o << std::get<0>(x) << "=" << b.*(std::get<1>(x));
	});
	return o << "]";
}


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

//hmm, fields enable_if detect_v works for Eigen struct below, but not for this ... why so?
template<typename real>
std::ostream& operator<<(std::ostream& o, const Prim_<real>& b) {
	o << "[";
	Common::TupleForEach(Prim_<real>::fields, [&o, &b](auto x, size_t i) constexpr {
		if (i > 0) o << ", ";
		o << std::get<0>(x) << "=" << b.*(std::get<1>(x));
	});
	return o << "]";
}


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
		using This = Eigen;
		
		Prim WL, WR;
		real3 v;
		real vSq;
		real hTotal;
		real Cs;
		real CsSq;
	
		static constexpr auto fields = std::make_tuple(
			std::make_pair("WL", &This::WL),
			std::make_pair("WR", &This::WR),
			std::make_pair("v", &This::v),
			std::make_pair("vSq", &This::vSq),
			std::make_pair("hTotal", &This::hTotal),
			std::make_pair("Cs", &This::Cs),
			std::make_pair("CsSq", &This::CsSq)
		);
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


	WaveVec applyEigL(Cons X, const Eigen& vars, real3 n) {
		const real& Cs = vars.Cs;
		const real3& v = vars.v;	//v^i ... upper
		const auto& vL = v;			//v_i ... lower (not left)
		const real& vSq = real3::dot(v, vL);

		real CsSq = Cs * Cs;
		real heatRatioMinusOne = heatCapacityRatio - 1;

		real denom = 2. * CsSq;
		real invDenom = 1. / denom;

		real3 nSq = n * n;	//per-element multiply

		const real gU_xx = 1;
		const real gU_yy = 1;
		const real gU_zz = 1;
		const real gU_xy = 0;
		const real gU_xz = 0;
		const real gU_yz = 0;
		const real sqrt_gUxx = 1.;//sqrt_gUxx;
		const real sqrt_gUyy = 1.;//sqrt_gUjj;
		const real sqrt_gUzz = 1.;//sqrt_gUjj;

		WaveVec Y;

		//x dir
		Y.ptr[0] += nSq(0) * (
			(
				X.ptr[0] * (.5 * heatRatioMinusOne * vSq + Cs * v(0) / sqrt_gUxx)
				+ X.ptr[1] * (-heatRatioMinusOne * vL(0) - Cs / sqrt_gUxx)
				+ X.ptr[2] * -heatRatioMinusOne * vL(1)
				+ X.ptr[3] * -heatRatioMinusOne * vL(2)
				+ X.ptr[4] * heatRatioMinusOne
			) * invDenom
		);
		Y.ptr[1] += nSq(0) * (
			(
				X.ptr[0] * (denom - heatRatioMinusOne * vSq)
				+ X.ptr[1] * 2. * heatRatioMinusOne * vL(0)
				+ X.ptr[2] * 2. * heatRatioMinusOne * vL(1)
				+ X.ptr[3] * 2. * heatRatioMinusOne * vL(2)
				+ X.ptr[4] * -2. * heatRatioMinusOne
			) * invDenom
		);
		Y.ptr[2] += nSq(0) * (
			X.ptr[0] * (v(0) * gU_xy / gU_xx - v(1))
			+ X.ptr[1] * -gU_xy / gU_xx
			+ X.ptr[2]
		);
		Y.ptr[3] += nSq(0) * (
			X.ptr[0] * (v(0) * gU_xz / gU_xx - v(2))
			+ X.ptr[1] * -gU_xz / gU_xx
			+ X.ptr[3]
		);
		Y.ptr[4] += nSq(0) * (
			(
				X.ptr[0] * (.5 * heatRatioMinusOne * vSq - Cs * v(0) / sqrt_gUxx)
				+ X.ptr[1] * (-heatRatioMinusOne * vL(0) + Cs / sqrt_gUxx)
				+ X.ptr[2] * -heatRatioMinusOne * vL(1)
				+ X.ptr[3] * -heatRatioMinusOne * vL(2)
				+ X.ptr[4] * heatRatioMinusOne
			) * invDenom
		);

		//y dir
		Y.ptr[0] += nSq(1) * (
			(
				X.ptr[0] * (.5 * heatRatioMinusOne * vSq + Cs * v(1) / sqrt_gUyy)
				+ X.ptr[1] * -heatRatioMinusOne * vL(0)
				+ X.ptr[2] * (-heatRatioMinusOne * vL(1) - Cs / sqrt_gUyy)
				+ X.ptr[3] * -heatRatioMinusOne * vL(2)
				+ X.ptr[4] * heatRatioMinusOne
			) * invDenom
		);
		Y.ptr[1] += nSq(1) * (
			X.ptr[0] * (v(1) * gU_xy / gU_yy - v(0))
			+ X.ptr[1]
			+ X.ptr[2] * -gU_xy / gU_yy
		);
		Y.ptr[2] += nSq(1) * (	
			(
				X.ptr[0] * (denom - heatRatioMinusOne * vSq)
				+ X.ptr[1] * 2. * heatRatioMinusOne * vL(0)
				+ X.ptr[2] * 2. * heatRatioMinusOne * vL(1)
				+ X.ptr[3] * 2. * heatRatioMinusOne * vL(2)
				+ X.ptr[4] * -2. * heatRatioMinusOne
			) * invDenom
		);
		Y.ptr[3] += nSq(1) * (
			X.ptr[0] * (v(1) * gU_yz / gU_yy - v(2))
			+ X.ptr[2] * -gU_yz / gU_yy
			+ X.ptr[3]
		);
		Y.ptr[4] += nSq(1) * (
			(
				X.ptr[0] * (.5 * heatRatioMinusOne * vSq - Cs * v(1) / sqrt_gUyy)
				+ X.ptr[1] * -heatRatioMinusOne * vL(0)
				+ X.ptr[2] * (-heatRatioMinusOne * vL(1) + Cs / sqrt_gUyy)
				+ X.ptr[3] * -heatRatioMinusOne * vL(2)
				+ X.ptr[4] * heatRatioMinusOne
			) * invDenom
		);

		//z dir
		Y.ptr[0] += nSq(2) * (
			(
				X.ptr[0] * (.5 * heatRatioMinusOne * vSq + Cs * v(2) / sqrt_gUzz)
				+ X.ptr[1] * -heatRatioMinusOne * vL(0)
				+ X.ptr[2] * -heatRatioMinusOne * vL(1)
				+ X.ptr[3] * (-heatRatioMinusOne * vL(2) - Cs / sqrt_gUzz)
				+ X.ptr[4] * heatRatioMinusOne
			) * invDenom
		);
		Y.ptr[1] += nSq(2) * (
			X.ptr[0] * (v(2) * gU_xz / gU_zz - v(0))
			+ X.ptr[1]
			+ X.ptr[3] * -gU_xz / gU_zz
		);
		Y.ptr[2] += nSq(2) * (
			X.ptr[0] * (v(2) * gU_yz / gU_zz - v(1))
			+ X.ptr[2]
			+ X.ptr[3] * -gU_yz / gU_zz
		);
		Y.ptr[3] += nSq(2) * (	
			(
				X.ptr[0] * (denom - heatRatioMinusOne * vSq)
				+ X.ptr[1] * 2. * heatRatioMinusOne * vL(0)
				+ X.ptr[2] * 2. * heatRatioMinusOne * vL(1)
				+ X.ptr[3] * 2. * heatRatioMinusOne * vL(2)
				+ X.ptr[4] * -2. * heatRatioMinusOne
			) * invDenom
		);
		Y.ptr[4] += nSq(2) * (
			(
				X.ptr[0] * (.5 * heatRatioMinusOne * vSq - Cs * v(2) / sqrt_gUzz)
				+ X.ptr[1] * -heatRatioMinusOne * vL(0)
				+ X.ptr[2] * -heatRatioMinusOne * vL(1)
				+ X.ptr[3] * (-heatRatioMinusOne * vL(2) + Cs / sqrt_gUzz)
				+ X.ptr[4] * heatRatioMinusOne
			) * invDenom
		);

		return Y;
	}

	Cons applyEigR(Cons X, const Eigen& vars, real3 n) {
		const real& Cs = vars.Cs;
		const real3& v = vars.v;	//v^i ... upper
		const auto& vL = v;			//v_i ... lower (not left)
		const real& vSq = real3::dot(v, vL);
		const real& hTotal = vars.hTotal;

		real3 nSq = n * n;	//per-component multiply
		
//		const real gU_xx = 1;
//		const real gU_yy = 1;
//		const real gU_zz = 1;
		const real gU_xy = 0;
		const real gU_xz = 0;
		const real gU_yz = 0;
		const real sqrt_gUxx = 1.;//sqrt_gUxx;
		const real sqrt_gUyy = 1.;//sqrt_gUjj;
		const real sqrt_gUzz = 1.;//sqrt_gUjj;

		Cons Y;
	
		//x dir
		Y.ptr[0] += nSq(0) * (
			X.ptr[0] + X.ptr[1] + X.ptr[4]
		);
		Y.ptr[1] += nSq(0) * (
			X.ptr[0] * (v(0) - Cs * sqrt_gUxx)
			+ X.ptr[1] * v(0)
			+ X.ptr[4] * (v(0) + Cs * sqrt_gUxx)
		);
		Y.ptr[2] += nSq(0) * (
			X.ptr[0] * (v(1) - Cs * gU_xy / sqrt_gUxx)
			+ X.ptr[1] * v(1)
			+ X.ptr[2]
			+ X.ptr[4] * (v(1) + Cs * gU_xy / sqrt_gUxx)
		);
		Y.ptr[3] += nSq(0) * (
			X.ptr[0] * (v(2) - Cs * gU_xz / sqrt_gUxx)
			+ X.ptr[1] * v(2)
			+ X.ptr[3]
			+ X.ptr[4] * (v(2) + Cs * gU_xz / sqrt_gUxx)
		);
		Y.ptr[4] += nSq(0) * (
			X.ptr[0] * (hTotal - Cs * v(0) / sqrt_gUxx)
			+ X.ptr[1] * vSq / 2.
			+ X.ptr[2] * vL(1)
			+ X.ptr[3] * vL(2)
			+ X.ptr[4] * (hTotal + Cs * v(0) / sqrt_gUxx)
		);

		//y dir
		Y.ptr[0] += nSq(1) * (
			X.ptr[0] + X.ptr[2] + X.ptr[4]
		);
		Y.ptr[1] += nSq(1) * (
			X.ptr[0] * (v(0) - Cs * gU_xy / sqrt_gUyy)
			+ X.ptr[1]
			+ X.ptr[2] * v(0)
			+ X.ptr[4] * (v(0) + Cs * gU_xy / sqrt_gUyy)
		);
		Y.ptr[2] += nSq(1) * (
			X.ptr[0] * (v(1) - Cs * sqrt_gUyy)
			+ X.ptr[2] * v(1)
			+ X.ptr[4] * (v(1) + Cs * sqrt_gUyy)
		);
		Y.ptr[3] += nSq(1) * (
			X.ptr[0] * (v(2) - Cs * gU_yz / sqrt_gUyy)
			+ X.ptr[2] * v(2)
			+ X.ptr[3]
			+ X.ptr[4] * (v(2) + Cs * gU_yz / sqrt_gUyy)
		);
		Y.ptr[4] += nSq(1) * (
			X.ptr[0] * (hTotal - Cs * v(1) / sqrt_gUyy)
			+ X.ptr[1] * vL(0)
			+ X.ptr[2] * vSq / 2.
			+ X.ptr[3] * vL(2)
			+ X.ptr[4] * (hTotal + Cs * v(1) / sqrt_gUyy)
		);

		//z dir
		Y.ptr[0] += nSq(2) * (
			X.ptr[0] + X.ptr[3] + X.ptr[4]
		);
		Y.ptr[1] += nSq(2) * (
			X.ptr[0] * (v(0) - Cs * gU_xz / sqrt_gUzz)
			+ X.ptr[1]
			+ X.ptr[3] * v(0)
			+ X.ptr[4] * (v(0) + Cs * gU_xz / sqrt_gUzz)
		);
		Y.ptr[2] += nSq(2) * (
			X.ptr[0] * (v(1) - Cs * gU_yz / sqrt_gUzz)
			+ X.ptr[2]
			+ X.ptr[3] * v(1)
			+ X.ptr[4] * (v(1) + Cs * gU_yz / sqrt_gUzz)
		);
		Y.ptr[3] += nSq(2) * (
			X.ptr[0] * (v(2) - Cs * sqrt_gUzz)
			+ X.ptr[3] * v(2)
			+ X.ptr[4] * (v(2) + Cs * sqrt_gUzz)
		);
		Y.ptr[4] += nSq(2) * (
			X.ptr[0] * (hTotal - Cs * v(2) / sqrt_gUzz)
			+ X.ptr[1] * vL(0)
			+ X.ptr[2] * vL(1)
			+ X.ptr[3] * vSq / 2.
			+ X.ptr[4] * (hTotal + Cs * v(2) / sqrt_gUzz)
		);

		return Y;
	}

#if 0
	//notice applyFlux(vars, U) == applyEigL(vars, lambdas * applyEigR(vars, U)) when vars are derived from U alone (instead of from a pair across an interface)
	Cons applyFlux(WaveVec x, const Eigen& vars, real3 n) {
		const real3& v = vars.v;	//v^i ... upper
		const auto& vL = v;			//v_i ... lower (not left)
		
		real heatRatioMinusOne = heatCapacityRatio - 1;

		Cons Y;
		
		//x dir
		//TODO replace n with (1,0,0)
		Y.ptr[0] =
			X.ptr[1] * n(0) 
			+ X.ptr[2] * n(1) 
			+ X.ptr[3] * n(2);
		Y.ptr[1] =
			X.ptr[0] * (-v_n * v(0) + heatRatioMinusOne * .5 * vSq * gUj.x)
			+ X.ptr[1] * (v(0) * n(0) - heatRatioMinusOne * gUj.x * vL(0) + v_n)
			+ X.ptr[2] * (v(0) * n(1) - heatRatioMinusOne * gUj.x * vL(1))
			+ X.ptr[3] * (v(0) * n(2) - heatRatioMinusOne * gUj.x * vL(2))
			+ X.ptr[4] * heatRatioMinusOne * n(0);
		Y.ptr[2] = 
			X.ptr[0] * (-v_n * v(1) + heatRatioMinusOne * .5 * vSq * gUj.y)
			+ X.ptr[1] * (v(1) * n(0) - heatRatioMinusOne * gUj.y * vL(0))
			+ X.ptr[2] * (v(1) * n(1) - heatRatioMinusOne * gUj.y * vL(1) + v_n)
			+ X.ptr[3] * (v(1) * n(2) - heatRatioMinusOne * gUj.y * vL(2))
			+ X.ptr[4] * heatRatioMinusOne * n(1);
		Y.ptr[3] = 
			X.ptr[0] * (-v_n * v(2) + heatRatioMinusOne * .5 * vSq * gUj.z)
			+ X.ptr[1] * (v(2) * n(0) - heatRatioMinusOne * gUj.z * vL(0))
			+ X.ptr[2] * (v(2) * n(1) - heatRatioMinusOne * gUj.z * vL(1))
			+ X.ptr[3] * (v(2) * n(2) - heatRatioMinusOne * gUj.z * vL(2) + v_n)
			+ X.ptr[4] * heatRatioMinusOne * n(2);
		Y.ptr[4] = 
			X.ptr[0] * v_n * (heatRatioMinusOne * .5 * vSq - hTotal)
			+ X.ptr[1] * (-heatRatioMinusOne * v_n * vL(0) + n(0) * hTotal)
			+ X.ptr[2] * (-heatRatioMinusOne * v_n * vL(1) + n(1) * hTotal)
			+ X.ptr[3] * (-heatRatioMinusOne * v_n * vL(2) + n(2) * hTotal)
			+ X.ptr[4] * heatCapacityRatio * v_n:

		//y dir
		//TODO replace n with (0,1,0)
		Y.ptr[0] =
			X.ptr[1] * n(0) 
			+ X.ptr[2] * n(1) 
			+ X.ptr[3] * n(2);
		Y.ptr[1] =
			X.ptr[0] * (-v_n * v(0) + heatRatioMinusOne * .5 * vSq * gUj.x)
			+ X.ptr[1] * (v(0) * n(0) - heatRatioMinusOne * gUj.x * vL(0) + v_n)
			+ X.ptr[2] * (v(0) * n(1) - heatRatioMinusOne * gUj.x * vL(1))
			+ X.ptr[3] * (v(0) * n(2) - heatRatioMinusOne * gUj.x * vL(2))
			+ X.ptr[4] * heatRatioMinusOne * n(0);
		Y.ptr[2] = 
			X.ptr[0] * (-v_n * v(1) + heatRatioMinusOne * .5 * vSq * gUj.y)
			+ X.ptr[1] * (v(1) * n(0) - heatRatioMinusOne * gUj.y * vL(0))
			+ X.ptr[2] * (v(1) * n(1) - heatRatioMinusOne * gUj.y * vL(1) + v_n)
			+ X.ptr[3] * (v(1) * n(2) - heatRatioMinusOne * gUj.y * vL(2))
			+ X.ptr[4] * heatRatioMinusOne * n(1);
		Y.ptr[3] = 
			X.ptr[0] * (-v_n * v(2) + heatRatioMinusOne * .5 * vSq * gUj.z)
			+ X.ptr[1] * (v(2) * n(0) - heatRatioMinusOne * gUj.z * vL(0))
			+ X.ptr[2] * (v(2) * n(1) - heatRatioMinusOne * gUj.z * vL(1))
			+ X.ptr[3] * (v(2) * n(2) - heatRatioMinusOne * gUj.z * vL(2) + v_n)
			+ X.ptr[4] * heatRatioMinusOne * n(2);
		Y.ptr[4] = 
			X.ptr[0] * v_n * (heatRatioMinusOne * .5 * vSq - hTotal)
			+ X.ptr[1] * (-heatRatioMinusOne * v_n * vL(0) + n(0) * hTotal)
			+ X.ptr[2] * (-heatRatioMinusOne * v_n * vL(1) + n(1) * hTotal)
			+ X.ptr[3] * (-heatRatioMinusOne * v_n * vL(2) + n(2) * hTotal)
			+ X.ptr[4] * heatCapacityRatio * v_n:

		//z dir
		//TODO replace n with (0,0,1)
		Y.ptr[0] =
			X.ptr[1] * n(0) 
			+ X.ptr[2] * n(1) 
			+ X.ptr[3] * n(2);
		Y.ptr[1] =
			X.ptr[0] * (-v_n * v(0) + heatRatioMinusOne * .5 * vSq * gUj.x)
			+ X.ptr[1] * (v(0) * n(0) - heatRatioMinusOne * gUj.x * vL(0) + v_n)
			+ X.ptr[2] * (v(0) * n(1) - heatRatioMinusOne * gUj.x * vL(1))
			+ X.ptr[3] * (v(0) * n(2) - heatRatioMinusOne * gUj.x * vL(2))
			+ X.ptr[4] * heatRatioMinusOne * n(0);
		Y.ptr[2] = 
			X.ptr[0] * (-v_n * v(1) + heatRatioMinusOne * .5 * vSq * gUj.y)
			+ X.ptr[1] * (v(1) * n(0) - heatRatioMinusOne * gUj.y * vL(0))
			+ X.ptr[2] * (v(1) * n(1) - heatRatioMinusOne * gUj.y * vL(1) + v_n)
			+ X.ptr[3] * (v(1) * n(2) - heatRatioMinusOne * gUj.y * vL(2))
			+ X.ptr[4] * heatRatioMinusOne * n(1);
		Y.ptr[3] = 
			X.ptr[0] * (-v_n * v(2) + heatRatioMinusOne * .5 * vSq * gUj.z)
			+ X.ptr[1] * (v(2) * n(0) - heatRatioMinusOne * gUj.z * vL(0))
			+ X.ptr[2] * (v(2) * n(1) - heatRatioMinusOne * gUj.z * vL(1))
			+ X.ptr[3] * (v(2) * n(2) - heatRatioMinusOne * gUj.z * vL(2) + v_n)
			+ X.ptr[4] * heatRatioMinusOne * n(2);
		Y.ptr[4] = 
			X.ptr[0] * v_n * (heatRatioMinusOne * .5 * vSq - hTotal)
			+ X.ptr[1] * (-heatRatioMinusOne * v_n * vL(0) + n(0) * hTotal)
			+ X.ptr[2] * (-heatRatioMinusOne * v_n * vL(1) + n(1) * hTotal)
			+ X.ptr[3] * (-heatRatioMinusOne * v_n * vL(2) + n(2) * hTotal)
			+ X.ptr[4] * heatCapacityRatio * v_n:

		return Y;
	}
#endif

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
