#pragma once

#include "Tensor/Vector.h"
#include "Tensor/Quat.h"
#include <cassert>

using bool2 = Tensor::bool2;
using bool3 = Tensor::bool3;
using bool4 = Tensor::bool4;
using uchar2 = Tensor::uchar2;
using uchar3 = Tensor::uchar3;
using uchar4 = Tensor::uchar4;
using int2 = Tensor::int2;
using int3 = Tensor::int3;
using int4 = Tensor::int4;
using float2 = Tensor::float2;
using float3 = Tensor::float3;
using float4 = Tensor::float4;
using double2 = Tensor::double2;
using double3 = Tensor::double3;
using double4 = Tensor::double4;


template<typename real3>
inline real3 cross(real3 a, real3 b) {
	return real3(
		a(1) * b(2) - a(2) * b(1),
		a(2) * b(0) - a(0) * b(2),
		a(0) * b(1) - a(1) * b(0));
}


//for giving operators to the Cons and Prim vector classes
//how can you add correctly-typed ops via crtp to a union?
//unions can't inherit.
//until then...

#define ADD_VECTOR_OP(classname, op)\
	classname operator op(const classname& b) const {\
		classname c;\
		for (int i = 0; i < size; ++i) {\
			c.ptr[i] = ptr[i] op b.ptr[i];\
		}\
		return c;\
	}

#define ADD_SCALAR_OP(classname, op)\
	classname operator op(real b) const {\
		classname c;\
		for (int i = 0; i < size; ++i) {\
			c.ptr[i] = ptr[i] op b;\
		}\
		return c;\
	}
#define ADD_VECTOR_OP_EQ(classname, op)\
	classname& operator op(const classname& b) {\
		for (int i = 0; i < size; ++i) {\
			ptr[i] op b.ptr[i];\
		}\
		return *this;\
	}

#define ADD_CAST_OP(classname)\
	template<typename T>\
	operator classname<T>() const {\
		classname<T> res;\
		for (int i = 0; i < size; ++i) {\
			res.ptr[i] = (T)ptr[i];\
		}\
		return res;\
	}


#define ADD_OPS(classname)\
	real& operator()(int i) { return ptr[i]; }\
	const real& operator()(int i) const { return ptr[i]; }\
\
	ADD_VECTOR_OP(classname, +)\
	ADD_VECTOR_OP(classname, -)\
	ADD_VECTOR_OP(classname, *)\
	ADD_SCALAR_OP(classname, *)\
	ADD_VECTOR_OP_EQ(classname, +=)\
	ADD_VECTOR_OP_EQ(classname, -=)\
	ADD_CAST_OP(classname)

#if 0	//hmm, this isn't working when it is run
	classname& operator=(const classname& o) {\
		for (int i = 0; i < size; ++i) {\
			ptr[i] = o.ptr[i];\
		}\
		return *this;\
	}
#endif


namespace CFDMesh {

template<typename V>
struct Rotate {};

template<typename real>
struct Rotate<Tensor::Vector<real, 2>> {
	using real2 = Tensor::Vector<real, 2>;
	
	static real2 from(real2 v, real2 n) {
		return real2(
			v(0) * n(0) + v(1) * n(1),
			v(1) * n(0) - v(0) * n(1));
	}

	static real2 to(real2 v, real2 n) {
		return real2(
			v(0) * n(0) - v(1) * n(1),
			v(1) * n(0) + v(0) * n(1));
	}
};

template<typename real>
struct Rotate<Tensor::Vector<real, 3>> {
	using real3 = Tensor::Vector<real, 3>;
	using Quat = Tensor::Quat<real>;

	// rotate vx,vy,vz such that n now points along the x dir
	static real3 to(real3 v, real3 n) {
		/*
		axis is n cross x-axis
		[ 1  0  0] x [nx ny nz] = [0, -nz, ny] / (ny^2 + nz^2)	
		angle = acos(n(0))
		cos angle = n(0)
		sin angle = sqrt(1 - nx^2)
		cos (angle/2) = 

		cos^2 theta + sin^2 theta = 1
		sin^2 theta = 1 - cos^2
		cos^2 theta - sin^2 theta = cos(2 theta) <-> 
		cos(theta/2) = sqrt((1 + cos(theta))/2) <-> 
		2 sin theta cos theta = sin(2 theta)
		*/
		real cosTheta = n(0);
		if (cosTheta >= 1 - 1e-7) {
			return v;
		} else if (cosTheta <= -1 + 1e-7) {
			//rotating 180' on any axis will do, as long as it is consistent between 'from' and 'to'
			return real3(-v(0), -v(1), v(2));
		}	
		real cosHalfTheta = sqrt(std::clamp<real>(.5 * (1. + cosTheta), 0, 1));
		real sinHalfTheta = sqrt(std::clamp<real>(1. - cosHalfTheta * cosHalfTheta, 0, 1));
		real n2 = sqrt(n(1) * n(1) + n(2) * n(2));
		real ax = 0.;
		real ay = -n(2) / n2;
		real az = n(1) / n2;
		Quat q(
			ax * sinHalfTheta,
			ay * sinHalfTheta,
			az * sinHalfTheta,
			cosHalfTheta);
		Quat qInv = q.unitConj();
		Quat _v(v(0), v(1), v(2), 0);
		Quat vres = Quat::mul(Quat::mul(q, _v), qInv);
		real3 vp = real3(vres(0), vres(1), vres(2));
		return vp;
	}

	// rotate vx,vy,vz such that the x dir now points along n 
	static real3 from(real3 v, real3 n) {
		real cosTheta = n(0);
		if (cosTheta >= 1 - 1e-7) {
			return v;
		} else if (cosTheta <= -1 + 1e-7) {
			//rotating 180' on any axis will do, as long as it is consistent between 'from' and 'to'
			return real3(-v(0), -v(1), v(2));
		}
assert(std::isfinite(cosTheta));
		real cosHalfTheta = sqrt(std::clamp<real>(.5 * (1. + cosTheta), 0, 1));
assert(std::isfinite(cosHalfTheta));
		real sinHalfTheta = sqrt(std::clamp<real>(1. - cosHalfTheta * cosHalfTheta, 0, 1));
assert(std::isfinite(sinHalfTheta));
		real n2 = sqrt(n(1) * n(1) + n(2) * n(2));
assert(std::isfinite(n2));
		real ax = 0.;
assert(std::isfinite(ax));
		real ay = n(2) / n2;
assert(std::isfinite(ay));
		real az = -n(1) / n2;
assert(std::isfinite(az));
		Quat q(
			ax * sinHalfTheta,
			ay * sinHalfTheta,
			az * sinHalfTheta,
			cosHalfTheta);
		Quat qInv = q.unitConj();
		Quat _v(v(0), v(1), v(2), 0);
		Quat vres = Quat::mul(Quat::mul(q, _v), qInv);
assert(std::isfinite(vres(0)));
assert(std::isfinite(vres(1)));
assert(std::isfinite(vres(2)));
		return real3(vres(0), vres(1), vres(2));

	}
};

template<typename T>
inline T rotateTo(T v, T n) {
	return Rotate<T>::to(v, n);
}

template<typename T>
inline T rotateFrom(T v, T n) {
	return Rotate<T>::from(v, n);
}

}
