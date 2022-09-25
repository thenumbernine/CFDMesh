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
