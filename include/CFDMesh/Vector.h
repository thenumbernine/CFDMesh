#pragma once

#include "Tensor/Vector.h"

using uchar2 = Tensor::Vector<unsigned char, 2>;
using uchar3 = Tensor::Vector<unsigned char, 3>;
using uchar4 = Tensor::Vector<unsigned char, 4>;

using int2 = Tensor::Vector<int, 2>;
using int3 = Tensor::Vector<int, 3>;
using int4 = Tensor::Vector<int, 4>;

using float2 = Tensor::Vector<float, 2>;
using float3 = Tensor::Vector<float, 3>;
using float4 = Tensor::Vector<float, 4>;

using double2 = Tensor::Vector<double, 2>;
using double3 = Tensor::Vector<double, 3>;
using double4 = Tensor::Vector<double, 4>;


//for giving operators to the Cons and Prim vector classes
//how can you add correctly-typed ops via crtp to a union?  unions can't inherit.
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

#define ADD_OPS(classname)	\
	real& operator()(int i) { return ptr[i]; }\
	const real& operator()(int i) const { return ptr[i]; }\
\
	ADD_VECTOR_OP(classname, +)\
	ADD_VECTOR_OP(classname, -)\
	ADD_VECTOR_OP(classname, *)\
	ADD_SCALAR_OP(classname, *)\
	ADD_VECTOR_OP_EQ(classname, +=)\
	ADD_VECTOR_OP_EQ(classname, -=)

#if 0	//hmm, why isn't this working
	classname& operator=(const classname& o) {\
		for (int i = 0; i < size; ++i) {\
			ptr[i] = o.ptr[i];\
		}\
		return *this;\
	}
#endif
