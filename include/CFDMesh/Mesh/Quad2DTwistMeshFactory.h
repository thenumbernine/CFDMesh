#pragma once

#include "CFDMesh/Mesh/Quad2DMeshFactory.h"

namespace CFDMesh {
namespace Mesh {

template<typename real, int dim, typename Cons>
struct Quad2DTwistMeshFactory : public Quad2DMeshFactory<real, dim, Cons> {
	using Super = Quad2DMeshFactory<real, dim, Cons>;
	using real2 = typename Super::real2;

	Quad2DTwistMeshFactory() : Super("unit square of quads, twist in the middle") {}
	virtual real2 coordChart(real2 v) const {
		real r = v.length();
		//real theta = std::max(0., 1. - r);
		real sigma = 3.;	//almost 0 at r=1
		real const rotationAmplitude = 3.;
		real theta = rotationAmplitude*sigma*r*exp(-sigma*sigma*r*r);
		real costh = cos(theta), sinth = sin(theta);
		return real2(
			costh * v(0) - sinth * v(1),
			sinth * v(0) + costh * v(1));	
	}
};

}
}
