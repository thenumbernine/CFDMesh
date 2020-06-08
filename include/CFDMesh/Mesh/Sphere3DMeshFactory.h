#pragma once

#include "CFDMesh/Mesh/Cube3DMeshFactory.h"
#include "CFDMesh/Vector.h"

namespace CFDMesh {
namespace Mesh {

template<typename real, int dim, typename Cons>
struct Sphere3DMeshFactory : public Cube3DMeshFactory<real, dim, Cons> {
	using Super = Cube3DMeshFactory<real, dim, Cons>;
	using real3 = typename Super::real3;

	Sphere3DMeshFactory() : Super("sphere") {
		Super::size = int3(10, 10, 10);
		Super::mins = real3(.5, .5, 0);
		Super::maxs = real3(1, 1., 1);
		Super::repeat = bool3(false, false, true);
		//Super::capmin = bool3(true, false, false);	//TODO
	}
	
	virtual real3 coordChart(real3 x) const {
		real r = x(0);
		real theta = x(1) * M_PI;
		real phi = x(2) * 2 * M_PI;
		real sinth = sin(theta);
		return real3(
			r * cos(phi) * sinth,
			r * sin(phi) * sinth,
			r * cos(theta));
	}
};

}
}
