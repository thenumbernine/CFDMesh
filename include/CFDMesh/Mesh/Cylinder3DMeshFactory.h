#pragma once

#include "CFDMesh/Mesh/Cube3DMeshFactory.h"
#include "CFDMesh/Vector.h"

namespace CFDMesh {
namespace Mesh {

template<typename real, int dim, typename Cons>
struct Cylinder3DMeshFactory : public Cube3DMeshFactory<real, dim, Cons> {
	using Super = Cube3DMeshFactory<real, dim, Cons>;
	using real3 = typename Super::real3;

	Cylinder3DMeshFactory() : Super("cylinder") {
		Super::size = int3(10, 10, 10);
		Super::mins = real3(0, 0, 0);
		Super::maxs = real3(1, 1, 1);
		Super::wrap = bool3(false, true, false);
		Super::capmin = bool3(true, false, false);
	}
	
	virtual real3 coordChart(real3 x) const {
		real theta = x(1) * 2 * M_PI;
		return real3(
			x(0) * cos(theta),
			x(0) * sin(theta),
			x(2));
	}
};

}
}
