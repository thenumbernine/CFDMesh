#pragma once

#include "CFDMesh/Mesh/Quad2DMeshFactory.h"
#include "CFDMesh/Vector.h"

namespace CFDMesh {
namespace Mesh {

template<typename real, int dim, typename Cons>
struct PolarMeshFactory : public Quad2DMeshFactory<real, dim, Cons> {
	using Super = Quad2DMeshFactory<real, dim, Cons>;
	using real2 = typename Super::real2;

	PolarMeshFactory() : Super("polar") {
		Super::size = int2(20, 50);
		Super::mins = real2(.1, 0);
		Super::maxs = real2(1, 1);
		Super::repeat = bool2(false, true);
		Super::capmin = bool2(false, false);
	}
	virtual real2 coordChart(real2 v) const {
		real theta = 2. * M_PI * v(1);
		return real2(cos(theta), sin(theta)) * v(0);
	}
};

}
}
