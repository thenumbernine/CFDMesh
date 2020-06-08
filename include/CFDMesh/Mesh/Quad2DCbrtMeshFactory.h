#pragma once

#include "CFDMesh/Mesh/Quad2DMeshFactory.h"

namespace CFDMesh {
namespace Mesh {

template<typename real, int dim, typename Cons>
struct Quad2DCbrtMeshFactory : public Quad2DMeshFactory<real, dim, Cons> {
	using Super = Quad2DMeshFactory<real, dim, Cons>;
	using real2 = typename Super::real2;
	
	Quad2DCbrtMeshFactory() : Super("unit square of quads, cbrt mapping") {}
	virtual real2 coordChart(real2 v) const {
		return real2(cbrt(v(0)), cbrt(v(1)));
	}
};

}
}
