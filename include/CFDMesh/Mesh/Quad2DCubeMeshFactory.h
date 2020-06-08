#pragma once

#include "CFDMesh/Mesh/Quad2DMeshFactory.h"

namespace CFDMesh {
namespace Mesh {

template<typename real, int dim, typename Cons>
struct Quad2DCubeMeshFactory : public Quad2DMeshFactory<real, dim, Cons> {
	using Super = Quad2DMeshFactory<real, dim, Cons>;
	using real2 = typename Super::real2;

	template<typename T> 
	static T cubed(const T& t) { return t * t * t; }
	
	Quad2DCubeMeshFactory() : Super("unit square of quads, cubed mapping") {}
	virtual real2 coordChart(real2 v) const {
		return real2(cubed(v(0)), cubed(v(1)));
	}
};

}
}
