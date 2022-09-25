#pragma once

#include "CFDMesh/Mesh/Quad2DMeshFactory.h"
#include "ImGuiCommon/Reflect.h"

namespace CFDMesh {
namespace Mesh {

template<typename real, int dim, typename Cons>
struct Quad2DRotateMeshFactory : public Quad2DMeshFactory<real, dim, Cons> {
	using This = Quad2DRotateMeshFactory;
	using Super = Quad2DMeshFactory<real, dim, Cons>;
	using real2 = typename Super::real2;

	real thetaOver2Pi = {};
	
	static constexpr auto fields = std::tuple_cat(
		Super::fields,
		std::make_tuple(
			std::make_pair("theta / 2pi", &This::thetaOver2Pi)
		)
	);

	Quad2DRotateMeshFactory() : Super("quad grid with fixed rotation") {}
	
	virtual real2 coordChart(real2 v) const {
		real theta = 2 * M_PI * thetaOver2Pi;
		real costh = cos(theta), sinth = sin(theta);
		return real2(
			costh * v(0) - sinth * v(1),
			sinth * v(0) + costh * v(1));
	}

	virtual void updateGUI() {
		ImGuiCommon::updateGUI(this);
	}
};

}
}
