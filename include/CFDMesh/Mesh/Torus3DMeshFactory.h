#pragma once

#include "CFDMesh/Mesh/Cube3DMeshFactory.h"
#include "ImGuiCommon/Reflect.h"

namespace CFDMesh {
namespace Mesh {

template<typename real, int dim, typename Cons>
struct Torus3DMeshFactory : public Cube3DMeshFactory<real, dim, Cons> {
	using This = Torus3DMeshFactory;
	using Super = Cube3DMeshFactory<real, dim, Cons>;
	using real3 = typename Super::real3;

	real R = 2.;

	static constexpr auto fields = std::tuple_cat(
		Super::fields,
		std::make_tuple(
			std::make_pair("R", &This::R)
		)
	);

	Torus3DMeshFactory() : Super("torus") {
		Super::size = int3(1, 4, 4);
		Super::mins = real3(0, 0, 0);
		Super::maxs = real3(1, 1, 1);
		Super::wrap = bool3(false, true, true);
	}

	virtual real3 coordChart(real3 x) const {
		real r = x(0);
		real theta = x(1) * 2 * M_PI;
		real phi = x(2) * 2 * M_PI;
		return real3(
			(r * cos(theta) + R) * cos(phi), 
			(r * cos(theta) + R) * sin(phi), 
			-r * sin(theta)
		);
	}

	virtual void updateGUI() {
		ImGuiCommon::updateGUI(this);
	}
};

}
}
