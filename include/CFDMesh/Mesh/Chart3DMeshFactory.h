#pragma once

#include "CFDMesh/Mesh/MeshFactory.h"
#include "CFDMesh/Vector.h"
#include "CFDMesh/GUI.h"

namespace CFDMesh {
namespace Mesh {

template<typename real, int dim, typename Cons>
struct Chart3DMeshFactory : public MeshFactory<real, dim, Cons> {
	using This = Chart3DMeshFactory;
	using Super = MeshFactory<real, dim, Cons>;
	using real3 = typename Super::real3;

	int3 size = int3(10,10,10);
	float3 mins = float3(-1, -1, -1);
	float3 maxs = float3(1, 1, 1);
	bool3 repeat;
	bool3 capmin;

	//TODO support for inheritence and reflection
	static constexpr auto fields = std::make_tuple(
		std::make_pair("size", &This::size),
		std::make_pair("mins", &This::mins),
		std::make_pair("maxs", &This::maxs),
		std::make_pair("repeat", &This::repeat),
		std::make_pair("capmin", &This::capmin)
	);

	Chart3DMeshFactory(const char* name_ = "3D chart mesh") : ::CFDMesh::Mesh::MeshFactory<real, dim, Cons>(name_) {}

	virtual real3 coordChart(real3 x) const { return x; }

	virtual void updateGUI() {
		CFDMesh::updateGUI(this);
	}
};

}
}
