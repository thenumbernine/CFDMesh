#pragma once

#include "CFDMesh/Mesh/MeshFactory.h"
#include "CFDMesh/Vector.h"
#include "CFDMesh/GUI.h"

namespace CFDMesh {
namespace Mesh {

template<typename real, int dim, typename Cons>
struct Chart2DMeshFactory : public MeshFactory<real, dim, Cons> {
	using This = Chart2DMeshFactory;
	using Super = MeshFactory<real, dim, Cons>;
	using real2 = typename Super::real2;
	
	//int2 size = int2(100, 100);
	int2 size = int2(20, 20);
	float2 mins = real2(-1, -1);
	float2 maxs = real2(1, 1);
	bool2 repeat = bool2(false, false);
	bool2 capmin = bool2(false, false);
	
	static constexpr auto fields = std::make_tuple(
		std::make_pair("size", &This::size),
		std::make_pair("mins", &This::mins),
		std::make_pair("maxs", &This::maxs),
		std::make_pair("repeat", &This::repeat),
		std::make_pair("capmin", &This::capmin)
	);

	Chart2DMeshFactory(const char* name_) : Super(name_) {}

	virtual real2 coordChart(real2 x) const { return x; }
	
	virtual void updateGUI() {
		CFDMesh::updateGUI(this);
	}
};

}
}
