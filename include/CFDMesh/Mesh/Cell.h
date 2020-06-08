#pragma once

#include "CFDMesh/Vector.h"

namespace CFDMesh {
namespace Mesh {

//n-forms
template<typename real, typename Cons>
struct Cell_ {
	using This = Cell_;
	using real3 = Tensor::Vector<real, 3>;

	real3 pos;
	real volume = 0;
	
	//required by the finite volume algorithm
	//std::vector<int> faces;
	int faceOffset = 0;
	int faceCount = 0;

	//not required by finite volume algorithm
	//however the cell volume is required, and is calculated using vtxs
	//also the renderer requires the vertexes
	//std::vector<int> vtxs;
	int vtxOffset = 0;
	int vtxCount = 0;

	float displayValue = 0;
	Cons U;
	
	static constexpr auto fields = std::make_tuple(
		std::make_pair("pos", &This::pos),
		std::make_pair("volume", &This::volume),
		std::make_pair("faceOffset", &This::faceOffset),
		std::make_pair("faceCount", &This::faceCount),
		std::make_pair("vtxOffset", &This::vtxOffset),
		std::make_pair("vtxCount", &This::vtxCount),
		std::make_pair("displayValue", &This::displayValue),
		std::make_pair("U", &This::U)
	);
};

}
}
