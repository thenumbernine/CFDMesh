#pragma once

#include "CFDMesh/Vector.h"

namespace CFDMesh {

//zero-forms
template<typename real>
struct Vertex_ {	//not required by finite volume algorithm
	using This = Vertex_;
	using real3 = Tensor::Vector<real, 3>;
	
	real3 pos;
	
	//keeping track of Vertex_::faces isn't used by the renderer, mesh generation, or finite-volume integration
	//std::vector<int> faces;
	
	static constexpr auto fields = std::make_tuple(
		std::make_pair("pos", &This::pos)
	);
	
	Vertex_() {}
	Vertex_(real3 pos_) : pos(pos_) {}
};

}
