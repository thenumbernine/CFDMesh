#pragma once

#include "CFDMesh/Vector.h"

namespace CFDMesh {
namespace Mesh {

//zero-forms
template<typename real>
struct Vertex {	//not required by finite volume algorithm
	using This = Vertex;
	using real3 = Tensor::Vector<real, 3>;
	
	real3 pos;
	
	//keeping track of Vertex::faces isn't used by the renderer, mesh generation, or finite-volume integration
	//std::vector<int> faces;
	
	static constexpr auto fields = std::make_tuple(
		std::make_pair("pos", &This::pos)
	);
	
	Vertex() {}
	Vertex(real3 pos_) : pos(pos_) {}
};

}
}
