#pragma once

#include "CFDMesh/Vector.h"

namespace CFDMesh {
namespace Mesh {

//(n-1)-forms
template<typename real, typename Cons>
struct Face {
	using This = Face;
	using real3 = Tensor::Vector<real, 3>;
	
	real3 pos;
	real3 normal;
	real area = 0;	//space taken up by the face
	real cellDist = 0;	//dist between adjacent cell centers
	
	int2 cells = int2(-1, -1);	//there are always only 2 n-forms on either side of a (n-1)-form
	
	//the vertexes of a (n-1)-form are not required by finite volume algorithm
	//they are required for determining the interface area which is then used by the finite volume algorithm
	//also used by the renderer
	//however, note, in n>2 dimensions, a (n-1)-form will be defined by an arbitrary number of vertices (more than just two)
	//int vtxs[2];
	int vtxOffset = 0;
	int vtxCount = 0;
	
	Cons flux;
	
	static constexpr auto fields = std::make_tuple(
		std::make_pair("pos", &This::pos),
		std::make_pair("normal", &This::normal),
		std::make_pair("area", &This::area),
		std::make_pair("cellDist", &This::cellDist),
		std::make_pair("cells", &This::cells),
		std::make_pair("vtxOffset", &This::vtxOffset),
		std::make_pair("vtxCount", &This::vtxCount),
		std::make_pair("flux", &This::flux)
	);

	bool removeCell(int cellIndex) {
		if (cells(1) == cellIndex) cells(1) = -1;
		if (cells(0) == cellIndex) {
			cells(0) = cells(1);
			cells(1) = -1;
		}
		return cells(0) == -1 && cells(1) == -1;
	}
};

}
}
