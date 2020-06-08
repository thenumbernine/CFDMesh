#pragma once

#include "CFDMesh/Mesh/Vertex.h"
#include "CFDMesh/Mesh/Face.h"
#include "CFDMesh/Mesh/Cell.h"
#include <memory>

namespace CFDMesh {

template<typename, int, typename>
struct MeshNamespace;

template<typename real, int dim, typename Cons>
struct MeshFactory {
	using real2 = Tensor::Vector<real, 2>;
	using real3 = Tensor::Vector<real, 3>;
		
	using Vertex = Vertex_<real>;
	using Face = Face_<real, Cons>;
	using Cell = Cell_<real, Cons>;

	using Mesh = typename MeshNamespace<real, dim, Cons>::Mesh;


	const char* name = nullptr;
	
	MeshFactory(const char* name_) : name(name_) {}
	virtual ~MeshFactory() {}
	virtual void updateGUI() {}
	virtual std::shared_ptr<Mesh> createMesh() = 0;
protected:
	virtual std::shared_ptr<Mesh> createMeshSuper() const {
		return Mesh::create();
	}
};

}
