#pragma once

#include "CFDMesh/Mesh/Vertex.h"
#include "CFDMesh/Mesh/Face.h"
#include "CFDMesh/Mesh/Cell.h"
#include "CFDMesh/Vector.h"
#include <memory>

namespace CFDMesh {
namespace Mesh {

template<typename, int, typename>
struct Mesh;

template<typename real, int dim, typename Cons>
struct MeshFactory {
	using Mesh = typename ::CFDMesh::Mesh::Mesh<real, dim, Cons>;
	using real2 = typename Mesh::real2;
	using real3 = typename Mesh::real3;
	using Vertex = typename Mesh::Vertex;
	using Face = typename Mesh::Face;
	using Cell = typename Mesh::Cell;

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
}
