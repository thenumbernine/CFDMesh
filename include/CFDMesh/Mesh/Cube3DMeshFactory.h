#pragma once

#include "CFDMesh/Mesh/MeshFactory.h"
#include "CFDMesh/Vector.h"
#include "CFDMesh/GUI.h"

namespace CFDMesh {
namespace Mesh {

template<typename real, int dim, typename Cons>
struct Cube3DMeshFactory : public MeshFactory<real, dim, Cons> {
	using This = Cube3DMeshFactory;
	using Super = MeshFactory<real, dim, Cons>;
	using Mesh = typename Super::Mesh;
	using real3 = typename Super::real3;

	int3 size = int3(8,8,8);
	float3 mins = float3(-1, -1, -1);
	float3 maxs = float3(1, 1, 1);
	bool3 wrap;
	bool3 capmin;

	static constexpr auto fields = std::make_tuple(
		std::make_pair("size", &This::size),
		std::make_pair("mins", &This::mins),
		std::make_pair("maxs", &This::maxs),
		std::make_pair("wrap", &This::wrap),
		std::make_pair("capmin", &This::capmin)
	);

	Cube3DMeshFactory(const char* name_ = "cube mesh") : Super(name_) {}

	virtual void updateGUI() {
		CFDMesh::updateGUI(this);
	}

	virtual real3 coordChart(real3 x) const { return x; }

	virtual std::shared_ptr<Mesh> createMesh() {
		std::shared_ptr<Mesh> mesh = MeshFactory<real, dim, Cons>::createMeshSuper();

		int3 n = size + 1;
		int3 step(1, n(0), n(0) * n(1));
		
		int vtxsize = n.volume();
		mesh->vtxs.resize(vtxsize);
		
		int3 vtxmax = size;
		for (int j = 0; j < 3; ++j) {
			if (wrap(j)) ++vtxmax(j);
		}
		
		int3 i;
		for (i(2) = 0; i(2) < n(2); ++i(2)) {
			for (i(1) = 0; i(1) < n(1); ++i(1)) {
				for (i(0) = 0; i(0) < n(0); ++i(0)) {
					real3 x = (real3)i / (real3)vtxmax * (maxs - mins) + mins;
					mesh->vtxs[int3::dot(i, step)].pos = this->coordChart(x);
				}
			}
		}
		
		int3 imax;
		for (int j = 0; j < 3; ++j) {
			imax(j) = wrap(j) ? n(j) : n(j)-1;
		}
			
		int3 in;
		for (i(2) = 0; i(2) < imax(2); ++i(2)) {
			in(2) = (i(2) + 1) % n(2);
			for (i(1) = 0; i(1) < imax(1); ++i(1)) {
				in(1) = (i(1) + 1) % n(1);
				for (i(0) = 0; i(0) < imax(0); ++i(0)) {
					in(0) = (i(0) + 1) % n(0);
					mesh->addCell(std::vector<int>{
						//using z-order
						i(0) + n(0) * (i(1) + n(1) * i(2)),
						in(0) + n(0) * (i(1) + n(1) * i(2)),
						i(0) + n(0) * (in(1) + n(1) * i(2)),
						in(0) + n(0) * (in(1) + n(1) * i(2)),
						
						i(0) + n(0) * (i(1) + n(1) * in(2)),
						in(0) + n(0) * (i(1) + n(1) * in(2)),
						i(0) + n(0) * (in(1) + n(1) * in(2)),
						in(0) + n(0) * (in(1) + n(1) * in(2)),
					});
				}
			}
		}

		mesh->calcAux();
		return mesh;
	}
};

}
}
