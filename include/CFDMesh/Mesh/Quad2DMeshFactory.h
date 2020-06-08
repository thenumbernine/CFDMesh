#pragma once

#include "CFDMesh/Mesh/Chart2DMeshFactory.h"

namespace CFDMesh {
namespace Mesh {

template<typename real, int dim, typename Cons>
struct Quad2DMeshFactory : public Chart2DMeshFactory<real, dim, Cons> {
	using Super = Chart2DMeshFactory<real, dim, Cons>;
	using Mesh = typename Super::Mesh;
	using real2 = typename Super::real2;
	using real3 = typename Super::real3;

	Quad2DMeshFactory(const char* name_ = "unit square of quads") : Super(name_) {}
	
	virtual std::shared_ptr<Mesh> createMesh() {
		std::shared_ptr<Mesh> mesh = MeshFactory<real, dim, Cons>::createMeshSuper();

		int2 n = Super::size + 1;
		int2 step(1, n(0));
		int vtxsize = n.volume();
		if (Super::capmin(0)) vtxsize++;
		mesh->vtxs.resize(vtxsize);
		
		int2 coordRangeMax = Super::size;
		if (Super::repeat(0) || Super::capmin(0)) ++coordRangeMax(0);
		if (Super::repeat(1) || Super::capmin(1)) ++coordRangeMax(1);

		int2 iofs;
		if (Super::capmin(0)) iofs(0) = 1;
		if (Super::capmin(1)) iofs(1) = 1;
		

		int2 i;
		for (i(1) = 0; i(1) < n(1); ++i(1)) {
			for (i(0) = 0; i(0) < n(0); ++i(0)) {
				real2 x = (real2)(i + iofs) / (real2)coordRangeMax * (Super::maxs - Super::mins) + Super::mins;
				real2 u = this->coordChart(x);
				mesh->vtxs[int2::dot(i, step)].pos = real3([&u](int i) -> real { 
					return i < real2::size ? u(i) : 0.;
				});
			}
		}
		
		int capindex = n.volume();
		if (Super::capmin(0)) {
			real3 sum;
			for (int j = 0; j < n(1); ++j) {
				sum += mesh->vtxs[0 + n(0) * j].pos;
			}
			mesh->vtxs[capindex].pos = sum / (real)n(1);
		}

		int2 imax;
		for (int j = 0; j < 2; ++j) {
			imax(j) = Super::repeat(j) ? n(j) : n(j)-1;
		}
		int2 in;
		for (i(1) = 0; i(1) < imax(1); ++i(1)) {
			in(1) = (i(1) + 1) % n(1);
			for (i(0) = 0; i(0) < imax(0); ++i(0)) {
				in(0) = (i(0) + 1) % n(0);
				mesh->addCell(std::vector<int>{
					i(0) + n(0) * i(1),
					in(0) + n(0) * i(1),
					in(0) + n(0) * in(1),
					i(0) + n(0) * in(1)
				});
			}
		}

		if (Super::capmin(0)) {
			for (int j = 0; j < imax(1); ++j) {
				int jn = (j + 1) % n(1);
				mesh->addCell(std::vector<int>{ 0 + n(0) * j, 0 + n(0) * jn, capindex });
			}
		}

		mesh->calcAux();
		return mesh;
	}
};

}
}
