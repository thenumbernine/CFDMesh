#pragma once

#include "CFDMesh/Mesh/MeshFactory.h"
#include "CFDMesh/Vector.h"
#include "ImGuiCommon/Reflect.h"
#include <vector>
#include <memory>

namespace CFDMesh {
namespace Mesh {

template<typename real, int dim, typename Cons>
struct Quad2DMeshFactory : public MeshFactory<real, dim, Cons> {
	using This = Quad2DMeshFactory;
	using Super = MeshFactory<real, dim, Cons>;
	using Mesh = typename Super::Mesh;
	using real2 = typename Super::real2;
	using real3 = typename Super::real3;
	
	//int2 size = int2(100, 100);
	int2 size = int2(20, 20);
	float2 mins = real2(-1, -1);
	float2 maxs = real2(1, 1);
	bool2 wrap = bool2(false, false);
	bool2 capmin = bool2(false, false);
	bool triangulate = false;
	
	static constexpr auto fields = std::make_tuple(
		std::make_pair("size", &This::size),
		std::make_pair("mins", &This::mins),
		std::make_pair("maxs", &This::maxs),
		std::make_pair("wrap", &This::wrap),
		std::make_pair("capmin", &This::capmin),
		std::make_pair("triangulate", &This::triangulate)
	);

	Quad2DMeshFactory(const char* name_ = "unit square of quads") : Super(name_) {}

	virtual void updateGUI() {
		ImGuiCommon::updateGUI(this);
	}

	virtual real2 coordChart(real2 x) const { 
		//TODO enum and gui dropdown
		return x;
	}
	
	virtual void addPoly(std::shared_ptr<Mesh> mesh, std::vector<int> poly) {
		//TODO behavior?  GUI can't handle that well.  dropdown enum works best with GUI.
		if (!triangulate) {
			mesh->addCell(poly);
		} else {
			int va = poly[0];
			int vb = poly[1];
			for (int i = 2; i < (int)poly.size(); ++i) {
				int vc = poly[i];
				mesh->addCell({va, vb, vc});
				vb = vc;
			}
		}
	}

	virtual bool testMakeCell(int2 i) {
		//TODO enum between this and picking a background image?
		return true;
	}

	virtual std::shared_ptr<Mesh> createMesh() {
		std::shared_ptr<Mesh> mesh = MeshFactory<real, dim, Cons>::createMeshSuper();

		int2 n = size + 1;
		int2 step(1, n(0));
		
		int vtxsize = n.volume();
		if (capmin(0)) vtxsize++;
		
		mesh->vtxs.resize(vtxsize);
		
		int2 vtxmax = size;
		for (int j = 0; j < 2; ++j) {
			if (wrap(j) || capmin(j)) ++vtxmax(j);
		}

		int2 iofs;
		for (int j = 0; j < 2; ++j) {
			if (capmin(j)) iofs(j) = 1;
		}

		int2 i;
		for (i(1) = 0; i(1) < n(1); ++i(1)) {
			for (i(0) = 0; i(0) < n(0); ++i(0)) {
				real2 x = ((real2)(i + iofs) / (real2)vtxmax).elemMul((real2)(maxs - mins)) + (real2)mins;
				real2 u = this->coordChart(x);
				mesh->vtxs[i.dot(step)].pos = (real3)u;
			}
		}
		
		int capindex = n.volume();
		if (capmin(0)) {
			real3 sum;
			for (int j = 0; j < n(1); ++j) {
				sum += mesh->vtxs[0 + n(0) * j].pos;
			}
			mesh->vtxs[capindex].pos = sum / (real)n(1);
		}

		int2 imax;
		for (int j = 0; j < 2; ++j) {
			imax(j) = wrap(j) ? n(j) : n(j)-1;
		}
		
		int2 in;
		for (i(1) = 0; i(1) < imax(1); ++i(1)) {
			in(1) = (i(1) + 1) % n(1);
			for (i(0) = 0; i(0) < imax(0); ++i(0)) {
				if (testMakeCell(i)) {
					in(0) = (i(0) + 1) % n(0);
					addPoly(mesh, {
						i(0) + n(0) * i(1),
						in(0) + n(0) * i(1),
						in(0) + n(0) * in(1),
						i(0) + n(0) * in(1)
					});
				}
			}
		}

		if (capmin(0)) {
			for (int j = 0; j < imax(1); ++j) {
				int jn = (j + 1) % n(1);
				addPoly(mesh, { 
					0 + n(0) * j,
					0 + n(0) * jn,
					capindex
				});
			}
		}

		mesh->calcAux();
		return mesh;
	}
};

}
}
