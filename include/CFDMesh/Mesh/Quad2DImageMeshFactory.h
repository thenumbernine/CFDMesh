#pragma once

#include "CFDMesh/Mesh/Chart2DMeshFactory.h"
#include "CFDMesh/GUI.h"
#include "Image/Image.h"

namespace CFDMesh {
namespace Mesh {

//not inheriting from Quad2DMeshFactory because it has variable size and we want fixed size (based on image size)
template<typename real, int dim, typename Cons>
struct Quad2DImageMeshFactory : public Chart2DMeshFactory<real, dim, Cons> {
	using This = Quad2DImageMeshFactory;
	using Super = Chart2DMeshFactory<real, dim, Cons>;
	using Mesh = typename Super::Mesh;
	using real2 = typename Super::real2;
	using real3 = typename Super::real3;

	Quad2DImageMeshFactory() : Super("unit based on image") {}

	std::string imageFilename = "layout.png";
	//std::string imageFilename = "layout.bmp";
	//std::string imageFilename = "layout.tiff";

	virtual std::shared_ptr<Mesh> createMesh() {
		std::shared_ptr<Mesh> mesh = MeshFactory<real, dim, Cons>::createMeshSuper();
		
		//auto iimg = ::Image::system->read(imageFilename);	//TODO fixme, it's not working
		auto iimg = ::Image::pngIO->read(imageFilename);
		
		auto img = std::dynamic_pointer_cast<Image::Image>(iimg);

		Super::size = img->getSize();
		int2 n = Super::size + 1;
		int2 step(1, n(0));	
		mesh->vtxs.resize(n.volume());
		int2 i;
		for (i(1) = 0; i(1) < n(1); ++i(1)) {
			for (i(0) = 0; i(0) < n(0); ++i(0)) {
				real2 x = (real2)i / (real2)Super::size * (Super::maxs - Super::mins) + Super::mins;
				
				real2 u = Super::coordChart(x);
				mesh->vtxs[int2::dot(i, step)].pos = real3([&u](int i) -> real { return i < real2::size ? u(i) : 0.; });
			}
		}

		int2 imax = n - 1;
		int2 in;
		for (i(1) = 0; i(1) < imax(1); ++i(1)) {
			in(1) = (i(1) + 1) % n(1);
			for (i(0) = 0; i(0) < imax(0); ++i(0)) {
				in(0) = (i(0) + 1) % n(0);
				if ((*img)(i(0), imax(1)-1-i(1))) {
					mesh->addCell(std::vector<int>{
						i(0) + n(0) * i(1),
						in(0) + n(0) * i(1),
						in(0) + n(0) * in(1),
						i(0) + n(0) * in(1)
					});
				}
			}
		}

		mesh->calcAux();
		return mesh;
	}

	//override Chart2DMeshFactory and get rid of size
	//TODO add in imageFilename
	//TODO TODO if you use fields for anything else, consider a readonly modifier
	static constexpr auto fields = std::make_tuple(
		std::make_pair("mins", &Super::mins),
		std::make_pair("maxs", &Super::maxs),
		std::make_pair("repeat", &Super::repeat),
		std::make_pair("capmin", &Super::capmin),
		std::make_pair("imageFilename", &This::imageFilename)
	);

	virtual void updateGUI() {
		CFDMesh::updateGUI(this);
	}
};

}
}
