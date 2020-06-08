#pragma once

#include "CFDMesh/Mesh/Quad2DMeshFactory.h"
#include "CFDMesh/GUI.h"
#include "Image/Image.h"

namespace CFDMesh {
namespace Mesh {

//not inheriting from Quad2DMeshFactory because it has variable size and we want fixed size (based on image size)
template<typename real, int dim, typename Cons>
struct Quad2DImageMeshFactory : public Quad2DMeshFactory<real, dim, Cons> {
	using This = Quad2DImageMeshFactory;
	using Super = Quad2DMeshFactory<real, dim, Cons>;
	using Mesh = typename Super::Mesh;
	using real2 = typename Super::real2;
	using real3 = typename Super::real3;

	std::string imageFilename;
	std::shared_ptr<Image::Image> img;

	Quad2DImageMeshFactory(
		std::string imageFilename_ = "layout.png"
		//std::string imageFilename_ = "layout.bmp"
		//std::string imageFilename_ = "layout.tiff"
	) : Super("unit based on image"), imageFilename(imageFilename_) {}

	virtual bool testMakeCell(int2 i) {
		return (*img)(i(0), Super::size(1)-1-i(1));
	}

	virtual std::shared_ptr<Mesh> createMesh() {
		//auto iimg = ::Image::system->read(imageFilename);	//TODO fixme, it's not working
		auto iimg = ::Image::pngIO->read(imageFilename);
		
		img = std::dynamic_pointer_cast<Image::Image>(iimg);
		Super::size = img->getSize();
		
		auto result = Super::createMesh();
		
		img = nullptr;
		
		return result;
	}

	static constexpr auto fields = std::tuple_cat(
		Super::fields,		//TODO remove Super::size
		std::make_tuple(
			std::make_pair("imageFilename", &This::imageFilename)
		)
	);

	virtual void updateGUI() {
		CFDMesh::updateGUI(this);
	}
};

}
}
