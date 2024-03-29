#pragma once

#include "CFDMesh/Vector.h"
#include "Common/File.h"
#include "Common/Meta.h"
#include "Common/String.h"	//Common::split
#include <list>
#include <vector>

namespace CFDMesh {
namespace Mesh {

template<typename real, int dim, typename Cons>
struct P2DFMTMeshFactory : public MeshFactory<real, dim, Cons> {
	using Super = MeshFactory<real, dim, Cons>;
	using Mesh = typename Super::Mesh;
	using real3 = typename Super::real3;

	std::string filename = {"grids/n0012_113-33.p2dfmt"};
	
	P2DFMTMeshFactory() : MeshFactory<real, dim, Cons>("p2dfmt mesh") {}

	virtual std::shared_ptr<Mesh> createMesh() {
		std::shared_ptr<Mesh> mesh = MeshFactory<real, dim, Cons>::createMeshSuper();
		
		std::list<std::string> ls = Common::split<std::list<std::string>>(Common::File::read(filename), "\n");
	
		std::string first = ls.front();
		ls.pop_front();
		std::vector<std::string> m_n = Common::split<std::vector<std::string>>(ls.front(), "\\s+");
		ls.pop_front();
		int m = std::stoi(m_n[0]);
		int n = std::stoi(m_n[1]);
		std::list<std::string> x2 = Common::split<std::list<std::string>>(Common::concat<std::list<std::string>>(ls, " "), "\\s+");
		if (x2.front() == "") x2.pop_front();
		if (x2.back() == "") x2.pop_back();
		auto x = Common::mapElems<std::vector<real>>(x2, [](std::string const & s) -> real { return std::stod(s); });
		assert(x.size() == (size_t)(2 * m * n));
	
		auto us = std::vector(x.begin(), x.begin() + m*n);
		auto vs = std::vector(x.begin() + m*n, x.end());
		assert(us.size() == vs.size());

		mesh->vtxs.resize(m*n);
		for (int i = 0; i < (int)us.size(); ++i) {
			mesh->vtxs[i].pos = real3(us[i], vs[i], 0);
		}
	
		for (int j = 0; j < n-1; ++j) {
			for (int i = 0; i < m-1; ++i) {
				mesh->addCell(std::vector<int>{i + m * j, i + m * (j+1), i+1 + m * (j+1), i+1 + m * j});
			}
		}
	
		mesh->calcAux();
		return mesh;
	}

	virtual void updateGUI() {
		//TODO filename popup
	}
};

}
}
