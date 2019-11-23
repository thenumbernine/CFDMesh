#pragma once

#include "CFDMesh/Vector.h"
#include "CFDMesh/Util.h"
#include "GLApp/gl.h"
#include "Image/Image.h"
#include "Common/File.h"
#include "Common/Macros.h"
#include "Common/Exception.h"
#include <vector>
#include <list>
#include <string>
#include <cassert>


template<typename T> void glVertex2v(const T* v);
template<> void glVertex2v<float>(const float* v) { glVertex2fv(v); }
template<> void glVertex2v<double>(const double* v) { glVertex2dv(v); }

template<typename T> void glVertex3v(const T* v);
template<> void glVertex3v<float>(const float* v) { glVertex3fv(v); }
template<> void glVertex3v<double>(const double* v) { glVertex3dv(v); }


template<typename real3>
inline real3 cross(real3 a, real3 b) {
	return real3(
		a(1) * b(2) - a(2) * b(1),
		a(2) * b(0) - a(0) * b(2),
		a(0) * b(1) - a(1) * b(0));
}


namespace CFDMesh {

//zero-forms
template<typename real>
struct Vertex_ {	//not required by finite volume algorithm
	using This = Vertex_;
	using real3 = Tensor::Vector<real, 3>;
	
	real3 pos;
	
	//keeping track of Vertex_::faces isn't used by the renderer, mesh generation, or finite-volume integration
	//std::vector<int> faces;
	
	static constexpr auto fields = std::make_tuple(
		std::make_pair("pos", &This::pos)
	);
	
	Vertex_() {}
	Vertex_(real3 pos_) : pos(pos_) {}
};

//(n-1)-forms
template<typename real, typename Cons>
struct Face_ {
	using This = Face_;
	using real3 = Tensor::Vector<real, 3>;
	
	real3 pos;
	real3 delta;
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
		std::make_pair("delta", &This::delta),
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

//n-forms
template<typename real, typename Cons>
struct Cell_ {
	using This = Cell_;
	using real3 = Tensor::Vector<real, 3>;

	real3 pos;
	real volume = 0;
	
	//required by the finite volume algorithm
	//std::vector<int> faces;
	int faceOffset = 0;
	int faceCount = 0;

	//not required by finite volume algorithm
	//however the cell volume is required, and is calculated using vtxs
	//also the renderer requires the vertexes
	//std::vector<int> vtxs;
	int vtxOffset = 0;
	int vtxCount = 0;

	float displayValue = 0;
	Cons U;
	
	static constexpr auto fields = std::make_tuple(
		std::make_pair("pos", &This::pos),
		std::make_pair("volume", &This::volume),
		std::make_pair("faceOffset", &This::faceOffset),
		std::make_pair("faceCount", &This::faceCount),
		std::make_pair("vtxOffset", &This::vtxOffset),
		std::make_pair("vtxCount", &This::vtxCount),
		std::make_pair("displayValue", &This::displayValue),
		std::make_pair("U", &This::U)
	);
};

}


template<typename T>
std::ostream& operator<<(std::ostream& a, const CFDMesh::Vertex_<T>& b) {
	return ostreamForFields(a, b);
}

template<typename T, typename C>
std::ostream& operator<<(std::ostream& a, const CFDMesh::Face_<T, C>& b) {
	return ostreamForFields(a, b);
}

template<typename T, typename C>
std::ostream& operator<<(std::ostream& a, const CFDMesh::Cell_<T, C>& b) {
	return ostreamForFields(a, b);
}

namespace std {

template<typename T>
std::string to_string(const CFDMesh::Vertex_<T>& x) {
	return objectStringFromOStream(x);
}

template<typename T, typename C>
std::string to_string(const CFDMesh::Face_<T, C>& x) {
	return objectStringFromOStream(x);
}

template<typename T, typename C>
std::string to_string(const CFDMesh::Cell_<T, C>& x) {
	return objectStringFromOStream(x);
}

}


namespace CFDMesh {

//dim is the dimension of the manifold, not of the vectors (which are all 3D atm)
template<typename real, int dim, typename Cons>
struct MeshNamespace {
using real2 = Tensor::Vector<real, 2>;
using real3 = Tensor::Vector<real, 3>;
	
using Vertex = Vertex_<real>;
using Face = Face_<real, Cons>;
using Cell = Cell_<real, Cons>;


struct MeshFactory;

struct Mesh {
protected:
	//https://stackoverflow.com/questions/8147027/how-do-i-call-stdmake-shared-on-a-class-with-only-protected-or-private-const
	struct ctorkey {
		explicit ctorkey(int) {}
	};

	Mesh(const Mesh&) = delete;
	const Mesh& operator=(const Mesh&) = delete;

	template <typename... T>
	static ::std::shared_ptr<Mesh> create(T &&...args) {
		return ::std::make_shared<Mesh>(ctorkey{0}, ::std::forward<T>(args)...);
	}

public:
	std::vector<Vertex> vtxs;	//0-forms, which construct n and n-1 forms
	std::vector<Face> faces;	//n-1-forms, hold flux information between cells
	std::vector<Cell> cells;	//n-forms, hold the info at each cell
	
	std::vector<int> cellFaceIndexes;
	std::vector<int> cellVtxIndexes;
	std::vector<int> faceVtxIndexes;

	friend struct MeshFactory;
	
	explicit Mesh(const ctorkey&) {}
	virtual ~Mesh() {}


	int addFaceForVtxs(const int* vs, int n) {
		int fi = 0;
		for (; fi < (int)faces.size(); ++fi) {
			Face* f = &faces[fi];
			if constexpr (dim == 2) {
				int va = vs[0];
				int vb = vs[1];
				if (f->vtxCount == 2) {
					if ((faceVtxIndexes[f->vtxOffset+0] == va && faceVtxIndexes[f->vtxOffset+1] == vb) ||
						(faceVtxIndexes[f->vtxOffset+0] == vb && faceVtxIndexes[f->vtxOffset+1] == va)) 
					{
						return fi;
					}
				}
			} else if constexpr (dim == 3) {
				for (int j = 0; j < n; ++j) {
					//check in one direction
					bool matches = true;
					for (int i = 0; i < n; ++i) {
						if (faceVtxIndexes[f->vtxOffset+i] != vs[(j+i)%n]) {
							matches = false;
							break;
						}
					}
					if (matches) return fi;
				
					//check in the other direction
					matches = true;
					for (int i = 0; i < n; ++i) {
						if (faceVtxIndexes[f->vtxOffset+i] != vs[(j+n-i)%n]) {
							matches = false;
							break;
						}
					}
					if (matches) return fi;
				}
			}
		}
		
		assert(fi == (int)faces.size());
		faces.push_back(Face());
		Face& f = faces.back();
		
		f.vtxOffset = faceVtxIndexes.size();
		for (int i = 0; i < n; ++i) {
			faceVtxIndexes.push_back(vs[i]);
		}
		f.vtxCount = faceVtxIndexes.size() - f.vtxOffset;

		//TODO calc these for n=3
		if constexpr (dim == 2) {
			auto& a = vtxs[vs[0]];
			auto& b = vtxs[vs[1]];
			f.pos = (a.pos + b.pos) * .5;
			f.delta = a.pos - b.pos;
			f.area = real3::length(f.delta);
			f.normal = real3(f.delta(1), -f.delta(0));
			f.normal *= 1. / real3::length(f.normal);
		} else if constexpr (dim == 3) {
			std::vector<real3> polyVtxs(n);
			for (int i = 0; i < n; ++i) {
				polyVtxs[i] = vtxs[vs[i]].pos;
			}
			f.normal = cross(polyVtxs[1] - polyVtxs[0], polyVtxs[2] - polyVtxs[1]);
			f.normal *= 1. / real3::length(f.normal);
			f.area = polygon3DVolume(polyVtxs, f.normal);
			f.pos = polygon3DCOM(polyVtxs, f.area, f.normal);
		} else {
			throw Common::Exception() << "here";
		}
		return fi;
	}

	
	int addFace(const int* vs, int n, int ci) {
		int fi = addFaceForVtxs(vs, n);
		Face& f = faces[fi];	
		if (f.cells(0) == -1) {
			f.cells(0) = ci;
		} else if (f.cells(1) == -1) {
			f.cells(1) = ci;
		} else {
			throw Common::Exception() << "tried to add too many cells to an edge";
		}
		return fi;
	}


	void addCell(std::vector<int> vis) {
		int ci = (int)cells.size();
		cells.push_back(Cell());
		Cell& c = cells.back();
		
		size_t n = vis.size();
	
		//c.vtxs = vis;
		c.vtxOffset = (int)cellVtxIndexes.size();
		cellVtxIndexes.insert(cellVtxIndexes.end(), vis.begin(), vis.end());
		c.vtxCount = (int)cellVtxIndexes.size() - c.vtxOffset;

		c.faceOffset = (int)cellFaceIndexes.size();
		if constexpr (dim == 2) {
			//face is a 1-form
			for (size_t i = 0; i < n; ++i) {
				int vtxs[2] = {vis[i], vis[(i+1)%n]}; 
				cellFaceIndexes.push_back(addFace(vtxs, numberof(vtxs), ci));
			}
			std::vector<real2> polyVtxs = map<
				std::vector<int>,
				std::vector<real2>
			>(vis, [this](int vi) -> real2 {
				return real2([this, vi](int i) -> real { return vtxs[vi].pos(i); });
			});
			
			c.volume = polygonVolume(polyVtxs);
			real2 com = polygonCOM(polyVtxs, c.volume);
			c.pos = real3([&com](int i) -> real { return i < 2 ? com(i) : 0; });
		
		} else if constexpr (dim == 3) {
			//face is a 2-form
			assert(vis.size() == 8);	//only adding cubes at the moment
/*
  6----7
 /|   /|
4----5 |
| |  | |  z
| 2--|-3  ^ y
|/   |/   |/
0----1    *-> x
*/
			std::vector<std::vector<int>> cubeSides = {
				{0,4,6,2},	//x-
				{1,3,7,5},	//x+
				{0,1,5,4},	//y-
				{2,6,7,3},	//y+
				{0,2,3,1},	//z-
				{4,5,7,6},	//z+
			};
		
			for (auto side : cubeSides) {
				std::vector<int> vtxIndexes = map<
					std::vector<int>,
					std::vector<int>
				>(side, [&vis](int side_i) -> int {
					return vis[side_i];
				});
				int f = addFace(vtxIndexes.data(), vtxIndexes.size(), ci);
				cellFaceIndexes.push_back(f);
			}

			std::vector<real3> polyVtxs = map<
				std::vector<int>,
				std::vector<real3>
			>(vis, [this](int vi) -> real3 {
				return vtxs[vi].pos;
			});

			c.volume = polyhedronVolume(polyVtxs, cubeSides);
			c.pos = polyhedronCOM(polyVtxs, cubeSides, c.volume);

		} else {
			throw Common::Exception() << "you are here";
		}
		c.faceCount = (int)cellFaceIndexes.size() - c.faceOffset;

		assert(c.volume > 0);
	}

	//calculate edge info
	//calculate cell volume info
	void calcAux() {
		for (auto& e : faces) {
			int a = e.cells(0);
			int b = e.cells(1);
			if (a != -1 && b != -1) {
				if (real3::dot(cells[a].pos - cells[b].pos, e.normal) < 0) {
					std::swap(a, b);
					e.cells(0) = a;
					e.cells(1) = b;
					e.normal *= -1;
				}
				//distance between cell centers
				//e.cellDist = real3::length(cells[b].pos - cells[a].pos);
				//distance projected to edge normal
				e.cellDist = fabs(real3::dot(e.normal, cells[b].pos - cells[a].pos));
			} else if (a != -1) {
				e.cellDist = real3::length(cells[a].pos - e.pos) * 2.;
			} else {
				throw Common::Exception() << "looks like you created a face that isn't touching any cells...";
			}
		}
	}

	struct DrawArgs {
		using This = DrawArgs;
		
		int gradientTex = 0;
		std::pair<float, float> displayValueRange = {};
		int selectedCellIndex = -1;
		bool showCells = true;
		bool showVtxs = false;
		bool showFaces = false;
		bool showFaceCenters = false;
		bool showCellCenters = false;
		float cellScale = 1;

		static constexpr auto fields = std::make_tuple(
			std::make_pair("display value range", &This::displayValueRange),
			std::make_pair("show cells", &This::showCells),
			std::make_pair("show vertexes", &This::showVtxs),
			std::make_pair("show faces", &This::showFaces),
			std::make_pair("show face centers", &This::showFaceCenters),
			std::make_pair("show cell centers", &This::showCellCenters),
			std::make_pair("cell scale", &This::cellScale)
		);
	};

	void draw(DrawArgs args) {
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		if (args.selectedCellIndex != -1) {
			Cell& c = cells[args.selectedCellIndex];
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
			glLineWidth(3);
			glColor3f(1,1,0);
			glBegin(GL_POLYGON);
			for (int vi = 0; vi < c.vtxCount; ++vi) {
				const auto& v = vtxs[cellVtxIndexes[vi + c.vtxOffset]].pos;
				glVertex3v(((v - c.pos) * args.cellScale + c.pos).v);
			}
			glEnd();
			glColor3f(1,1,1);
			glLineWidth(1);
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		}

		glEnable(GL_TEXTURE_1D);
		glBindTexture(GL_TEXTURE_1D, args.gradientTex);

		if (args.showCells) {
			for (const auto& c : cells) {
				real f = (c.displayValue - args.displayValueRange.first) / (args.displayValueRange.second - args.displayValueRange.first);
				glColor3f(1,1,1);
				glTexCoord1f(f);
				
				if constexpr (dim == 2) {
					glBegin(GL_POLYGON);
					for (int vi = 0; vi < c.vtxCount; ++vi) {
						const auto& v = vtxs[cellVtxIndexes[vi + c.vtxOffset]].pos;
						glVertex3v(((v - c.pos) * args.cellScale + c.pos).v);
					}
					glEnd();
				} else if constexpr (dim == 3) {
					for (int i = 0; i < c.faceCount; ++i) {
						Face& f = faces[cellFaceIndexes[i + c.faceOffset]];
						glBegin(GL_POLYGON);
						for (int vi = 0; vi < f.vtxCount; ++vi) {
							const auto& v = vtxs[faceVtxIndexes[vi + f.vtxOffset]].pos;
							glVertex3v(((v - c.pos) * args.cellScale + c.pos).v);
						}
						glEnd();
					}
				}
			}
		}

		glBindTexture(GL_TEXTURE_1D, 0);
		glDisable(GL_TEXTURE_1D);
		
		if (args.showVtxs || args.showCellCenters || args.showFaceCenters) {
			glPointSize(3);
			glBegin(GL_POINTS);
			if (args.showVtxs) {
				glColor3f(1,1,0);
				for (const auto& v : vtxs) {
					glVertex3v(v.pos.v);
				}
			}
			if (args.showFaceCenters) {
				glColor3f(1,0,1);
				for (const auto& f : faces) {
					glVertex3v(f.pos.v);
				}
			}
			if (args.showCellCenters) {
				glColor3f(0,1,1);
				for (const auto& c : cells) {
					glVertex3v(c.pos.v);
				}
			}
			glEnd();
			glColor3f(1,1,1);
			glPointSize(1);
		}
		
		if (args.showFaces) {
			for (const auto& f : faces) {
				glBegin(GL_LINE_LOOP);
				for (int vi = 0; vi < f.vtxCount; ++vi) {
					const auto& v = vtxs[faceVtxIndexes[vi + f.vtxOffset]].pos;
					glVertex3v(((v - f.pos) * args.cellScale + f.pos).v);
				}
				glEnd();
			}
		}

	}
};


//	3D - polygon


static real polygon3DVolume(const std::vector<real3>& vs, real3 normal) {
	size_t n = vs.size();
	real area = 0;
	for (size_t i = 0; i < n; ++i) {
		const real3 &a = vs[i];
		const real3 &b = vs[(i+1)%n];
		area += parallelepipedVolume(a, b, normal);
	}
	return .5 * area;
}

static real3 polygon3DCOM(const std::vector<real3>& vs, real area, real3 normal) {
	int n = (int)vs.size();
	real3 com;
#if 0
	for (int i = 0; i < n; ++i) {
		const real3& a = vs[i];
		const real3& b = vs[(i+1)%n];
		com += (a + b) * parallelepipedVolume(a, b, normal);
	}
	return com * (1. / (6. * area));
#else
	const real3& a = vs[0];
	for (int i = 2; i < n; ++i) {
		const real3& b = vs[i-1];
		const real3& c = vs[i];
		com += (a + b + c) * real3::dot(cross(c - a, c - b), normal);
	}
	return com * (1. / (6. * area));
#endif
}



//	2D - parallelogram


static real parallelogramVolume(real2 a, real2 b) {
	//epsilon_ij a^i b^j
	return a(0) * b(1)
		- a(1) * b(0);
}


//	2D - polygon


static real polygonVolume(const std::vector<real2>& vs) {
	size_t n = vs.size();
	real volume = 0;
	for (size_t i = 0; i < n; ++i) {
		const real2 &a = vs[i];
		const real2 &b = vs[(i+1)%n];
		volume += parallelogramVolume(a, b);
	}
	return .5 * volume;
}

static real2 polygonCOM(const std::vector<real2>& vs, real volume) {
	size_t n = vs.size();
	real2 com;
	for (int i = 0; i < (int)n; ++i) {
		const real2& a = vs[i];
		const real2& b = vs[(i+1)%n];
		com += (a + b) * parallelogramVolume(a, b);
	}
	return com * (1. / (6. * volume));
}

template<typename I, typename F>
static bool polygonContains(real2 pos, I begin, I end, F f) {
	auto check = [pos](real2 prev, real2 cur) -> bool {
		return parallelogramVolume(prev - pos, cur - pos) > 0;
	};
	
	I i = begin;
	real2 first = f(*i);
	real2 prev = first;
	if (i != end) {
		for (++i; i != end; ++i) {
			real2 cur = f(*i);
			if (!check(prev, cur)) return false;
			prev = cur;
		}
	}
	if (!check(prev, first)) return false;
	return true;
}


//	3D - parallelepiped


static real parallelepipedVolume(real3 a, real3 b, real3 c) {
	//epsilon_ijk a^i b^j c^k
	return a(0) * b(1) * c(2)
		+ a(1) * b(2) * c(0)
		+ a(2) * b(0) * c(1)
		- c(0) * b(1) * a(2)
		- c(1) * b(2) * a(0)
		- c(2) * b(0) * a(1);
}


//	3D - polyhedron


static real polyhedronVolume(const std::vector<real3>& vtxs, const std::vector<std::vector<int>>& faces) {
	real volume = 0;
	for (const auto& face : faces) {
		for (int i = 2; i < (int)face.size(); ++i) {
			//tri from 0, i-1, 1
			const real3& a = vtxs[face[0]];
			const real3& b = vtxs[face[i-1]];
			const real3& c = vtxs[face[i]];
			volume += parallelepipedVolume(a, b, c);
		}
	}
	//volume of n-sided pyramid in nD is volume of parallelogram divided by 1/n!
	return volume / 6.;
}

static real3 polyhedronCOM(const std::vector<real3>& vtxs, const std::vector<std::vector<int>>& faces, real volume) {
	real3 com;
	for (const auto& face : faces) {
		for (int i = 2; i < (int)face.size(); ++i) {
			//tri from 0, i-1, 1
			const real3& a = vtxs[face[0]];
			const real3& b = vtxs[face[i-1]];
			const real3& c = vtxs[face[i]];
			com += ((a + b) * (a + b) + (b + c) * (b + c) + (c + a) * (c + a)) * cross(b - a, c - a) / 48.;
		}
	}
	return com / volume;
}


struct MeshFactory {
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

struct FileMeshFactory : public MeshFactory {
	std::string filename = {"grids/n0012_113-33.p2dfmt"};
	
	FileMeshFactory() : MeshFactory("p2dfmt mesh") {}

	virtual std::shared_ptr<Mesh> createMesh() {
		std::shared_ptr<Mesh> mesh = MeshFactory::createMeshSuper();
		
		std::list<std::string> ls = split<std::list<std::string>>(Common::File::read(filename), "\n");
	
		std::string first = ls.front();
		ls.pop_front();
		std::vector<std::string> m_n = split<std::vector<std::string>>(ls.front(), "\\s+");
		ls.pop_front();
		int m = std::stoi(m_n[0]);
		int n = std::stoi(m_n[1]);
		std::list<std::string> _x = split<std::list<std::string>>(concat<std::list<std::string>>(ls, " "), "\\s+");
		if (_x.front() == "") _x.pop_front();
		if (_x.back() == "") _x.pop_back();
		std::vector<real> x = map<std::list<std::string>, std::vector<real>>(_x, [](const std::string& s) -> real { return std::stod(s); });
		assert(x.size() == (size_t)(2 * m * n));
	
		auto us = std::vector(x.begin(), x.begin() + m*n);
		auto vs = std::vector(x.begin() + m*n, x.end());
		assert(us.size() == vs.size());

		mesh->vtxs.resize(m*n);
		for (int i = 0; i < (int)us.size(); ++i) {
			mesh->vtxs[i].pos = real3(us[i], vs[i]);
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

struct Chart2DMeshFactory : public MeshFactory {
	using This = Chart2DMeshFactory; 
	
	int2 size = int2(100, 100);
	float2 mins = real2(-1, -1);
	float2 maxs = real2(1, 1);
	bool2 repeat = bool2(false, false);
	bool2 capmin = bool2(false, false);
	
	static constexpr auto fields = std::make_tuple(
		std::make_pair("size", &This::size),
		std::make_pair("mins", &This::mins),
		std::make_pair("maxs", &This::maxs),
		std::make_pair("repeat", &This::repeat),
		std::make_pair("capmin", &This::capmin)
	);

	Chart2DMeshFactory(const char* name_) : MeshFactory(name_) {}

	virtual real2 coordChart(real2 x) const { return x; }
	
	virtual void updateGUI() {
		CFDMesh::updateGUI(this);
	}
};

struct TriUnitMeshFactory : public Chart2DMeshFactory {
	using Super = Chart2DMeshFactory;
	TriUnitMeshFactory() : Super("unit square of triangles") {}
	
	virtual std::shared_ptr<Mesh> createMesh() {
		std::shared_ptr<Mesh> mesh = MeshFactory::createMeshSuper();
		
		int m = Super::size(0) + 1;
		int n = Super::size(1) + 1;
	
		mesh->vtxs.resize(m * n);
		for (int j = 0; j < n; ++j) {
			for (int i = 0; i < m; ++i) {
				real2 x = real2(
					(real)i / (real)Super::size(0) * (Super::maxs(0) - Super::mins(0)) + Super::mins(0),
					(real)j / (real)Super::size(1) * (Super::maxs(1) - Super::mins(1)) + Super::mins(1));
				
				real2 u = Super::coordChart(x);
				mesh->vtxs[i + m * j].pos = real3([&u](int i) -> real { return i < real2::size ? u(i) : 0.; });
			}
		}
		
		int imax = Super::repeat(0) ? m : m-1;
		int jmax = Super::repeat(1) ? n : n-1;
		for (int j = 0; j < jmax; ++j) {
			int jn = (j + 1) % n;
			for (int i = 0; i < imax; ++i) {
				int in = (i + 1) % m;
				mesh->addCell(std::vector<int>{i + m * j, in + m * j, in + m * jn});
				mesh->addCell(std::vector<int>{in + m * jn, i + m * jn, i + m * j});
			}
		}
	
		mesh->calcAux();
		return mesh;
	}
};

struct QuadUnitMeshFactory : public Chart2DMeshFactory {
	using Super = Chart2DMeshFactory;
	QuadUnitMeshFactory(const char* name_ = "unit square of quads") : Super(name_) {}
	
	virtual std::shared_ptr<Mesh> createMesh() {
		std::shared_ptr<Mesh> mesh = MeshFactory::createMeshSuper();

		int2 n = Super::size + 1;
		int2 step(1, n(0));
		int vtxsize = n.volume();
		if (Super::capmin(0)) vtxsize++;
		mesh->vtxs.resize(vtxsize);
		int2 i;
		for (i(1) = 0; i(1) < n(1); ++i(1)) {
			for (i(0) = 0; i(0) < n(0); ++i(0)) {
				real2 x = (real2)i / (real2)Super::size * (Super::maxs - Super::mins) + Super::mins;
				
				//if I do not put Super:: then I get: "error: there are no arguments to ‘coordChart’ that depend on a template parameter, so a declaration of ‘coordChart’ must be available"
				//real2 u = coordChart(x);
			
				//if I say 'Super::coordChart' then it will explicitly call that class' coordChart method (and not the vtable entry)
				//real2 u = Super::coordChart(x);
				
				real2 u = this->coordChart(x);
				
				mesh->vtxs[int2::dot(i, step)].pos = real3([&u](int i) -> real { return i < real2::size ? u(i) : 0.; });
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

struct QuadUnitCbrtMeshFactory : public QuadUnitMeshFactory {
	QuadUnitCbrtMeshFactory() : QuadUnitMeshFactory("unit square of quads, cbrt mapping") {}
	virtual real2 coordChart(real2 v) const {
		return real2(cbrt(v(0)), cbrt(v(1)));
	}
};

template<typename T> 
static T cubed(const T& t) { return t * t * t; }

struct QuadUnitCubedMeshFactory : public QuadUnitMeshFactory {
	QuadUnitCubedMeshFactory() : QuadUnitMeshFactory("unit square of quads, cubed mapping") {}
	virtual real2 coordChart(real2 v) const {
		return real2(cubed(v(0)), cubed(v(1)));
	}
};

struct TwistQuadUnitMeshFactory : public QuadUnitMeshFactory {
	TwistQuadUnitMeshFactory() : QuadUnitMeshFactory("unit square of quads, twist in the middle") {}
	virtual real2 coordChart(real2 v) const {
		real r = real2::length(v);
		//real theta = std::max(0., 1. - r);
		real sigma = 3.;	//almost 0 at r=1
		const real rotationAmplitude = 3.;
		real theta = rotationAmplitude*sigma*r*exp(-sigma*sigma*r*r);
		real costh = cos(theta), sinth = sin(theta);
		return real2(
			costh * v(0) - sinth * v(1),
			sinth * v(0) + costh * v(1));	
	}
};

struct DonutQuadUnitMeshFactory : public QuadUnitMeshFactory {
	using Super = QuadUnitMeshFactory;
	DonutQuadUnitMeshFactory() : Super("polar") {
		Super::size = int2(50, 200);
		Super::mins = real2(.5, 0);
		Super::maxs = real2(1, 2*M_PI);
		Super::repeat = bool2(false, true);
		//Super::capmin = int2(true, false);
	}
	virtual real2 coordChart(real2 v) const {
		return real2(cos(v(1)), sin(v(1))) * v(0);
	}
};

//not inheriting from QuadUnitMeshFactory because it has variable size and we want fixed size (based on image size)
struct QuadUnitBasedOnImageMeshFactory : public Chart2DMeshFactory {
	using This = QuadUnitBasedOnImageMeshFactory;
	using Super = Chart2DMeshFactory;
	QuadUnitBasedOnImageMeshFactory() : Super("unit based on image") {}

	std::string imageFilename = "layout.png";
	//std::string imageFilename = "layout.bmp";
	//std::string imageFilename = "layout.tiff";

	virtual std::shared_ptr<Mesh> createMesh() {
		std::shared_ptr<Mesh> mesh = MeshFactory::createMeshSuper();
		
		//auto iimg = ::Image::system->read(imageFilename);	//TODO fixme, it's not working
		auto iimg = ::Image::pngIO->read(imageFilename);
		
		auto img = std::dynamic_pointer_cast<Image::Image>(iimg);

		Super::size = img->getSize();
		
		int m = Super::size(0) + 1;
		int n = Super::size(1) + 1;
	
		mesh->vtxs.resize(m * n);
		for (int j = 0; j < n; ++j) {
			for (int i = 0; i < m; ++i) {
				real2 x = real2(
					(real)i / (real)Super::size(0) * (Super::maxs(0) - Super::mins(0)) + Super::mins(0),
					(real)j / (real)Super::size(1) * (Super::maxs(1) - Super::mins(1)) + Super::mins(1));
				
				real2 u = Super::coordChart(x);
				mesh->vtxs[i + m * j].pos = real3([&u](int i) -> real { return i < real2::size ? u(i) : 0.; });
			}
		}

		int imax = m-1;
		int jmax = n-1;
		for (int j = 0; j < jmax; ++j) {
			int jn = (j + 1) % n;
			for (int i = 0; i < imax; ++i) {
				int in = (i + 1) % m;
				if ((*img)(i, jmax-1-j)) {
					mesh->addCell(std::vector<int>{i + m * j, in + m * j, in + m * jn, i + m * jn});
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

struct Chart3DMeshFactory : public MeshFactory {
	using This = Chart3DMeshFactory;
	
	int3 size = int3(10,10,10);
	float3 mins = float3(-1, -1, -1);
	float3 maxs = float3(1, 1, 1);
	bool3 repeat;

	//TODO support for inheritence and reflection
	static constexpr auto fields = std::make_tuple(
		std::make_pair("size", &This::size),
		std::make_pair("mins", &This::mins),
		std::make_pair("maxs", &This::maxs),
		std::make_pair("repeat", &This::maxs)
	);

	Chart3DMeshFactory(const char* name_ = "3D chart mesh") : MeshFactory(name_) {}

	virtual real3 coordChart(real3 x) const { return x; }

	virtual void updateGUI() {
		CFDMesh::updateGUI(this);
	}
};

struct CubeUnitMeshFactory : public Chart3DMeshFactory {
	using Super = Chart3DMeshFactory;
	CubeUnitMeshFactory(const char* name_ = "unit cube of cubes") : Super(name_) {}
	
	virtual std::shared_ptr<Mesh> createMesh() {
		std::shared_ptr<Mesh> mesh = MeshFactory::createMeshSuper();
	
		int3 n = Super::size + 1;
		int vtxsize = n.volume();
		mesh->vtxs.resize(vtxsize);
		int3 i;
		for (i(2) = 0; i(2) < n(2); ++i(2)) {
			for (i(1) = 0; i(1) < n(1); ++i(1)) {
				for (i(0) = 0; i(0) < n(0); ++i(0)) {
					real3 x = (real3)i / (real3)Super::size * (Super::maxs - Super::mins) + Super::mins;
					mesh->vtxs[i(0) + n(0) * (i(1) + n(1) * i(2))].pos = Super::coordChart(x);
				}
			}
		}
		
		int3 imax;
		for (int j = 0; j < 3; ++j) {
			imax(j) = Super::repeat(j) ? n(j) : n(j)-1;
		}
		int3 in;
		for (i(2) = 0; i(2) < imax(2); ++i(2)) {
			in(2) = (i(2) + 1) % n(1);
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

static std::vector<std::shared_ptr<MeshFactory>> getGens() {
	if constexpr (dim == 2) {
		return std::vector<std::shared_ptr<MeshFactory>>{
			std::make_shared<QuadUnitMeshFactory>(),
			std::make_shared<TriUnitMeshFactory>(),
			std::make_shared<QuadUnitCbrtMeshFactory>(),
			std::make_shared<QuadUnitCubedMeshFactory>(),
			std::make_shared<TwistQuadUnitMeshFactory>(),
			std::make_shared<DonutQuadUnitMeshFactory>(),
			std::make_shared<QuadUnitBasedOnImageMeshFactory>(),
			std::make_shared<FileMeshFactory>(),
		};
	} else if constexpr (dim == 3) {
		return std::vector<std::shared_ptr<MeshFactory>>{
			std::make_shared<CubeUnitMeshFactory>(),
		};
	}
	throw Common::Exception() << "here";
}

};	//MeshNamespace

}
