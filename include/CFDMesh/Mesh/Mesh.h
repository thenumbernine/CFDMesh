#pragma once

#include "CFDMesh/Mesh/MeshFactory.h"
#include "CFDMesh/Mesh/P2DFMTMeshFactory.h"
#include "CFDMesh/Mesh/Quad2DMeshFactory.h"
#include "CFDMesh/Mesh/Quad2DCbrtMeshFactory.h"
#include "CFDMesh/Mesh/Quad2DCubeMeshFactory.h"
#include "CFDMesh/Mesh/Quad2DRotateMeshFactory.h"
#include "CFDMesh/Mesh/Quad2DTwistMeshFactory.h"
#include "CFDMesh/Mesh/PolarMeshFactory.h"
#include "CFDMesh/Mesh/Quad2DImageMeshFactory.h"
#include "CFDMesh/Mesh/Cube3DMeshFactory.h"
#include "CFDMesh/Mesh/Sphere3DMeshFactory.h"
#include "CFDMesh/Mesh/Cylinder3DMeshFactory.h"
#include "CFDMesh/Mesh/Torus3DMeshFactory.h"
#include "CFDMesh/Mesh/Vertex.h"
#include "CFDMesh/Mesh/Face.h"
#include "CFDMesh/Mesh/Cell.h"
#include "CFDMesh/Vector.h"
#include "CFDMesh/Util.h"
#include "CFDMesh/GUI.h"
#include "GLApp/gl.h"
#include "Common/Macros.h"
#include "Common/Exception.h"
#include <vector>
#include <list>
#include <string>
#include <cassert>


template<typename T> void glVertex2v(T const * v);
template<> void glVertex2v<float>(float const * v) { glVertex2fv(v); }
template<> void glVertex2v<double>(double const * v) { glVertex2dv(v); }

template<typename T> void glVertex3v(T const * v);
template<> void glVertex3v<float>(float const * v) { glVertex3fv(v); }
template<> void glVertex3v<double>(double const * v) { glVertex3dv(v); }

template<typename T> void glVertex3(T x, T y, T z);
template<> void glVertex3<float>(float x, float y, float z) { glVertex3f(x,y,z); }
template<> void glVertex3<double>(double x, double y, double z) { glVertex3d(x,y,z); }


namespace CFDMesh {
namespace Mesh {

//dim is the dimension of the manifold, not of the vectors (which are all 3D atm)
template<typename real, int dim, typename Cons>
struct Mesh {
	using real2 = Tensor::Vector<real, 2>;
	using real3 = Tensor::Vector<real, 3>;
		
	using Vertex = CFDMesh::Mesh::Vertex<real>;
	using Face = CFDMesh::Mesh::Face<real, Cons>;
	using Cell = CFDMesh::Mesh::Cell<real, Cons>;

protected:
	//https://stackoverflow.com/questions/8147027/how-do-i-call-stdmake-shared-on-a-class-with-only-protected-or-private-const
	struct ctorkey {
		explicit ctorkey(int) {}
	};

	Mesh(Mesh const &) = delete;
	Mesh const & operator=(Mesh const &) = delete;

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

	template<typename, int, typename> 
	friend struct MeshFactory;
	
	explicit Mesh(ctorkey const &) {}
	virtual ~Mesh() {}


	int addFaceForVtxs(int const * vs, int n) {
#if 0
		std::vector<int> uniquevs(vs, vs + n);
std::cout << "requesting to build a face from " << uniquevs << std::endl;
		if constexpr (dim == 3) {
			//remove any degenerate vertexes
			for (int i = uniquevs.size()-2; i >= 0; --i) {
std::cout << "comparing vtx " << vtxs[uniquevs[i]].pos << " and " << vtxs[uniquevs[i+1]].pos << std::endl;
				if (i < uniquevs.size()-1 && (vtxs[uniquevs[i]].pos - vtxs[uniquevs[i+1]].pos).lenLInf() < 1e-7) {
std::cout << " ... is less than threshold, removing" << std::endl;					
					uniquevs.erase(uniquevs.begin() + i);
					++i;
				}
			}
#if 0			
			if (uniquevs.size() < 3) {
std::cout << "removed too many vtxs, bailing on face" << std::endl;				
				return -1;
			}
#endif		
		}
//		n = vs.size();
//std::cout << "after pruning redundant vertexes, building a face from " << vs << std::endl;
#endif

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
				if (f->vtxCount == n) {
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
			real3 delta = a.pos - b.pos;
			f.area = delta.length();
			real3 normal = real3(delta(1), -delta(0)).unit();
			f.normal(0,0) = normal(0);
			f.normal(0,1) = normal(1);
			f.normal(0,2) = normal(2);
		} else if constexpr (dim == 3) {
			std::vector<real3> polyVtxs(n);
			for (int i = 0; i < n; ++i) {
				polyVtxs[i] = vtxs[vs[i]].pos;
			}
			real3 normal;
			for (int i = 0; i < n; ++i) {
				int i2 = (i+1)%n;
				int i3 = (i2+1)%n;
				normal += cross(polyVtxs[i2] - polyVtxs[i], polyVtxs[i3] - polyVtxs[i2]);
			}
			normal = normal.unit();
			f.normal(0,0) = normal(0);
			f.normal(0,1) = normal(1);
			f.normal(0,2) = normal(2);
#if 0
std::cout << "building face with vertexes " << polyVtxs << std::endl;
std::cout << "assigning normal to " 
	<< (polyVtxs[1] - polyVtxs[0]) << " x " 
	<< (polyVtxs[2] - polyVtxs[1]) << " = "
	<< normal << std::endl;
#endif
			f.area = polygon3DVolume(polyVtxs, normal);
			f.pos = polygon3DCOM(polyVtxs, f.area, normal);
#if 0
std::cout << "getting area " << f.area << " pos " << f.pos << std::endl;
#endif
		} else {
			throw Common::Exception() << "here";
		}
		return fi;
	}

	
	int addFace(int const * vs, int n, int ci) {
		int fi = addFaceForVtxs(vs, n);
//		if (fi != -1) return fi;
		Face& f = faces[fi];	
		if (f.cells(0) == -1) {
			f.cells(0) = ci;
		} else if (f.cells(1) == -1) {
			f.cells(1) = ci;
		} else {
			throw Common::Exception() << "tried to add too many cells to an edge";
		}
#if 0
std::cout << "adding face " << fi << " as " << std::to_string(f) << " with vtx indexes " << std::vector<int>(vs, vs+n) << std::endl;
if (f.cells(0) == -1 && f.cells(1) == -1) throw Common::Exception() << "here " << __FILE__ << ":" << __LINE__;
#endif
		return fi;
	}


	void addCell(std::vector<int> vis) {
#if 0
std::cout << "adding cell " << vis << std::endl;		
#endif
		int ci = (int)cells.size();
		cells.push_back(Cell());
		Cell& c = cells.back();
		
		size_t n = vis.size();
	
		//c.vtxs = vis;
		c.vtxOffset = (int)cellVtxIndexes.size();
		cellVtxIndexes.insert(cellVtxIndexes.end(), vis.begin(), vis.end());
		c.vtxCount = (int)cellVtxIndexes.size() - c.vtxOffset;

		int lastFaceSize = (int)faces.size();
		c.faceOffset = (int)cellFaceIndexes.size();
		if constexpr (dim == 2) {
			//face is a 1-form
			for (size_t i = 0; i < n; ++i) {
				int vtxs[2] = {vis[i], vis[(i+1)%n]}; 
				int fi = addFace(vtxs, numberof(vtxs), ci);
//				if (fi != -1) 
				{
					cellFaceIndexes.push_back(fi);
				}
			}
			std::vector<real2> polyVtxs = map<
				decltype(vis),
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
			//identity cube sides
			std::vector<std::vector<int>> identityCubeSides = {
				{0,4,6,2},	//x-
				{1,3,7,5},	//x+
				{0,1,5,4},	//y-
				{2,6,7,3},	//y+
				{0,2,3,1},	//z-
				{4,5,7,6},	//z+
			};

			//vector of per-face vector-of-vertexes
			//to pass to the polyhedron functions
			std::vector<std::vector<real3>> cubeVtxs;

			for (auto const & side : identityCubeSides) {
				
				std::vector<int> thisFaceVtxIndexes = map<
					std::vector<int>,
					std::vector<int>
				>(side, [&vis](int side_i) -> int {
					return vis[side_i];
				});
				
				int fi = addFace(thisFaceVtxIndexes.data(), thisFaceVtxIndexes.size(), ci);
//				if (fi != -1) 
				{
					Face const & f = faces[fi];
					
					//if face area is zero then don't add it to cell's faces
					// and don't use it later for determining cell's volume
					if (f.area <= 1e-7) {
#if 0
std::cout << "disregarding face with area " << f.area << std::endl;
#endif
					} else {
						cellFaceIndexes.push_back(fi);

						cubeVtxs.push_back(map<
							decltype(thisFaceVtxIndexes),
							std::vector<real3>
						>(
							thisFaceVtxIndexes,
							[this](int i) -> real3 { return vtxs[i].pos; }
						));
					}
				}
			}

			if (cellFaceIndexes.size() - lastFaceSize < 4) {	//not enough sides to even form a tetrahedron
#if 0
std::cout << "cell didn't have enough faces, setting volume to zero" << std::endl;				
#endif				
				c.volume = 0;
			} else {
				c.volume = polyhedronVolume(cubeVtxs);
				c.pos = polyhedronCOM(cubeVtxs, c.volume);
			}
		} else {
			throw Common::Exception() << "you are here";
		}
		
		c.faceCount = (int)cellFaceIndexes.size() - c.faceOffset;

#if 0
std::cout << "cell before testing validity " << c << std::endl;
#endif

#if 0
		//remove the cell
		if (c.volume <= 1e-7) {
std::cout << "cell had zero volume -- removing." << std::endl;			
			for (int i = 0; i < c.faceCount; ++i) {
				int fi = cellFaceIndexes[c.faceOffset + i];
				if (fi < lastFaceSize) {
					auto& f = faces[fi];
int2 before = f.cells;
					f.removeCell(ci);
if (f.cells(0) == -1 && f.cells(1) == -1) {
	throw Common::Exception() << "before " << before << " after = (-1, -1), here " << __FILE__ << ":" << __LINE__;
}
				}
			}
			faces.resize(lastFaceSize);
			cellFaceIndexes.resize(c.faceOffset);
			cellVtxIndexes.resize(c.vtxOffset);
			cells.resize(ci);
		}
#endif

	}

	static std::pair<real3, real3> getPerpendicularBasis(real3 n) {
		real3 n_x_x = cross(n, real3(1,0,0));
		real3 n_x_y = cross(n, real3(0,1,0));
		real3 n_x_z = cross(n, real3(0,0,1));
		real n_x_xSq = n_x_x.lenSq();
		real n_x_ySq = n_x_y.lenSq();
		real n_x_zSq = n_x_z.lenSq();
		real3 n2;
		if (n_x_xSq > n_x_ySq) {
			if (n_x_xSq > n_x_zSq) {
				n2 = n_x_x;	//use x
			} else {
				n2 = n_x_z;	//use z
			}
		} else {
			if (n_x_ySq > n_x_zSq) {
				n2 = n_x_y;	//use y
			} else {
				n2 = n_x_z;	//use z
			}
		}
		n2 = n2.unit();
		real3 n3 = cross(n, n2);
		return std::make_pair(n2, n3);
	}

	//calculate edge info
	//calculate cell volume info
	void calcAux() {

#if 0
		//while we're here ...
		//cycle through all vertexes, collapse mathcing vertexes, build a map from src to collapsed
		// then apply that map to the cellVtxIndexes and faceVtxIndexes
		{
			size_t oldsize = vtxs.size();
			
			//collapseVtxs[from] = to
			std::map<int, int> collapseVtxs;
			for (int j = (int)vtxs.size() - 1; j > 1; --j) {
				//as long as we cycle first-to-last, our collapseVtxs map old entries will not need to be updated as we add new entries into it
				for (int i = 0; i < j; ++i) {
					if ((vtxs[i].pos - vtxs[j].pos).lenLInf() < 1e-7) {
						collapseVtxs[j] = i;			//replace all 'j's with 'i'
						vtxs.erase(vtxs.begin() + j);	//erase 'j'
						break;
					}
				}
			}

			size_t newsize = vtxs.size();
			std::cout << "collapsing #vtxs from " << oldsize << " to " << newsize << std::endl;

			for (auto v : std::vector<std::vector<int>*>{ 
				&cellVtxIndexes,
				&faceVtxIndexes
			}) {
				for (auto& i : *v) {
					auto iter = collapseVtxs.find(i);
					if (iter != collapseVtxs.end()) {
						i = iter->second;
					}
				}
			}
		}
#endif

#if 0
		std::cout << "vtxs" << std::endl;
		for (auto const & v : vtxs) {
			std::cout << v << std::endl;
		}
		std::cout << "faces" << std::endl;
		for (auto const & f : faces) {
			std::cout << f << std::endl;
			for (int i = 0; i < f.vtxCount; ++i) {
				std::cout << " local vtx[" << i << "] = global vtx[" << faceVtxIndexes[i + f.vtxOffset] << "] = " << vtxs[faceVtxIndexes[i + f.vtxOffset]] << std::endl;;
			}
		}
		std::cout << "faceVtxIndexes " << faceVtxIndexes << std::endl;
		std::cout << "cells" << std::endl;
		for (auto const & c : cells) {
			std::cout << c << std::endl;
		}
		std::cout << "cellVtxIndexes " << cellVtxIndexes << std::endl;
		std::cout << "cellFaceIndexes " << cellFaceIndexes << std::endl;
#endif
		
		for (auto& f : faces) {
			int a = f.cells(0);
			int b = f.cells(1);
			if (a != -1 && b != -1) {
				if (real3::dot(cells[a].pos - cells[b].pos, real3(f.normal(0,0), f.normal(0,1), f.normal(0,2))) < 0) {
#if 0
std::cout << "swapping cells so normal points to b" << std::endl;					
#endif					
					std::swap(a, b);
					f.cells(0) = a;
					f.cells(1) = b;
					f.normal(0,0) *= -1;
					f.normal(0,1) *= -1;
					f.normal(0,2) *= -1;
				}
				//distance between cell centers
				//f.cellDist = (cells[b].pos - cells[a].pos).length();
				//distance projected to edge normal
				f.cellDist = fabs(real3::dot(real3(f.normal(0,0), f.normal(0,1), f.normal(0,2)), cells[b].pos - cells[a].pos));
			} else if (a != -1) {
				f.cellDist = (cells[a].pos - f.pos).length() * 2.;
			} else {
				throw Common::Exception() << "looks like you created a face that isn't touching any cells...";
			}
	
#if 0	//DEBUG	//if we allow degenerate cells then we can have non-positive cell distances
if (f.cellDist <= 1e-7) throw Common::Exception() << "got non-positive cell distance " << std::to_string(f);
#endif
		}

#if 0
		//validity check
		for (int ci = 0; ci < (int)cells.size(); ++ci) {
			auto& c = cells[ci];
			std::cout << "cell " << ci << " has faces " << std::vector(cellFaceIndexes.begin() + c.faceOffset, cellFaceIndexes.begin() + c.faceOffset + c.faceCount) << std::endl;
		}
		for (int fi = 0; fi < (int)faces.size(); ++fi) {
			auto& f = faces[fi];
			std::cout << "face " << fi << " has cells " << f.cells << std::endl;
		}

		for (int ci = 0; ci < (int)cells.size(); ++ci) {
			auto& c = cells[ci];
			for (int i = 0; i < c.faceCount; ++i) {
				int fi = cellFaceIndexes[c.faceOffset + i];
				auto& f = faces[fi];

				bool found = false;
				for (int j = 0; j < f.cells.size; ++j) {
					if (f.cells(j) == ci) {
						found = true;
						break;
					}
				}
				if (!found) {
					throw Common::Exception() << "cell " << ci << " attached to face " << fi << " that doesn't have that cell\n"
						<< "face=" << f << "\n"
						<< "cell=" << c;
				}
			}
		}
#endif
	
		//build basis from normal
		for (auto& f : faces) {
			auto [n2, n3] = getPerpendicularBasis(real3(f.normal(0,0), f.normal(0,1), f.normal(0,2)));
			for (int i = 0; i < 3; ++i) {
				f.normal(1,i) = n2(i);
				f.normal(2,i) = n3(i);
			}
		}
	}

	struct DrawArgs {
		using This = DrawArgs;
	
		using ValueRange = std::pair<float, float>; //min, max

		int gradientTex = 0;
		ValueRange displayValueRange = {};
		int selectedCellIndex = -1;
		bool showCells = true;
		bool showVtxs = false;
		bool showFaces = false;
		bool showFaceNormals = false;
		bool showFaceCenters = false;
		bool showCellCenters = false;
		float cellScale = 1;

		static constexpr auto fields = std::make_tuple(
			std::make_pair("display value range", &This::displayValueRange),
			GUISeparator(),
			
			std::make_pair("show cells", &This::showCells),
			std::make_pair("show vertexes", &This::showVtxs),
			std::make_pair("show faces", &This::showFaces),
			std::make_pair("show face normals", &This::showFaceNormals),
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
				auto const & v = vtxs[cellVtxIndexes[vi + c.vtxOffset]].pos;
				glVertex3v(((v - c.pos) * args.cellScale + c.pos).v);
			}
			glEnd();
			glColor3f(1,1,1);
			glLineWidth(1);
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		}

		
		if (args.showVtxs) {
			glPointSize(1);
			glColor3f(1,1,1);
			glEnable(GL_TEXTURE_1D);
			glBindTexture(GL_TEXTURE_1D, args.gradientTex);
			glBegin(GL_POINTS);
			for (size_t i = 0; i < vtxs.size(); ++i) {
				float f = ((float)i + .5) / (float)vtxs.size();
				glTexCoord1f(f);
				glVertex3v(vtxs[i].pos.v);
			}
			glEnd();
			glBindTexture(GL_TEXTURE_1D, 0);
			glDisable(GL_TEXTURE_1D);
		
			glPointSize(3);
			glColor3f(0,0,0);
			glBegin(GL_POINTS);
			for (size_t i = 0; i < vtxs.size(); ++i) {
				float f = ((float)i + .5) / (float)vtxs.size();
				glTexCoord1f(f);
				glVertex3v(vtxs[i].pos.v);
			}
			glEnd();
		}
		
		glPointSize(3);
		glBegin(GL_POINTS);
		if (args.showFaceCenters) {
			glColor3f(1,0,1);
			for (auto const & f : faces) {
				glVertex3v(f.pos.v);
			}
		}
		if (args.showCellCenters) {
			glColor3f(0,1,1);
			for (auto const & c : cells) {
				glVertex3v(c.pos.v);
			}
		}
		glEnd();
		
		glColor3f(1,1,1);
		glPointSize(1);
		
		if (args.showFaces) {
			for (auto const & f : faces) {
				glBegin(GL_LINE_LOOP);
				for (int vi = 0; vi < f.vtxCount; ++vi) {
					auto const & v = vtxs[faceVtxIndexes[vi + f.vtxOffset]].pos;
					glVertex3v(((v - f.pos) * args.cellScale + f.pos).v);
				}
				glEnd();
			}
		}
		
		if (args.showFaceNormals) {
			glColor3f(1,0,1);
			glBegin(GL_LINES);
			for (auto const & f : faces) {
				for (int i = 0; i < 3; ++i) {
					glVertex3v(f.pos.v);
					glVertex3(
						(f.pos(0) + f.normal(i,0) * f.area * .25),
						(f.pos(1) + f.normal(i,1) * f.area * .25),
						(f.pos(2) + f.normal(i,2) * f.area * .25));
				}
			}	
			glEnd();
		}


		glEnable(GL_TEXTURE_1D);
		glBindTexture(GL_TEXTURE_1D, args.gradientTex);

		if (args.showCells) {
			for (auto const & c : cells) {
				real f = (c.displayValue - args.displayValueRange.first) / (args.displayValueRange.second - args.displayValueRange.first);
				glColor3f(1,1,1);
				glTexCoord1f(f);
				
				if constexpr (dim == 2) {
					glBegin(GL_POLYGON);
					for (int vi = 0; vi < c.vtxCount; ++vi) {
						auto const & v = vtxs[cellVtxIndexes[vi + c.vtxOffset]].pos;
						glVertex3v(((v - c.pos) * args.cellScale + c.pos).v);
					}
					glEnd();
				} else if constexpr (dim == 3) {
					for (int i = 0; i < c.faceCount; ++i) {
						Face& f = faces[cellFaceIndexes[i + c.faceOffset]];
						glBegin(GL_POLYGON);
						for (int vi = 0; vi < f.vtxCount; ++vi) {
							auto const & v = vtxs[faceVtxIndexes[vi + f.vtxOffset]].pos;
							glVertex3v(((v - c.pos) * args.cellScale + c.pos).v);
						}
						glEnd();
					}
				}
			}
		}

		glBindTexture(GL_TEXTURE_1D, 0);
		glDisable(GL_TEXTURE_1D);
		
	}


	//	3D - polygon


	static real polygon3DVolume(std::vector<real3> const & vs, real3 normal) {
		size_t n = vs.size();
		real volume = 0;
		for (size_t i = 0; i < n; ++i) {
			real3 const & a = vs[i];
			real3 const & b = vs[(i+1)%n];
			volume += parallelepipedVolume(a, b, normal);
		}
		return .5 * volume;
	}

	static real3 polygon3DCOM(std::vector<real3> const & vs, real area, real3 normal) {
		int n = (int)vs.size();
#if 0
		for (int i = 0; i < n; ++i) {
			real3 const & a = vs[i];
			real3 const & b = vs[(i+1)%n];
			com += (a + b) * parallelepipedVolume(a, b, normal);
		}
		return com * (1. / (6. * area));
#else
		if (area == 0) {
			if (n == 0) {
				throw Common::Exception() << "you can't get a COM without any vertexes";
			}
			return std::accumulate(vs.begin()+1, vs.end(), vs.front()) / (real)n;
		}
		real3 com;
		real3 const & a = vs[0];
		for (int i = 2; i < n; ++i) {
			real3 const & b = vs[i-1];
			real3 const & c = vs[i];
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


	static real polygonVolume(std::vector<real2> const & vs) {
		size_t n = vs.size();
		real volume = 0;
		for (size_t i = 0; i < n; ++i) {
			real2 const & a = vs[i];
			real2 const & b = vs[(i+1)%n];
			volume += parallelogramVolume(a, b);
		}
		return .5 * volume;
	}

	static real2 polygonCOM(std::vector<real2> const & vs, real volume) {
		size_t n = vs.size();
		if (volume == 0) {
			if (n == 0) {
				throw Common::Exception() << "you can't get a COM without any vertexes";
			}
			return std::accumulate(vs.begin()+1, vs.end(), vs.front()) / (real)n;
		}
		real2 com;
		for (int i = 0; i < (int)n; ++i) {
			real2 const & a = vs[i];
			real2 const & b = vs[(i+1)%n];
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


	static real polyhedronVolume(std::vector<std::vector<real3>> const & faces) {
		real volume = 0;
		for (auto const & face : faces) {
			for (int i = 2; i < (int)face.size(); ++i) {
				//tri from 0, i-1, 1
				real3 const & a = face[0];
				real3 const & b = face[i-1];
				real3 const & c = face[i];
				volume += parallelepipedVolume(a, b, c);
			}
		}
		//volume of n-sided pyramid in nD is volume of parallelogram divided by 1/n!
		return volume / 6.;
	}

	static real3 polyhedronCOM(std::vector<std::vector<real3>> const & faces, real volume) {
		if (volume == 0) {
			if (faces.size() == 0) {
				throw Common::Exception() << "you can't get a COM without any vertexes";
			}
			real3 sum;
			real total = {};
			for (auto const & face : faces) {
				for (auto const & vtx : face) {
					sum += vtx;
					++total;
				}
			}
			return sum * (1. / (real)total);
		}
		real3 com;
		for (auto const & face : faces) {
			for (int i = 2; i < (int)face.size(); ++i) {
				//tri from 0, i-1, 1
				real3 const & a = face[0];
				real3 const & b = face[i-1];
				real3 const & c = face[i];
				com += ((a + b) * (a + b) + (b + c) * (b + c) + (c + a) * (c + a)) * cross(b - a, c - a) / 48.;
			}
		}
		return com / volume;
	}


	static std::vector<std::shared_ptr<MeshFactory<real, dim, Cons>>> getGens() {
		if constexpr (dim == 2) {
			return std::vector<std::shared_ptr<MeshFactory<real, dim, Cons>>>{
				std::make_shared<PolarMeshFactory<real, dim, Cons>>(),
				std::make_shared<Quad2DMeshFactory<real, dim, Cons>>(),
				std::make_shared<Quad2DCbrtMeshFactory<real, dim, Cons>>(),
				std::make_shared<Quad2DCubeMeshFactory<real, dim, Cons>>(),
				std::make_shared<Quad2DRotateMeshFactory<real, dim, Cons>>(),
				std::make_shared<Quad2DTwistMeshFactory<real, dim, Cons>>(),
				std::make_shared<Quad2DImageMeshFactory<real, dim, Cons>>(),
				std::make_shared<P2DFMTMeshFactory<real, dim, Cons>>(),
			};
		} else if constexpr (dim == 3) {
			return std::vector<std::shared_ptr<MeshFactory<real, dim, Cons>>>{
				std::make_shared<Cube3DMeshFactory<real, dim, Cons>>(),
				std::make_shared<Cylinder3DMeshFactory<real, dim, Cons>>(),
				std::make_shared<Sphere3DMeshFactory<real, dim, Cons>>(),
				std::make_shared<Torus3DMeshFactory<real, dim, Cons>>(),
			};
		}
		throw Common::Exception() << "here";
	}
};

}
}
