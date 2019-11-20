#pragma once

#include "CFDMesh/Vector.h"
#include "CFDMesh/Util.h"
#include "GLApp/gl.h"
#include "Common/File.h"
#include "Common/Exception.h"
#include <vector>
#include <list>
#include <string>
#include <cassert>


template<typename T> void glVertex2v(const T* v);
template<> void glVertex2v<float>(const float* v) { glVertex2fv(v); }
template<> void glVertex2v<double>(const double* v) { glVertex2dv(v); }


namespace CFDMesh {

template<typename real, typename Cons>
struct MeshNamespace {
using real2 = Tensor::Vector<real, 2>;
using real3 = Tensor::Vector<real, 3>;


struct Vertex;
struct Face;
struct Cell;

//zero-forms
struct Vertex {	//not required by finite volume algorithm
	real3 pos;
	
	//keeping track of Vertex::faces isn't used by the renderer, mesh generation, or finite-volume integration
	//std::vector<int> faces;
	
	Vertex() {}
	Vertex(real3 pos_) : pos(pos_) {}
};

//(n-1)-forms
struct Face {
	real3 pos;
	real3 delta;
	real3 normal;
	real length;
	real cellDist;
	Cons flux;
	
	int cells[2];	//there are always only 2 n-forms on either side of a (n-1)-form
	
	//the vertexes of a (n-1)-form are not required by finite volume algorithm
	//they are required for determining the interface area which is then used by the finite volume algorithm
	//also used by the renderer
	//however, note, in n>2 dimensions, a (n-1)-form will be defined by an arbitrary number of vertices (more than just two)
	int vtxs[2];
	
	Face(int va, int vb) : length(0), cellDist(0) {
		vtxs[0] = va;
		vtxs[1] = vb;
		cells[0] = -1;
		cells[1] = -1;
	}

	bool removeCell(int cellIndex) {
		if (cells[1] == cellIndex) cells[1] = -1;
		if (cells[0] == cellIndex) {
			cells[0] = cells[1];
			cells[1] = -1;
		}
		return cells[0] == -1 && cells[1] == -1;
	}
};

//n-forms
struct Cell {
	real3 pos;
	real volume;
	Cons U;
	
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

	Cell() : volume(0) {}
};

struct MeshFactory;
struct Mesh {
//liu kang wins.  friendship.
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


	friend struct MeshFactory;
	
	explicit Mesh(const ctorkey&) {}
	virtual ~Mesh() {}


	int addEdge(int va, int vb) {
		int ei = 0;
		for (; ei < (int)faces.size(); ++ei) {
			Face* e = &faces[ei];
			if ((e->vtxs[0] == va && e->vtxs[1] == vb) ||
				(e->vtxs[0] == vb && e->vtxs[1] == va)) 
			{
				return ei;
			}
		}
		assert(ei == (int)faces.size());
		faces.push_back(Face(va, vb));
		//vtxs[va].faces.push_back(ei);
		//vtxs[vb].faces.push_back(ei);
		return ei;
	}

	void addCell(std::vector<int> vis) {
		Cell c;
		
		size_t n = vis.size();
	
		//c.vtxs = vis;
		c.vtxOffset = (int)cellVtxIndexes.size();
		cellVtxIndexes.insert(cellVtxIndexes.end(), vis.begin(), vis.end());
		c.vtxCount = (int)cellVtxIndexes.size() - c.vtxOffset;
	
		c.faceOffset = (int)cellFaceIndexes.size();
		for (size_t i = 0; i < n; ++i) {
			cellFaceIndexes.push_back(addEdge(vis[i], vis[(i+1)%n]));
		}
		c.faceCount = (int)cellFaceIndexes.size() - c.faceOffset;

		c.volume = polyVol(map<std::vector<int>, std::vector<real3>>(vis, [this](int vi) -> real3 { return vtxs[vi].pos; }));
		assert(c.volume > 0);

#if 1	//vertex average
		//TODO use COM
		c.pos = sum(map<
			std::vector<int>,
			std::vector<real3>
		>(vis, [this](int vi) {
			return vtxs[vi].pos;
		})) * (1. / (real)n);
#else	//COM
		for (int j = 0; j < 3; ++j) {
			for (int i = 0; i < n; ++i) {
				real3 x1 = vtxs[vis[i]].pos;
				real3 x2 = vtxs[vis[(i+1)%n]].pos;
				c.pos(j) += (x1(j) + x2(j)) / (x1(0) * x2(1) - x1(1) * x2(0));
			}
			c.pos(j) /= (6 * c.volume);
		}
#endif
		int ci = (int)cells.size();
		cells.push_back(c);
	
		//for (int ei : c.faces) {
		for (int ei = 0; ei < c.faceCount; ++ei) {
			Face* e = &faces[cellFaceIndexes[ei+c.faceOffset]];
			if (e->cells[0] == -1) {
				e->cells[0] = ci;
			} else if (e->cells[1] == -1) {
				e->cells[1] = ci;
			} else {
				throw Common::Exception() << "tried to add too many cells to an edge";
			}
		}
	}

	//calculate edge info
	//calculate cell volume info
	void calcAux() {		
		for (auto& e : faces) {
			{
				auto& a = vtxs[e.vtxs[0]];
				auto& b = vtxs[e.vtxs[1]];
				e.pos = (a.pos + b.pos) * .5;
				e.delta = a.pos - b.pos;
				e.length = real3::length(e.delta);
				e.normal = real3(e.delta(1), -e.delta(0));
				e.normal *= 1. / real3::length(e.normal);
			}

			{
				int a = e.cells[0];
				int b = e.cells[1];
				if (a != -1 && b != -1) {
					if (real3::dot(cells[a].pos - cells[b].pos, e.normal) < 0) {
						std::swap(a, b);
						e.cells[0] = a;
						e.cells[1] = b;
						e.normal *= -1;
					}
					//distance between cell centers
					//e.cellDist = real3::length(cells[b].pos - cells[a].pos);
					//distance projected to edge normal
					e.cellDist = fabs(real3::dot(e.normal, cells[b].pos - cells[a].pos));
				} else if (a != -1) {
					e.cellDist = real3::length(cells[a].pos - e.pos) * 2.;
				} else {
					throw Common::Exception() << "you are here";
				}
			}
		}
	}

	void draw(
		int gradientTex,
		float displayValueMin,
		float displayValueMax,
		int selectedCellIndex,
		bool showVtxs,
		bool showEdges,
		bool showCellCenters
	) {
		glClear(GL_COLOR_BUFFER_BIT);

		glEnable(GL_TEXTURE_1D);
		glBindTexture(GL_TEXTURE_1D, gradientTex);

		for (const auto& c : cells) {
			real f = (c.displayValue - displayValueMin) / (displayValueMax - displayValueMin);
			glTexCoord1f(f);
			glBegin(GL_POLYGON);
			for (int vi = 0; vi < c.vtxCount; ++vi) {
				glVertex2v(vtxs[cellVtxIndexes[vi + c.vtxOffset]].pos.v);
			}
			glEnd();
		}
		
		glBindTexture(GL_TEXTURE_1D, 0);
		glDisable(GL_TEXTURE_1D);
		
		if (selectedCellIndex != -1) {
			Cell* selectedCell = &cells[selectedCellIndex];
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
			glLineWidth(3);
			glColor3f(1,1,0);
			glBegin(GL_POLYGON);
			for (int vi = 0; vi < selectedCell->vtxCount; ++vi) {
				glVertex2v(vtxs[cellVtxIndexes[vi + selectedCell->vtxOffset]].pos.v);
			}
			glEnd();
			glLineWidth(1);
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		}

		if (showVtxs || showCellCenters) {
			glColor3f(1,1,1);
			glPointSize(3);
			glBegin(GL_POINTS);
			if (showVtxs) {
				for (const auto& v : vtxs) {
					glVertex2v(v.pos.v);
				}
			}
			if (showCellCenters) {
				for (const auto& c : cells) {
					glVertex2v(c.pos.v);
				}
			}
			glEnd();
			glPointSize(1);
		}
		if (showEdges) {
			glColor3f(1,1,1);
			glBegin(GL_LINES);
			for (const auto& e : faces) {
				for (int vi : e.vtxs) {
					glVertex2v(vtxs[vi].pos.v);
				}
			}
			glEnd();
		}

	}
};


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
	

//2D polygon volume
static real polyVol(const std::vector<real3>& vs) {
	size_t n = vs.size();
	real v = 0;
	for (size_t i = 0; i < n; ++i) {
		real3 pi = vs[i];
		real3 pj = vs[(i+1)%n];
		v += pi(0) * pj(1) - pi(1) * pj(0);
	}
	return .5 * v;
}

template<typename I, typename F>
static bool contains(
	real3 pos,
	I begin,
	I end,
	F f
) {
	I i = begin;
	using V = decltype(f(*i));
	V first = f(*i);
	V prev = first;

	auto check = [pos](V prev, V cur) -> bool {
		real3 pi = prev - pos;
		real3 pj = cur - pos;
		real vz = pi(0) * pj(1) - pi(1) * pj(0);
		return vz < 0;
	};
	if (i != end) {
		for (++i; i != end; ++i) {
			V cur = f(*i);
			if (!check(prev, cur)) return false;
			prev = cur;
		}
	}
	check(prev, first);
	return true;

}

};

}
