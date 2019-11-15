#pragma once

#include "CFDMesh/Vector.h"
#include "CFDMesh/Util.h"
#include "Common/File.h"
#include "Common/Exception.h"
#include <vector>
#include <list>
#include <string>
#include <cassert>

namespace CFDMesh {

template<typename MeshConfig>
struct MeshNamespace {
using real = typename MeshConfig::real;
using real2 = typename MeshConfig::real2;
using vec = typename MeshConfig::vec;
using StateVec = typename MeshConfig::StateVec;
using Cons = typename MeshConfig::Cons;


struct Vertex;
struct Face;
struct Cell;

//zero-forms
struct Vertex {	//not required by finite volume algorithm
	vec pos;
	
	//keeping track of Vertex::faces isn't used by the renderer, mesh generation, or finite-volume integration
	//std::vector<int> faces;
	
	Vertex() {}
	Vertex(vec pos_) : pos(pos_) {}
};

//(n-1)-forms
struct Face {
	vec pos;
	vec delta;
	vec normal;
	real length;
	real cellDist;
	StateVec flux;
	
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
	vec pos;
	real volume;
	Cons U;
	
	//required by the finite volume algorithm
	std::vector<int> faces;
	
	//not required by finite volume algorithm
	//however the cell volume is required, and is calculated using vtxs
	//also the renderer requires the vertexes
	std::vector<int> vtxs;
	
	Cell() : volume(0) {}
};

struct MeshArgs {
	int2 size_ = int2(101, 101);
	MeshArgs& size(int2 x) { size_ = x; return *this; }
	
	real2 mins_ = real2(-1, -1);
	MeshArgs& mins(real2 x) { mins_ = x; return *this; }

	real2 maxs_ = real2(1, 1);
	MeshArgs& maxs(real2 x) { maxs_ = x; return *this; }

	using GridFunc = std::function<real2(real2)>;
	GridFunc grid_ = [](real2 x) -> real2 { return x; };
	MeshArgs& grid(GridFunc x) { grid_ = x; return *this; }
	
	int2 repeat_ = int2(0, 0);
	MeshArgs& repeat(int2 x) { repeat_ = x; return *this; }

	int2 capmin_ = int2(0, 0);
	MeshArgs& capmin(int2 x) { capmin_ = x; return *this; }
};

struct Mesh {
	std::vector<Vertex> vtxs;	//0-forms, which construct n and n-1 forms
	std::vector<Face> faces;	//n-1-forms, hold flux information between cells
	std::vector<Cell> cells;	//n-forms, hold the info at each cell

//liu kang wins.  friendship.
private:
	Mesh() {}

public:

	static Mesh buildQuadChart(
		MeshArgs args = MeshArgs()
	) {
		int2 size = args.size_;
		real2 mins = args.mins_;
		real2 maxs = args.maxs_;
		auto grid = args.grid_;
		int2 repeat = args.repeat_;
		int2 capmin = args.capmin_;

		Mesh mesh;
		int m = size(0);
		int n = size(1);

		int vtxsize = m * n;
		if (capmin(0)) vtxsize++;
		mesh.vtxs.resize(vtxsize);
		for (int j = 0; j < n; ++j) {
			for (int i = 0; i < m; ++i) {
				real2 x = real2(
					((real)i + .5) / (real)size(0) * (maxs(0) - mins(0)) + mins(0),
					((real)j + .5) / (real)size(1) * (maxs(1) - mins(1)) + mins(1));
				
				real2 u = grid(x);
				std::function<real(int)> f = [&u](int i) -> real { return i < real2::size ? u(i) : 0.; };
				mesh.vtxs[i + m * j].pos = vec(f);
			}
		}
		
		int capindex = m * n;
		if (capmin(0)) {
			vec sum;
			for (int j = 0; j < n; ++j) {
				sum += mesh.vtxs[0 + m * j].pos;
			}
			mesh.vtxs[capindex].pos = sum / (real)n;
		}

		int imax = repeat(0) ? m : m-1;
		int jmax = repeat(1) ? n : n-1;
		for (int j = 0; j < jmax; ++j) {
			int jn = (j + 1) % n;
			for (int i = 0; i < imax; ++i) {
				int in = (i + 1) % m;
				mesh.addCell(std::vector<int>{i + m * j, in + m * j, in + m * jn, i + m * jn});
			}
		}

		if (capmin(0)) {
			for (int j = 0; j < jmax; ++j) {
				int jn = (j + 1) % n;
				mesh.addCell(std::vector<int>{ 0 + m * j, 0 + m * jn, capindex });
			}
		}


		mesh.calcAux();
		
		return mesh;
	}

	static Mesh buildTriChart(
		MeshArgs args = MeshArgs()
	) {
		int2 size = args.size_;
		real2 mins = args.mins_;
		real2 maxs = args.maxs_;
		auto grid = args.grid_;
		int2 repeat = args.repeat_;
		//int2 capmin = args.capmin_;
	
		Mesh mesh;
		int m = size(0);
		int n = size(1);
	
		mesh.vtxs.resize(m * n);
		for (int j = 0; j < n; ++j) {
			for (int i = 0; i < m; ++i) {
				real2 x = real2(
					((real)i + .5) / (real)size(0) * (maxs(0) - mins(0)) + mins(0),
					((real)j + .5) / (real)size(1) * (maxs(1) - mins(1)) + mins(1));
				
				real2 u = grid(x);
				std::function<real(int)> f = [&u](int i) -> real { return i < real2::size ? u(i) : 0.; };
				mesh.vtxs[i + m * j].pos = vec(f);
			}
		}
		
		int imax = repeat(0) ? m : m-1;
		int jmax = repeat(1) ? n : n-1;
		for (int j = 0; j < jmax; ++j) {
			int jn = (j + 1) % n;
			for (int i = 0; i < imax; ++i) {
				int in = (i + 1) % m;
				mesh.addCell(std::vector<int>{i + m * j, in + m * j, in + m * jn});
				mesh.addCell(std::vector<int>{in + m * jn, i + m * jn, i + m * j});
			}
		}
	
		mesh.calcAux();
		
		return mesh;
	}


	static Mesh buildFromFile(const std::string& fn) {
		Mesh mesh;	
		
		std::list<std::string> ls = split<std::list<std::string>>(Common::File::read(fn), "\n");
	
		std::string first = ls.front();
		ls.pop_front();
		std::vector<std::string> m_n = split<std::vector<std::string>>(ls.front(), "\\s+");
		ls.pop_front();
		int m = std::stoi(m_n[0]);
		int n = std::stoi(m_n[1]);
		std::list<std::string> _x = split<std::list<std::string>>(concat<std::list<std::string>>(ls, " "), "\\s+");
		if (_x.front() == "") _x.pop_front();
		if (_x.back() == "") _x.pop_back();
		std::function<real(const std::string&)> f = [](const std::string& s) -> real { return std::stod(s); };
		std::vector<real> x = map<std::list<std::string>, std::vector<real>>(_x, f);
		assert(x.size() == (size_t)(2 * m * n));
	
		auto us = std::vector(x.begin(), x.begin() + m*n);
		auto vs = std::vector(x.begin() + m*n, x.end());
		assert(us.size() == vs.size());

		mesh.vtxs.resize(m*n);
		for (int i = 0; i < (int)us.size(); ++i) {
			mesh.vtxs[i].pos = vec(us[i], vs[i]);
		}
	
		for (int j = 0; j < n-1; ++j) {
			for (int i = 0; i < m-1; ++i) {
				mesh.addCell(std::vector<int>{i + m * j, i + m * (j+1), i+1 + m * (j+1), i+1 + m * j});
			}
		}
	
		mesh.calcAux();
	
		return mesh;
	}

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
		c.vtxs = vis;
		size_t n = vis.size();
		for (size_t i = 0; i < n; ++i) {
			c.faces.push_back(addEdge(vis[i], vis[(i+1)%n]));
		}

		c.volume = polyVol(map<std::vector<int>, std::vector<vec>>(c.vtxs, [this](int vi) -> vec { return vtxs[vi].pos; }));
		assert(c.volume > 0);

#if 1	//vertex average
		//TODO use COM
		c.pos = sum(map<
			std::vector<int>,
			std::vector<vec>
		>(c.vtxs, [this](int vi) {
			return vtxs[vi].pos;
		})) * (1. / (real)n);
#else	//COM
		for (int j = 0; j < 3; ++j) {
			for (int i = 0; i < n; ++i) {
				vec x1 = vtxs[vis[i]].pos;
				vec x2 = vtxs[vis[(i+1)%n]].pos;
				c.pos(j) += (x1(j) + x2(j)) / (x1(0) * x2(1) - x1(1) * x2(0));
			}
			c.pos(j) /= (6 * c.volume);
		}
#endif
		int ci = (int)cells.size();
		cells.push_back(c);
	
		for (int ei : c.faces) {
			Face* e = &faces[ei];
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
				e.length = vec::length(e.delta);
				e.normal = vec(e.delta(1), -e.delta(0));
				e.normal *= 1. / vec::length(e.normal);
			}

			{
				int a = e.cells[0];
				int b = e.cells[1];
				if (a != -1 && b != -1) {
					if (vec::dot(cells[a].pos - cells[b].pos, e.normal) < 0) {
						std::swap(a, b);
						e.cells[0] = a;
						e.cells[1] = b;
						e.normal *= -1;
					}
					//distance between cell centers
					//e.cellDist = vec::length(cells[b].pos - cells[a].pos);
					//distance projected to edge normal
					e.cellDist = fabs(vec::dot(e.normal, cells[b].pos - cells[a].pos));
				} else if (a != -1) {
					e.cellDist = vec::length(cells[a].pos - e.pos) * 2.;
				} else {
					throw Common::Exception() << "you are here";
				}
			}
		}
	}
};

//2D polygon volume
static real polyVol(const std::vector<vec>& vs) {
	size_t n = vs.size();
	real v = 0;
	for (size_t i = 0; i < n; ++i) {
		vec pi = vs[i];
		vec pj = vs[(i+1)%n];
		v += pi(0) * pj(1) - pi(1) * pj(0);
	}
	return .5 * v;
}

static bool contains(const vec pos, std::vector<vec> vtxs) {
	size_t n = vtxs.size();
	for (size_t i = 0; i < n; ++i) {
		vec pi = vtxs[i] - pos;
		vec pj = vtxs[(i+1)%n] - pos;
		real vz = pi(0) * pj(1) - pi(1) * pj(0);
		if (vz > 0) return false;
	}
	return true;
}

};
}
