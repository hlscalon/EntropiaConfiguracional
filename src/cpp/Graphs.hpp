#ifndef GRAPHS_HPP
#define GRAPHS_HPP

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <vector>
#include <random>
#include <tuple>

#include "Graph.hpp"
#include "SearchTree.hpp"

namespace py = pybind11;

using Point = std::tuple<double, double, double>;

template <typename T>
using Vector = std::vector<T>;

template <typename T>
using Vector2D = Vector<Vector<T>>;

struct Graphs {
	Graphs(int m, int n, const Graph & completeGraph, int nDimensions, int nParticles, bool pbcX, bool pbcY, bool pbcZ, const py::array_t<double> & dMin, const py::array_t<double> & dMax) :
		_m(m), _n(n), _graphs(Vector<Graph>(m)), _completeGraph(completeGraph),
		_searchTree(nDimensions, nParticles, pbcX, pbcY, pbcZ) {
			auto _dMin = dMin.unchecked<1>();
			auto _dMax = dMax.unchecked<1>();
			_generator = std::mt19937(_randDevice());
			_distrX = std::uniform_real_distribution<>(_dMin(0), _dMax(0));
			_distrY = std::uniform_real_distribution<>(_dMin(1), _dMax(1));
			_distrZ = std::uniform_real_distribution<>(_dMin(2), _dMax(2));
		}

	inline void insert(int pos, const Graph & _graph) { _graphs[pos] = _graph; }
	py::tuple check_isomorfism(double c);
	py::tuple calculate_configurational_entropy(double c);
	void generate_subgraphs(const Vector2D<int> & closestNeighbors);
	const Point generate_random_point();
	const Vector2D<int> get_closest_neighbors();
	inline void add_positions(const py::array_t<double> & positions) { _searchTree.add_positions(positions); }
	inline void init_search(double xMin, double xMax, double yMin, double yMax, double zMin, double zMax) { _searchTree.init_search(xMin, xMax, yMin, yMax, zMin, zMax); }
	inline Vector<int> search_nearest_neighbors(double x, double y, double z, unsigned int n) { return _searchTree.search_nearest_neighbors(x, y, z, n); }
private:
	int _m;
	int _n;
	Vector<Graph> _graphs;
	const Graph & _completeGraph;
	SearchTree _searchTree;
	std::random_device _randDevice;
	std::mt19937 _generator;
	std::uniform_real_distribution<> _distrX;
	std::uniform_real_distribution<> _distrY;
	std::uniform_real_distribution<> _distrZ;
};

#endif /* GRAPHS_HPP */
