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

using Point = std::tuple<float, float, float>;

template <typename T>
using Vector = std::vector<T>;

template <typename T>
using Vector2D = Vector<Vector<T>>;

struct ConfigurationalEntropy {
	ConfigurationalEntropy(const Graph & completeGraph, int nDimensions, int nParticles, bool pbcX, bool pbcY, bool pbcZ, const py::array_t<double> & dMin, const py::array_t<double> & dMax) :
		_completeGraph(completeGraph), _searchTree(nDimensions, nParticles, pbcX, pbcY, pbcZ) {
			auto _dMin = dMin.unchecked<1>();
			auto _dMax = dMax.unchecked<1>();
			_generator = std::mt19937(_randDevice());
			_distrX = std::uniform_real_distribution<float>(_dMin(0), _dMax(0));
			_distrY = std::uniform_real_distribution<float>(_dMin(1), _dMax(1));
			_distrZ = std::uniform_real_distribution<float>(_dMin(2), _dMax(2));
	}

	py::tuple check_isomorfism(Vector<Graph> & graphs, double c, int m, int n);

	py::tuple calculate(int m, int n, double c);

	const Point generate_random_point(int precision);

	const Vector<Graph> get_subgraphs(int m, int n);

	Graph generate_subgraph(const Vector<int> & closestNeighbors);

	inline void add_positions(const py::array_t<double> & positions) {
		_searchTree.add_positions(positions);
	}

	inline void init_search(double xMin, double xMax, double yMin, double yMax, double zMin, double zMax) {
		_searchTree.init_search(xMin, xMax, yMin, yMax, zMin, zMax);
	}

	inline Vector<int> search_nearest_neighbors(float x, float y, float z, unsigned int n) {
		return _searchTree.search_nearest_neighbors(x, y, z, n);
	}

	double calc_shannon_entropy(double Hn, double fi, double m);
private:
	const Graph & _completeGraph;
	SearchTree _searchTree;
	std::random_device _randDevice;
	std::mt19937 _generator;
	std::uniform_real_distribution<float> _distrX;
	std::uniform_real_distribution<float> _distrY;
	std::uniform_real_distribution<float> _distrZ;
};

#endif /* GRAPHS_HPP */
