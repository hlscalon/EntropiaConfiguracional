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

template <typename T>
using Vector = std::vector<T>;

template <typename T>
using Vector2D = Vector<Vector<T>>;

template <typename T>
using PyArray = pybind11::detail::unchecked_reference<T, 1>;

struct ConfigurationalEntropy {
	ConfigurationalEntropy(const Graph & completeGraph, int nDimensions, int nParticles, bool pbcX, bool pbcY, bool pbcZ, int maxM) :
		_completeGraph(completeGraph), _searchTree(nDimensions, nParticles, pbcX, pbcY, pbcZ) {
			_generator = std::mt19937(_randDevice());
			_distr = std::uniform_int_distribution<unsigned int>(0, maxM - 1);
	}

	void check_isomorfism(Vector<Graph> & graphs, Graph & graph1, int & iso_label, Vector<int> & label_total);

	py::tuple calculate(int m, int n, double c);

	const int generate_random_point();

	const std::tuple<int, Vector<int>> generate_subgraphs(int m, int n);

	void generate_subgraph(Graph & graph, const Vector<int> & closestNeighbors);

	inline void add_positions(const py::array_t<double> & positions) {
		_searchTree.add_positions(positions);
	}

	void init_search(double xMin, double xMax, double yMin, double yMax, double zMin, double zMax, int maxN, int maxM, const py::array_t<double> & dMin, const py::array_t<double> & dMax);

	inline Vector<int> get_nearest_neighbors(int randomPosition, int n) {
		return _searchTree.get_nearest_neighbors(randomPosition, n);
	}

	double calc_shannon_entropy(double Hn, double fi, double m);
private:
	const Graph & _completeGraph;
	SearchTree _searchTree;
	std::random_device _randDevice;
	std::mt19937 _generator;
	std::uniform_int_distribution<unsigned int> _distr;
};

#endif /* GRAPHS_HPP */
