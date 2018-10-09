#ifndef SEARCH_TREE_HPP
#define SEARCH_TREE_HPP

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <Aboria.h>

#include "Graph.hpp"

#include <vector>
#include <random>
#include <tuple>

namespace py = pybind11;

using Particle_t = Aboria::Particles<std::tuple<>, 3, std::vector, Aboria::Kdtree>;
using Particle_position = Particle_t::position;

template <typename T>
using Vector = std::vector<T>;

template <typename T>
using Vector2D = Vector<Vector<T>>;

template <typename T>
using PyArray = pybind11::detail::unchecked_reference<T, 1>;

struct SearchTree {
	SearchTree(int nDimensions, int nParticles, bool pbcX, bool pbcY, bool pbcZ) :
		_nDimensions(nDimensions), _nParticles(nParticles), _pbcX(pbcX), _pbcY(pbcY), _pbcZ(pbcZ), _largestRadius(0), _particles(nParticles) {}

	void add_positions(const py::array_t<double> & positions);

	void init_search(double xMin, double xMax, double yMin, double yMax, double zMin, double zMax, const PyArray<double> & dMin, const PyArray<double> & dMax);

	void generate_random_positions(int maxN, int maxM);

	Vector<int> get_nearest_neighbors(int randomPosition, int n);

	Vector<int> search_nearest_neighbors(float x, float y, float z, unsigned int n);

private:
	int _nDimensions;
	int _nParticles;
	bool _pbcX;
	bool _pbcY;
	bool _pbcZ;
	int _largestRadius;
	std::random_device _randDevice;
	std::mt19937 _generator;
	std::uniform_real_distribution<float> _distrX;
	std::uniform_real_distribution<float> _distrY;
	std::uniform_real_distribution<float> _distrZ;
	Particle_t _particles;
	Vector2D<int> _randomPositions;
};

#endif /* SEARCH_TREE_HPP */
