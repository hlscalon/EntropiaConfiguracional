#ifndef SEARCH_TREE_HPP
#define SEARCH_TREE_HPP

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <Aboria.h>

#include <vector>
#include <random>
#include <tuple>

namespace py = pybind11;

using Particle_t = Aboria::Particles<std::tuple<>, 3, std::vector, Aboria::Kdtree>;
using Particle_position = Particle_t::position;
using Point = std::tuple<double, double, double>;

template <typename T>
using Vector = std::vector<T>;

struct SearchTree {
	SearchTree(int nDimensions, int nParticles, bool pbcX, bool pbcY, bool pbcZ) :
		_nDimensions(nDimensions), _nParticles(nParticles), _pbcX(pbcX), _pbcY(pbcY), _pbcZ(pbcZ), _largestRadius(0), _particles(nParticles) {}

	void add_positions(const py::array_t<double> & positions);
	void init_search(double xMin, double xMax, double yMin, double yMax, double zMin, double zMax);
	Vector<int> search_nearest_neighbors(double x, double y, double z, unsigned int n);

private:
	int _nDimensions;
	int _nParticles;
	bool _pbcX;
	bool _pbcY;
	bool _pbcZ;
	int _largestRadius;
	Particle_t _particles;
};

#endif /* SEARCH_TREE_HPP */