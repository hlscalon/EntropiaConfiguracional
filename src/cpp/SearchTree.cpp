#include "SearchTree.hpp"

#include <set>
#include <algorithm>
#include <cassert>

void SearchTree::add_positions(const py::array_t<double> & positions) {
	auto rows = positions.unchecked<2>();

	for (int i = 0; i < _nParticles; ++i) {
		//if (_nDimensions == 3) {}
		Aboria::get<Particle_position>(_particles)[i] = Aboria::vdouble3(rows(i, 0), rows(i, 1), rows(i, 2));
	}
}

void SearchTree::init_search(double xMin, double xMax, double yMin, double yMax, double zMin, double zMax, const PyArray<double> & dMin, const PyArray<double> & dMax) {
	Aboria::vdouble3 min = Aboria::vdouble3(xMin, yMin, zMin);
	Aboria::vdouble3 max = Aboria::vdouble3(xMax, yMax, zMax);
	Aboria::vbool3 periodic = Aboria::vbool3(_pbcX, _pbcY, _pbcZ);

	_largestRadius = std::ceil((min - max).norm());
	_particles.init_neighbour_search(min, max, periodic);
	_generator = std::mt19937(_randDevice());
	_distrX = std::uniform_real_distribution<float>(dMin(0), dMax(0));
	_distrY = std::uniform_real_distribution<float>(dMin(1), dMax(1));
	_distrZ = std::uniform_real_distribution<float>(dMin(2), dMax(2));
}

void SearchTree::generate_random_positions(int maxN, int maxM) {
	// inicializa com M linhas, e cada linha um vetor com N elementos (0)
	_randomPositions = Vector2D<int>(maxM, Vector<int>(maxN));

	#pragma omp parallel for
	for (int i = 0; i < maxM; ++i) {
		float x = _distrX(_generator); float y = _distrY(_generator); float z = _distrZ(_generator);
		_randomPositions[i] = this->search_nearest_neighbors(x, y, z, maxN);
	}
}

Vector<int> SearchTree::get_nearest_neighbors(int randPos, int n) {
	assert(randPos >= 0 && randPos < static_cast<int>(_randomPositions.size()));

	return {_randomPositions[randPos].begin(), _randomPositions[randPos].begin() + n};
}

Vector<int> SearchTree::search_nearest_neighbors(float x, float y, float z, unsigned int n) {
	using Pair = std::pair<double, int>;
	std::set<Pair> distNeighbors;
	std::unordered_set<int> uniqueNeighbors;

	double radius = 3; // optimize this value?
	int maxRadius = 5;
	do {
		for (auto i = Aboria::euclidean_search(_particles.get_query(), Aboria::vdouble3(x, y, z), radius); i != false; ++i) {
			int node = Aboria::get<Aboria::id>(*i);
			auto inserted = distNeighbors.emplace(Pair{(i.dx()).norm(), node});
			if (inserted.second) {
				uniqueNeighbors.emplace(inserted.first->second);
			}
		}

		if (radius >= maxRadius && uniqueNeighbors.size() < n) {
			maxRadius += 1;
		}

		radius += 1;
	} while (radius <= maxRadius && radius <= _largestRadius);

	Vector<int> neighbors(n); int idxNeighbors = 0;
	std::set<int> alreadyReturned;
	for (const auto & neighbor : distNeighbors) {
		auto inserted = alreadyReturned.emplace(neighbor.second);
		if (inserted.second) {
			neighbors[idxNeighbors++] = neighbor.second;
			if (alreadyReturned.size() >= n) {
				break;
			}
		}
	}

	return neighbors;
}
