#include "SearchTree.hpp"

#include <set>

void SearchTree::add_positions(const py::array_t<double> & positions) {
	auto rows = positions.unchecked<2>();

	for (int i = 0; i < _nParticles; ++i) {
		//if (_nDimensions == 3) {}
		Aboria::get<Particle_position>(_particles)[i] = Aboria::vdouble3(rows(i, 0), rows(i, 1), rows(i, 2));
	}
}

void SearchTree::init_search(double xMin, double xMax, double yMin, double yMax, double zMin, double zMax) {
	Aboria::vdouble3 min = Aboria::vdouble3(xMin, yMin, zMin);
	Aboria::vdouble3 max = Aboria::vdouble3(xMax, yMax, zMax);
	Aboria::vbool3 periodic = Aboria::vbool3(_pbcX, _pbcY, _pbcZ);
	_largestRadius = std::ceil((min - max).norm());
	_particles.init_neighbour_search(min, max, periodic);
}

Vector<int> SearchTree::search_nearest_neighbors(float x, float y, float z, unsigned int n, const Graph & completeGraph) {
	using Pair = std::pair<double, int>;
	std::set<Pair> distNeighbors;
	std::unordered_set<int> uniqueNeighbors;

	auto checkIsNeighbor = [&completeGraph](const std::unordered_set<int> & uN, int node) {
		if (uN.size() == 0) return true;

		for (const auto & n : uN) {
			if (n != node && completeGraph.has_neighbor(n, node)) {
				return true;
			}
		}

		return false;
	};

	double radius = 3; // optimize this value?
	int maxRadius = 5;
	do {
		for (auto i = Aboria::euclidean_search(_particles.get_query(), Aboria::vdouble3(x, y, z), radius); i != false; ++i) {
			int node = Aboria::get<Aboria::id>(*i);
			if (!checkIsNeighbor(uniqueNeighbors, node)) {
				continue;
			}

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

	#ifdef DEBUG
	std::cout << "x = " << x << " \ty = " << y << " \tz = " << z << " \tN = ";
	#endif

	Vector<int> neighbors(n); int idxNeighbors = 0;
	std::set<int> alreadyReturned;
	for (const auto & neighbor : distNeighbors) {
		auto inserted = alreadyReturned.emplace(neighbor.second);
		if (inserted.second) {
			neighbors[idxNeighbors++] = neighbor.second;

			#ifdef DEBUG
			std::cout << neighbor.second << " ";
			#endif

			if (alreadyReturned.size() >= n) {
				break;
			}
		}
	}

	#ifdef DEBUG
	std::cout << "\n";
	#endif

	return neighbors;
}
