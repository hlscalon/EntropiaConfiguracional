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

	#ifdef DEBUG
	std::cout << "distNeighbors = ";
	for (const auto & d : distNeighbors) std::cout << d.second << " ";
	std::cout.precision(2);
	std::cout << "\nx = " << std::fixed << x << " \t\ty = " << y << " \t\tz = " << z << " \t\tN = \n";
	#endif

	auto checkIsNeighbor = [&](const Vector<int> & neighbors, int node) {
		if (neighbors.size() == 0) return true;

		for (const auto & n : neighbors) {
			if (n != node && completeGraph.has_neighbor(n, node)) {
				return true;
			}
		}

		return false;
	};

	Vector<int> neighbors;
	neighbors.reserve(n);

	bool achouAlgum = false, pararBusca = false; // cada passada tem que achar 1 pelo menos (loop infinito)
	int i = 0;
	comecar_denovo:
	do {
		achouAlgum = false;
		// for (const auto & neighbor : distNeighbors) {
		auto itneighbors = distNeighbors.begin();
		std::advance(itneighbors, i);

		for (; itneighbors != distNeighbors.end(); ++itneighbors) {
			auto neighbor = *itneighbors;
			if (!checkIsNeighbor(neighbors, neighbor.second)) {
				continue;
			}

			if (std::find(neighbors.begin(), neighbors.end(), neighbor.second) == neighbors.end()) {
				achouAlgum = true;
				neighbors.push_back(neighbor.second);

				if (neighbors.size() >= n) {
					pararBusca = true;
					break;
				}
			}
		}
	} while (achouAlgum && !pararBusca);

	if (i < static_cast<int>(distNeighbors.size() - n) && neighbors.size() < n) {
		achouAlgum = false; pararBusca = false;
		i++;

		#ifdef DEBUG
		std::cout << "N Intermediario = ";
		for (const auto & n : neighbors) std::cout << n << ";";
		std::cout << "\n\n";
		#endif

		neighbors.clear();
		goto comecar_denovo;
	}

	#ifdef DEBUG
	std::cout << "N Final = ";
	for (const auto & n : neighbors) std::cout << n << ";";
	std::cout << "\n\n";
	#endif

	return neighbors;
}
