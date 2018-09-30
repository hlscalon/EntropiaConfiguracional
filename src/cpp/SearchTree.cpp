#include "SearchTree.hpp"

#include <set>
#include <algorithm>

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

bool SearchTree::check_is_connected(const Graph & g, const Vector<int> & verticesNewGraph) {
	int total = verticesNewGraph.size();
	for (int i = 0; i < total; ++i) {
		int j = i + 1;
		bool hasSomeNeighbor = j >= total; // ultima iteracao = sempre true
		for (; j < total; ++j) {
			// std::cout << "neighbor: " << verticesNewGraph[i] << " : " << verticesNewGraph[j] << " = ";
			if (g.has_neighbor(verticesNewGraph[i], verticesNewGraph[j])) {
				// std::cout << " yes\n";
				hasSomeNeighbor = true;
				break;
			}
			// std::cout << " no\n";
		}

		if (!hasSomeNeighbor) {
			return false;
		}
	}

	return true;
}

std::pair<bool, Vector<int>> SearchTree::generate_all_combinations(const Graph & g, const Vector<int> & vertices, int size, unsigned int n) {
	std::string bitmask(n, 1); // K leading 1's
	bitmask.resize(size, 0); // N-K trailing 0's

	Vector<int> verticesNewGraph;
	// print integers and permute bitmask
	bool isConnected = false;
	do {
		// std::cout << "bitmask: "; for (int i = 0; i < size; i++) std::cout << (bitmask[i] ? "1" : "0"); std::cout << "\n";

		verticesNewGraph.clear();
		verticesNewGraph.reserve(n);
		for (int i = 0; i < size; ++i) { // [0..N-1] integers
			if (bitmask[i]) verticesNewGraph.push_back(vertices[i]);
		}

		// for (const auto & v: verticesNewGraph) std::cout << v << " ";
		// std::cout << "\n";

		isConnected = this->check_is_connected(g, verticesNewGraph);

		// std::cout << "is_connected: " << isConnected << "\n\n";
	} while (std::prev_permutation(bitmask.begin(), bitmask.end()) && !isConnected);

	return {isConnected, verticesNewGraph};
}

const Vector<int> SearchTree::get_neighbors_connected(const Graph & g, const Vector<int> & vertices, unsigned int n) {
	int i = n; // comeca com n primeiros
	int N = vertices.size();
	for (; i < N; ++i) {
		Vector<int> verticesCheck(vertices.begin(), vertices.begin() + i);
		auto ret = this->generate_all_combinations(g, verticesCheck, i, n);
		if (ret.first) {
			// std::cout << "Achou: "; for (auto n : ret.second) std::cout << n << " "; std::cout << "\n";
			return ret.second;
		}
	}

	return {vertices.begin(), vertices.begin() + n}; // pega os n primeiros
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

	/*
	Vector<int> neighbors;
	neighbors.reserve(n);

	for (auto itneighbors = distNeighbors.begin(); itneighbors != distNeighbors.end(); ++itneighbors) {
		auto neighbor = *itneighbors;

		if (std::find(neighbors.begin(), neighbors.end(), neighbor.second) == neighbors.end()) {
			neighbors.push_back(neighbor.second);

			if (neighbors.size() >= n) {
				break;
			}
		}
	}
	*/

	Vector<int> neighborsSorted;
	neighborsSorted.reserve(distNeighbors.size());

	for (const auto & itneighbors : distNeighbors) {
		if (std::find(neighborsSorted.begin(), neighborsSorted.end(), itneighbors.second) == neighborsSorted.end()) {
			neighborsSorted.push_back(itneighbors.second);
		}
	}

	Vector<int> neighbors = this->get_neighbors_connected(completeGraph, neighborsSorted, n);

	return neighbors;
}
