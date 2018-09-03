#include "Graphs.hpp"

double calc_shannon_entropy(double Hn, double fi, double m) {
	double pi = fi / m;
	Hn -= (pi * log(pi));
	return Hn;
}

const Point Graphs::generate_random_point() {
	Point rp;
	std::get<0>(rp) = _distrX(_generator);
	std::get<1>(rp) = _distrY(_generator);
	std::get<2>(rp) = _distrZ(_generator);
	return rp;
}

py::tuple Graphs::calculate_configurational_entropy(double c) {
	Vector2D<int> closestNeighbors = this->get_closest_neighbors();
	this->generate_subgraphs(closestNeighbors);
	return this->check_isomorfism(c);
}

const Vector2D<int> Graphs::get_closest_neighbors() {
	Vector2D<int> closestNeighbors(_m);
	for (int i = 0; i < _m; ++i) {
		Point rp = this->generate_random_point();
		double x = std::get<0>(rp);
		double y = std::get<1>(rp);
		double z = std::get<2>(rp);
		closestNeighbors[i] = this->search_nearest_neighbors(x, y, z, _n);
	}
	return closestNeighbors;
}

void Graphs::generate_subgraphs(const Vector2D<int> & closestNeighbors) {
	int pos = 0;
	for (const auto & nClosestNeighbors : closestNeighbors) {
		Graph graph = Graph();
		for (const auto & node : nClosestNeighbors) {
			if (_completeGraph.has_node(node)) {
				if (!graph.has_node(node)) {
					graph.add_node(node);
				}

				const Vector<int> & neighbors = _completeGraph.get_neighbors(node);
				for (const auto & neighbor : neighbors) {
					if (std::find(nClosestNeighbors.begin(), nClosestNeighbors.end(), neighbor) != nClosestNeighbors.end()) {
						graph.add_edge(node, neighbor);
					}
				}
			}
		}

		this->insert(pos++, graph);
	}
}

py::tuple Graphs::check_isomorfism(double c) {
	std::map<int, int> label_total;
	int iso_label = 1;

	int size = _graphs.size();
	for (int i = 0; i < size; ++i) {
		for (int j = i + 1; j < size; ++j) {
			int iso_label_i = _graphs[i].get_iso_label();
			int iso_label_j = _graphs[j].get_iso_label();

			if (iso_label_i == 0 || iso_label_j == 0) {
				if (is_isomorphic(*_graphs[i].getGraph(), *_graphs[j].getGraph())) {
					if (iso_label_i == 0 && iso_label_j == 0) {
						_graphs[i].set_iso_label(iso_label);
						_graphs[j].set_iso_label(iso_label);
						label_total[iso_label] = 2;
						iso_label += 1; // label already used
					} else if (iso_label_i > 0 && iso_label_j == 0) {
						_graphs[j].set_iso_label(iso_label_i);
						label_total[iso_label_i] += 1;
					} else if (iso_label_j > 0 && iso_label_i == 0) {
						_graphs[i].set_iso_label(iso_label_j);
						label_total[iso_label_j] += 1;
					} else if (iso_label_i != iso_label_j) {
						std::cout << "Error while checking isomorphism:\nlabelGi " << iso_label_i << " : labelGj " << iso_label_j << "\n";
					}
				}
			}
		}
	}

	// get all graphs that are not isomorphic with any other
	for (const auto & g : _graphs) {
		if (g.get_iso_label() == 0) {
			label_total[iso_label] = 1;
			iso_label += 1;
		}
	}

	double H_n = 0.0, H1n = 0.0;
	for (int i = 1; i < iso_label; ++i) {
		double fi = double(label_total[i]);
		H_n = calc_shannon_entropy(H_n, fi, _m);
		if (fi == 1.0) {
			H1n = calc_shannon_entropy(H1n, fi, _m);
		}
	}

	double H1nDiv = 0.0;
	if (H_n > 0) {
		H1nDiv = (H1n / H_n);
	}

	double H_n_extrapolated = H_n + (c * H1nDiv);
	double g_n = 2 * log(_n); // (spatial_dimensions - 1)
	double Hc_n = H_n_extrapolated - g_n;

	bool valid = true;
	if (H1n > (H_n / 100)) {
		valid = false;
	}

	return py::make_tuple(Hc_n, valid);
}
