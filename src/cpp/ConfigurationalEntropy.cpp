#include "ConfigurationalEntropy.hpp"

#include <omp.h>

const Point ConfigurationalEntropy::generate_random_point(int precision) {
	auto round = [](float number, int precision) {
		int decimals = std::pow(10, precision);
		return (std::round(number * decimals)) / decimals;
	};

	Point rp;
	std::get<0>(rp) = round(_distrX(_generator), precision);
	std::get<1>(rp) = round(_distrY(_generator), precision);
	std::get<2>(rp) = round(_distrZ(_generator), precision);
	return rp;
}

py::tuple ConfigurationalEntropy::calculate(int m, int n, double c) {
	Vector<Graph> graphs = this->get_subgraphs(m, n);
	return this->check_isomorfism(graphs, c, m, n);
}

const Vector<Graph> ConfigurationalEntropy::get_subgraphs(int m, int n) {
	std::map<Vector<int>, int> differentGraphs;
	std::map<Point, Vector<int>> nearestNeighborsFromPoint;

	for (int i = 0; i < m; ++i) {
		Point rp = this->generate_random_point(1);

		Vector<int> nearestNeighbors;

		auto itNN = nearestNeighborsFromPoint.find(rp);
		if (itNN != nearestNeighborsFromPoint.end()) {
			nearestNeighbors = itNN->second;
		} else {
			float x = std::get<0>(rp); float y = std::get<1>(rp); float z = std::get<2>(rp);
			nearestNeighbors = this->search_nearest_neighbors(x, y, z, n);
			std::sort(nearestNeighbors.begin(), nearestNeighbors.end());
			nearestNeighborsFromPoint[rp] = nearestNeighbors;
		}

		differentGraphs[nearestNeighbors]++;
	}

	Vector<Graph> graphs(differentGraphs.size());
	auto it = differentGraphs.begin();
	for (int i = 0; it != differentGraphs.end(); ++it, ++i) {
		this->generate_subgraph(graphs[i], it->first);
		graphs[i].add_qty(it->second - 1);
	}

	return graphs;
}

Graph ConfigurationalEntropy::generate_subgraph(Graph & graph, const Vector<int> & closestNeighbors) {
	for (const auto & node : closestNeighbors) {
		if (_completeGraph.has_node(node)) {
			if (!graph.has_node(node)) {
				graph.add_node(node);
			}

			const Vector<int> & neighbors = _completeGraph.get_neighbors(node); // fazer uma cache disso
			for (const auto & neighbor : neighbors) {
				if (std::find(closestNeighbors.begin(), closestNeighbors.end(), neighbor) != closestNeighbors.end()) {
					graph.add_edge(node, neighbor);
				}
			}
		}
	}

	return graph;
}

py::tuple ConfigurationalEntropy::check_isomorfism(Vector<Graph> & graphs, double c, int m, int n) {
	std::map<int, int> label_total;
	int iso_label = 1;

	int size = graphs.size();
	for (int i = 0; i < size; ++i) {
		for (int j = i + 1; j < size; ++j) {
			int iso_label_i = graphs[i].get_iso_label();
			int iso_label_j = graphs[j].get_iso_label();

			if (iso_label_i == 0 || iso_label_j == 0) {
				if (is_isomorphic(*graphs[i].getGraph(), *graphs[j].getGraph())) {
					if (iso_label_i == 0 && iso_label_j == 0) {
						graphs[i].set_iso_label(iso_label);
						graphs[j].set_iso_label(iso_label);
						label_total[iso_label] = graphs[i].get_qty() + graphs[j].get_qty();
						iso_label += 1; // label already used
					} else if (iso_label_i > 0 && iso_label_j == 0) {
						graphs[j].set_iso_label(iso_label_i);
						label_total[iso_label_i] += graphs[j].get_qty();
					} else if (iso_label_j > 0 && iso_label_i == 0) {
						graphs[i].set_iso_label(iso_label_j);
						label_total[iso_label_j] += graphs[i].get_qty();
					} else if (iso_label_i != iso_label_j) {
						std::cout << "Error while checking isomorphism:\nlabelGi " << iso_label_i << " : labelGj " << iso_label_j << "\n";
					}
				}
			}
		}
	}

	// get all graphs that are not isomorphic with any other
	for (const auto & g : graphs) {
		if (g.get_iso_label() == 0) {
			label_total[iso_label] = 1;
			iso_label += 1;
		}
	}

	double H_n = 0.0, H1n = 0.0;
	for (int i = 1; i < iso_label; ++i) {
		double fi = double(label_total[i]);
		H_n = this->calc_shannon_entropy(H_n, fi, m);
		if (fi == 1.0) {
			H1n = this->calc_shannon_entropy(H1n, fi, m);
		}
	}

	double H1nDiv = 0.0;
	if (H_n > 0) {
		H1nDiv = (H1n / H_n);
	}

	double H_n_extrapolated = H_n + (c * H1nDiv);
	double g_n = 2 * log(n); // (spatial_dimensions - 1)
	double Hc_n = H_n_extrapolated - g_n;

	bool valid = true;
	if (H1n > (H_n / 100)) {
		valid = false;
	}

	return py::make_tuple(Hc_n, valid);
}

double ConfigurationalEntropy::calc_shannon_entropy(double Hn, double fi, double m) {
	double pi = fi / m;
	Hn -= (pi * log(pi));
	return Hn;
}
