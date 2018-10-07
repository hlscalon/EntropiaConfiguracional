#include "ConfigurationalEntropy.hpp"

#include <parallel/algorithm>

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
	std::tuple<int, Vector<int>> ret = this->generate_subgraphs(m, n);
	int iso_label = std::get<0>(ret);
	Vector<int> v_label_total = std::get<1>(ret);

	double H_n = 0.0, H1n = 0.0;
	for (int i = 1; i < iso_label; ++i) {
		double fi = double(v_label_total[i]);

		double H_n_antes = H_n;
		H_n = this->calc_shannon_entropy(H_n, fi, m);
		if (fi == 1.0) {
			H1n += (H_n - H_n_antes);
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
	if (H1n > (H_n_extrapolated / 100)) {
		valid = false;
	}

	#ifdef LOG
	std::cout << "H1nDiv = " << H1nDiv << "; H_n_extrapolated = " << H_n_extrapolated << "; Hc_n = " << Hc_n << "; H1n = " << H1n << "; valid = " << valid << "\n";
	#endif

	return py::make_tuple(Hc_n, valid);
}

const std::tuple<int, Vector<int>> ConfigurationalEntropy::generate_subgraphs(int m, int n) {
	std::map<Vector<int>, Graph> differentGraphs;
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

		// testa antes de colocar
		if (differentGraphs.find(nearestNeighbors) == differentGraphs.end()) {
			Graph graph;
			this->generate_subgraph(graph, nearestNeighbors);
			if (!graph.is_connected()) {
				i--; // continua na mesma iteracao, gera outro ponto aleatorio
			} else {
				// graph.set_cannonical_label();
				differentGraphs[nearestNeighbors] = std::move(graph);
			}
		} else {
			differentGraphs[nearestNeighbors].add_qty(1);
		}
	}

	int iso_label = 1;
	Vector<int> label_total(m); // inicia todos no zero (maximo)
	Vector<Graph> graphs;
	graphs.reserve(differentGraphs.size());

	auto it = differentGraphs.begin();
	for (int i = 0; it != differentGraphs.end(); ++it, ++i) {
		Graph graphTmp(std::move(it->second));
		this->check_isomorfism(graphs, graphTmp, iso_label, label_total, i);

		#ifdef DEBUG
		graphTmp.print_graph();
		#endif

		if (graphTmp.get_iso_label() == 0) { // nao achou iso
			graphs.push_back(std::move(graphTmp));
		}
	}

	// get all graphs that are not isomorphic with any other
	for (const auto & g : graphs) {
		if (g.get_iso_label() == 0) {
			label_total[iso_label] = g.get_qty();
			iso_label++;
		}
	}

	return {iso_label, label_total};
}

void ConfigurationalEntropy::generate_subgraph(Graph & graph, const Vector<int> & closestNeighbors) {
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
}

void ConfigurationalEntropy::check_isomorfism(Vector<Graph> & graphs, Graph & graph, int & iso_label, Vector<int> & label_total, int size) {
	auto itIdx = __gnu_parallel::find_if(graphs.begin(), graphs.end(), [&graph](const Graph & g) {
		//return is_isomorphic(g, graph);
		return is_isomorphic(*g.get_ugraph(), *graph.get_ugraph());
	});

	int idx = itIdx == graphs.end() ? -1 : itIdx - graphs.begin();
	if (idx >= 0) {
		int iso_label_idx = graphs[idx].get_iso_label();

		if (iso_label_idx == 0) {
			graph.set_iso_label(iso_label);
			graphs[idx].set_iso_label(iso_label);

			label_total[iso_label] = graphs[idx].get_qty() + graph.get_qty();
			iso_label += 1; // label already used
		} else if (iso_label_idx > 0) {
			graph.set_iso_label(iso_label_idx);
			label_total[iso_label_idx] += graph.get_qty();
		}
	}
}

double ConfigurationalEntropy::calc_shannon_entropy(double Hn, double fi, double m) {
	double pi = fi / m;
	Hn -= (pi * log(pi));
	return Hn;
}
