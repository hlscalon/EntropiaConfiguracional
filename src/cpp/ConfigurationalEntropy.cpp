#include "ConfigurationalEntropy.hpp"
#include "FindIsomorphicIndex.hpp"
#include "GraphIsomorphism.hpp"

#include <parallel/algorithm>

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

const int ConfigurationalEntropy::generate_random_point() {
	return _distr(_generator);
}

void ConfigurationalEntropy::init_search(double xMin, double xMax, double yMin, double yMax, double zMin, double zMax, int maxN, int maxM, const py::array_t<double> & dMin, const py::array_t<double> & dMax) {
	PyArray<double> _dMin = dMin.unchecked<1>();
	PyArray<double> _dMax = dMax.unchecked<1>();

	_searchTree.init_search(xMin, xMax, yMin, yMax, zMin, zMax, _dMin, _dMax);
	_searchTree.generate_random_positions(maxN, maxM);
}

const std::tuple<int, Vector<int>> ConfigurationalEntropy::generate_subgraphs(int m, int n) {
	std::map<Vector<int>, GraphIsomorphism> cannonicalLabels;

	for (int i = 0; i < m; ++i) {
		const int rp = this->generate_random_point();
		Vector<int> nearestNeighbors(this->get_nearest_neighbors(rp, n));
		std::sort(nearestNeighbors.begin(), nearestNeighbors.end());

		// testa antes de colocar
		if (cannonicalLabels.find(nearestNeighbors) == cannonicalLabels.end()) {
			Graph graph;
			this->generate_subgraph(graph, nearestNeighbors);
			if (!graph.is_connected()) {
				i--; // continua na mesma iteracao, gera outro ponto aleatorio
			} else {
				cannonicalLabels[nearestNeighbors] = GraphIsomorphism(graph.get_cannonical_label());
			}
		} else {
			cannonicalLabels[nearestNeighbors].add_qty(1);
		}
	}

	int isoLabel = 1;
	Vector<int> labelTotal(m); // inicia todos no zero (maximo)
	Vector<GraphIsomorphism> graphsIso;
	graphsIso.reserve(cannonicalLabels.size());

	auto it = cannonicalLabels.begin();
	for (int i = 0; it != cannonicalLabels.end(); ++it, ++i) {
		GraphIsomorphism graphIsoTmp(std::move(it->second));
		this->check_isomorfism(graphsIso, graphIsoTmp, isoLabel, labelTotal, n);

		#ifdef DEBUG
		graphIsoTmp.print_graph();
		#endif

		if (graphIsoTmp.get_iso_label() == 0) { // nao achou iso
			graphsIso.push_back(std::move(graphIsoTmp));
		}
	}

	// get all graphsIso that are not isomorphic with any other
	for (const auto & g : graphsIso) {
		if (g.get_iso_label() == 0) {
			labelTotal[isoLabel] = g.get_qty();
			isoLabel++;
		}
	}

	return {isoLabel, labelTotal};
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

void ConfigurationalEntropy::check_isomorfism(Vector<GraphIsomorphism> & graphsIso, GraphIsomorphism & graphIsoTmp, int & isoLabel, Vector<int> & labelTotal, int totalNodes) {
	graph * clgraph = graphIsoTmp.get_cannonical_label();
	int idx = find_isomorphic_index(graphsIso, clgraph, totalNodes);

	if (idx >= 0) {
		int isoLabelIdx = graphsIso[idx].get_iso_label();

		if (isoLabelIdx == 0) {
			graphIsoTmp.set_iso_label(isoLabel);
			graphsIso[idx].set_iso_label(isoLabel);

			labelTotal[isoLabel] = graphsIso[idx].get_qty() + graphIsoTmp.get_qty();
			isoLabel += 1; // label already used
		} else if (isoLabelIdx > 0) {
			graphIsoTmp.set_iso_label(isoLabelIdx);
			labelTotal[isoLabelIdx] += graphIsoTmp.get_qty();
		}
	}
}

double ConfigurationalEntropy::calc_shannon_entropy(double Hn, double fi, double m) {
	double pi = fi / m;
	Hn -= (pi * log(pi));
	return Hn;
}
