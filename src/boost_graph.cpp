#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <iostream>
#include <memory>
#include <map>

#include "boost/graph/graph_traits.hpp"
#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/isomorphism.hpp"
#include "boost/graph/vf2_sub_graph_iso.hpp"

namespace py = pybind11;

/*
g++ -O3 -Wall -std=c++14 boost_graph.cpp -o boost_graph
g++ -O3 -Wall -shared -std=c++14 -Ilibs/pybind11/include -fPIC `python-config --includes` boost_graph.cpp -o boost_graph.so

python2 main.py ../arquivos_xyz/fcc.xyz 1.12 0 3 10 Y
python2 main.py graph_files/fcc.xyz 1.12 0 3 10
python2 -m cProfile -s time main.py ../../../../python/graph_files/fcc.xyz 1.12 0 3 10
*/

using Edge = std::pair<int, int>;
using Vertex = int;
using UndirectedGraph = boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, Vertex, Edge>;
using edge_iterator = boost::graph_traits<UndirectedGraph>::edge_iterator;
using vertex_iterator = boost::graph_traits<UndirectedGraph>::vertex_iterator;
using vertex_descriptor = boost::graph_traits<UndirectedGraph>::vertex_descriptor;
using adjacency_iterator = boost::graph_traits<UndirectedGraph>::adjacency_iterator;

struct Graph {
	Graph() : isoLabel(0), graph(std::make_shared<UndirectedGraph>()) {}

	inline std::shared_ptr<UndirectedGraph> getGraph() const { return graph; }
	inline int get_iso_label() const { return isoLabel; }
	inline void set_iso_label(int _isoLabel) { isoLabel = _isoLabel; }
	vertex_descriptor add_node(int node);
	void add_edge(int e1, int e2);
	bool has_node(int node) const;
	//py::list get_neighbors(int node);
	std::vector<int> get_neighbors(int node) const;
	bool has_neighbor(int node, int neighbor);
	int get_total_nodes() const;
	int get_total_edges() const;
	void print_graph() const;

private:
	int isoLabel;
	std::shared_ptr<UndirectedGraph> graph;
	std::map<int, vertex_descriptor> mVertexDesc;
};

vertex_descriptor Graph::add_node(int node) {
	// std::cout << "add_node\n" << node << "\n";

	vertex_descriptor v = boost::add_vertex(node, *this->graph);
	mVertexDesc[node] = v;

	return v;
}

void Graph::add_edge(int e1, int e2) {
	// if vertices are not present, add them
	auto itE1 = mVertexDesc.find(e1);
	vertex_descriptor ve1;
	if (itE1 == mVertexDesc.end()) {
		ve1 = this->add_node(e1);
	} else {
		ve1 = itE1->second;
	}

	auto itE2 = mVertexDesc.find(e2);
	vertex_descriptor ve2;
	if (itE2 == mVertexDesc.end()) {
		ve2 = this->add_node(e2);
	} else {
		ve2 = itE2->second;
	}

	// std::cout << "add_edge\n" << e1 << " " << e2 << "\n";
	boost::add_edge(ve1, ve2, *this->graph);
}

bool Graph::has_node(int node) const {
	// std::cout << "has_node\n" << node << "\n";

	return mVertexDesc.find(node) != mVertexDesc.end();
}

bool Graph::has_neighbor(int node, int neighbor) {
	// std::cout << "has_neighbor\n" << node << " " << neighbor << "\n";

	if (mVertexDesc.find(node) == mVertexDesc.end()) { return false; }
	if (mVertexDesc.find(neighbor) == mVertexDesc.end()) { return false; }

	return boost::edge(mVertexDesc[node], mVertexDesc[neighbor], *this->graph).second;
}

std::vector<int> Graph::get_neighbors(int node) const {
	std::vector<int> neighbors;

	// std::cout << "get_neighbors\n" << node << "\n";

	auto it = mVertexDesc.find(node); // int, vertex_descriptor
	if (it == mVertexDesc.end()) { return neighbors; }

	adjacency_iterator neighbor, neighbor_end;
	for (tie(neighbor, neighbor_end) = boost::adjacent_vertices(it->second, *this->graph); neighbor != neighbor_end; ++neighbor) {
		neighbors.push_back((*this->graph)[*neighbor]);
	}

	return neighbors;
}

int Graph::get_total_nodes() const {
	return boost::num_vertices(*this->graph);
}

int Graph::get_total_edges() const {
	return boost::num_edges(*this->graph);
}

void Graph::print_graph() const {
	// std::cout << "print2\n";
	std::cout << "-----------------------------\n";
	std::cout << "vertices:\n";
	std::cout << boost::num_vertices(*this->graph) << "\n";
	std::pair<vertex_iterator, vertex_iterator> vi = boost::vertices(*this->graph);
	for (vertex_iterator vertex_iter = vi.first; vertex_iter != vi.second; ++vertex_iter) {
		std::cout << "(" << (*this->graph)[*vertex_iter] << ")\n";
	}

	std::cout << "edges:\n";
	std::cout << boost::num_edges(*this->graph) << "\n";

	std::pair<edge_iterator, edge_iterator> ei = boost::edges(*this->graph);
	for (edge_iterator edge_iter = ei.first; edge_iter != ei.second; ++edge_iter) {
		vertex_descriptor vs = boost::source(*edge_iter, *this->graph);
		vertex_descriptor vt = boost::target(*edge_iter, *this->graph);
		std::cout << "(" << (*this->graph)[vs] << ", " << (*this->graph)[vt] << ")\n";
	}

	std::cout << "-----------------------------\n";
}

template <typename Graph1, typename Graph2>
struct vf2_callback {

	vf2_callback(const Graph1 & graph1, const Graph2 & graph2) {}

	template <typename CorrespondenceMap1To2, typename CorrespondenceMap2To1>
	bool operator()(CorrespondenceMap1To2, CorrespondenceMap2To1) const {
		// return on the first mapping found
		return false;
	}
};

bool is_isomorphic(const UndirectedGraph & uGraph1, const UndirectedGraph & uGraph2) {
	vf2_callback<UndirectedGraph, UndirectedGraph> callback(uGraph1, uGraph2);

	return vf2_graph_iso(uGraph1, uGraph2, callback);
}

double calc_shannon_entropy(double Hn, double fi, double m) {
	double pi = fi / m;
	Hn -= (pi * log(pi));
	return Hn;
}

struct Graphs {
	Graphs() {}
	Graphs(int qty) : graphs(std::vector<Graph>(qty)) {}

	inline void insert(int pos, const Graph & _graph) { /*graphs.push_back(_graph);*/ graphs[pos] = _graph; }
	py::tuple check_isomorfism(double n, double m, double c);
	void generate_subgraphs(const Graph & G, int n, const py::list & closestNeighbors);
private:
	std::vector<Graph> graphs;
};

void Graphs::generate_subgraphs(const Graph & G, int n, const py::list & closestNeighbors) {
	std::vector<std::vector<int>> vecClosestNeighbors = closestNeighbors.cast<std::vector<std::vector<int>>>();

	int pos = 0;
	for (const auto & nClosestNeighbors : vecClosestNeighbors) {
		Graph graph = Graph();
		for (const auto & node : nClosestNeighbors) {
			//int node = pynode.cast<int>();
			if (G.has_node(node)) {
				if (!graph.has_node(node)) {
					graph.add_node(node);
				}

				const std::vector<int> & neighbors = G.get_neighbors(node);
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

py::tuple Graphs::check_isomorfism(double n, double m, double c) {
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
						label_total[iso_label] = 2;
						iso_label += 1; // label already used
					} else if (iso_label_i > 0 && iso_label_j == 0) {
						graphs[j].set_iso_label(iso_label_i);
						label_total[iso_label_i] += 1;
					} else if (iso_label_j > 0 && iso_label_i == 0) {
						graphs[i].set_iso_label(iso_label_j);
						label_total[iso_label_j] += 1;
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
		H_n = calc_shannon_entropy(H_n, fi, m);
		if (fi == 1.0) {
			H1n = calc_shannon_entropy(H1n, fi, m);
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

PYBIND11_MODULE(boost_graph, m) {
	m.doc() = "pybind11 boost_graph plugin"; // optional module docstring

	py::class_<UndirectedGraph, std::shared_ptr<UndirectedGraph>>(m, "UndirectedGraph");

	py::class_<Graph, std::shared_ptr<Graph>>(m, "Graph")
		.def(py::init<>())
		.def("add_node", &Graph::add_node)
		.def("add_edge", &Graph::add_edge)
		.def("has_node", &Graph::has_node)
		.def("has_neighbor", &Graph::has_neighbor)
		.def("get_neighbors", &Graph::get_neighbors)
		.def("get_total_nodes", &Graph::get_total_nodes)
		.def("get_total_edges", &Graph::get_total_edges)
		.def("print_graph", &Graph::print_graph)
	;

	py::class_<Graphs>(m, "Graphs")
		.def(py::init<>())
		.def(py::init<int>())
		.def("insert", &Graphs::insert)
		.def("generate_subgraphs", &Graphs::generate_subgraphs)
		.def("check_isomorfism", &Graphs::check_isomorfism)
	;
}
