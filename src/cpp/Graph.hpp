#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <nauty.h>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <utility>
#include <memory>
#include <vector>
#include <map>

using Edge = std::pair<int, int>;
using Vertex = int;
using UndirectedGraph = boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, Vertex, Edge>;
using EdgeIterator = boost::graph_traits<UndirectedGraph>::edge_iterator;
using VertexIterator = boost::graph_traits<UndirectedGraph>::vertex_iterator;
using VertexDescriptor = boost::graph_traits<UndirectedGraph>::vertex_descriptor;
using AdjacencyIterator = boost::graph_traits<UndirectedGraph>::adjacency_iterator;

template <typename T>
using Vector = std::vector<T>;

template <typename T>
using Vector2D = Vector<Vector<T>>;

struct Graph {
	Graph() : isoLabel(0), qty(1), ugraph(std::make_unique<UndirectedGraph>()) {}

	inline UndirectedGraph * get_ugraph() const {
		return ugraph.get();
	}

	inline int get_iso_label() const {
		return isoLabel;
	}

	inline void set_iso_label(int _isoLabel) {
		isoLabel = _isoLabel;
	}

	inline void add_qty(int _qty) {
		qty += _qty;
	}

	inline int get_qty() const {
		return qty;
	}

	inline graph * get_cannonical_label() const {
		return cannonicalLabel.get();
	}

	void set_cannonical_label();

	VertexDescriptor add_node(int node);

	VertexDescriptor get_node_descriptor(int node) const;

	bool is_connected() const;

	void add_edge(int e1, int e2);

	bool has_node(int node) const;

	Vector<int> get_neighbors(int node) const;

	bool has_neighbor(int node, int neighbor) const;

	int get_total_nodes() const;

	int get_total_edges() const;

	void print_graph() const;

private:
	int isoLabel;
	int qty;
	std::unique_ptr<UndirectedGraph> ugraph;
	std::map<int, VertexDescriptor> mVertexDesc;
	std::unique_ptr<graph[]> cannonicalLabel;
};

bool is_isomorphic(const Graph & graph1, const Graph & graph2);

#endif /* GRAPH_HPP */
