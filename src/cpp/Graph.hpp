#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/isomorphism.hpp>
#include <boost/graph/vf2_sub_graph_iso.hpp>

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
	Graph() : isoLabel(0), graph(std::make_shared<UndirectedGraph>()) {}

	inline std::shared_ptr<UndirectedGraph> getGraph() const { return graph; }
	inline int get_iso_label() const { return isoLabel; }
	inline void set_iso_label(int _isoLabel) { isoLabel = _isoLabel; }
	VertexDescriptor add_node(int node);
	void add_edge(int e1, int e2);
	bool has_node(int node) const;
	Vector<int> get_neighbors(int node) const;
	bool has_neighbor(int node, int neighbor);
	int get_total_nodes() const;
	int get_total_edges() const;
	void print_graph() const;

private:
	int isoLabel;
	std::shared_ptr<UndirectedGraph> graph;
	std::map<int, VertexDescriptor> mVertexDesc;
};

template <typename Graph1, typename Graph2>
struct vf2_callback {

	vf2_callback(const Graph1 & graph1, const Graph2 & graph2) {}

	template <typename CorrespondenceMap1To2, typename CorrespondenceMap2To1>
	bool operator()(CorrespondenceMap1To2, CorrespondenceMap2To1) const {
		// return on the first mapping found
		return false;
	}
};

bool is_isomorphic(const UndirectedGraph & uGraph1, const UndirectedGraph & uGraph2);

#endif /* GRAPH_HPP */
