#include "Graph.hpp"

#include <iostream>
#include <boost/graph/connected_components.hpp>

std::unique_ptr<graph[]> Graph::get_cannonical_label() const {
	static DEFAULTOPTIONS_GRAPH(options);
	options.getcanon = TRUE;

	int n = this->get_total_nodes();
	int m = SETWORDSNEEDED(n);
	nauty_check(WORDSIZE, m, n, NAUTYVERSIONID);

	std::unique_ptr<graph[]> ngraph(new graph[n * m]);

	EMPTYGRAPH(ngraph.get(), m, n);  // start with no edges

	auto ei = boost::edges(*this->get_ugraph());
	for (EdgeIterator edge_iter = ei.first; edge_iter != ei.second; ++edge_iter) {
		VertexDescriptor vs = boost::source(*edge_iter, *this->get_ugraph());
		VertexDescriptor vt = boost::target(*edge_iter, *this->get_ugraph());

		ADDONEEDGE(ngraph.get(), vs, vt, m);
	}

	std::unique_ptr<int[]> lab(new int[n]);
	std::unique_ptr<int[]> ptn(new int[n]);
	std::unique_ptr<int[]> orbits(new int[n]);
	std::unique_ptr<graph[]> cannonicalLabel(new graph[n * m]);

	statsblk stats;
	densenauty(ngraph.get(), lab.get(), ptn.get(), orbits.get(), &options, &stats, m, n, cannonicalLabel.get());

	return cannonicalLabel;
}

VertexDescriptor Graph::add_node(int node) {
	VertexDescriptor v = boost::add_vertex(node, *this->get_ugraph());
	mVertexDesc[node] = v;

	return v;
}

VertexDescriptor Graph::get_node_descriptor(int node) const {
	auto itNode = mVertexDesc.find(node);
	if (itNode != mVertexDesc.end()) {
		return itNode->second;
	}

	return -1;
}

bool Graph::is_connected() const {
    Vector<int> component(this->get_total_nodes());
    return boost::connected_components(*this->get_ugraph(), &component[0]) == 1;
}

void Graph::add_edge(int e1, int e2) {
	if (e1 == e2) {
		return;
	}

	// if vertices are not present, add them
	auto itE1 = mVertexDesc.find(e1);
	VertexDescriptor ve1;
	if (itE1 == mVertexDesc.end()) {
		ve1 = this->add_node(e1);
	} else {
		ve1 = itE1->second;
	}

	auto itE2 = mVertexDesc.find(e2);
	VertexDescriptor ve2;
	if (itE2 == mVertexDesc.end()) {
		ve2 = this->add_node(e2);
	} else {
		ve2 = itE2->second;
	}

	if (!boost::edge(ve1, ve2, *this->get_ugraph()).second) {
		boost::add_edge(ve1, ve2, *this->get_ugraph());
	}
}

bool Graph::has_node(int node) const {
	return mVertexDesc.find(node) != mVertexDesc.end();
}

bool Graph::has_neighbor(int node, int neighbor) const {
	auto itNode = mVertexDesc.find(node);
	VertexDescriptor vNode;
	if (itNode == mVertexDesc.end()) {
		return false;
	} else {
		vNode = itNode->second;
	}

	auto itNeighbor = mVertexDesc.find(neighbor);
	VertexDescriptor vNeighbor;
	if (itNeighbor == mVertexDesc.end()) {
		return false;
	} else {
		vNeighbor = itNeighbor->second;
	}

	return boost::edge(vNode, vNeighbor, *this->get_ugraph()).second;
}

Vector<int> Graph::get_neighbors(int node) const {
	Vector<int> neighbors;
	neighbors.reserve(5); // provavelmente nao vai ter mais vizinhos que isso

	auto it = mVertexDesc.find(node); // int, VertexDescriptor
	if (it == mVertexDesc.end()) { return neighbors; }

	AdjacencyIterator neighbor, neighbor_end;
	for (tie(neighbor, neighbor_end) = boost::adjacent_vertices(it->second, *this->get_ugraph()); neighbor != neighbor_end; ++neighbor) {
		neighbors.push_back((*this->get_ugraph())[*neighbor]);
	}

	return neighbors;
}

int Graph::get_total_nodes() const {
	return boost::num_vertices(*this->get_ugraph());
}

int Graph::get_total_edges() const {
	return boost::num_edges(*this->get_ugraph());
}

void Graph::print_graph() const {
	std::cout << "-----------------------------\n";
	std::cout << "vertices:\n";
	std::cout << boost::num_vertices(*this->get_ugraph()) << "\n";
	std::pair<VertexIterator, VertexIterator> vi = boost::vertices(*this->get_ugraph());
	for (VertexIterator vertex_iter = vi.first; vertex_iter != vi.second; ++vertex_iter) {
		std::cout << "(" << (*this->get_ugraph())[*vertex_iter] << ") ";
	}

	std::cout << "\nedges:\n";
	std::cout << boost::num_edges(*this->get_ugraph()) << "\n";

	std::pair<EdgeIterator, EdgeIterator> ei = boost::edges(*this->get_ugraph());
	for (EdgeIterator edge_iter = ei.first; edge_iter != ei.second; ++edge_iter) {
		VertexDescriptor vs = boost::source(*edge_iter, *this->get_ugraph());
		VertexDescriptor vt = boost::target(*edge_iter, *this->get_ugraph());
		std::cout << "(" << (*this->get_ugraph())[vs] << ", " << (*this->get_ugraph())[vt] << ") ";
	}

	std::cout << "\n-----------------------------\n";
}
