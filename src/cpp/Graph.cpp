#include "Graph.hpp"

#include <iostream>

VertexDescriptor Graph::add_node(int node) {
	VertexDescriptor v = boost::add_vertex(node, *this->graph);
	mVertexDesc[node] = v;

	return v;
}

void Graph::add_edge(int e1, int e2) {
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

	if (!boost::edge(ve1, ve2, *this->graph).second) {
		boost::add_edge(ve1, ve2, *this->graph);
	}
}

bool Graph::has_node(int node) const {
	return mVertexDesc.find(node) != mVertexDesc.end();
}

bool Graph::has_neighbor(int node, int neighbor) {
	if (mVertexDesc.find(node) == mVertexDesc.end()) { return false; }
	if (mVertexDesc.find(neighbor) == mVertexDesc.end()) { return false; }

	return boost::edge(mVertexDesc[node], mVertexDesc[neighbor], *this->graph).second;
}

Vector<int> Graph::get_neighbors(int node) const {
	Vector<int> neighbors;
	neighbors.reserve(5); // provavelmente nao vai ter mais vizinhos que isso

	auto it = mVertexDesc.find(node); // int, VertexDescriptor
	if (it == mVertexDesc.end()) { return neighbors; }

	AdjacencyIterator neighbor, neighbor_end;
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
	std::cout << "-----------------------------\n";
	std::cout << "vertices:\n";
	std::cout << boost::num_vertices(*this->graph) << "\n";
	std::pair<VertexIterator, VertexIterator> vi = boost::vertices(*this->graph);
	for (VertexIterator vertex_iter = vi.first; vertex_iter != vi.second; ++vertex_iter) {
		std::cout << "(" << (*this->graph)[*vertex_iter] << ")\n";
	}

	std::cout << "edges:\n";
	std::cout << boost::num_edges(*this->graph) << "\n";

	std::pair<EdgeIterator, EdgeIterator> ei = boost::edges(*this->graph);
	for (EdgeIterator edge_iter = ei.first; edge_iter != ei.second; ++edge_iter) {
		VertexDescriptor vs = boost::source(*edge_iter, *this->graph);
		VertexDescriptor vt = boost::target(*edge_iter, *this->graph);
		std::cout << "(" << (*this->graph)[vs] << ", " << (*this->graph)[vt] << ")\n";
	}

	std::cout << "-----------------------------\n";
}
