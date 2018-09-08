#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <memory>

#include "Graph.hpp"
#include "ConfigurationalEntropy.hpp"

namespace py = pybind11;

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

	py::class_<ConfigurationalEntropy>(m, "ConfigurationalEntropy")
		.def(py::init<const Graph &, int, int, bool, bool, bool, const py::array_t<double> &, const py::array_t<double> &>())
		.def("add_positions", &ConfigurationalEntropy::add_positions)
		.def("init_search", &ConfigurationalEntropy::init_search)
		.def("calculate", &ConfigurationalEntropy::calculate)
	;
}
