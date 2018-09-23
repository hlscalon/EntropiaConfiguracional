#ifndef IS_ISOMORPHIC_HPP
#define IS_ISOMORPHIC_HPP

#include "Graph.hpp"
#include "CudaInclude.hpp"

#include <boost/graph/isomorphism.hpp>
#include <boost/graph/vf2_sub_graph_iso.hpp>

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

#endif // IS_ISOMORPHIC_HPP