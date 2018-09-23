#include "IsIsomorphic.hpp"

bool is_isomorphic(const UndirectedGraph & uGraph1, const UndirectedGraph & uGraph2) {
	vf2_callback<UndirectedGraph, UndirectedGraph> callback(uGraph1, uGraph2);

	// return vf2_graph_iso(uGraph1, uGraph2, callback);
	return true;
}
