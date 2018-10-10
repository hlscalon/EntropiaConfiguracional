#include "FindIsomorphicIndex.hpp"

#include <parallel/algorithm>

bool is_isomorphic(const graph * clgraph1, const graph * clgraph2, int totalNodes) {
	int m = SETWORDSNEEDED(totalNodes);
	return memcmp(clgraph1, clgraph2, m * sizeof(graph) * totalNodes) == 0;
}

int find_isomorphic_index(const Vector<GraphIsomorphism> & graphsIso, const graph * clgraph, int totalNodes) {
	auto itIdx = __gnu_parallel::find_if(graphsIso.begin(), graphsIso.end(), [&](const GraphIsomorphism & g) {
		return is_isomorphic(g.get_cannonical_label(), clgraph, totalNodes);
	});

	return itIdx == graphsIso.end() ? -1 : itIdx - graphsIso.begin();
}
