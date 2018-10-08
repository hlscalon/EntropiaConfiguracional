#include "FindIsomorphicIndex.hpp"

#include <parallel/algorithm>

bool is_isomorphic(const graph * clgraph1, const graph * clgraph2, int totalNodes) {
	int m = SETWORDSNEEDED(totalNodes);
	return memcmp(clgraph1, clgraph2, m * sizeof(graph) * totalNodes) == 0;
}

int find_isomorphic_index(const Vector<Graph> & graphs, const graph * clgraph, int totalNodes) {
	auto itIdx = __gnu_parallel::find_if(graphs.begin(), graphs.end(), [&](const Graph & g) {
		return is_isomorphic(g.get_cannonical_label(), clgraph, totalNodes);
	});

	return itIdx == graphs.end() ? -1 : itIdx - graphs.begin();
}
