#include "FindIsomorphicIndex.hpp"

#include <parallel/algorithm>

bool is_isomorphic(const graph * clgraph1, const graph * clgraph2, int totalNodes) {
	int m = SETWORDSNEEDED(totalNodes);
	return memcmp(clgraph1, clgraph2, m * sizeof(graph) * totalNodes) == 0;
}

int find_isomorphic_index(const Vector<graph*> & ngraphs, const graph * clgraph, int totalNodes) {
	auto itIdx = __gnu_parallel::find_if(ngraphs.begin(), ngraphs.end(), [&](const graph * clg) {
		return is_isomorphic(clg, clgraph, totalNodes);
	});

	return itIdx == ngraphs.end() ? -1 : itIdx - ngraphs.begin();
}
