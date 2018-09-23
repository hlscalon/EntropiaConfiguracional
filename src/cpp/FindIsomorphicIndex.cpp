#include "FindIsomorphicIndex.hpp"
#include "IsIsomorphic.hpp"

#include <parallel/algorithm>

int findIsomorphicIndex(const Vector<UndirectedGraph> & uGraphs, const UndirectedGraph & uGraph, int size) {
	auto itIdx = __gnu_parallel::find_if(uGraphs.begin(), uGraphs.begin() + size, [&uGraph](const UndirectedGraph & ug) {
		return is_isomorphic(ug, uGraph);
	});

	return itIdx == uGraphs.begin() + size ? -1 : itIdx - uGraphs.begin();
}
