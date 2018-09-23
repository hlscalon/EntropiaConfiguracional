#ifndef FIND_ISOMORPHIC_INDEX_HPP
#define FIND_ISOMORPHIC_INDEX_HPP

#include "Graph.hpp"
#include "CudaInclude.hpp"

int findIsomorphicIndex(const Vector<UndirectedGraph> & uGraphs, const UndirectedGraph & uGraph, int size);

#endif // FIND_ISOMORPHIC_INDEX_HPP