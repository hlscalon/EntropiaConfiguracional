#include "FindIsomorphicIndex.hpp"
#include "IsIsomorphic.hpp"

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/find.h>
#include <thrust/execution_policy.h>

int findIsomorphicIndex(const Vector<UndirectedGraph> & uGraphs, const UndirectedGraph & uGraph, int size) {
	thrust::device_vector<UndirectedGraph> dev_vec(uGraphs.begin(), uGraphs.begin() + size);

	// thrust::device_vector<UndirectedGraph> dev_vec = host_vec;

	// auto itIdx = thrust::find_if(thrust::device, dev_vec.begin(), dev_vec.end(), [uGraph] CUDA_HOST_DEVICE (const UndirectedGraph & ug) {
	// 	// return is_isomorphic(ug, uGraph);
	// 	return true;
	// });

	auto itIdx = dev_vec.begin();
	return itIdx == dev_vec.end() ? -1 : itIdx - dev_vec.begin();
}
