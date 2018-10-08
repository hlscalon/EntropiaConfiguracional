#ifndef CUDA_INCLUDE_HPP
#define CUDA_INCLUDE_HPP

#if defined(__CUDACC__)

#define CUDA_HOST_DEVICE __host__ __device__
#define CUDA_DEVICE __device__

#else

#define CUDA_HOST_DEVICE
#define CUDA_DEVICE

#endif

#endif // CUDA_INCLUDE_HPP