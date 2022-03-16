// #define ENABLE_CUDA_MEMORY_CHECK
#define SIM_CANNOT_COPY(class_name)                                            \
    class_name(const class_name &) = delete;                                   \
    class_name &operator=(const class_name &) = delete;

// ! -----------------CUDA compatibility-----------------
#if defined(__CUDACC__)

#define SIM_INLINE __inline__
#define SIM_ALIGN(n) __align__(n)
#define SIM_CUDA_CALLABLE_DEVICE __device__
#define SIM_CUDA_CALLABLE __host__ __device__
#define SIM_CUDA_CALLABLE_INLINE __host__ __device__ __forceinline__
#define SIM_CUDA_CALLABLE_UNROLL __pragma(unroll)
#define SIM_ASSERT(expression)
#else
#define SIM_INLINE __inline
#define SIM_ALIGN(n) __declspec(align(n)) //!	Only for windows platform.
#define SIM_CUDA_CALLABLE_DEVICE
#define SIM_CUDA_CALLABLE
#define SIM_CUDA_CALLABLE_INLINE __forceinline
#define SIM_CUDA_CALLABLE_UNROLL
// #define SIM_ASSERT(expression)		assert(expression)

#endif

// ! ----------------CUDA handy macro--------------------
#ifndef SIM_MAX
#define SIM_MAX(a, b) ((a > b) ? a : b)
#endif

#ifndef SIM_MIN
#define SIM_MIN(a, b) ((a > b) ? b : a)
#endif

#define CUDA_at(total, pre_block)                                              \
<<<SIM_MAX(1, ((int)total + pre_block - 1) / pre_block), pre_block>>>

#define CUDA_for(thread_index, total)                                          \
    CUDA_function;                                                             \
    const auto thread_index = blockDim.x * blockIdx.x + threadIdx.x;           \
    if (thread_index >= total)                                                 \
        return;

#define CUDA_ERR(func)                                                         \
    {                                                                          \
        cudaDeviceSynchronize();                                               \
        cudaError_t err = cudaGetLastError();                                  \
        if (err != cudaSuccess)                                                \
        {                                                                      \
            printf("CUDA_Error: func [%s][%d]: %s.\n", func, __LINE__,         \
                   cudaGetErrorString(err));                                   \
            throw(err);                                                        \
        }                                                                      \
    }