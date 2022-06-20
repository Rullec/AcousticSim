#ifndef CUDA_INTRINSI_H_
#define CUDA_INTRINSI_H_

#ifndef __CUDACC__
#include <atomic>
#include <intrin.h>
#endif
#ifdef __CUDACC__
#include <cuda_fp16.h>
#endif
#include <cuda.h>

namespace cCudaIntrinsic
{

#ifdef __CUDACC__
template <typename Type> __device__ Type AtomicOr(Type *Address, Type Val)
{
    return atomicOr(Address, Val);
}
template <typename Type> __device__ Type AtomicXor(Type *Address, Type Val)
{
    return atomicXor(Address, Val);
}
template <typename Type> __device__ Type AtomicAnd(Type *Address, Type Val)
{
    return atomicAnd(Address, Val);
}
template <typename Type> __device__ Type AtomicAdd(Type *Address, Type Val)
{
    return atomicAdd(Address, Val);
}

template <typename Type> __device__ Type AtomicSub(Type *Address, Type Val)
{
    return atomicSub(Address, Val);
}
template <typename Type> __device__ Type AtomicExch(Type *Address, Type Val)
{
    return atomicExch(Address, Val);
}
template <typename Type>
__device__ Type AtomicCAS(Type *Address, Type Exp, Type Val)
{
    return atomicCAS(Address, Exp, Val);
}
__inline__ __device__ void AtomicAdd(tCudaMatrix3f *address,
                                     tCudaMatrix3f value)
{
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            atomicAdd(&(*address)(i, j), value(i, j));
}
__inline__ __device__ void AtomicAdd(float *address, float value)
{
    atomicAdd(address, value);
}
__inline__ __device__ void AtomicAdd(tCudaVector3f *address,
                                     tCudaVector3f value)
{
    atomicAdd(&((*address)[0]), value[0]);
    atomicAdd(&((*address)[1]), value[1]);
    atomicAdd(&((*address)[2]), value[2]);
}
#else
template <typename Type> inline Type AtomicOr(Type *Address, Type Val)
{
    return reinterpret_cast<std::atomic<Type> *>(Address)->fetch_or(Val);
}
template <typename Type> inline Type AtomicXor(Type *Address, Type Val)
{
    return reinterpret_cast<std::atomic<Type> *>(Address)->fetch_xor(Val);
}
template <typename Type> inline Type AtomicAnd(Type *Address, Type Val)
{
    return reinterpret_cast<std::atomic<Type> *>(Address)->fetch_and(Val);
}
template <typename Type> inline Type AtomicAdd(Type *Address, Type Val)
{
    return reinterpret_cast<std::atomic<Type> *>(Address)->fetch_add(Val);
}
template <typename Type> inline Type AtomicSub(Type *Address, Type Val)
{
    return reinterpret_cast<std::atomic<Type> *>(Address)->fetch_sub(Val);
}
template <typename Type> inline Type AtomicExch(Type *Address, Type Val)
{
    return reinterpret_cast<std::atomic<Type> *>(Address)->exchange(Val);
}
template <typename Type>
inline Type AtomicCAS(Type *Address, Type Exp, Type Val)
{
    reinterpret_cast<std::atomic<Type> *>(Address)->compare_exchange_weak(Exp,
                                                                          Val);
    return Exp;
}
#endif
} // namespace cCudaIntrinsic
#endif