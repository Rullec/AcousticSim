/**
 * 公司：凌笛数码
 * 版权信息：凌笛数码所有
 * 描述：用核函数进行数据初始化和GPU-GPU内存拷贝
 * 作者：黄文超
 * 日期：2019/06/26
 */

#define DISABLE_NAN_CHECK

#include "CudaAsync.h"
// #include "SeCudaRuntime.h"

// SE_USING_NAMESPACE

namespace CuKernel
{
template <typename Type>
__global__ void MemsetAsync(devPtr<Type> pData, Type value,
                            CudaAsync::size_type count)
{
    CUDA_for(i, count);
    pData[i] = value;
}
template <typename Type>
__global__ void MemcpyAsnyc(devPtr<Type> pDst, devPtr<const Type> pSrc,
                            CudaAsync::size_type count)
{
    CUDA_for(i, count);
    pDst[i] = pSrc[i];
}
template <typename Type>
__global__ void AddAsnyc(devPtr<Type> pResult, devPtr<const Type> pA,
                         devPtr<const Type> pB, CudaAsync::size_type count)
{
    CUDA_for(i, count);
    pResult[i] = pA[i] + pB[i];
}
template <typename Type>
__global__ void SubAsnyc(devPtr<Type> pResult, devPtr<const Type> pMinuend,
                         devPtr<const Type> pSubtrahend,
                         CudaAsync::size_type count)
{
    CUDA_for(i, count);
    pResult[i] = pMinuend[i] - pSubtrahend[i];
}
} // namespace CuKernel

/*************************************************************************
******************************    Memset    ******************************
*************************************************************************/

// void CudaAsync::Memset(devPtr<unsigned int> pDev, unsigned int value,
// size_type count)
// {
// 	CuKernel::MemsetAsync CUDA_at(count, 256)(pDev, value, count);
// CUDA_ERR("memset");
// }
void CudaAsync::Memset(devPtr<int> pDev, int value, size_type count)
{
    CuKernel::MemsetAsync CUDA_at(count, 256)(pDev, value, count);
    CUDA_ERR("memset int");
}
// void CudaAsync::Memset(devPtr<Int2> pDev, Int2 Value, size_type Count)
// {
// 	CuKernel::MemsetAsync CUDA_at(Count, 256)(pDev, Value, Count);
// CUDA_ERR("memset");
// }
// void CudaAsync::Memset(devPtr<Int3> pDev, Int3 Value, size_type Count)
// {
// 	CuKernel::MemsetAsync CUDA_at(Count, 256)(pDev, Value, Count);
// CUDA_ERR("memset");
// }
// void CudaAsync::Memset(devPtr<bool> pDev, bool value, size_type count)
// {
// 	CuKernel::MemsetAsync CUDA_at(count, 256)(pDev, value, count);
// CUDA_ERR("memset");
// }
void CudaAsync::Memset(devPtr<float> pDev, float value, size_type count)
{
    CuKernel::MemsetAsync CUDA_at(count, 256)(pDev, value, count);
    CUDA_ERR("memset float");
}
// void CudaAsync::Met floatmset(devPtr<double> pDev, double value, size_type
// count)
// {
// 	CuKernel::MemsetAsync CUDA_at(count, 256)(pDev, value, count);
// CUDA_ERR("memset");
// }
// void CudaAsync::Memset(devPtr<Float2> pDev, Float2 Value, size_type Count)
// {
// 	CuKernel::MemsetAsync CUDA_at(Count, 256)(pDev, Value, Count);
// CUDA_ERR("memset");
// }
// void CudaAsync::Memset(devPtr<Float4> pDev, Float4 Value, size_type Count)
// {
// 	CuKernel::MemsetAsync CUDA_at(Count, 256)(pDev, Value, Count);
// CUDA_ERR("memset");
// }
// void CudaAsync::Memset(devPtr<SeBool> pDev, SeBool value, size_type count)
// {
// 	CuKernel::MemsetAsync CUDA_at(count, 256)(pDev, value, count);
// CUDA_ERR("memset");
// }
// void CudaAsync::Memset(devPtr<CuFloat3> pDev, CuFloat3 value, size_type
// count)
// {
// 	CuKernel::MemsetAsync CUDA_at(count, 256)(pDev, value, count);
// CUDA_ERR("memset");
// }
// void CudaAsync::Memset(devPtr<IntFloat> pDev, IntFloat value, size_type
// count)
// {
// 	CuKernel::MemsetAsync CUDA_at(count, 256)(pDev, value, count);
// CUDA_ERR("memset");
// }
// void CudaAsync::Memset(devPtr<IntFloat3> pDev, IntFloat3 value, size_type
// count)
// {
// 	CuKernel::MemsetAsync CUDA_at(count, 256)(pDev, value, count);
// CUDA_ERR("memset");
// }
void CudaAsync::Memset(devPtr<tCudaMatrix9f> pDev, tCudaMatrix9f value,
                       size_type count)
{
    CuKernel::MemsetAsync CUDA_at(count, 256)(pDev, value, count);
    CUDA_ERR("memset mat9f");
}

void CudaAsync::Memset(devPtr<tCudaMatrix3f> pDev, tCudaMatrix3f value,
                       size_type count)
{
    CuKernel::MemsetAsync CUDA_at(count, 256)(pDev, value, count);
    CUDA_ERR("memset mat3f");
}
// void CudaAsync::Memset(devPtr<SeFrictionInfo> pDev, SeFrictionInfo value,
// size_type count)
// {
// 	CuKernel::MemsetAsync CUDA_at(count, 256)(pDev, value, count);
// CUDA_ERR("memset");
// }

// /*************************************************************************
// ******************************    Memcpy    ******************************
// *************************************************************************/

// void CudaAsync::Memcpy(devPtr<int> pDst, devPtr<const int> pSrc, size_type
// count)
// {
// 	CuKernel::MemcpyAsnyc CUDA_at(count, 256)(pDst, pSrc, count);
// CUDA_ERR("memset");
// }
// void CudaAsync::Memcpy(devPtr<float> pDst, devPtr<const float> pSrc,
// size_type count)
// {
// 	CuKernel::MemcpyAsnyc CUDA_at(count, 256)(pDst, pSrc, count);
// CUDA_ERR("memset");
// }
// void CudaAsync::Memcpy(devPtr<SePin> pDst, devPtr<const SePin> pSrc,
// size_type count)
// {
// 	CuKernel::MemcpyAsnyc CUDA_at(count, 256)(pDst, pSrc, count);
// CUDA_ERR("memset");
// }
// void CudaAsync::Memcpy(devPtr<CuFloat3> pDst, devPtr<const CuFloat3> pSrc,
// size_type count)
// {
// 	CuKernel::MemcpyAsnyc CUDA_at(count, 256)(pDst, pSrc, count);
// CUDA_ERR("memset");
// }
// void CudaAsync::Memcpy(devPtr<Float2> pDst, devPtr<const Float2> pSrc,
// size_type Count)
// {
// 	CuKernel::MemcpyAsnyc CUDA_at(Count, 256)(pDst, pSrc, Count);
// CUDA_ERR("memset");
// }
// void CudaAsync::Memcpy(devPtr<Int4> pDst, devPtr<const Int4> pSrc, size_type
// count)
// {
// 	CuKernel::MemcpyAsnyc CUDA_at(count, 256)(pDst, pSrc, count);
// CUDA_ERR("memset");
// }
// void CudaAsync::Memcpy(devPtr<Int2> pDst, devPtr<const Int2> pSrc, size_type
// count)
// {
// 	CuKernel::MemcpyAsnyc CUDA_at(count, 256)(pDst, pSrc, count);
// CUDA_ERR("memset");
// }
// void CudaAsync::Memcpy(devPtr<SeBool> pDst, devPtr<const SeBool> pSrc,
// size_type count)
// {
// 	CuKernel::MemcpyAsnyc CUDA_at(count, 256)(pDst, pSrc, count);
// CUDA_ERR("memset");
// }

// /*************************************************************************
// *******************************    Add    ********************************
// *************************************************************************/

// void CudaAsync::Add(devPtr<CuFloat3> pSum, devPtr<const CuFloat3> pA,
// devPtr<const CuFloat3> pB, size_type count)
// {
// 	CuKernel::AddAsnyc CUDA_at(count, 256)(pSum, pA, pB, count);
// CUDA_ERR("memset");
// }
// void CudaAsync::Add(devPtr<float> pSum, devPtr<const float> pA, devPtr<const
// float> pB, size_type count)
// {
// 	CuKernel::AddAsnyc CUDA_at(count, 256)(pSum, pA, pB, count);
// CUDA_ERR("memset");
// }

// /*************************************************************************
// *******************************    Sub    ********************************
// *************************************************************************/

// void CudaAsync::Sub(devPtr<CuFloat3> pResult, devPtr<const CuFloat3>
// pMinuend, devPtr<const CuFloat3> pSubtrahend, size_type count)
// {
// 	CuKernel::SubAsnyc CUDA_at(count, 256)(pResult, pMinuend, pSubtrahend,
// count); CUDA_ERR("memset");
// }
// void CudaAsync::Sub(devPtr<float> pResult, devPtr<const float> pMinuend,
// devPtr<const float> pSubtrahend, size_type count)
// {
// 	CuKernel::SubAsnyc CUDA_at(count, 256)(pResult, pMinuend, pSubtrahend,
// count); CUDA_ERR("memset");
// }