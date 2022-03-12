#pragma once
#include "sim/gpu_utils/CudaDevPtr.h"
#include "sim/gpu_utils/CudaMatrix.h"
class CudaAsync
{
public:
    using size_type = unsigned int;

    static void Memset(devPtr<int> pDev, int Value, size_type count);
    // static void Memset(devPtr<Int2> pDev, Int2 Value, size_type count);
    // static void Memset(devPtr<Int3> pDev, Int3 Value, size_type count);
    // static void Memset(devPtr<bool> pDev, bool Value, size_type count);
    static void Memset(devPtr<float> pDev, float Value, size_type count);
	static void Memset(devPtr<tCudaMatrix9f> pDev, tCudaMatrix9f value, size_type count);
	static void Memset(devPtr<tCudaVector9f> pDev, tCudaVector9f value, size_type count);
	static void Memset(devPtr<tCudaMatrix3f> pDev, tCudaMatrix3f value, size_type count);
	static void Memset(devPtr<tCudaVector3f> pDev, tCudaVector3f value, size_type count);
	
    // static void Memset(devPtr<double> pDev, double Value, size_type count);
    // static void Memset(devPtr<Float2> pDev, Float2 Value, size_type count);
    // static void Memset(devPtr<Float4> pDev, Float4 Value, size_type count);
    // static void Memset(devPtr<SeBool> pDev, SeBool Value, size_type count);
    // static void Memset(devPtr<CuFloat3> pDev, CuFloat3 Value, size_type count);
    // static void Memset(devPtr<IntFloat> pDev, IntFloat Value, size_type count);
    // static void Memset(devPtr<IntFloat3> pDev, IntFloat3 Value,
    //                    size_type count);
    // static void Memset(devPtr<SeMatrix3f> pDev, SeMatrix3f Value,
    //                    size_type count);
    // static void Memset(devPtr<unsigned int> pDev, unsigned int Value,
    //                    size_type count);
    // static void Memset(devPtr<SeFrictionInfo> pDev, SeFrictionInfo Value,
    //                    size_type count);

    // static void Memcpy(devPtr<int> pDst, devPtr<const int> pSrc,
    //                    size_type count);
    // static void Memcpy(devPtr<Int4> pDst, devPtr<const Int4> pSrc,
    //                    size_type count);
    // static void Memcpy(devPtr<Int2> pDst, devPtr<const Int2> pSrc,
    //                    size_type count);
    // static void Memcpy(devPtr<float> pDst, devPtr<const float> pSrc,
    //                    size_type count);
    // static void Memcpy(devPtr<SePin> pDst, devPtr<const SePin> pSrc,
    //                    size_type count);
    // static void Memcpy(devPtr<Float2> pDst, devPtr<const Float2> pSrc,
    //                    size_type count);
    // static void Memcpy(devPtr<SeBool> pDst, devPtr<const SeBool> pSrc,
    //                    size_type count);
    // static void Memcpy(devPtr<CuFloat3> pDst, devPtr<const CuFloat3> pSrc,
    //                    size_type count);

    // static void Add(devPtr<float> pResult, devPtr<const float> pA,
    //                 devPtr<const float> pB, size_type count);
    // static void Add(devPtr<CuFloat3> pResult, devPtr<const CuFloat3> pA,
    //                 devPtr<const CuFloat3> pB, size_type count);

    // static void Sub(devPtr<float> pResult, devPtr<const float> pMinuend,
    //                 devPtr<const float> pSubtrahend, size_type count);
    // static void Sub(devPtr<CuFloat3> pResult, devPtr<const CuFloat3> pMinuend,
    //                 devPtr<const CuFloat3> pSubtrahend, size_type count);
};