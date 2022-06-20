#include "CudaMemory.h"
#include "utils/DefUtil.h"
#include "utils/LogUtil.h"
#include <cuda_runtime_api.h>
size_t cCudaMemory::mTotalNumOfBytes = 0;
cCudaMemory::PFN_FreeMemory cCudaMemory::mTotalPFNFree =
    nullptr; // used to inform upper application, we will free memory
cCudaMemory::PFN_RequestMemory cCudaMemory::mTotalPFNRequeset =
    nullptr; // used to inform upper app, we will request memory
cCudaMemory::cCudaMemory()
{
    mpDevData = nullptr;
    mNumOfBytes = 0;
}

cCudaMemory::~cCudaMemory() {}

size_t cCudaMemory::CalcTotalBytes() { return mTotalNumOfBytes; }
void cCudaMemory::RegisterSignals(PFN_RequestMemory pfn_request,
                                  PFN_FreeMemory pfn_free)
{
    mTotalPFNRequeset = pfn_request;
    mTotalPFNFree = pfn_free;
}

void cCudaMemory::Allocate(size_t num_of_bytes)
{
    // if the requested size is different
    if (num_of_bytes != this->mNumOfBytes)
    {
        void *dev_mem = nullptr;
        // info
        if (mTotalPFNRequeset != nullptr)
        {
            cCudaMemory::mTotalPFNRequeset(num_of_bytes); //
        }

        // request new memory
        cudaError_t err = cudaMalloc(&dev_mem, num_of_bytes);
        if (err != cudaSuccess)
        {
            cudaGetLastError(); // pop last error
            SIM_ERROR(cudaGetErrorString(err));
        }

#ifdef ENABLE_CUDA_MEMORY_CHECK
        err = cudaMemset(dev_mem, -1, num_of_bytes);
        if (err != cudaSuccess)
        {
            cudaGetLastError(); // pop last error
            SIM_ERROR(cudaGetErrorString(err));
        }
#endif

        // release old memory
        this->Free();

        mNumOfBytes = num_of_bytes;
        mpDevData = dev_mem;
        mTotalNumOfBytes += num_of_bytes;
    }
}

/**
 * \brief           read device memory to host
 */
bool cCudaMemory::ReadFromDevice(void *p_host, size_t num_of_bytes,
                                 size_t offset_bytes) const
{
    cudaError_t err =
        cudaMemcpy(p_host, (void *)((size_t)mpDevData + offset_bytes),
                   num_of_bytes, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        SIM_ERROR(cudaGetErrorString(err));
        cudaGetLastError();
        return false;
    }
    else
    {
        return true;
    }
}

/**
 * \brief           write host memory to device.
 */
bool cCudaMemory::WriteToDevice(const void *p_host, size_t num_of_bytes,
                                size_t offset_bytes)
{
    cudaError_t err = cudaMemcpy((void *)((size_t)mpDevData + offset_bytes),
                                 p_host, num_of_bytes, cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        SIM_ERROR(cudaGetErrorString(err));
        cudaGetLastError();
        return false;
    }
    else
    {
        return true;
    }
}

bool cCudaMemory::IsEmpty() const { return mpDevData == nullptr; }
const void *cCudaMemory::Ptr() const { return this->mpDevData; }
void *cCudaMemory::Ptr() { return this->mpDevData; }

/**
 * \brief           swap two cuda memory, just their:
 *      1. memory pointer
 *      2. counting
 */
void cCudaMemory::Swap(cCudaMemory &rhs)
{
    SIM_SWAP(mpDevData, rhs.mpDevData);
    SIM_SWAP(mNumOfBytes, rhs.mNumOfBytes);
}

void cCudaMemory::Free()
{
    if (mpDevData != nullptr)
    {
        cudaError_t err = cudaSuccess;

        if (mTotalPFNFree != nullptr)
        {
            mTotalPFNFree(mNumOfBytes);
        }

        err = cudaFree(mpDevData);

        if (err != cudaSuccess)
        {
            SIM_ERROR(cudaGetErrorString(err));
            cudaGetLastError(); // pop last error
        }

        mTotalNumOfBytes -= mNumOfBytes;
        mpDevData = nullptr;
        mNumOfBytes = 0;
    }
}

size_t cCudaMemory::Bytes() const { return mNumOfBytes; }

bool cCudaMemory::Read(void *pHost, size_t SizeBytes, size_t OffsetBytes) const
{
    cudaError_t err =
        cudaMemcpy(pHost, (void *)((size_t)mpDevData + OffsetBytes), SizeBytes,
                   cudaMemcpyDeviceToHost);

    if (err != cudaSuccess)
    {
        printf("SDCudaMemory::Read() => %s.", cudaGetErrorString(err));

        cudaGetLastError(); //!	pop last error

        return false;
    }

    return true;
}

bool cCudaMemory::Write(const void *pHost, size_t SizeBytes, size_t OffsetBytes)
{
    cudaError_t err = cudaMemcpy((void *)((size_t)mpDevData + OffsetBytes),
                                 pHost, SizeBytes, cudaMemcpyHostToDevice);

    if (err != cudaSuccess)
    {
        printf("SDCudaMemory::Write() => %s.", cudaGetErrorString(err));

        cudaGetLastError(); //!	pop last error

        return false;
    }

    return true;
}
