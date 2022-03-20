#ifndef CUDA_ARRAY_H_
#pragma once
#include "CudaDevPtr.h"
#include "CudaMath.h"
#include "CudaMemory.h"
#include "CudaAsync.h"
/**
 * \brief       CUDA 1D array object
 */
template <typename Type> class cCudaArray
{
public:
    cCudaArray(){};
    ~cCudaArray(){};
    using size_type = size_t;

    //!	@brief	Resize array. If count == 0, memory will dealocated automaticly.
    void Resize(size_type Count) { mMemory.Allocate(Count * sizeof(Type)); }

    //!	@brief	Memory copy form device to host.
    void Download(Type *pHost, size_type Count) const
    {
        mMemory.ReadFromDevice(pHost, Count * sizeof(Type), 0);
    }

    //!	@brief	Memory copy form device to host vector.
    void Download(std::vector<Type> &HostVec) const
    {
        size_t num_of_ele = this->Size();
        HostVec.resize(num_of_ele);

        mMemory.ReadFromDevice(HostVec.data(), this->Size() * sizeof(Type), 0);
    }

    //!	@brief	Memory copy form host to device.
    void Upload(const Type *pHost, size_type Count)
    {
        if (this->Size() < Count)
        {
            this->Resize(Count);
#ifdef ENABLE_CUDA_MEMORY_CHECK
            printf("Device Memory Not Yet Allocated!\n");
#endif
        }

#ifdef ENABLE_CUDA_MEMORY_CHECK
        for (size_type i = 0; i < Count; i++)
        {
            if (cCudaMath::IsNan(pHost[i]))
            {
                printf("Uploading invalid value!");

                break;
            }
        }
#endif

        mMemory.WriteToDevice(pHost, Count * sizeof(Type), 0);
    }

    //!	@brief	Memory copy form host vector to device.
    void Upload(const std::vector<Type> &HostVec)
    {
        mMemory.Allocate(HostVec.size() * sizeof(Type));

        this->Upload(HostVec.data(), (size_type)HostVec.size());
    }

    //!	@brief	Fill memory with Value.
    // void MemsetAsync(Type Value, size_type Offset, size_type Count)
    // {
    //     if (Offset + Count <= this->Size())
    //     {
    //         CuAsync::Memset(this->Ptr() + Offset, Value, Count);
    //     }
    // }

    // //!	@brief	Fill memory with Value.
    void MemsetAsync(Type Value)
    {
        if (!this->IsEmpty())
        {
            CudaAsync::Memset(this->Ptr(), Value, this->Size());
        }
    }

    //!	@brief	Return count of elements.
    size_type Size() const
    {
        return static_cast<size_type>(mMemory.Bytes()) / sizeof(Type);
    }

    //!	@brief	Return count of bytes.
    size_type Bytes() const { return static_cast<size_type>(mMemory.Bytes()); }

    //!	@brief	If device array is empty.
    bool IsEmpty() const { return mMemory.IsEmpty(); }

    //!	@brief	Erase all.
    void Clear() { mMemory.Free(); }

#ifdef ENABLE_CUDA_MEMORY_CHECK

    //!	@brief	Return constant address of device memory.
    devPtr<const Type> Ptr() const
    {
        // static_assert(false);
        return devPtr<const Type>((const Type *)mMemory.Ptr(), this->Size());
    }

    //!	@brief	Return address of device memory.
    devPtr<Type> Ptr()
    {
        // static_assert(false);
        return devPtr<Type>((Type *)mMemory.Ptr(), this->Size());
    }

#else

    //!	@brief	Return constant address of device memory.
    devPtr<const Type> Ptr() const { return devPtr<const Type>(mMemory.Ptr()); }

    //!	@brief	Return address of device memory.
    devPtr<Type> Ptr() { return devPtr<Type>(mMemory.Ptr()); }
#endif

protected:
    cCudaMemory mMemory;
};
#endif