#pragma once
#include "CudaDef.h"
#include <functional>
class cCudaMemory
{
    SIM_CANNOT_COPY(cCudaMemory);

public:
    using PFN_FreeMemory = std::function<void(size_t)>;
    using PFN_RequestMemory = std::function<void(size_t)>;

    cCudaMemory();
    ~cCudaMemory();

    static size_t CalcTotalBytes(); // total cuda memory bytes
    static void RegisterSignals(PFN_RequestMemory pfg_request,
                                PFN_FreeMemory pfn_free); //

    void Allocate(size_t num_of_bytes);
    bool ReadFromDevice(void *p_host, size_t num_of_bytes,
                        size_t offset_bytes) const;
    bool WriteToDevice(const void *p_host, size_t num_of_bytes,
                       size_t offset_bytes);
    /**
     *	@brief		Memory copy from device to host.
     *	@param[in]	Pointer to host data.
     *	@param[in]	Size of data in bytes.
     *	@param[in]	Device address offset.
     *	@return		Success -> true, else ->false.
     */
    bool Read(void *pHost, size_t SizeBytes, size_t OffsetBytes) const;

    /**
     *	@brief		Memory copy from host to device.
     *	@param[in]	Pointer to host data.
     *	@param[in]	Size of data in bytes.
     *	@param[in]	Device address offset.
     *	@return		Success -> true, else ->false.
     */
    bool Write(const void *pHost, size_t SizeBytes, size_t OffsetBytes);
    bool IsEmpty() const;
    size_t Bytes() const;
    const void *Ptr() const;
    void *Ptr();
    void Swap(cCudaMemory &rhs);
    void Free();

protected:
    void *mpDevData;    // pointer to device memory
    size_t mNumOfBytes; // num of bytes
private:
    static size_t mTotalNumOfBytes;
    static PFN_FreeMemory
        mTotalPFNFree; // used to inform upper application, we will free memory
    static PFN_RequestMemory
        mTotalPFNRequeset; // used to inform upper app, we will request memory
};