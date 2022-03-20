#ifndef CUDA_2D_ARRAY_H_
#define CUDA_2D_ARRAY_H_
#include "CudaMemory.h"
#include "sim/gpu_utils/CudaAsync.h"
#include "sim/gpu_utils/CudaDevPtr.h"
template <typename Type> class cCuda2DArray
{

public:
    cCuda2DArray() : m_Rows(0), m_Columns(0) {}

    ~cCuda2DArray() {}

public:
    using size_type = unsigned int;

    void Resize(size_type _Rows, size_type _Columns)
    {
        m_DevMemory.Allocate(_Rows * _Columns * sizeof(Type));

        m_Columns = _Columns;

        m_Rows = _Rows;
    }

    void Download(cCuda2DArray<Type> &hostArray2D) const
    {
        hostArray2D.Resize(m_Rows, m_Columns);

        m_DevMemory.Read(hostArray2D.Ptr(), m_Rows * m_Columns * sizeof(Type),
                         0);
    }

    void Download(std::vector<Type> &hostArray) const
    {
        hostArray.resize(Size());
        m_DevMemory.Read(hostArray.data(), Size() * sizeof(Type), 0);
    }
    void Upload(const cCuda2DArray<Type> &hostArray2D)
    {
        this->Resize((size_t)(hostArray2D.Rows()),
                     (size_t)(hostArray2D.Columns()));

        m_DevMemory.Write(hostArray2D.Ptr(), m_Rows * m_Columns * sizeof(Type),
                          0);
    }

    void MemsetAsync(Type Value, size_type Offset, size_type Count)
    {
        if (Offset + Count <= this->Size())
        {
            CudaAsync::Memset(this->Ptr()[0] + Offset, Value, Count);
        }
    }

    void MemsetAsync(Type Value)
    {
        CudaAsync::Memset(this->Ptr(), Value, this->Size());
    }

    devPtr2<Type> Ptr()
    {
        return devPtr2<Type>(static_cast<Type *>(m_DevMemory.Ptr()), m_Rows,
                             m_Columns);
    }

    devPtr2<const Type> Ptr() const
    {
        return devPtr2<const Type>(static_cast<const Type *>(m_DevMemory.Ptr()),
                                   m_Rows, m_Columns);
    }

    size_type Bytes() const
    {
        return static_cast<size_type>(m_DevMemory.Bytes());
    }

    bool IsEmpty() const { return m_DevMemory.IsEmpty(); }

    size_type Size() const { return m_Rows * m_Columns; }

    size_type Columns() const { return m_Columns; }

    size_type Rows() const { return m_Rows; }

    void Clear() noexcept
    {
        m_DevMemory.Free();

        m_Columns = 0;

        m_Rows = 0;
    }

private:
    size_type m_Rows;

    size_type m_Columns;

    cCudaMemory m_DevMemory;
};

#endif