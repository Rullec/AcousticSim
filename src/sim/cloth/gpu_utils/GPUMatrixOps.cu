#include "sim/cloth/gpu_utils/Cuda2DArray.h"
#include "sim/cloth/gpu_utils/CudaArray.h"
#include "sim/cloth/gpu_utils/CudaMatrix.h"
namespace GPUMatrixOps
{
__global__ void ELLMatrixAddKernel(int num_of_items,
                                   devPtr2<const tCudaMatrix3f> mat0,
                                   devPtr2<const tCudaMatrix3f> mat1,
                                   devPtr2<tCudaMatrix3f> mat2)
{
    CUDA_function;
    int item_idx = threadIdx.x + blockDim.x * blockIdx.x;
    if (item_idx >= num_of_items)
    {
        return;
    }
    // printf("[error] matrix add is ignored\n");
    // return;
    for (int i = 0; i < mat2.Columns(); i++)
    {
        // printf("try to handle %d %d\n", item_idx, i);
        mat2[item_idx][i] = mat0[item_idx][i] + mat1[item_idx][i];
    }
}

__global__ void VectorAddKernel(int num_of_items,
                                devPtr<const tCudaVector3f> vec0,
                                devPtr<const tCudaVector3f> vec1,
                                devPtr<tCudaVector3f> vec2)
{
    CUDA_function;
    int item_idx = threadIdx.x + blockDim.x * blockIdx.x;
    if (item_idx >= num_of_items)
    {
        return;
    }
    // printf("[vec] begin to handle %d\n", item_idx);
    vec2[item_idx] = vec0[item_idx] + vec1[item_idx];
    // vec0[item_idx];
    // vec1[item_idx];
    // vec2[item_idx];
}

void ELLMatrixAdd(const cCuda2DArray<tCudaMatrix3f> &mat0,
                  const cCuda2DArray<tCudaMatrix3f> &mat1,
                  cCuda2DArray<tCudaMatrix3f> &mat2)
{
    if (mat0.Size() != mat1.Size() || mat0.Size() != mat2.Size())
    {
        printf(
            "[error] mat0 size {} != mat1 size {} != mat2 size {}, cannot be "
            "added together",
            mat0.Size(), mat1.Size(), mat2.Size());
        exit(1);
    }
    int num_of_items = mat0.Rows();
    ELLMatrixAddKernel CUDA_at(num_of_items, 128)(num_of_items, mat0.Ptr(),
                                                  mat1.Ptr(), mat2.Ptr());
    CUDA_ERR("ell matrix add");
}

void VectorAdd(const cCudaArray<tCudaVector3f> &vec0,
               const cCudaArray<tCudaVector3f> &vec1,
               cCudaArray<tCudaVector3f> &vec2)
{
    if (vec0.Size() != vec1.Size() || vec0.Size() != vec2.Size())
    {
        printf(
            "[error] vec0 size {} != vec1 size {} != vec2 size {}, cannot be "
            "added together",
            vec0.Size(), vec1.Size(), vec2.Size());
        exit(1);
    }
    // std::cout << "[vec_add] vec0 size " << vec0.Size() << " vec1 size "
    //           << vec1.Size() << " vec2 size  " << vec2.Size() << std::endl;
    int num_of_items = vec0.Size();
    VectorAddKernel CUDA_at(num_of_items, 128)(num_of_items, vec0.Ptr(),
                                               vec1.Ptr(), vec2.Ptr());
    CUDA_ERR("vector add");
}

}; // namespace GPUMatrixOps