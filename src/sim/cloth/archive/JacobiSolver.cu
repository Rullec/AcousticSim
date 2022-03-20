#include "sim/gpu_utils/Cuda2DArray.h"
#include "sim/gpu_utils/CudaArray.h"
#include "sim/gpu_utils/CudaDevPtr.h"
#include "sim/gpu_utils/CudaIntrinsic.h"
#include "sim/gpu_utils/CudaMatrix.h"
#include <cassert>
#include <iostream>
#include <map>
__global__ void CalcResidual(
    int num_of_v,
    devPtr<const tCudaVector32i> ELL_local_vertex_id_to_global_vertex_id,
    devPtr2<const tCudaMatrix3f> A, devPtr<const tCudaVector3f> b,
    devPtr<const tCudaVector3f> x, devPtr<tCudaVector3f> r)
{
    CUDA_function;
    int v_id = threadIdx.x + blockDim.x * blockIdx.x;
    if (v_id >= num_of_v)
    {
        return;
    }
    tCudaVector32i ell_local_to_global_id =
        ELL_local_vertex_id_to_global_vertex_id[v_id];
    tCudaVector3f cur_sum_vec = tCudaVector3f::Zero();
    for (int i = 0; i < ell_local_to_global_id.size(); i++)
    {
        int v_column_global_id = ell_local_to_global_id[i];
        if (v_column_global_id == -1)
            break;
        // if (v_id != v_column_global_id)
        cur_sum_vec += A[v_id][i] * x[v_column_global_id];
    }

    r[v_id] = b[v_id] - cur_sum_vec;
    // printf("[calc_r] v %d, b %.1e, %.1e, %.1e, r %.1e, %.1e, %.1e\n", v_id,
    //        b[v_id][0], b[v_id][1], b[v_id][2], r[v_id][0], r[v_id][1],
    //        r[v_id][2]);
}

__global__ void CalcUpdateVec(
    int num_of_v,
    devPtr<const tCudaVector32i> ELL_local_vertex_id_to_global_vertex_id,
    devPtr2<const tCudaMatrix3f> A, devPtr<const tCudaVector3f> b,
    devPtr<const tCudaVector3f> x, devPtr<tCudaVector3f> update_dir)
{
    // 1. calculate
    CUDA_function;
    int v_id = threadIdx.x + blockDim.x * blockIdx.x;
    if (v_id >= num_of_v)
    {
        return;
    }
    tCudaVector32i ell_local_to_global_id =
        ELL_local_vertex_id_to_global_vertex_id[v_id];
    tCudaVector3f b_minu_Rx = tCudaVector3f::Zero();
    tCudaVector3f cur_diag_inv = tCudaVector3f::Zero();
    for (int i = 0; i < ell_local_to_global_id.size(); i++)
    {
        int v_column_global_id = ell_local_to_global_id[i];
        if (v_column_global_id == -1)
            break;
        tCudaMatrix3f cur_block = A[v_id][i];
        if (v_column_global_id == v_id)
        {
            cur_diag_inv[0] = 1.0 / cur_block(0, 0);
            cur_diag_inv[1] = 1.0 / cur_block(1, 1);
            cur_diag_inv[2] = 1.0 / cur_block(2, 2);
            cur_block(0, 0) = 0;
            cur_block(1, 1) = 0;
            cur_block(2, 2) = 0;
        }

        b_minu_Rx -= cur_block * x[v_column_global_id];
    }
    b_minu_Rx += b[v_id];
    b_minu_Rx[0] *= cur_diag_inv[0];
    b_minu_Rx[1] *= cur_diag_inv[1];
    b_minu_Rx[2] *= cur_diag_inv[2];
    update_dir[v_id] = b_minu_Rx;
}

__global__ void Update(int num_of_v, devPtr<const tCudaVector3f> dir,
                       devPtr<tCudaVector3f> x)
{
    CUDA_function;
    int v_id = threadIdx.x + blockDim.x * blockIdx.x;
    if (v_id >= num_of_v)
    {
        return;
    }
    float w = 0.67;
    x[v_id] = w * dir[v_id] + (1 - w) * x[v_id];
}

namespace JacobSolver
{
void Solve(
    const cCudaArray<tCudaVector32i> &ELL_local_vertex_id_to_global_vertex_id,
    const cCuda2DArray<tCudaMatrix3f> &A, const cCudaArray<tCudaVector3f> &b,
    cCudaArray<tCudaVector3f> &x, cCudaArray<tCudaVector3f> &r_buf)
{
    int max_iter = 100;
    int num_of_v = x.Size();
    std::vector<tCudaVector3f> r_cpu(x.Size());
    for (int i = 0; i < max_iter; i++)
    {

        CalcResidual CUDA_at(num_of_v, 128)(
            num_of_v, ELL_local_vertex_id_to_global_vertex_id.Ptr(), A.Ptr(),
            b.Ptr(), x.Ptr(), r_buf.Ptr());

        if (i % 10 == 0)
        {
            r_buf.Download(r_cpu);
            float max_r = -1;
            for (int i = 0; i < x.Size(); i++)
            {
                max_r = std::max(r_cpu[i].cwiseAbs().maxCoeff(), max_r);
            }
            printf("[log] cur iter %d max residual %.1e\n", i, max_r);
        }
        CalcUpdateVec CUDA_at(num_of_v, 128)(
            num_of_v, ELL_local_vertex_id_to_global_vertex_id.Ptr(), A.Ptr(),
            b.Ptr(), x.Ptr(), r_buf.Ptr());
        Update CUDA_at(num_of_v, 128)(num_of_v, r_buf.Ptr(), x.Ptr());
    }
}
} // namespace JacobSolver