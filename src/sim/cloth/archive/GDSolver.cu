#include "sim/gpu_utils/Cuda2DArray.h"
#include "sim/gpu_utils/CudaArray.h"
#include "sim/gpu_utils/CudaDevPtr.h"
#include "sim/gpu_utils/CudaIntrinsic.h"
#include "sim/gpu_utils/CudaMatrix.h"
#include <cassert>
#include <iostream>
#include <map>

// r = b - Ax
extern __global__ void CalcResidual(
    int num_of_v,
    devPtr<const tCudaVector32i> ELL_local_vertex_id_to_global_vertex_id,
    devPtr2<const tCudaMatrix3f> A, devPtr<const tCudaVector3f> b,
    devPtr<const tCudaVector3f> x, devPtr<tCudaVector3f> r);
// {
//     CUDA_function;
//     int v_id = threadIdx.x + blockDim.x * blockIdx.x;
//     if (v_id >= num_of_v)
//     {
//         return;
//     }
//     tCudaVector32i ell_local_to_global_id =
//         ELL_local_vertex_id_to_global_vertex_id[v_id];
//     tCudaVector3f cur_sum_vec = tCudaVector3f::Zero();
//     for (int i = 0; i < ell_local_to_global_id.size(); i++)
//     {
//         int v_column_global_id = ell_local_to_global_id[i];
//         if (v_column_global_id == -1)
//             break;

//         cur_sum_vec += A[v_id][i] * x[v_column_global_id];
//     }

//     r[v_id] = b[v_id] - cur_sum_vec;
//     // printf("[calc_r] v %d, b %.1e, %.1e, %.1e, r %.1e, %.1e, %.1e\n", v_id,
//     //        b[v_id][0], b[v_id][1], b[v_id][2], r[v_id][0], r[v_id][1],
//     //        r[v_id][2]);
// }

// d = 2 * AT * r = 2 * A * r, A is symmetric
__global__ void CalcDirection(
    int num_of_v, devPtr<const float> alpha,
    devPtr<const tCudaVector32i> ELL_local_vertex_id_to_global_vertex_id,
    devPtr2<const tCudaMatrix3f> A, devPtr<const tCudaVector3f> r,
    devPtr<tCudaVector3f> x)
{
    CUDA_function;
    int v_id = threadIdx.x + blockDim.x * blockIdx.x;
    if (v_id >= num_of_v)
    {
        return;
    }
    // find the v_id th column
    tCudaVector32i ell_local_to_global_id =
        ELL_local_vertex_id_to_global_vertex_id[v_id];
    tCudaVector3f dx = tCudaVector3f::Zero();
    for (int i = 0; i < ell_local_to_global_id.size(); i++)
    {
        int v_column_global_id = ell_local_to_global_id[i];
        if (v_column_global_id != -1)
            dx += A[v_id][i] * r[v_column_global_id];
    }

    dx *= alpha[0];
    // printf("v%d, cur x %.1e,%.1e,%.1e dx %.1e,%.1e,%.1e\n", v_id, x[v_id][0],
    //        x[v_id][1], x[v_id][2], dx[0], dx[1], dx[2]);
    x[v_id] += dx;
}

namespace GDSolver
{
void Solve(
    const cCudaArray<tCudaVector32i> &ELL_local_vertex_id_to_global_vertex_id,
    const cCuda2DArray<tCudaMatrix3f> &A, const cCudaArray<tCudaVector3f> &b,
    cCudaArray<tCudaVector3f> &x, cCudaArray<tCudaVector3f> &r_buf)
{
    cCudaArray<float> lr;
    // lr.Upload({1e-6});
    lr.Upload({10});
    // x.MemsetAsync(tCudaVector3f::Zero());

    int max_iters = 50;
    int num_of_v = b.Size();
    std::vector<tCudaVector3f> r_cpu(x.Size());
    for (int i = 0; i < max_iters; i++)
    {
        // 1. calculate r = b - Ax
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
            printf("[log_gd] cur iter %d max residual %.1e\n", i, max_r);
        }
        // check the residual, judge convergence
        CalcDirection CUDA_at(num_of_v, 128)(
            num_of_v, lr.Ptr(), ELL_local_vertex_id_to_global_vertex_id.Ptr(),
            A.Ptr(), r_buf.Ptr(), x.Ptr());
        /*
            2. calcualte d = 2 * AT * r = 2 * rT * A,
                then x = x + alpha * d
        */
    }
}
} // namespace GDSolver