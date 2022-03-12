#include "sim/gpu_utils/Cuda2DArray.h"
#include "sim/gpu_utils/CudaArray.h"
#include "sim/gpu_utils/CudaDevPtr.h"
#include "sim/gpu_utils/CudaIntrinsic.h"
#include "sim/gpu_utils/CudaMatrix.h"
#include <cassert>
#include <iostream>
#include <map>

__global__ void AtomicAdd(int num_of_v, devPtr<const float> data_array,
                          devPtr<float> sum)
{
    CUDA_function;
    int v_id = threadIdx.x + blockDim.x * blockIdx.x;
    if (v_id >= num_of_v)
    {
        return;
    }
    float *ptr = sum;
    cCudaIntrinsic::AtomicAdd(ptr, data_array[v_id]);
}
//

__global__ void BlockJacobPreconditioner(
    int num_of_v,
    devPtr<const tCudaVector32i> ELL_local_vertex_id_to_global_vertex_id,
    devPtr2<const tCudaMatrix3f> W, devPtr<tCudaMatrix3f> PreconditionerWinv)
{
    CUDA_function;
    int v_id = threadIdx.x + blockDim.x * blockIdx.x;
    if (v_id >= num_of_v)
    {
        return;
    }

    // find the diagnoal entry in ELL matrix
    tCudaVector32i ELL_local2global_map =
        ELL_local_vertex_id_to_global_vertex_id[v_id];
    for (int i = 0; i < ELL_local2global_map.size(); i++)
    {
        if (ELL_local2global_map[i] == v_id)
        {
            PreconditionerWinv[v_id] = W[v_id][i].inverse();
            return;
        }
    }
    assert(false && "cannot find diag entry in ELL matrix");
}

__global__ void JacobPreconditioner(
    int num_of_v,
    devPtr<const tCudaVector32i> ELL_local_vertex_id_to_global_vertex_id,
    devPtr2<const tCudaMatrix3f> W, devPtr<tCudaMatrix3f> PreconditionerWinv)
{
    CUDA_function;
    int v_id = threadIdx.x + blockDim.x * blockIdx.x;
    if (v_id >= num_of_v)
    {
        return;
    }

    // find the diagnoal entry in ELL matrix
    tCudaVector32i ELL_local2global_map =
        ELL_local_vertex_id_to_global_vertex_id[v_id];
    for (int i = 0; i < ELL_local2global_map.size(); i++)
    {
        if (ELL_local2global_map[i] == v_id)
        {
            PreconditionerWinv[v_id].setZero();
            for (int k = 0; k < 3; k++)
            {
                PreconditionerWinv[v_id](k, k) = 1.0 / W[v_id][i](k, k);
            }
            return;
        }
    }
    assert(false && "cannot find diag entry in ELL matrix");
}
__global__ void NoPreconditioner(
    int num_of_v,
    devPtr<const tCudaVector32i> ELL_local_vertex_id_to_global_vertex_id,
    devPtr2<const tCudaMatrix3f> W, devPtr<tCudaMatrix3f> PreconditionerWinv)
{
    CUDA_function;
    int v_id = threadIdx.x + blockDim.x * blockIdx.x;
    if (v_id >= num_of_v)
    {
        return;
    }

    // find the diagnoal entry in ELL matrix
    tCudaVector32i ELL_local2global_map =
        ELL_local_vertex_id_to_global_vertex_id[v_id];
    for (int i = 0; i < ELL_local2global_map.size(); i++)
    {
        if (ELL_local2global_map[i] == v_id)
        {
            PreconditionerWinv[v_id].setIdentity();
            return;
        }
    }
    assert(false && "cannot find diag entry in ELL matrix");
}

/**
 * \brief       init PCG:
 *  1. calculate r0 = b - A x0
 *  2. calculate r0Minvr0
 *  3. calculate d0 = Minv * r0
 *
 */
__global__ void Calc_r0_d0_r0Minvr0(
    int num_of_v, devPtr<const tCudaVector3f> x0, devPtr<const tCudaVector3f> b,
    devPtr<tCudaVector3f> r0, devPtr<tCudaVector3f> d0, devPtr<float> r0Minvr0,
    devPtr<const tCudaMatrix3f> Minv,
    devPtr<const tCudaVector32i> ELL_local_vertex_id_to_global_vertex_id,
    devPtr2<const tCudaMatrix3f> A)
{
    CUDA_function;
    int v_id = threadIdx.x + blockIdx.x * blockDim.x;
    if (v_id >= num_of_v)
    {
        return;
    }

    // 1. calculate r0[v_id]

    const tCudaVector32i &ELL_local_to_global_id =
        ELL_local_vertex_id_to_global_vertex_id[v_id];
    tCudaVector3f r0_entry = tCudaVector3f::Zero();
    for (int i = 0; i < ELL_local_to_global_id.size(); i++)
    {
        int cur_column_v_global_id = ELL_local_to_global_id[i];
        if (cur_column_v_global_id == -1)
        {
            break;
        }

        r0_entry += A[v_id][i] * x0[cur_column_v_global_id];
    }
    r0_entry = b[v_id] - r0_entry;
    r0[v_id] = r0_entry;
    // 2. calculate d0
    d0[v_id] = Minv[v_id] * r0_entry;
    // 3. calculate r0Minvr0

    cCudaIntrinsic::AtomicAdd((float *)(r0Minvr0), r0_entry.dot(d0[v_id]));
}

/**
 * \brief           calculate:
 *      1. di.T * A * di
 *      2. zi = A * di
 */
__global__ void Calc_dTAd_zi(
    int num_of_v, devPtr<const tCudaVector3f> di,
    devPtr<const tCudaVector32i> ELL_local_vertex_id_to_global_vertex_id,
    devPtr2<const tCudaMatrix3f> A, devPtr<tCudaVector3f> zi,
    devPtr<float> diAdi)
{
    CUDA_function;
    int v_id = threadIdx.x + blockDim.x * blockIdx.x;
    if (v_id >= num_of_v)
        return;

    /*
        2. di.T * A * di

        3. zi = A * di
    */
    {
        const tCudaVector32i &ell_local_to_global_id =
            ELL_local_vertex_id_to_global_vertex_id[v_id];
        tCudaVector3f cur_sum_vec = tCudaVector3f::Zero();
        for (int i = 0; i < ell_local_to_global_id.size(); i++)
        {
            int v_column_global_id = ell_local_to_global_id[i];
            if (v_column_global_id != -1)
            {
                cur_sum_vec += A[v_id][i] * di[v_column_global_id];
            }
        }
        zi[v_id] = cur_sum_vec;
        float *sum = diAdi;
        cCudaIntrinsic::AtomicAdd(sum, di[v_id].dot(cur_sum_vec));
    }
}

/**
 * \brief
 *      alpha = riMinvri / diAdi
 *
 *      1. r_{i+1} = r_i - alpha * z_i
 *
 *      2. r_{i+1} Minv * r_{i+1}
 */
__global__ void Calc_rnext_rnextMinvrnext_xnext(
    int num_of_v, devPtr<const float> riMinvri, devPtr<const float> diAdi,
    devPtr<const tCudaVector3f> zi, devPtr<tCudaVector3f> r,
    devPtr<const tCudaMatrix3f> Minv, devPtr<float> rnextMinvrnext,
    devPtr<const tCudaVector3f> d, devPtr<tCudaVector3f> x)
{
    CUDA_function;
    int v_id = threadIdx.x + blockDim.x * blockIdx.x;
    if (v_id >= num_of_v)
    {
        return;
    }

    float alpha = riMinvri[0] / diAdi[0];
    r[v_id] -= alpha * zi[v_id];

    cCudaIntrinsic::AtomicAdd((float *)(rnextMinvrnext),
                              r[v_id].dot(Minv[v_id] * r[v_id]));
    x[v_id] += alpha * d[v_id];
}

__global__ void UpdateDirection(int num_of_v,
                                devPtr<const float> rnextMinvrnext,
                                devPtr<const float> rcurMinvrcur,
                                devPtr<const tCudaMatrix3f> Minv,
                                devPtr<const tCudaVector3f> rnext,
                                devPtr<tCudaVector3f> d)
{
    CUDA_function;
    int v_id = threadIdx.x + blockDim.x * blockIdx.x;
    if (v_id >= num_of_v)
    {
        return;
    }

    float beta = rnextMinvrnext[0] / rcurMinvrcur[0];
    d[v_id] = Minv[v_id] * rnext[v_id] + beta * d[v_id];
}

namespace PCGSolver
{
void TestAtomicAdd(const cCudaArray<float> &data_array, cCudaArray<float> &sum)
{
    int size = data_array.Size();
    AtomicAdd CUDA_at(size, 128)(size, data_array.Ptr(), sum.Ptr());
}
/**
 * \brief       solve W x = b
 */
void Solve(
    const cCudaArray<tCudaMatrix3f> &Winv_preconditioner,
    const cCudaArray<tCudaVector32i> &ELL_local_vertex_id_to_global_vertex_id,
    const cCuda2DArray<tCudaMatrix3f> &W, const cCudaArray<tCudaVector3f> &b,
    cCudaArray<tCudaVector3f> &pcg_rbuf, cCudaArray<tCudaVector3f> &pcg_dbuf,
    cCudaArray<tCudaVector3f> &pcg_zbuf, cCudaArray<float> &pcg_rMinvr_array,
    cCudaArray<tCudaVector3f> &x)
{
    // ! try to group work together
    cCudaArray<float> pcg_dTAd;

    int max_iters = 200;
    pcg_dTAd.Resize(1);
    pcg_dTAd.MemsetAsync(0);

    pcg_rbuf.MemsetAsync(tCudaVector3f::Zero());
    pcg_dbuf.MemsetAsync(tCudaVector3f::Zero());
    pcg_rMinvr_array.Resize(max_iters + 1);
    pcg_rMinvr_array.MemsetAsync(0.f);
    x.MemsetAsync(tCudaVector3f::Zero());

    int num_of_v = ELL_local_vertex_id_to_global_vertex_id.Size();

    std::vector<tCudaVector3f> x_cpu(x.Size());
    std::vector<tCudaVector3f> r_cpu(x.Size());
    std::vector<tCudaVector3f> z_cpu(x.Size());
    std::vector<tCudaVector3f> d_cpu(x.Size());
    /*
        1. calculate r0, d0

        rcur = r0 = b - A x0
        dcur = d0 = Minv r0
        r_0^T M^{-1} r_0
    */
    Calc_r0_d0_r0Minvr0 CUDA_at(num_of_v, 128)(
        num_of_v, x.Ptr(), b.Ptr(), pcg_rbuf.Ptr(), pcg_dbuf.Ptr(),
        pcg_rMinvr_array.Ptr(), Winv_preconditioner.Ptr(),
        ELL_local_vertex_id_to_global_vertex_id.Ptr(), W.Ptr());

    // printf("---------------start-------------\n");
    // x.Download(x_cpu);
    // pcg_rbuf.Download(r_cpu);
    // pcg_dbuf.Download(d_cpu);
    // for (int i = 0; i < x.Size(); i++)
    // {
    //     std::cout << "x" << i << " = " << x_cpu[i].transpose() << std::endl;
    //     std::cout << "r" << i << " = " << r_cpu[i].transpose() << std::endl;
    //     std::cout << "d" << i << " = " << d_cpu[i].transpose() << std::endl;
    // }

    for (int cur_iter = 0; cur_iter < max_iters; cur_iter++)
    {
        /*
            2. part1:
                calculate di.T * A * di
                calculate zi = A * di
        */
        pcg_dTAd.MemsetAsync(0);

        Calc_dTAd_zi CUDA_at(num_of_v,
                             128)(num_of_v, pcg_dbuf.Ptr(),
                                  ELL_local_vertex_id_to_global_vertex_id.Ptr(),
                                  W.Ptr(), pcg_zbuf.Ptr(), pcg_dTAd.Ptr());
        // {
        //     std::vector<float> pcg_dtAd_cpu;
        //     pcg_dTAd.Download(pcg_dtAd_cpu);
        //     pcg_zbuf.Download(z_cpu);
        //     printf("---------------iter %d 1st-------------\n", cur_iter);
        //     for (int i = 0; i < x.Size(); i++)
        //     {
        //         std::cout << "z" << i << " = " << z_cpu[i].transpose()
        //                   << std::endl;
        //     }
        //     assert(pcg_dtAd_cpu.size() == 1);
        //     std::cout << "dTAd = " << pcg_dtAd_cpu[0] << std::endl;
        // }
        /*
            3. part2:
                Input: riTMinvri, diTAdi
                Output: r_{i+1}, r_{i+1}^T Minv r_{i+1}, x_{i+1}
        */

        Calc_rnext_rnextMinvrnext_xnext CUDA_at(num_of_v, 128)(
            num_of_v, pcg_rMinvr_array.Ptr() + cur_iter, pcg_dTAd.Ptr(),
            pcg_zbuf.Ptr(), pcg_rbuf.Ptr(), Winv_preconditioner.Ptr(),
            pcg_rMinvr_array.Ptr() + cur_iter + 1, pcg_dbuf.Ptr(), x.Ptr());

        /*
            4. part3:
                d_{i+1} = Minv * r_next  + beta * di
        */

        UpdateDirection CUDA_at(num_of_v, 128)(
            num_of_v, pcg_rMinvr_array.Ptr() + cur_iter + 1,
            pcg_rMinvr_array.Ptr() + cur_iter, Winv_preconditioner.Ptr(),
            pcg_rbuf.Ptr(), pcg_dbuf.Ptr());

        // x.Download(x_cpu);
        if (cur_iter % 50 == 0 && cur_iter != 0)
        {
            pcg_rbuf.Download(r_cpu);
            // printf("---------------iter %d-------------\n", cur_iter);
            float max_r = -1;
            for (int i = 0; i < x.Size(); i++)
            {
                max_r = std::max(r_cpu[i].cwiseAbs().maxCoeff(), max_r);
            }
            if (max_r < 1e-7)
            {
                printf("[log] max residual %.1e, break\n", max_r);
                break;
            }
        }
        // break;
    }
    // exit(1);
}

void CalcJacobiPreconditioner(
    const cCudaArray<tCudaVector32i> &ELL_local_vertex_id_to_global_vertex_id,
    const cCuda2DArray<tCudaMatrix3f> &W,
    cCudaArray<tCudaMatrix3f> &Winv_preconditioner)
{
    int num_of_v = ELL_local_vertex_id_to_global_vertex_id.Size();
    JacobPreconditioner CUDA_at(num_of_v, 128)(
        num_of_v, ELL_local_vertex_id_to_global_vertex_id.Ptr(), W.Ptr(),
        Winv_preconditioner.Ptr());
}
void CalcBlockJacobiPreconditioner(
    const cCudaArray<tCudaVector32i> &ELL_local_vertex_id_to_global_vertex_id,
    const cCuda2DArray<tCudaMatrix3f> &W,
    cCudaArray<tCudaMatrix3f> &Winv_preconditioner)
{
    int num_of_v = ELL_local_vertex_id_to_global_vertex_id.Size();
    BlockJacobPreconditioner CUDA_at(num_of_v, 128)(
        num_of_v, ELL_local_vertex_id_to_global_vertex_id.Ptr(), W.Ptr(),
        Winv_preconditioner.Ptr());
}

void CalcNoPreconditioner(
    const cCudaArray<tCudaVector32i> &ELL_local_vertex_id_to_global_vertex_id,
    const cCuda2DArray<tCudaMatrix3f> &W,
    cCudaArray<tCudaMatrix3f> &Winv_preconditioner)
{
    int num_of_v = ELL_local_vertex_id_to_global_vertex_id.Size();
    NoPreconditioner CUDA_at(num_of_v, 128)(
        num_of_v, ELL_local_vertex_id_to_global_vertex_id.Ptr(), W.Ptr(),
        Winv_preconditioner.Ptr());
}

} // namespace PCGSolver