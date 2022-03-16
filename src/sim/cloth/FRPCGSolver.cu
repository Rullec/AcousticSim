#include "sim/gpu_utils/Cuda2DArray.h"
#include "sim/gpu_utils/CudaArray.h"
#include "sim/gpu_utils/CudaDef.h"
#include "sim/gpu_utils/CudaDevPtr.h"
#include "sim/gpu_utils/CudaIntrinsic.h"
#include "sim/gpu_utils/CudaMatrix.h"

__global__ void ComputeResidual(
    int num_of_v,
    devPtr<const tCudaVector32i> ELL_local_vertex_id_to_global_vertex_id,
    devPtr2<const tCudaMatrix3f> A, devPtr<const tCudaVector3f> b,
    devPtr<const tCudaVector3f> x0, devPtr<tCudaVector3f> r0)
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
}

__global__ void ApplyPreconditioner(int num_of_v,
                                    devPtr<const tCudaMatrix3f> Minv,
                                    devPtr<const tCudaVector3f> raw_vector,
                                    devPtr<tCudaVector3f> new_vector)
{
    CUDA_function;
    int v_id = threadIdx.x + blockIdx.x * blockDim.x;
    if (v_id >= num_of_v)
    {
        return;
    }

    // 1. calculate r0[v_id]
    new_vector[v_id] = Minv[v_id] * raw_vector[v_id];
}
__global__ void ComputeRdotZ(int num_of_v,

                             devPtr<const tCudaVector3f> r,
                             devPtr<const tCudaVector3f> z, devPtr<float> sum)
{
    CUDA_function;
    int v_id = threadIdx.x + blockIdx.x * blockDim.x;
    if (v_id >= num_of_v)
    {
        return;
    }

    // 1. calculate r0[v_id]

    cCudaIntrinsic::AtomicAdd((float *)(sum), r[v_id].dot(z[v_id]));
}

__global__ void UpdateDirection(int num_of_v, devPtr<tCudaVector3f> z,
                                devPtr<const float> rdotz,
                                devPtr<tCudaVector3f> p)
{
    CUDA_function;
    int v_id = threadIdx.x + blockIdx.x * blockDim.x;
    if (v_id >= num_of_v)
    {
        return;
    }
    if (rdotz[0] == 0.f)
    {

        p[v_id] = z[v_id];
    }
    else
    {

        p[v_id] = z[v_id] + rdotz[1] / rdotz[0] * p[v_id];
    }
}

__global__ void
Ap_and_pAp(int num_of_v,
           devPtr<const tCudaVector32i> ELL_local_vertex_id_to_global_vertex_id,
           devPtr2<const tCudaMatrix3f> A, devPtr<const tCudaVector3f> p,
           devPtr<tCudaVector3f> Ap, devPtr<float> pAp)
{
    CUDA_function;
    int v_id = threadIdx.x + blockIdx.x * blockDim.x;
    if (v_id >= num_of_v)
    {
        return;
    }

    // 1. calculate Ap
    tCudaVector3f cur_Ap = tCudaVector3f::Zero();
    const tCudaVector32i &local_to_global =
        ELL_local_vertex_id_to_global_vertex_id[v_id];
    for (int i = 0; i < local_to_global.size(); i++)
    {
        int cur_column_v_global_id = local_to_global[i];
        if (cur_column_v_global_id == -1)
            break;

        cur_Ap += A[v_id][i] * p[cur_column_v_global_id];
    }
    Ap[v_id] = cur_Ap;
    cCudaIntrinsic::AtomicAdd((float *)pAp, p[v_id].dot(cur_Ap));
}

// alpha = rdotz / pAp; x = x + alpha * p; r = r - alpha * Ap
__global__ void StepCG(int num_of_v, devPtr<const float> rdotz,
                       devPtr<const float> pAp, devPtr<const tCudaVector3f> p,
                       devPtr<const tCudaVector3f> Ap, devPtr<tCudaVector3f> x,
                       devPtr<tCudaVector3f> r)
{
    CUDA_function;
    int v_id = threadIdx.x + blockIdx.x * blockDim.x;
    if (v_id >= num_of_v)
    {
        return;
    }
    float alpha = rdotz[0] / pAp[0];
    x[v_id] = x[v_id] + alpha * p[v_id];
    r[v_id] = r[v_id] - alpha * Ap[v_id];
}

// E(x) = 0.5 xT * A * x - b^T * x
__global__ void CollectCGEnergy(
    int num_of_v,
    devPtr<const tCudaVector32i> ELL_local_vertex_id_to_global_vertex_id,
    devPtr2<const tCudaMatrix3f> A, devPtr<const tCudaVector3f> b,
    devPtr<const tCudaVector3f> x, devPtr<float> energy)
{
    CUDA_function;
    int v_id = threadIdx.x + blockIdx.x * blockDim.x;
    if (v_id >= num_of_v)
    {
        return;
    }

    // 1. calculate A x
    tCudaVector3f Ax = tCudaVector3f::Zero();
    tCudaVector32i local_to_global_id =
        ELL_local_vertex_id_to_global_vertex_id[v_id];
    for (int i = 0; i < local_to_global_id.size(); i++)
    {
        int v_column_global_id = local_to_global_id[i];
        if (v_column_global_id == -1)
            break;
        Ax += A[v_id][i] * x[v_column_global_id];
    }
    float sum = 0.5 * x[v_id].dot(Ax);
    sum -= b[v_id].dot(x[v_id]);
    // 2. add b.dot(x)
    cCudaIntrinsic::AtomicAdd((float *)energy, sum);
}
namespace FRPCGSolver
{

void Solve(
    const cCudaArray<tCudaMatrix3f> &Winv_preconditioner,
    const cCudaArray<tCudaVector32i> &ELL_local_vertex_id_to_global_vertex_id,
    const cCuda2DArray<tCudaMatrix3f> &A, const cCudaArray<tCudaVector3f> &b,
    cCudaArray<tCudaVector3f> &x, cCudaArray<float> &debugEnergy,
    cCudaArray<float> &pap, cCudaArray<float> &rdotZ,
    cCudaArray<tCudaVector3f> &p, cCudaArray<tCudaVector3f> &Ap,
    cCudaArray<tCudaVector3f> &residual, cCudaArray<tCudaVector3f> &z)
{
    int max_iters = 50;
    int num_of_v = b.Size();
    debugEnergy.Resize(max_iters);
    debugEnergy.MemsetAsync(0);
    pap.Resize(max_iters + 1);
    pap.MemsetAsync(0);
    rdotZ.Resize(max_iters + 1);
    rdotZ.MemsetAsync(0);
    residual.Resize(num_of_v);
    z.Resize(num_of_v);

    // 1. p = 0
    p.MemsetAsync(tCudaVector3f::Zero());

    // 2. r = b - A x
    ComputeResidual CUDA_at(num_of_v, 128)(
        num_of_v, ELL_local_vertex_id_to_global_vertex_id.Ptr(), A.Ptr(),
        b.Ptr(), x.Ptr(), residual.Ptr());
    // CUDA_ERR("ComputeResidual");
    std::vector<float> energy_cpu;
    CollectCGEnergy CUDA_at(num_of_v, 128)(
        num_of_v, ELL_local_vertex_id_to_global_vertex_id.Ptr(), A.Ptr(),
        b.Ptr(), x.Ptr(), debugEnergy.Ptr());
    // CUDA_ERR("CG energy");
    debugEnergy.Download(energy_cpu);
    std::cout << "init energy = " << energy_cpu[0] << std::endl;
    for (int iter = 0; iter < max_iters; iter++)
    {
        // 3. z = Minv * r
        // std::cout << "Winv size = " << Winv_preconditioner.Size() <<
        // std::endl; std::cout << "residual size = " << residual.Size() <<
        // std::endl; std::cout << "z size = " << z.Size() << std::endl;
        ApplyPreconditioner CUDA_at(num_of_v, 128)(
            num_of_v, Winv_preconditioner.Ptr(), residual.Ptr(), z.Ptr());
        // CUDA_ERR("ApplyPreconditioner");
        // 4. rdotz
        ComputeRdotZ CUDA_at(num_of_v, 128)(num_of_v, z.Ptr(), residual.Ptr(),
                                            rdotZ.Ptr() + iter + 1);
        // CUDA_ERR("ComputeRdotZ");

        // 5. p = z + rdotz_{k-1} / rdotz_{k-2} * p
        UpdateDirection CUDA_at(num_of_v, 128)(num_of_v, z.Ptr(),
                                               rdotZ.Ptr() + iter, p.Ptr());
        // CUDA_ERR("UpdateDirection");
        // 6. pAp, Ap
        Ap_and_pAp CUDA_at(num_of_v, 128)(
            num_of_v, ELL_local_vertex_id_to_global_vertex_id.Ptr(), A.Ptr(),
            p.Ptr(), Ap.Ptr(), pap.Ptr() + iter + 1);
        // CUDA_ERR("Ap_and_pAp");
        // 7. alpha = rdotz / pAp; x = x + alpha * p; r = r - alpha * Ap
        StepCG CUDA_at(num_of_v, 128)(num_of_v, rdotZ.Ptr() + iter + 1,
                                      pap.Ptr() + iter + 1, p.Ptr(), Ap.Ptr(),
                                      x.Ptr(), residual.Ptr());
        // CUDA_ERR("StepCG");

        CollectCGEnergy CUDA_at(num_of_v, 128)(
            num_of_v, ELL_local_vertex_id_to_global_vertex_id.Ptr(), A.Ptr(),
            b.Ptr(), x.Ptr(), debugEnergy.Ptr() + iter + 1);
        // CUDA_ERR("CG energy");
        debugEnergy.Download(energy_cpu);

        if (iter % 10 == 0)
            std::cout << "iter " << iter << " energy = " << energy_cpu[iter + 1]
                      << std::endl;
        if (std::isnan(energy_cpu[iter + 1]))
        {
            exit(1);
        }
    }
}
}; // namespace FRPCGSolver