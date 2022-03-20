#include "sim/cloth/BaraffClothGpu.h"
#include "sim/gpu_utils/CudaDef.h"
#include "sim/gpu_utils/CudaIntrinsic.h"
namespace BaraffClothGpu
{

__device__ void CalcFuAndFv(const tCudaVector3f &v0, const tCudaVector3f &v1,
                            const tCudaVector3f &v2,
                            const tCudaVector3f &Fu_coef,
                            const tCudaVector3f &Fv_coef, tCudaVector3f &Fu,
                            tCudaVector3f &Fv)
{
    CUDA_function;
    Fu = Fu_coef[0] * v0 + Fu_coef[1] * v1 + Fu_coef[2] * v2;
    Fv = Fv_coef[0] * v0 + Fv_coef[1] * v1 + Fv_coef[2] * v2;
}

__device__ tCudaMatrix3f CalcFint(float tri_area, float Ku, float Kv, float Cu,
                                  float Cv, const tCudaVector3f &coef_u,
                                  const tCudaVector3f &coef_v,
                                  const tCudaVector3f &Fu,
                                  const tCudaVector3f &Fv)
{
    CUDA_function;
    tCudaVector3f u_base = -tri_area * Ku * Cu * Fu.normalized();
    tCudaVector3f v_base = -tri_area * Kv * Cv * Fv.normalized();
    tCudaMatrix3f fint = tCudaMatrix3f::Zero();
    // printf("[calc_fint] tri area %.1e\n", tri_area);
    // printf("[calc_fint] u_base %.1e %.1e %.1e, Ku %.1e, Cu %.1e, Fu_normed "
    //        "%.1e, %.1e, %.1e\n",
    //        u_base[0], u_base[1], u_base[2], Ku, Cu, Fu.normalized()[0],
    //        Fu.normalized()[1], Fu.normalized()[2]);
    // printf("[calc_fint] v_base %.1e %.1e %.1e\n", v_base[0], v_base[1],
    //        v_base[2]);

    for (int i = 0; i < 3; i++)
    {
        // printf("coef u %d = %.1e, coef v = %.1e\n", i, coef_u[i], coef_v[i]);
        tCudaVector3f f = u_base * coef_u[i] + v_base * coef_v[i];
        // printf("[calc_fint] fint comp %.1e %.1e %.1e\n", f[0], f[1], f[2]);
        fint.setcol(f, i);
    }
    return fint;
}

__device__ void CalcHessian(float tri_area, float Ku, float Kv,
                            const tCudaVector3f &coef_u,
                            const tCudaVector3f &coef_v,
                            const tCudaVector3f &Fu, const tCudaVector3f &Fv,
                            tCudaMatrix3f *array)
{
    CUDA_function;
    // 1. calculate Du and Dv
    float Fu_norm = Fu.norm();
    float Fv_norm = Fu.norm();
    SIM_ASSERT(Fu_norm > 1e-6);
    SIM_ASSERT(Fv_norm > 1e-6);

    tCudaVector3f Fu_normalized = Fu.normalized();
    tCudaVector3f Fv_normalized = Fv.normalized();

    tCudaMatrix3f Du =
        1.0 / Fu_norm * Fu_normalized * Fu_normalized.transpose();
    tCudaMatrix3f Dv =
        1.0 / Fv_norm * Fv_normalized * Fv_normalized.transpose();

    if (Fu_norm > 1)
        Du += (1 - 1.0 / Fu_norm) * tCudaMatrix3f::Identity();
    if (Fv_norm > 1)
        Dv += (1 - 1.0 / Fv_norm) * tCudaMatrix3f::Identity();

    // 2. calculate array
    // H = -S * (Ku * ciu * cju * Du + Kv * civ * cjv * Dv)
    // row major storage in a list

    for (int i = 0; i < 3; i++)
    {
        float ciu = coef_u[i], civ = coef_v[i];
        for (int j = 0; j < 3; j++)
        {
            float cju = coef_u[j], cjv = coef_v[j];
            array[3 * i + j] +=
                (-tri_area * (Ku * ciu * cju * Du + Kv * civ * cjv * Dv));
        }
    }
}

__device__ void
UpdateComponent(float Ku, float Kv, const tCudaVector3f &coef_Fu,
                const tCudaVector3f &coef_Fv, const tCudaVector3f &v0,
                const tCudaVector3f &v1, const tCudaVector3f &v2,
                float tri_area, tCudaMatrix3f *Fint_total,
                tCudaMatrix3f *Hessian)
{
    CUDA_function;

    tCudaVector3f Fu, Fv;

    CalcFuAndFv(v0, v1, v2, coef_Fu, coef_Fv, Fu, Fv);

    float Cu = Fu.norm() - 1;
    float Cv = Fv.norm() - 1;
    // printf("Cu %.1e Cv %.1e, Fu %.1e %.1e %.1e Fv %.1e %.1e %.1e\n", Cu, Cv,
        //    Fu[0], Fu[1], Fu[2], Fv[0], Fv[1], Fv[2]);
    Fint_total[0] +=
        CalcFint(tri_area, Ku, Kv, Cu, Cv, coef_Fu, coef_Fv, Fu, Fv);

    CalcHessian(tri_area, Ku, Kv, coef_Fu, coef_Fv, Fu, Fv, Hessian);
}

template <int N>
__device__ int
global_to_local_same(const tCudaMatrix<int, N, 1> &local_to_global_id,
                     int global_id)
{
    CUDA_function;
    int res = -1;
    for (int i = 0; i < N; i++)
    {
        if (local_to_global_id[i] == global_id)
        {
            res = i;
        }
    }
    return res;
}
__global__ void UpdateStretch_K_fint(
    int num_of_tri,
    devPtr<const tCudaVector3i> tri_vertices_id,  // vertices id in triangles
    devPtr<const tCudaVector3f> x_cur,            // current vertices pos
    devPtr<const float> tri_area_lst,             // triangle area
    float Ku, float Kv, float Kb,                 // stiffness
    devPtr<const tCudaVector3f> coef_Fu_warpweft, // Fu coef, warp weft
    devPtr<const tCudaVector3f> coef_Fv_warpweft, // Fv coef, warp weft
    devPtr<const tCudaVector3f> coef_Fu_diag,     // Fu coef, diag antidiag
    devPtr<const tCudaVector3f> coef_Fv_diag,     // Fu coef, diag antidiag,
    devPtr<tCudaVector3f> fint_total,
    devPtr<const tCudaVector32i> ELL_local_vertex_id_to_global_vertex_id,
    devPtr2<tCudaMatrix3f> H_total)
{
    CUDA_function;
    int tri_id = threadIdx.x + blockIdx.x * blockDim.x;
    if (tri_id >= num_of_tri)
        return;
    tCudaMatrix3f Fint_total[1];
    tCudaMatrix3f Hessian[9];
    for (int i = 0; i < 9; i++)
    {
        Hessian[i].setZero();
    }
    float tri_area = tri_area_lst[tri_id];
    tCudaVector3i v_id = tri_vertices_id[tri_id];
    const tCudaVector3f v0 = x_cur[v_id[0]];
    const tCudaVector3f v1 = x_cur[v_id[1]];
    const tCudaVector3f v2 = x_cur[v_id[2]];

    // printf(
    //     "v0(%d) %.1e %.1e %.1e v1(%d) %.1e %.1e %.1e v2(%d) %.1e %.1e %.1e\n",
    //     v_id[0], v0[0], v0[1], v0[2], v_id[1], v1[0], v1[1], v1[2], v_id[2],
    //     v2[0], v2[1], v2[2]);

    UpdateComponent(Ku, Kv, coef_Fu_warpweft[tri_id], coef_Fv_warpweft[tri_id],
                    v0, v1, v2, tri_area, Fint_total, Hessian);
    UpdateComponent(Kb, Kb, coef_Fu_diag[tri_id], coef_Fv_diag[tri_id], v0, v1,
                    v2, tri_area, Fint_total, Hessian);

    // begin to dispatch fint
    for (int i = 0; i < 3; i++)
    {
        int v_global_id = v_id[i];
        tCudaVector3f fint = Fint_total[0].col(i);
        // printf("[fint] add fint%d %.1e, %.1e, %.1e\n", v_global_id, fint[0],
        //        fint[1], fint[2]);
        cCudaIntrinsic::AtomicAdd(&fint_total[v_global_id], fint);
    }

    // begin to dispatch hessian

    for (int i = 0; i < 3; i++)
    {
        int vi_global_id = v_id[i];
        tCudaVector32i local_to_global_id_map =
            ELL_local_vertex_id_to_global_vertex_id[vi_global_id];
        for (int j = 0; j < 3; j++)
        {
            int _idx = 3 * i + j;
            int vj_global_id = v_id[j];
            int ELL_local_id =
                global_to_local_same(local_to_global_id_map, vj_global_id);
            cCudaIntrinsic::AtomicAdd(&H_total[vi_global_id][ELL_local_id],
                                      Hessian[_idx]);
        }
    }
}

void UpdateStiffnessMatrixAndFint(
    const cCudaArray<tCudaVector3i> &mTriangleVertexIdCuda,
    const cCudaArray<float> &mTriangleInitAreaCuda,
    const cCudaArray<tCudaVector3f> &mXcurCuda, const tCudaVector3f &mStretchK,
    const cCudaArray<tCudaVector3f> &mCoefFu_warp_weft,
    const cCudaArray<tCudaVector3f> &mCoefFv_warp_weft,
    const cCudaArray<tCudaVector3f> &mCoefFu_diag,
    const cCudaArray<tCudaVector3f> &mCoefFv_diag,
    const cCudaArray<tCudaVector32i> &mELLVidToGlobalVid,
    cCudaArray<tCudaVector3f> &mIntForceCuda,
    cCuda2DArray<tCudaMatrix3f> &mStiffnessMatrixCuda)
{
    // printf("mStiffnessMatrixCuda = %d %d\n", mStiffnessMatrixCuda.Rows(),
    //        mStiffnessMatrixCuda.Columns());
    mStiffnessMatrixCuda.MemsetAsync(tCudaMatrix3f::Zero());
    CUDA_ERR("clear K");
    // printf("fint size last = %d\n", mIntForceCuda.Size());
    mIntForceCuda.MemsetAsync(tCudaVector3f::Zero());
    CUDA_ERR("clear f");

    int num_of_tri = mTriangleVertexIdCuda.Size();
    // printf("mTriangleVertexIdCuda size = %d\n",
    // mTriangleVertexIdCuda.Size()); printf("mXcurCuda size = %d\n",
    // mXcurCuda.Size()); printf("mTriangleInitAreaCuda size = %d\n",
    // mTriangleInitAreaCuda.Size());
    // 1. update stretch
    UpdateStretch_K_fint CUDA_at(num_of_tri, 128)(
        num_of_tri, mTriangleVertexIdCuda.Ptr(), mXcurCuda.Ptr(),
        mTriangleInitAreaCuda.Ptr(), mStretchK[0], mStretchK[1], mStretchK[2],
        mCoefFu_warp_weft.Ptr(), mCoefFv_warp_weft.Ptr(), mCoefFu_diag.Ptr(),
        mCoefFv_diag.Ptr(), mIntForceCuda.Ptr(), mELLVidToGlobalVid.Ptr(),
        mStiffnessMatrixCuda.Ptr());
    CUDA_ERR("update K");
}

/**
 * \brief
 *  dt2K = dt * dt * K
 *  W = M - dt2K
 *  b = dt * dt * (mGravityForce + mUserForce + mIntForce)
 *      + dt * M * V_cur
 *  V_cur = (mXcur - mXpre) / dt
 */
__global__ void AssembleSystemMatrixKernel(
    int num_of_v, float dt, const float rayleigh_damping_alpha,
    const float rayleigh_damping_beta,
    devPtr<const tCudaVector32i> ELL_local_vertex_id_to_global_vertex_id,
    devPtr2<const tCudaMatrix3f> K, devPtr<const float> vertices_mass,
    devPtr2<tCudaMatrix3f> W)
{
    CUDA_function;
    int v_id = threadIdx.x + blockDim.x * blockIdx.x;
    if (v_id >= num_of_v)
        return;

    float alpha = rayleigh_damping_alpha, beta = rayleigh_damping_beta;
    // handle W.row(3 * v_id + [0, 3])
    tCudaVector32i ell_local_id_to_global_id =
        ELL_local_vertex_id_to_global_vertex_id[v_id];
    for (int i = 0; i < ell_local_id_to_global_id.size(); i++)
    {
        int cur_column_global_id = ell_local_id_to_global_id[i];
        if (cur_column_global_id == -1)
            break;
        else
        {
            // copy this column block (3x3)
            W[v_id][i] = (dt * beta - dt * dt) * K[v_id][i];
            // if this ELL columdn represents the same vertex as "v_is" tells
            // the mass matrix is row-diagnozation lumped, only its diagnoal is
            // nonzero.
            if (cur_column_global_id == v_id)
            {
                tCudaVector3f m_diag = vertices_mass[v_id] * (1 + dt * alpha) *
                                       tCudaVector3f::Ones();
                W[v_id][i](0, 0) += m_diag[0];
                W[v_id][i](1, 1) += m_diag[1];
                W[v_id][i](2, 2) += m_diag[2];
                // printf("for v%d begin to add diagnoal entry %.3f, %.3f,
                // %.3f\n",
                //    v_id, m_diag[0], m_diag[1], m_diag[2]);
            }
        }
    }
}
void AssembleSystemMatrix(
    float dt, const float rayleigh_damping_a, const float rayleigh_damping_b,
    const cCudaArray<tCudaVector32i> &ELL_local_vertex_id_to_global_vertex_id,
    const cCuda2DArray<tCudaMatrix3f> &K,
    const cCudaArray<float> &vertices_mass, cCuda2DArray<tCudaMatrix3f> &A)
{
    A.MemsetAsync(tCudaMatrix3f::Zero());
    int num_of_v = ELL_local_vertex_id_to_global_vertex_id.Size();
    AssembleSystemMatrixKernel CUDA_at(num_of_v, 128)(
        num_of_v, dt, rayleigh_damping_a, rayleigh_damping_b,
        ELL_local_vertex_id_to_global_vertex_id.Ptr(), K.Ptr(),
        vertices_mass.Ptr(), A.Ptr());
    CUDA_ERR("assembly system matrix W");
}

/**
 * b = dt * dt * (mGravityForce + mUserForce + mIntForce)
 *      + dt * (W + dt2K) * V_cur
 */
__global__ void AssembleSystemRHSKernel(
    int num_of_v, float dt, devPtr<const tCudaVector3f> mGravity,
    devPtr<const tCudaVector3f> mIntForce,
    devPtr<const tCudaVector32i> ELL_local_vertex_id_to_global_vertex_id,
    devPtr<const float> v_mass, devPtr<const tCudaVector3f> xcur,
    devPtr<const tCudaVector3f> xpre, devPtr<tCudaVector3f> RHS)
{
    CUDA_function;
    /*
        part1: dt * dt * (mGravityForce + mUserForce + mIntForce)

    */
    int v_id = threadIdx.x + blockDim.x * blockIdx.x;
    if (v_id >= num_of_v)
        return;
    RHS[v_id] = dt * dt * (mGravity[v_id] + mIntForce[v_id]);
    // if (true == cCudaMath::IsNan(RHS[v_id]))
    // {
    //     printf("[error] RHS[%d] step1 is nan!");
    // }
    // part2: dt * (W + dt2K) * V_cur
    const tCudaVector32i &ell_local_id_to_global_id =
        ELL_local_vertex_id_to_global_vertex_id[v_id];
    // do row * vel, and get the sum
    tCudaVector3f sum = tCudaVector3f::Zero();

    for (int i = 0; i < ell_local_id_to_global_id.size(); i++)
    {
        int cur_column_v_global_id = ell_local_id_to_global_id[i];
        if (cur_column_v_global_id == v_id)
        {
            sum = xcur[cur_column_v_global_id] - xpre[cur_column_v_global_id];
            sum[0] *= v_mass[v_id];
            sum[1] *= v_mass[v_id];
            sum[2] *= v_mass[v_id];
        }
    }

    RHS[v_id] += sum;
    // if (true == cCudaMath::IsNan(RHS[v_id]))
    // {
    //     printf("[error] RHS[%d] step2 is nan!");
    // }
}

/**
 * \brief           update linear system RHS
 * b = dt * dt * (mGravityForce + mUserForce + mIntForce) + dt * (W + dt2K) *
 * V_cur
 */
void AssembleSystemRHS(
    float dt, const cCudaArray<tCudaVector3f> &Gravity,
    const cCudaArray<tCudaVector3f> &IntForce,
    const cCudaArray<tCudaVector32i> &ELL_local_vertex_id_to_global_vertex_id,
    const cCudaArray<float> &v_mass, const cCudaArray<tCudaVector3f> &x_cur,
    const cCudaArray<tCudaVector3f> &x_pre, cCudaArray<tCudaVector3f> &RHS)
{

    int num_of_v = Gravity.Size();
    // std::cout << "gravity size = " << Gravity.Size() << std::endl;
    // std::cout << "UserForce size = " << UserForce.Size() << std::endl;
    // std::cout << "IntForce size = " << IntForce.Size() << std::endl;
    // std::cout << "W size = " << W.Rows() << " " << W.Columns() << std::endl;
    // std::cout << "K size = " << K.Rows() << " " << K.Columns() << std::endl;
    // std::cout << "cur vel size = " << cur_v.Size() << std::endl;
    // std::cout << "RHS size = " << RHS.Size() << std::endl;
    RHS.MemsetAsync(tCudaVector3f::Zero());
    AssembleSystemRHSKernel CUDA_at(num_of_v, 128)(
        num_of_v, dt, Gravity.Ptr(), IntForce.Ptr(),
        ELL_local_vertex_id_to_global_vertex_id.Ptr(), v_mass.Ptr(),
        x_cur.Ptr(), x_pre.Ptr(), RHS.Ptr());

    CUDA_ERR("assemble system RHS");
    // printf("---begin to check RHS---\n");
    // std::vector<tCudaVector3f> rhs_cpu;
    // RHS.Download(rhs_cpu);
    // for (int i = 0; i < rhs_cpu.size(); i++)
    // {
    //     if (true == cCudaMath::IsNan(rhs_cpu[i]))
    //     {
    //         std::cout << "RHS for vertex " << i << " has Nan\n";
    //         exit(1);
    //     }
    // }
}

void UpdateLinearSystem() {}
void SolveLinearSystem() {}
} // namespace BaraffClothGpu