#include "sim/gpu_utils/Cuda2DArray.h"
#include "sim/gpu_utils/CudaArray.h"
#include "sim/gpu_utils/CudaDevPtr.h"
#include "sim/gpu_utils/CudaMatrix.h"
#include <cassert>
#include <iostream>
#include <map>

// __device__ void AddeleStiffComp(float K, const tCudaVector9f &g, float C,
//                                 float XNinormj, const tCudaVector3f &N,
//                                 const tCudaMatrix3f &P,
//                                 tCudaMatrix9f &eleStiff){
//     // eleStiff.Kronecker
//     /*
//         eleStiff +=
//             K * ggT - C/XNinorm * Kronecker(NNT, P)
//     */
// };
__global__ void UpdateEleStretchStiffAndFint_Kernel(
    int num_of_triangles, float k0, float k1,
    devPtr<const tCudaVector3f> pos_lst,
    devPtr<const tCudaVector3i> tri_vertices_id_lst,
    devPtr<const tCudaMatrix32f> N_lst, devPtr<tCudaMatrix32f> F_lst,
    devPtr<tCudaMatrix32f> n_lst,
    devPtr<tCudaVector2f> C_lst, // condition
    devPtr<tCudaMatrix92f> g_lst, devPtr<tCudaMatrix9f> K_lst,
    devPtr<tCudaVector9f> fint_lst_each_triangle)
{
    CUDA_function;
    // F = X * N
    int tri_id = blockIdx.x * blockDim.x + threadIdx.x;
    if (tri_id >= num_of_triangles)
        return;

    // for current triangle
    // 1. get X
    tCudaMatrix3f X;
    tCudaVector3i v_id = tri_vertices_id_lst[tri_id];
    // printf("triangle %d v0 %d v1 %d v2 %d\n", tri_id, v_id[0], v_id[1],
    //    v_id[2]);
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            X(j, i) = pos_lst[v_id[i]][j];
        }
    }

    // 2. calculate F
    tCudaMatrix32f F = X * N_lst[tri_id];

    // 3. calculate n
    tCudaMatrix32f n;
    {
        n.setcol(F.col(0).normalized(), 0);
        n.setcol(F.col(1).normalized(), 1);
        // printf("n = ")
    }
    // 4. calculate C
    tCudaVector2f C;
    {
        C[0] = F.col(0).norm() - 1;
        C[1] = F.col(1).norm() - 1;
    }
    // 5. calculate g
    tCudaMatrix92f g;
    {
        auto N = N_lst[tri_id];
        for (size_t i = 0; i < 2; i++)
        {
            auto Ni = N.col(i);
            auto ni = n.col(i);
            for (int j = 0; j < 3; j++)
            {
                float Nij = Ni[j];
                for (int k = 0; k < 3; k++)
                {
                    g(3 * j + k, i) = Nij * ni[k];
                }
            }
        }
    }
    F_lst[tri_id] = F;
    n_lst[tri_id] = n;
    C_lst[tri_id] = C;
    g_lst[tri_id] = g;

    // update K

    {
        // printf("----begin to update K----\n");
        // float Xni0 = ;
        // float Xni1 =;
        tCudaVector2f X_Ni_norm =
            tCudaVector2f({F.col(0).norm(), F.col(1).norm()});

        tCudaMatrix3f I3 = tCudaMatrix3f::Identity();
        // return;
        tCudaMatrix3f P[2];

        P[0] = I3 - n.col(0) * n.col(0).transpose();
        P[1] = I3 - n.col(1) * n.col(1).transpose();
        // return;
        for (int j = 0; j < 2; j++)
        {
            float cur_k = j == 0 ? k0 : k1;
            tCudaMatrix9f part1 = -g.col(j) * g.col(j).transpose();
            assert(false == cCudaMath::IsNan(part1));
            tCudaMatrix9f part2 =
                (N_lst[tri_id].col(j) * N_lst[tri_id].col(j).transpose())
                    .KroneckerProduct(P[j]);
            assert(false == cCudaMath::IsNan(part2));
            // return;
            K_lst[tri_id] += cur_k * (part1 - C[j] / X_Ni_norm[j] * part2

                                     );
            assert(false == cCudaMath::IsNan(K_lst[tri_id]));
        }
    }

    // update fint
    for (int i = 0; i < 2; i++)
    {
        float cur_k = i == 0 ? k0 : k1;
        fint_lst_each_triangle[tri_id] += -cur_k * C[i] * g.col(i);
    }
}

// __device__ int GetLocalId(const tCudaVector3i &v_id_lst, int v_id)
// {
//     for (int i = 0; i < 3; i++)
//     {
//         if (v_id_lst[i] == v_id)
//         {
//             return i;
//         }
//     }
//     printf("cannot find vertex %d in vertex id lst %d %d %d\n", v_id,
//            v_id_lst[0], v_id_lst[1], v_id_lst[2]);
//     assert(false);
// }

template <int N>
__device__ int global_to_local(const tCudaMatrix<int, N, 1> &local_to_global_id,
                               int global_id)
{
    for (int i = 0; i < N; i++)
    {
        if (local_to_global_id[i] == global_id)
        {
            return i;
        }
    }
    assert(false);
    return -1;
}
// template <int N>
// int global_to_local_host(const tCudaMatrix<int, N, 1> &local_to_global_id,
//                          int global_id)
// {
//     for (int i = 0; i < N; i++)
//     {
//         if (local_to_global_id[i] == global_id)
//         {
//             return i;
//         }
//     }
//     assert(false);
//     return -1;
// }
// void DispatchBendingKToGlobalK( int v_id,
__global__ void DispatchBendingKToGlobalK(
    int num_of_v, devPtr<const tCudaVector4i> bending_vertex_per_edge_lst,
    devPtr<const tCudaVector32i> vertex_connected_edges_lst,
    devPtr<const tCudaMatrix12f> ele_K_lst,
    devPtr<const tCudaVector32i> ELL_local_vertex_id_to_global_vertex_id,
    devPtr2<tCudaMatrix3f> global_K_ele)
{
    CUDA_function;
    // 1.
    int v_id = threadIdx.x + blockDim.x * blockIdx.x;
    if (v_id >= num_of_v)
    {
        return;
    }
    /*
    for each row (vertex)
    1. get all edges: e0, e1, .. ei
    2. for current edge e0:
        2.1 get four involved index from bending_vertex_lst: v0_global,
    v1_global, v2_global, v3_global

        2.2 get index in ELL representation, v0_local, .., v3_local

        2.3 get this constriant's (edge's) ele K from ele_K_lst: K_{ele}

        2.4 get current vertex's local id, cur_local

        2.5 append K_{ele}(cur_local, idx) -> K_{global}.row(v_id, idx_local)
    */
    const tCudaVector32i &connected_edge = vertex_connected_edges_lst[v_id];
    const tCudaVector32i &ell_local_to_global_id =
        ELL_local_vertex_id_to_global_vertex_id[v_id];
    for (int _idx = 0; _idx < connected_edge.size(); _idx++)
    {
        int cur_edge = connected_edge[_idx];
        if (cur_edge == -1)
            break;

        // four vertices global indices, involved in this constraint
        const tCudaVector4i &vertex_global_id_involved =
            bending_vertex_per_edge_lst[cur_edge];
        // printf("for vertex %d handle inner edge %d, included vertices: %d %d
        // %d %d\n", v_id, cur_edge, vertex_global_id_involved[0],
        // vertex_global_id_involved[1], vertex_global_id_involved[2],
        // vertex_global_id_involved[3]); return; 2.2 get index in ELL
        // representation
        tCudaVector4i vertex_ell_local_id_involed;
        for (int k = 0; k < 4; k++)
        {
            vertex_ell_local_id_involed[k] = global_to_local(
                ell_local_to_global_id, vertex_global_id_involved[k]);
        }
        // 2.3 get ele_K
        const tCudaMatrix12f &ele_K = ele_K_lst[cur_edge];
        // 2.4 get current local

        int cur_vertex_in_constraint_local_id =
            global_to_local(vertex_global_id_involved, v_id);

        // 2.5 append K_ele
        for (int k = 0; k < 4; k++)
        {
            // row major
            global_K_ele[v_id][vertex_ell_local_id_involed[k]] +=
                ele_K.block<3, 3>(3 * cur_vertex_in_constraint_local_id, 3 * k);
            // printf("take eleK part(%d %d) to global(%d, %d), means ELL
            // (%d,%d)\n",
            //        cur_vertex_in_constraint_local_id, k,
            //        v_id, vertex_global_id_involved[k],
            //        v_id,
            //        vertex_ell_local_id_involed[k]);
        }
    }
}

__global__ void DispatchStretch_StiffnessMatrix_and_fint(
    int num_of_v, devPtr<const tCudaMatrix9f> ele_K_lst_per_triangle,
    devPtr<const tCudaVector32i> vertex_connected_triangles_lst,
    devPtr<const tCudaVector3i> vertex_id_in_triangle_lst,
    devPtr<const tCudaVector32i> ELL_local_vertex_id_to_global_vertex_id,
    devPtr<const tCudaVector9f> ele_fint, const int num_of_fixed_vertices,
    devPtr<const int> fixed_vertex_indices_lst,
    devPtr<const tCudaVector3f> fixed_vertex_target_pos_lst,
    devPtr<const tCudaVector3f> vertex_pos_lst, devPtr2<tCudaMatrix3f> global_K,
    devPtr<tCudaVector3f> global_fint)
{
    CUDA_function;
    int v_id = blockIdx.x * blockDim.x + threadIdx.x;
    if (v_id >= num_of_v)
        return;
    const tCudaVector32i &connected_tri = vertex_connected_triangles_lst[v_id];
    // 1. dispatch common stretch stiffness matrix
    for (int i = 0; i < 32; i++)
    {
        int tri_id = connected_tri[i];
        if (tri_id == -1)
        {
            break;
        }
        else
        // if (v_id == 0 && i == 0)
        {
            tCudaVector3i triangle_vid_lst = vertex_id_in_triangle_lst[tri_id];
            // printf("[debug] triangle %d include v %d %d %d\n", tri_id,
            //        triangle_vid_lst[0], triangle_vid_lst[1],
            //        triangle_vid_lst[2]);
            int cur_v_local_id = global_to_local(triangle_vid_lst, v_id);

            // assemble stiffness matrix
            {
                const tCudaMatrix9f &ele_K = ele_K_lst_per_triangle[tri_id];
                // column major
                for (int j = 0; j < 3; j++)
                {
                    int v_global_id = triangle_vid_lst[j];
                    // int v_local_id = global_to_local[v_global_id];
                    int v_local_id = global_to_local(
                        ELL_local_vertex_id_to_global_vertex_id[v_id],
                        v_global_id);
                    // printf(
                    //     "handle vertex %d, triangle %d, vertex %d (local id
                    //     %d)\n", v_id, tri_id, v_global_id, v_local_id);
                    const tCudaMatrix3f &part =
                        ele_K.block<3, 3>(3 * cur_v_local_id, 3 * j);

                    // row major
                    global_K[v_id][v_local_id] += part;
                }
            }
            // assemble internal force
            {
                // printf("fint %d from triangle %d local %d\n", v_id, tri_id,
                //        cur_v_local_id);
                global_fint[v_id] += tCudaVector3f(
                    ele_fint[tri_id].segment<3>(3 * cur_v_local_id));
            }
        }
    }

    // 2. add fixed point stiffness (implicit spring)
    {
        float fixed_point_K = 1e1;
        for (int i = 0; i < num_of_fixed_vertices; i++)
        {
            if (v_id == fixed_vertex_indices_lst[i])
            {
                // printf("need to fix vertex %d\n", v_id);
                // add fixed!
                tCudaVector3f target_pos = fixed_vertex_target_pos_lst[i];
                // printf("vertex %d tar pos %.3f %.3f %.3f\n", v_id,
                //    target_pos[0], target_pos[1], target_pos[2]);
                tCudaVector3f cur_pos = vertex_pos_lst[v_id];
                // printf("vertex %d cur pos %.3f %.3f %.3f\n", v_id,
                // cur_pos[0],
                //    cur_pos[1], cur_pos[2]);
                // fint = K (p - x)
                // global_fint[v_id] += fixed_point_K * (target_pos - cur_pos);
                // H = -K, add to ELL
                int v_local_id = global_to_local(
                    ELL_local_vertex_id_to_global_vertex_id[v_id], v_id);
                global_K[v_id][v_local_id] += 1e12;
                // -fixed_point_K * tCudaMatrix3f::Identity();
            }
        }
    }
}

namespace BaraffClothGpu
{
void UpdateStretch_StiffnessMatrix_Fint(
    float K0, float K1, const cCudaArray<tCudaVector3f> &pos_lst,
    const cCudaArray<tCudaVector3i> &tri_vertices_lst,
    cCudaArray<tCudaMatrix32f> &N_lst, cCudaArray<tCudaMatrix32f> &F_lst,
    cCudaArray<tCudaMatrix32f> &n_lst,
    cCudaArray<tCudaVector2f> &C_lst, // condition
    cCudaArray<tCudaMatrix92f> &g_lst, cCudaArray<tCudaMatrix9f> &K_lst,
    cCudaArray<tCudaVector9f> &fint_lst_in_each_triangle)
{
    int num_of_tri = tri_vertices_lst.Size();
    UpdateEleStretchStiffAndFint_Kernel CUDA_at(num_of_tri, 128)(
        num_of_tri, K0, K1, pos_lst.Ptr(), tri_vertices_lst.Ptr(), N_lst.Ptr(),
        F_lst.Ptr(), n_lst.Ptr(), C_lst.Ptr(), g_lst.Ptr(), K_lst.Ptr(),
        fint_lst_in_each_triangle.Ptr());
    CUDA_ERR("UpdateStretch_StiffnessMatrix_Fint");
}

// one row, one line, one vertex
void AssembleStretch_StiffnessMatrix_Fint(
    const cCudaArray<tCudaMatrix9f> &ele_K_lst_per_triangle,
    const cCudaArray<tCudaVector32i> &vertex_connected_triangle_lst,
    const cCudaArray<tCudaVector3i> &vertex_id_of_triangle_lst,
    const cCudaArray<tCudaVector32i> &ELL_local_vertex_id_to_global_vertex_id,
    const cCudaArray<tCudaVector9f> &ele_fint,
    const cCudaArray<int> &fixed_vertex_indices_lst,
    const cCudaArray<tCudaVector3f> &fixed_vertex_target_pos_lst,
    const cCudaArray<tCudaVector3f> &vertex_pos_lst,

    cCuda2DArray<tCudaMatrix3f> &global_K,
    cCudaArray<tCudaVector3f> &global_fint)
{
    int num_of_v = vertex_connected_triangle_lst.Size();

    // std::cout << "[debug] ele K size = " << ele_K_lst_per_triangle.Size()
    //   << std::endl;
    DispatchStretch_StiffnessMatrix_and_fint CUDA_at(num_of_v, 128)(
        num_of_v, ele_K_lst_per_triangle.Ptr(),
        vertex_connected_triangle_lst.Ptr(), vertex_id_of_triangle_lst.Ptr(),
        ELL_local_vertex_id_to_global_vertex_id.Ptr(), ele_fint.Ptr(),
        fixed_vertex_indices_lst.Size(), fixed_vertex_indices_lst.Ptr(),
        fixed_vertex_target_pos_lst.Ptr(), vertex_pos_lst.Ptr(), global_K.Ptr(),
        global_fint.Ptr());

    CUDA_ERR("dispatch stretch stiffness matrix and fint");
}

void AssembleBendingStiffnessMatrix(
    const cCudaArray<tCudaVector4i> &bending_vertex_per_edge_lst,
    const cCudaArray<tCudaVector32i> &vertex_connected_edges_lst,
    const cCudaArray<tCudaMatrix12f> &ele_K_lst,
    const cCudaArray<tCudaVector32i> &ELL_local_vertex_id_to_global_vertex_id,
    cCuda2DArray<tCudaMatrix3f> &global_K)
{
    // std::cout << bending_vertex_per_edge_lst.Size() << std::endl;
    // std::cout << vertex_connected_edges_lst.Size() << std::endl;
    // std::cout << ele_K_lst.Size() << std::endl;
    // std::cout << ELL_local_vertex_id_to_global_vertex_id << std::endl;
    int num_of_v = vertex_connected_edges_lst.Size();
    DispatchBendingKToGlobalK CUDA_at(num_of_v, 128)(
        num_of_v, bending_vertex_per_edge_lst.Ptr(),
        vertex_connected_edges_lst.Ptr(), ele_K_lst.Ptr(),
        ELL_local_vertex_id_to_global_vertex_id.Ptr(), global_K.Ptr());
    // for (int i = 0; i < num_of_v; i++)
    // {
    //     DispatchBendingKToGlobalK(
    //         i, num_of_v, bending_vertex_per_edge_lst.Ptr(),
    //         vertex_connected_edges_lst.Ptr(), ele_K_lst.Ptr(),
    //         ELL_local_vertex_id_to_global_vertex_id.Ptr(), global_K.Ptr());
    // }
    CUDA_ERR("assemble bending");
    // exit(1);
}

/**
 * \brief           calc bending int force kernel
 *
 *
 *           f = H * x
 */
__global__ void CalcBendingIntForceKernel(
    int num_of_v,
    devPtr<const tCudaVector32i> ELL_local_vertex_id_to_global_vertex_id,
    devPtr2<const tCudaMatrix3f> bending_hessian,
    devPtr<const tCudaVector3f> vertex_pos_lst,
    devPtr<tCudaVector3f> bending_fint)
{
    CUDA_function;
    int v_id = threadIdx.x + blockDim.x * blockIdx.x;
    if (v_id >= num_of_v)
        return;

    tCudaVector32i ELL_vertex_local_id_to_global_id =
        ELL_local_vertex_id_to_global_vertex_id[v_id];
    for (int j = 0; j < ELL_vertex_local_id_to_global_id.size(); j++)
    {
        int cur_column_v = ELL_vertex_local_id_to_global_id[j];
        if (cur_column_v == -1)
            break;
        else
        {
            bending_fint[v_id] +=
                bending_hessian[v_id][j] * vertex_pos_lst[cur_column_v];
        }
    }
}
void UpdateBendingIntForce(
    const cCudaArray<tCudaVector32i> &ELL_local_vertex_id_to_global_vertex_id,
    const cCuda2DArray<tCudaMatrix3f> &bending_hessian,
    const cCudaArray<tCudaVector3f> &vertex_pos_lst,
    cCudaArray<tCudaVector3f> &bending_fint)
{
    int num_of_v = vertex_pos_lst.Size();
    CalcBendingIntForceKernel CUDA_at(num_of_v, 128)(
        num_of_v, ELL_local_vertex_id_to_global_vertex_id.Ptr(),
        bending_hessian.Ptr(), vertex_pos_lst.Ptr(), bending_fint.Ptr());
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
    int num_of_v, float dt,
    devPtr<const tCudaVector32i> ELL_local_vertex_id_to_global_vertex_id,
    devPtr2<const tCudaMatrix3f> K, devPtr<const tCudaVector3f> M,
    devPtr2<tCudaMatrix3f> W)
{
    CUDA_function;
    int v_id = threadIdx.x + blockDim.x * blockIdx.x;
    if (v_id >= num_of_v)
        return;

    float alpha = 0, beta = 0;
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
                tCudaVector3f m_diag = M[v_id] * (1 + dt * alpha);
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
    float dt,
    const cCudaArray<tCudaVector32i> &ELL_local_vertex_id_to_global_vertex_id,
    const cCuda2DArray<tCudaMatrix3f> &K, const cCudaArray<tCudaVector3f> &M,
    cCuda2DArray<tCudaMatrix3f> &W)
{
    int num_of_v = ELL_local_vertex_id_to_global_vertex_id.Size();
    AssembleSystemMatrixKernel CUDA_at(num_of_v, 128)(
        num_of_v, dt, ELL_local_vertex_id_to_global_vertex_id.Ptr(), K.Ptr(),
        M.Ptr(), W.Ptr());
    CUDA_ERR("assembly system matrix W");
}

/**
 * b = dt * dt * (mGravityForce + mUserForce + mIntForce)
 *      + dt * (W + dt2K) * V_cur
 */
__global__ void AssembleSystemRHSKernel(
    int num_of_v, float dt, devPtr<const tCudaVector3f> mGravity,
    devPtr<const tCudaVector3f> mUserForce,
    devPtr<const tCudaVector3f> mIntForce,
    devPtr<const tCudaVector32i> ELL_local_vertex_id_to_global_vertex_id,
    devPtr<const tCudaVector3f> M, devPtr<const tCudaVector3f> Vel,
    devPtr<tCudaVector3f> RHS)
{
    CUDA_function;
    /*
        part1: dt * dt * (mGravityForce + mUserForce + mIntForce)

    */
    int v_id = threadIdx.x + blockDim.x * blockIdx.x;
    if (v_id >= num_of_v)
        return;
    RHS[v_id] = dt * dt * (mGravity[v_id] + mUserForce[v_id] + mIntForce[v_id]);
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
            sum = Vel[cur_column_v_global_id] * dt;
            sum[0] *= M[v_id][0];
            sum[1] *= M[v_id][1];
            sum[2] *= M[v_id][2];
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
    const cCudaArray<tCudaVector3f> &UserForce,
    const cCudaArray<tCudaVector3f> &IntForce,
    const cCudaArray<tCudaVector32i> &ELL_local_vertex_id_to_global_vertex_id,
    const cCudaArray<tCudaVector3f> &M, const cCudaArray<tCudaVector3f> &cur_v,
    cCudaArray<tCudaVector3f> &RHS)
{

    int num_of_v = Gravity.Size();
    // std::cout << "gravity size = " << Gravity.Size() << std::endl;
    // std::cout << "UserForce size = " << UserForce.Size() << std::endl;
    // std::cout << "IntForce size = " << IntForce.Size() << std::endl;
    // std::cout << "W size = " << W.Rows() << " " << W.Columns() << std::endl;
    // std::cout << "K size = " << K.Rows() << " " << K.Columns() << std::endl;
    // std::cout << "cur vel size = " << cur_v.Size() << std::endl;
    // std::cout << "RHS size = " << RHS.Size() << std::endl;
    AssembleSystemRHSKernel CUDA_at(num_of_v, 128)(
        num_of_v, dt, Gravity.Ptr(), UserForce.Ptr(), IntForce.Ptr(),
        ELL_local_vertex_id_to_global_vertex_id.Ptr(), M.Ptr(), cur_v.Ptr(),
        RHS.Ptr());

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

/**
 * \brief           claculate the current velocity
 */
__global__ void CalcCurVelKernel(int num_of_v, float dt,
                                 devPtr<const tCudaVector3f> mXcur,
                                 devPtr<const tCudaVector3f> mXpre,
                                 devPtr<tCudaVector3f> mVel)
{
    CUDA_function;
    int v_id = threadIdx.x + blockIdx.x * blockDim.x;
    if (v_id >= num_of_v)
        return;

    mVel[v_id] = (mXcur[v_id] - mXpre[v_id]) / dt;
}
void CalcCurVelocity(float dt, const cCudaArray<tCudaVector3f> &mXcur,
                     const cCudaArray<tCudaVector3f> &mXpre,
                     cCudaArray<tCudaVector3f> &mVel)

{
    int num_of_v = mXpre.Size();

    CalcCurVelKernel CUDA_at(num_of_v, 128)(num_of_v, dt, mXcur.Ptr(),
                                            mXpre.Ptr(), mVel.Ptr());
    CUDA_ERR("calcualte current velocity");
}
} // namespace BaraffClothGpu
