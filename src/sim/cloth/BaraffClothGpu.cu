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
__global__ void
UpdateEleStiffKernel(int num_of_triangles, float k0, float k1,
                     devPtr<const tCudaVector3f> pos_lst,
                     devPtr<const tCudaVector3i> tri_vertices_id_lst,
                     devPtr<const tCudaMatrix32f> N_lst,
                     devPtr<tCudaMatrix32f> F_lst, devPtr<tCudaMatrix32f> n_lst,
                     devPtr<tCudaVector2f> C_lst, // condition
                     devPtr<tCudaMatrix92f> g_lst, devPtr<tCudaMatrix9f> K_lst)
{
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
        // printf("for vertex %d handle inner edge %d, included vertices: %d %d %d %d\n", v_id, cur_edge, vertex_global_id_involved[0], vertex_global_id_involved[1], vertex_global_id_involved[2], vertex_global_id_involved[3]);
        // return;
        // 2.2 get index in ELL representation
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
            // printf("take eleK part(%d %d) to global(%d, %d), means ELL (%d,%d)\n",
            //        cur_vertex_in_constraint_local_id, k, 
            //        v_id, vertex_global_id_involved[k],
            //        v_id,
            //        vertex_ell_local_id_involed[k]);
        }
    }
}
__global__ void DispatchStretchKToGlobalK(
    int num_of_v, devPtr<const tCudaMatrix9f> ele_K_lst_per_triangle,
    devPtr<const tCudaVector32i> vertex_connected_triangles_lst,
    devPtr<const tCudaVector3i> vertex_id_in_triangle_lst,
    devPtr<const tCudaVector32i> ELL_local_vertex_id_to_global_vertex_id,
    devPtr2<tCudaMatrix3f> global_K)
{

    int v_id = blockIdx.x * blockDim.x + threadIdx.x;
    if (v_id >= num_of_v)
        return;
    const tCudaVector32i &connected_tri = vertex_connected_triangles_lst[v_id];
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
            const tCudaMatrix9f &ele_K = ele_K_lst_per_triangle[tri_id];
            int cur_v_local_id = global_to_local(triangle_vid_lst, v_id);
            // column major
            for (int j = 0; j < 3; j++)
            {
                int v_global_id = triangle_vid_lst[j];
                // int v_local_id = global_to_local[v_global_id];
                int v_local_id = global_to_local(
                    ELL_local_vertex_id_to_global_vertex_id[v_id], v_global_id);
                // printf(
                //     "handle vertex %d, triangle %d, vertex %d (local id
                //     %d)\n", v_id, tri_id, v_global_id, v_local_id);
                const tCudaMatrix3f &part =
                    ele_K.block<3, 3>(3 * cur_v_local_id, 3 * j);

                // row major
                global_K[v_id][v_local_id] += part;
            }
        }
    }
}
namespace BaraffClothGpu
{
void UpdateK(float K0, float K1, const cCudaArray<tCudaVector3f> &pos_lst,
             const cCudaArray<tCudaVector3i> &tri_vertices_lst,
             cCudaArray<tCudaMatrix32f> &N_lst,
             cCudaArray<tCudaMatrix32f> &F_lst,
             cCudaArray<tCudaMatrix32f> &n_lst,
             cCudaArray<tCudaVector2f> &C_lst, // condition
             cCudaArray<tCudaMatrix92f> &g_lst,
             cCudaArray<tCudaMatrix9f> &K_lst)
{
    int num_of_tri = tri_vertices_lst.Size();
    UpdateEleStiffKernel CUDA_at(num_of_tri, 128)(
        num_of_tri, K0, K1, pos_lst.Ptr(), tri_vertices_lst.Ptr(), N_lst.Ptr(),
        F_lst.Ptr(), n_lst.Ptr(), C_lst.Ptr(), g_lst.Ptr(), K_lst.Ptr());
    CUDA_ERR("UpdateK");
}

// one row, one line, one vertex
void AssembleStretchStiffnessMatrix(
    const cCudaArray<tCudaMatrix9f> &ele_K_lst_per_triangle,
    const cCudaArray<tCudaVector32i> &vertex_connected_triangle_lst,
    const cCudaArray<tCudaVector3i> &vertex_id_of_triangle_lst,
    const cCudaArray<tCudaVector32i> &ELL_local_vertex_id_to_global_vertex_id,
    cCuda2DArray<tCudaMatrix3f> &global_K)
{
    int num_of_v = vertex_connected_triangle_lst.Size();

    std::cout << "[debug] ele K size = " << ele_K_lst_per_triangle.Size()
              << std::endl;
    DispatchStretchKToGlobalK CUDA_at(num_of_v, 128)(
        num_of_v, ele_K_lst_per_triangle.Ptr(),
        vertex_connected_triangle_lst.Ptr(), vertex_id_of_triangle_lst.Ptr(),
        ELL_local_vertex_id_to_global_vertex_id.Ptr(), global_K.Ptr());
}

template <typename T>
std::vector<T> DownloadFromGPUToCpu(const cCudaArray<T> &array)
{
    std::vector<T> res;
    array.Download(res);
    return res;
}
template <typename T>
std::vector<T> DownloadFromGPUToCpu2D(const cCudaArray<T> &array)
{
    std::vector<T> res;
    array.Download(res);
    return res;
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
} // namespace BaraffClothGpu
