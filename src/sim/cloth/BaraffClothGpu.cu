#include "sim/gpu_utils/CudaArray.h"
#include "sim/gpu_utils/CudaDevPtr.h"
#include "sim/gpu_utils/CudaMatrix.h"
#include <cassert>
#include <iostream>

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
__global__ void UpdateEleStiffKernel(
    int num_of_triangles, float k0, float k1, devPtr<tCudaVector3f> pos_lst,
    devPtr<tCudaVector3i> tri_vertices_id_lst, devPtr<tCudaMatrix32f> N_lst,
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
    printf("triangle %d v0 %d v1 %d v2 %d\n", tri_id, v_id[0], v_id[1],
           v_id[2]);
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
        printf("----begin to update K----\n");
        // float Xni0 = ;
        // float Xni1 =;
        tCudaVector2f X_Ni_norm =
            tCudaVector2f({F.col(0).norm(), F.col(1).norm()});

        tCudaMatrix3f I3;

        I3.setIdentity();
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
namespace BaraffClothGpu
{
void UpdateK(float K0, float K1, cCudaArray<tCudaVector3f> &pos_lst,
             cCudaArray<tCudaVector3i> &tri_vertices_lst,
             cCudaArray<tCudaMatrix32f> &N_lst,
             cCudaArray<tCudaMatrix32f> &F_lst,
             cCudaArray<tCudaMatrix32f> &n_lst,
             cCudaArray<tCudaVector2f> &C_lst, // condition
             cCudaArray<tCudaMatrix92f> &g_lst,
             cCudaArray<tCudaMatrix9f> &K_lst)
{
    // std::cout << pos_lst.Size() << std::endl;
    // std::cout << tri_vertices_lst.Size() << std::endl;
    // std::cout << N_lst.Size() << std::endl;
    // std::cout << F_lst.Size() << std::endl;
    // std::cout << n_lst.Size() << std::endl;
    // std::cout << C_lst.Size() << std::endl;
    // std::cout << g_lst.Size() << std::endl;
    int num_of_tri = tri_vertices_lst.Size();
    UpdateEleStiffKernel CUDA_at(num_of_tri, 128)(
        num_of_tri, K0, K1, pos_lst.Ptr(), tri_vertices_lst.Ptr(), N_lst.Ptr(),
        F_lst.Ptr(), n_lst.Ptr(), C_lst.Ptr(), g_lst.Ptr(), K_lst.Ptr());
    CUDA_ERR("UpdateK");
}

} // namespace BaraffClothGpu