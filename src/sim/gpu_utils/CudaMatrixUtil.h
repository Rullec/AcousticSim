#include "CudaMatrix.h"
#include "utils/EigenUtil.h"

namespace cCudaMatrixUtil
{
template <typename dtype, int N, int M>
Eigen::Matrix<dtype, N, M>
CudaMatrixToEigenMatrix(const tCudaMatrix<dtype, N, M> &cuda_mat)
{
    Eigen::Matrix<dtype, N, M> eigen;
    for (int i = 0; i < N; i++)
        for (int j = 0; j < M; j++)
        {
            eigen(i, j) = cuda_mat(i, j);
        }
    return eigen;
};

template <typename dtype, int N, int M>
tCudaMatrix<dtype, N, M>
EigenMatrixToCudaMatrix(const Eigen::Matrix<dtype, N, M> &eigen_mat)
{
    tCudaMatrix<dtype, N, M> mat;
    for (int i = 0; i < N; i++)
        for (int j = 0; j < M; j++)
        {
            mat(i, j) = eigen_mat(i, j);
        }
    return mat;
};

}; // namespace cCudaMatrixUtil
