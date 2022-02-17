#include "FourOrderTensor.h"
#include "utils/LogUtil.h"

cFourOrderTensor::cFourOrderTensor(size_t dimI, size_t dimJ, size_t dimK, size_t dimL) : mDimI(dimI), mDimJ(dimJ), mDimK(dimK), mDimL(dimL)
{
    Allocate();
}

cFourOrderTensor::cFourOrderTensor(const cFourOrderTensor &old_tensor)
{
    tVector4i shape = old_tensor.GetShape();
    mDimI = shape[0];
    mDimJ = shape[1];
    mDimK = shape[2];
    mDimL = shape[3];

    mInnerArray = old_tensor.GetInnerArray();
}

tEigenArr<tMatrixXd> cFourOrderTensor::GetInnerArray() const
{
    return mInnerArray;
}
void cFourOrderTensor::Allocate()
{
    mInnerArray.resize(mDimI * mDimJ);
    for (auto &x : mInnerArray)
    {
        x.noalias() = tMatrixXd::Zero(mDimK, mDimL);
    }
}
tMatrixXd cFourOrderTensor::operator()(size_t i, size_t j) const
{
    if (i >= mDimI || j >= mDimJ || i < 0 || j < 0)
    {
        SIM_ERROR("illegal access ({}, {}) for 4th order tensor shaped ({}, {}, {}, {})", i, j, mDimI, mDimJ, mDimK, mDimL);
        exit(1);
    }
    return mInnerArray[j * mDimI + i];
}
tMatrixXd &cFourOrderTensor::operator()(size_t i, size_t j)
{
    if (i >= mDimI || j >= mDimJ || i < 0 || j < 0)
    {
        SIM_ERROR("illegal access ({}, {}) for 4th order tensor shaped ({}, {}, {}, {})", i, j, mDimI, mDimJ, mDimK, mDimL);
        exit(1);
    }
    return mInnerArray[j * mDimI + i];
}

cFourOrderTensor cFourOrderTensor::operator+(const cFourOrderTensor &tensor)
{
    tVector4i shape = GetShape();
    if (shape != tensor.GetShape())
    {
        SIM_ERROR("illegal tensor add: shape {} and shape {}", shape.transpose(), tensor.GetShape().transpose());
        exit(1);
    }
    cFourOrderTensor new_tensor(*this);
    for (int i = 0; i < shape[0]; i++)
        for (int j = 0; j < shape[1]; j++)
        {
            new_tensor(i, j) += tensor(i, j);
        }
    return new_tensor;
}

cFourOrderTensor cFourOrderTensor::operator*(double scalar)
{
    cFourOrderTensor new_tensor(*this);
    for (int i = 0; i < mDimI; i++)
        for (int j = 0; j < mDimJ; j++)
        {
            new_tensor(i, j) *= scalar;
        }
    return new_tensor;
}
/**
 * \brief           expand a tensor T_{ijkl} to M_{pq}
 *      M_{i * prod(jkl) + j * prod(kl) + k * prod(l) + l} = T_{ijkl}
 *      For more details, please check "四阶张量转换为矩阵.md"
*/
tMatrixXd cFourOrderTensor::ExpandToMatrix()
{
    int rows = mDimI * mDimJ;
    int cols = mDimK * mDimL;
    // printf("expand to matrix (%d, %d)\n", rows, cols);
    tMatrixXd val = tMatrixXd::Zero(rows, cols);

    for (size_t i_idx = 0; i_idx < mDimI; i_idx++)
    {
        for (size_t j_idx = 0; j_idx < mDimJ; j_idx++)
        {
            const tMatrixXd &T_kl = (*this)(i_idx, j_idx);
            for (size_t k_idx = 0; k_idx < mDimK; k_idx++)
            {
                val.row(i_idx * mDimJ + j_idx).segment(k_idx * mDimL, mDimL).noalias() = T_kl.row(k_idx);
            }
        }
    }
    return val;
}

/**
 * \brief           convert a matrix to a tensor
*/
void cFourOrderTensor::LoadFromAMatrix(size_t dimI, size_t dimJ, size_t dimK, size_t dimL, const tMatrixXd &val)
{
    // shape assert
    if (dimI * dimJ != val.rows())
    {
        SIM_ERROR("dim i {} dim j {} is not compatible with val's row {}", dimI, dimJ, val.rows());
        exit(1);
    }
    if (dimK * dimL != val.cols())
    {
        SIM_ERROR("dim k {} dim l {} is not compatible with val's col {}", dimK, dimL, val.cols());
        exit(1);
    }

    // allocation
    Allocate();

    // fill in data
    for (size_t i = 0; i < mDimI; i++)
        for (size_t j = 0; j < mDimJ; j++)
        {
            tMatrixXd &compoenent = (*this)(i, j);
            const tVectorXd &val_row = val.row(i * mDimJ + j);
            for (size_t k = 0; k < dimK; k++)
            {
                compoenent.row(k) = val_row.segment(k * dimL, dimL);
            }
        }
}

tVector4i cFourOrderTensor::GetShape() const
{
    return tVector4i(mDimI, mDimJ, mDimK, mDimL);
}

tMatrixXd cFourOrderTensor::GetLastTwoComp(size_t k, size_t l)
{
    tMatrixXd ret = tMatrixXd::Zero(mDimI, mDimJ);
    for (size_t i = 0; i < mDimI; i++)
        for (size_t j = 0; j < mDimJ; j++)
        {
            ret(i, j) = (*this)(i, j)(k, l);
        }
    return ret;
}

void cFourOrderTensor::SetLastTwoComp(size_t k, size_t l, const tMatrixXd &val)
{
    for (size_t i = 0; i < mDimI; i++)
        for (size_t j = 0; j < mDimJ; j++)
        {
            (*this)(i, j)(k, l) = val(i, j);
        }
}

/**
 * \brief               Tensor * Matrix
*/
void cFourOrderTensor::TensorMatrixProductWithoutCopy(size_t target_dim0, size_t target_dim1, const tMatrixXd &matrix)
{
    size_t mat_rows = matrix.rows(),
           mat_cols = matrix.cols();
    if (target_dim0 == 0 && target_dim1 == 1)
    {
        // we need to get and set the last two component

        if (mat_rows != mDimJ)
        {
            SIM_ERROR("right matrix rows {} must equal to the second dim of tensor {}", mat_rows, mDimJ);
            exit(1);
        }

        for (size_t idx_k = 0; idx_k < mDimK; idx_k++)
            for (size_t idx_l = 0; idx_l < mDimL; idx_l++)
            {
                SetLastTwoComp(idx_k, idx_l, GetLastTwoComp(idx_k, idx_l) * matrix);
            }
    }
    else if (target_dim0 == 2 && target_dim1 == 3)
    {
        // multiplication on the last two dimensions
        for (auto &x : mInnerArray)
        {
            x = x * matrix;
        }
    }
    else
    {
        SIM_ERROR("unsupported target dim {} and {}", target_dim0, target_dim1);
        exit(1);
    }
    mDimJ = mat_cols;
}
void cFourOrderTensor::MatrixTensorProductWithoutCopy(size_t target_dim0, size_t target_dim1, const tMatrixXd &matrix)
{
    size_t mat_rows = matrix.rows(),
           mat_cols = matrix.cols();
    if (target_dim0 == 0 && target_dim1 == 1)
    {
        // we need to get and set the last two component
        if (mat_cols != mDimK)
        {
            SIM_ERROR("left matrix cols {} must equal to the third dim of tensor {}", mat_cols, mDimK);
            exit(1);
        }

        for (size_t idx_k = 0; idx_k < mDimK; idx_k++)
            for (size_t idx_l = 0; idx_l < mDimL; idx_l++)
            {
                SetLastTwoComp(idx_k, idx_l, matrix * GetLastTwoComp(idx_k, idx_l));
            }
    }
    else if (target_dim0 == 2 && target_dim1 == 3)
    {
        // multiplication on the last two dimensions
        for (auto &x : mInnerArray)
        {
            x = matrix * x;
        }
    }
    else
    {
        SIM_ERROR("unsupported target dim {} and {}", target_dim0, target_dim1);
        exit(1);
    }
    mDimK = mat_rows;
}

cFourOrderTensor TensorMatrixProductWithCopy(const cFourOrderTensor &tensor_, const tMatrixXd &matrix, size_t target_dim0, size_t target_dim1)
{
    cFourOrderTensor tensor = tensor_;
    tensor.TensorMatrixProductWithoutCopy(target_dim0, target_dim1, matrix);
    return tensor;
}

cFourOrderTensor MatrixTensorProductWithCopy(const tMatrixXd &matrix, const cFourOrderTensor &tensor_, size_t target_dim0, size_t target_dim1)
{

    cFourOrderTensor tensor = tensor_;
    tensor.MatrixTensorProductWithoutCopy(target_dim0, target_dim1, matrix);
    return tensor;
}
