#include "FourOrderTensor.h"
#include "utils/LogUtil.h"

cFourOrderTensor::cFourOrderTensor(size_t dimI, size_t dimJ, size_t dimK, size_t dimL) : mDimI(dimI), mDimJ(dimJ), mDimK(dimK), mDimL(dimL)
{
    Allocate();
}
void cFourOrderTensor::Allocate()
{
    mInnerArray.resize(mDimI * mDimJ);
    for (auto &x : mInnerArray)
    {
        x.noalias() = tMatrixXd::Zero(mDimK, mDimL);
    }
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

/**
 * \brief           expand a tensor T_{ijkl} to M_{pq}
 *      M_{i * prod(jkl) + j * prod(kl) + k * prod(l) + l} = T_{ijkl}
 *      For more details, please check "四阶张量转换为矩阵.md"
*/
tMatrixXd cFourOrderTensor::ExpandToMatrix()
{
    int rows = mDimI * mDimJ;
    int cols = mDimK * mDimL;
    printf("expand to matrix (%d, %d)\n", rows, cols);
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