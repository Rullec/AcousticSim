#include "ThreeOrderTensor.h"
#include "utils/LogUtil.h"
cThreeOrderTensor::cThreeOrderTensor(const tEigenArr<tMatrixXd> &value_lst)
{
    SetInnerArray(value_lst);
}
cThreeOrderTensor::cThreeOrderTensor(const tEigenArr<tMatrix3d> &value_lst)
{
    tEigenArr<tMatrixXd> value_lst_xd = {};
    for (auto &x : value_lst)
    {
        value_lst_xd.push_back(x);
    }
    SetInnerArray(value_lst_xd);
}
cThreeOrderTensor::cThreeOrderTensor(size_t dimI, size_t dimJ, size_t dimK) : mDimI(dimI), mDimJ(dimJ), mDimK(dimK)
{
    mInnerArray.resize(dimK, tMatrixXd::Zero(dimI, dimJ));
}

tMatrixXd &cThreeOrderTensor::operator()(size_t k)
{
    if (k >= mDimK || k < 0)
    {
        SIM_ERROR("illegal access {} for range [0, {}]", k, mDimK - 1);
        exit(1);
    }
    return mInnerArray[k];
}

tMatrixXd cThreeOrderTensor::operator()(size_t k) const
{
    if (k >= mDimK || k < 0)
    {
        SIM_ERROR("illegal access {} for range [0, {}]", k, mDimK - 1);
        exit(1);
    }
    return mInnerArray[k];
}

tVector3i cThreeOrderTensor::GetShape() const
{
    return tVector3i(mDimI, mDimJ, mDimK);
}

tEigenArr<tMatrixXd> cThreeOrderTensor::GetInnerArray() const
{
    return mInnerArray;
}

void cThreeOrderTensor::SetInnerArray(const tEigenArr<tMatrixXd> &value_lst)
{
    mInnerArray = value_lst;
    mDimI = mInnerArray[0].rows();
    mDimJ = mInnerArray[1].rows();
    mDimK = mInnerArray.size();
    for (auto &x : mInnerArray)
    {
        if (x.rows() != mDimI || x.cols() != mDimJ)
        {
            SIM_ERROR("the array is inconsistent\n");
            exit(1);
        }
    }
}

void cThreeOrderTensor::TensorMatrixProductWithoutCopy(size_t target_dim0, size_t target_dim1, const tMatrixXd &matrix)
{
    // matrix right product into the tensor
    if (target_dim0 != 0 || target_dim1 != 1)
    {
        SIM_ERROR("only support matrix contraction on the first two indices");
        exit(1);
    }
    for (auto &x : mInnerArray)
    {
        x = x * matrix;
    }
    mDimJ = matrix.cols();
}

void cThreeOrderTensor::MatrixTensorProductWithoutCopy(size_t target_dim0, size_t target_dim1, const tMatrixXd &matrix)
{
    // matrix left product into the tensor
    if (target_dim0 != 0 || target_dim1 != 1)
    {
        SIM_ERROR("only support matrix contraction on the first two indices");
        exit(1);
    }
    for (auto &x : mInnerArray)
    {
        x = matrix * x;
    }
    mDimI = matrix.rows();
}

cThreeOrderTensor cThreeOrderTensor::operator*(double scalar)
{
    cThreeOrderTensor new_tensor(*this);
    for (auto &x : mInnerArray)
        x *= scalar;

    return new_tensor;
}
