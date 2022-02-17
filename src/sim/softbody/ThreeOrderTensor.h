#pragma once
#include "utils/MathUtil.h"

struct cThreeOrderTensor
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    cThreeOrderTensor(const tEigenArr<tMatrixXd> &value_lst);
    cThreeOrderTensor(const tEigenArr<tMatrix3d> &value_lst);
    cThreeOrderTensor(size_t dimI, size_t dimJ, size_t dimK);
    tMatrixXd &operator()(size_t k);
    tMatrixXd operator()(size_t k) const;
    cThreeOrderTensor operator*(double scalar);
    tVector3i GetShape() const;

    tEigenArr<tMatrixXd> GetInnerArray() const;
    void SetInnerArray(const tEigenArr<tMatrixXd> &value_lst);
    void TensorMatrixProductWithoutCopy(size_t target_dim0, size_t target_dim1, const tMatrixXd &matrix); // matrix right product into the tensor
    void MatrixTensorProductWithoutCopy(size_t target_dim0, size_t target_dim1, const tMatrixXd &matrix); // matrix left product into the tensor

protected:
    size_t mDimI, mDimJ, mDimK;
    tEigenArr<tMatrixXd> mInnerArray;
};