#pragma once
#include "utils/MathUtil.h"

// column-major storaged tensor
struct cFourOrderTensor
{

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    cFourOrderTensor(size_t dimI, size_t dimJ, size_t dimK, size_t dimL);
    cFourOrderTensor(const cFourOrderTensor &old_tensor);

    tMatrixXd &operator()(size_t i, size_t j);
    tMatrixXd operator()(size_t i, size_t j) const;
    cFourOrderTensor operator+(const cFourOrderTensor &tensor);
    cFourOrderTensor operator-(const cFourOrderTensor &tensor);
    cFourOrderTensor operator*(double scalar) const;
    cFourOrderTensor operator/(double scalar) const;
    tMatrixXd GetLastTwoComp(size_t k, size_t l);
    void SetLastTwoComp(size_t k, size_t l, const tMatrixXd &val);

    tMatrixXd ExpandToMatrix() const;
    void LoadFromAMatrix(size_t dimI, size_t dimJ, size_t dimK, size_t dimL, const tMatrixXd &);
    tVector4i GetShape() const;
    tEigenArr<tMatrixXd> GetInnerArray() const;

    // multiplication with copy

    // multiplication with itself
    void TensorMatrixProductWithoutCopy(size_t target_dim0, size_t target_dim1, const tMatrixXd &matrix); // matrix right product into the tensor

    void MatrixTensorProductWithoutCopy(size_t target_dim0, size_t target_dim1, const tMatrixXd &matrix); // matrix left product into the tensor

protected:
    size_t mDimI, mDimJ, mDimK, mDimL;
    tEigenArr<tMatrixXd> mInnerArray;
    void Allocate();
};

cFourOrderTensor TensorMatrixProductWithCopy(const cFourOrderTensor &tensor, const tMatrixXd &matrix, size_t target_dim0, size_t target_dim1);
cFourOrderTensor MatrixTensorProductWithCopy(const tMatrixXd &matrix, const cFourOrderTensor &tensor, size_t target_dim0, size_t target_dim1);