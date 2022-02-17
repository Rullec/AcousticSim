#pragma once
#include "utils/MathUtil.h"

// column-major storaged tensor
struct cFourOrderTensor
{

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    cFourOrderTensor(size_t dimI, size_t dimJ, size_t dimK, size_t dimL);
    tMatrixXd &operator()(size_t i, size_t j);
    
    tMatrixXd ExpandToMatrix();
    void LoadFromAMatrix(size_t dimI, size_t dimJ, size_t dimK, size_t dimL, const tMatrixXd &);
    tVector4i GetShape() const;
protected:
    size_t mDimI, mDimJ, mDimK, mDimL;
    tEigenArr<tMatrixXd> mInnerArray;
    void Allocate();
};