#pragma once
#include "utils/MathUtil.h"

// column-major storaged tensor
struct cFourOrderTensor
{

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    cFourOrderTensor(size_t dimI, size_t dimJ, size_t dimK, size_t dimL) : mDimI(dimI), mDimJ(dimJ), mDimK(dimK), mDimL(dimL)
    {
        mInnerArray.resize(dimI * dimJ);
        for (auto &x : mInnerArray)
        {
            x.noalias() = tMatrixXd::Zero(dimK, dimL);
        }
    }
    tMatrixXd &operator()(size_t i, size_t j)
    {
        if (i >= mDimI || j >= mDimJ || i < 0 || j < 0)
        {
            SIM_ERROR("illegal access ({}, {}) for 4th order tensor shaped ({}, {}, {}, {})", i, j, mDimI, mDimJ, mDimK, mDimL);
            exit(1);
        }
        return mInnerArray[j * mDimI + i];
    }

    tMatrixXd ExpandToMatrix();
    void LoadFromAMatrix(size_t dimI, size_t dimJ, size_t dimK, size_t dimL, const tMatrixXd & )
protected:
    size_t mDimI, mDimJ, mDimK, mDimL;
    tEigenArray<tMatrixXd> mInnerArray;
}