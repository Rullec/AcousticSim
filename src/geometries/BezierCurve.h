#pragma once
#include "utils/MathUtil.h"
/**
 * \brief       discretize a cubic beizer curve (determined by 4 points)
 */
class cBezierCurve
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    explicit cBezierCurve(int num_of_div, const tVector2d &A,
                          const tVector2d &B, const tVector2d &C,
                          const tVector2d &D);
    virtual double GetTotalLength() const;
    virtual int GetNumOfDrawEdges() const;
    virtual const tVectorXf &GetDrawBuffer();
    tEigenArr<tVector2d> GetPointList();

protected:
    tVectorXf mDrawBuffer;
    int mNumOfDiv;
    tVector2d A, B, C, D;
    tMatrixXd mPointList;
    virtual void InitPointlist(tMatrixXd &point_lst);
};