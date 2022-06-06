#pragma once
#include "utils/MathUtil.h"

class cMonopole
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    explicit cMonopole(int _id);
    virtual void Init(double strength, const tVector3d &center_pos);
    virtual double EvaluatePressure(const tVector3d &pos);
    virtual tVector3d EvaluatePressureGrad(const tVector3d &pos);

protected:
    int mId;
    tVector3d mCenterPos;
    double mStrength;
};