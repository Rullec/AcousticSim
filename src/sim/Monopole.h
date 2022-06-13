#pragma once
#include "utils/MathUtil.h"
#include "utils/DefUtil.h"

class cMonopole
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    explicit cMonopole(int _id, double omega);
    virtual void Init(double strength, const tVector3d &center_pos);
    virtual double EvaluatePressure(const tVector3d &pos);
    virtual tVector3d EvaluatePressureGrad(const tVector3d &pos);
    virtual void CheckGrad();
protected:
    int mId;
    tVector3d mCenterPos;
    double mStrength;
    double mOmega; // 2 * pi * frequency
};

SIM_DECLARE_CLASS_AND_PTR(cMonopole);
