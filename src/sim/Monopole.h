#pragma once
#include "utils/DefUtil.h"
#include "utils/MathUtil.h"

class cMonopole
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    explicit cMonopole(int _id, double omega);
    virtual void Init(double strength, const tVector3d &center_pos);
    virtual tVector3d EvaluatePressure(const tVector3d &pos);
    virtual tMatrix3d EvaulatedPdcenter(const tVector3d &pos);
    virtual tVector3d EvaluatedPdcoef(const tVector3d &pos);
    virtual void CheckGrad_dPdo();
    virtual void CheckGrad_dPdc();

    int mId;
    tVector3d mCenterPos;
    double mStrength;
    double mOmega; // 2 * pi * frequency
    struct
    {
        int mId;
        tVector3d mCenterPos;
        double mStrength;
        double mOmega; // 2 * pi * frequency
    } tPreState;
    virtual void PushState();
    virtual void PopState();
};

SIM_DECLARE_CLASS_AND_PTR(cMonopole);
