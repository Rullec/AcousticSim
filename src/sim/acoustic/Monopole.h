#pragma once
#include "utils/DefUtil.h"
#include "utils/MathUtil.h"

/*
the sound pressure of monopole:

p(x) = A * exp(- i * k * r) / (4 * pi * r)

where r = |x - c_i|, c_i is the pole pos.

------------- sound pressure derivative on normal -------------
dp_j/dn_i = c_{ij} * A_j

c_{ij} can be calculated by "CalcCoef"
verified by "CheckdPdn"
*/
class cMonopole
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    explicit cMonopole(int _id, double omega);
    // virtual void Init(double strength, const tVector3d &center_pos);
    virtual double EvaluatePressure(const tVector3d &pos, double time);
    virtual tVector3d CalcdPdr(const tVector3d &pos);
    virtual void CheckdPdr();
    virtual double CalcdPdn(const tVector3d &pos, const tVector3d &normal);

    virtual double CalcCoef(const tVector3d &pos, const tVector3d &normal);
    virtual void CheckCoef();
    // virtual tMatrix3d EvaulatedPdcenter(const tVector3d &pos);
    // virtual tVector3d EvaluatedPdcoef(const tVector3d &pos);
    // virtual void CheckGrad_dPdo();
    // virtual void CheckGrad_dPdc();

    int mId;
    tVector3d mPos;
    double mA;     // amplitutde
    double mOmega; // 2 * pi * frequency
    struct
    {
        int mId;
        tVector3d mPos;
        double mA;
        double mOmega; // 2 * pi * frequency
    } tPreState;
    virtual void PushState();
    virtual void PopState();

protected:
    virtual double CalcBaseItem(const tVector3d &pos);
};

SIM_DECLARE_CLASS_AND_PTR(cMonopole);
