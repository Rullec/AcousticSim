#pragma once
#include "BasePole.h"

/*
the sound pressure of dipole:

p(x) = 

where r = |x - c_i|, c_i is the pole pos.

------------- sound pressure derivative on normal -------------
dp_j/dn_i = c_{ij}

c_{ij} can be calculated by "CalcCoef"
verified by "CheckdPdn"
*/
class cMonopole : public cBasePole
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    explicit cMonopole(int _id, double omega, const tVector3d &pos);
    virtual tVector3d CalcdPdx_Re(const tVector3d &pos) override;
    // virtual double CalcPressure(const tVector3d &pos) override;
    virtual double CalcPressureForSoundSynthesis(const tVector3d &pos, const tVectorXd & weight) override;
    virtual double CalcdPdn_Re(const tVector3d &pos,
                            const tVector3d &normal) override;

protected:
    virtual double CalcBaseItem(const tVector3d &pos) override;
};

SIM_DECLARE_CLASS_AND_PTR(cMonopole);
