#pragma once
#include "utils/DefUtil.h"
#include "utils/MathUtil.h"
#include <memory>

enum ePoleType
{
    eMonopole = 0,
    eDipole = 1,
    // there is no tripole
    eQuadraple = 2,
    NUM_OF_POLE_TYPES
};

class cBasePole : public std::enable_shared_from_this<cBasePole>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    explicit cBasePole(ePoleType type, int _id, double omega,
                       const tVector3d &pos);
    virtual tVector3d CalcdPdx_Re(const tVector3d &pos) = 0;
    virtual void CheckdPdx_Re(const tVector3d &pos);
    virtual double CalcdPdn_Re(const tVector3d &pos,
                               const tVector3d &normal) = 0;
    virtual double CalcPressureForSoundSynthesis(const tVector3d &pos,
                                                 const tVectorXd &weight) = 0;
    virtual int GetNumOfDof() const;
    ePoleType mPoleType;
    int mId;
    tVector3d mPos;
    double mOmega; // 2 * pi * frequency

protected:
    virtual double CalcBaseItem(const tVector3d &pos) = 0;
};