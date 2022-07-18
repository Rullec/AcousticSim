#include "BasePole.h"
#include "utils/LogUtil.h"

cBasePole::cBasePole(ePoleType type, int _id, double omega,
                     const tVector3d &pos)
    : mPoleType(type), mId(_id), mOmega(omega), mPos(pos)
{
}

void cBasePole::CheckdPdx_Re(const tVector3d &pos) {}

int cBasePole::GetNumOfDof() const
{
    int dof = 0;
    switch (mPoleType)
    {
    case ePoleType::eMonopole:
    {
        dof = 1;
        break;
    }
    case ePoleType::eDipole:
    {
        dof = 3;
        break;
    }
    case ePoleType::eQuadraple:
    {
        dof = 5;
        break;
    }
    default:
        SIM_ERROR("unsupported pole type {}", mPoleType);
        exit(1);
    }
    return dof;
}