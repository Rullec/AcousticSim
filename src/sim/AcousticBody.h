#pragma once
#include "sim/BaseObject.h"
#include "utils/MathUtil.h"

class cAcousticBody
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    explicit cAcousticBody();
    void SolveVibration(
        const tVectorXd &MassDiag,
        const tSparseMatd &StiffMat,
        const tVector2f &rayleigh_damping,
        const tVectorXd &xcur,
        const tVectorXd &xprev);
    
protected:
};