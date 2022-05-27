#pragma once
#include "sim/BaseObject.h"
#include "utils/MathUtil.h"
#include "utils/SparseUtil.h"
SIM_DECLARE_CLASS_AND_PTR(tDiscretedWave);

class cAcousticBody
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    explicit cAcousticBody();
    tDiscretedWavePtr SolveVibration(
        const tEigenArr<tVector> & normal_array,
        const tVectorXd &MassDiag,
        const tSparseMatd &StiffMat,
        const tVector2f &rayleigh_damping,
        const tVectorXd &xcur,
        const tVectorXd &xprev,
        int sampling_freq,
        float duration);
    
protected:
};