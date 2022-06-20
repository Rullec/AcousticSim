#pragma once
#include "utils/DefUtil.h"
#include "utils/EigenUtil.h"
/**
 * \brief       storage the solution of mode vibration equation
 */
/*
    modal vibration equation:
    qddot + (a  + b * eigen_val) qdot + eigen_val * q = UTf_init

    for mode j, solution:
    q_j(t)  = dt * f_j^U / (w_d) e^{- \xi * w * t} sin(w_d * t)
            = coef_j / w_d e^{- \xi * w * t} sin(w_d * t)
    
    a_j(t) = q_j''(t) \approx =  coef_j * w_d e^{- \xi * w * t} sin(w_d * t)
    f^U is the generalized force

    we need to calculate
        1. coef_j = dt * f_J^U
        3. w = \sqrt(eigen_val)
        4. xi = (a + b * eigen_val) / (2 * w)
        2. w_d = w * (1 - xi^2)
*/
struct tModeVibration
{
    tModeVibration(const tVector &coef_vec);
    tModeVibration(double coef, double w, double xi, double wd);
    tVector GetCoefVec();
    tVectorXd GetWave(double duration_sec, double sampling_freq);

    double mCoef, mW, mXi, mWd;
};

SIM_DECLARE_STRUCT_AND_PTR(tModeVibration);