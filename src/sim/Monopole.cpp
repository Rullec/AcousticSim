#include "sim/Monopole.h"
#include "utils/DefUtil.h"
#include <iostream>

cMonopole::cMonopole(int _id, double omega)
{
    mId = _id;
    mCenterPos.setZero();
    mStrength = 1.0;
    mOmega = omega;
}

void cMonopole::Init(double strength, const tVector3d &center_pos)
{

    mCenterPos = center_pos;
}

/*
 *
 */
tVector3d cMonopole::EvaluatePressure(const tVector3d &pos)
{
    tVector3d r = pos - this->mCenterPos;
    double r_norm = r.norm();
    r_norm = SIM_MAX(r_norm, 1e-10);
    tVector3d r_bar = r / r_norm;
    return mStrength / (4 * M_PI * r_norm) * r_bar;
}

/*
- cm / (4 * pi) * (I - 2 r_bar r_bar^T ) / (r^2)
*/
tMatrix3d cMonopole::EvaulatedPdcenter(const tVector3d &pos)
{
    tVector3d r = pos - this->mCenterPos;
    double r_norm = r.norm();
    r_norm = SIM_MAX(r_norm, 1e-10);
    tVector3d r_bar = r / r_norm;
    tMatrix3d grad = -this->mStrength / (4 * M_PI) *
                     (tMatrix3d::Identity() - 2 * r_bar * r_bar.transpose()) /
                     (r_norm * r_norm);
    return grad;
}

/*
r_bar / (4 * pi * r_norm)
*/
tVector3d cMonopole::EvaluatedPdcoef(const tVector3d &pos)
{
    tVector3d r = pos - this->mCenterPos;
    double r_norm = r.norm();
    r_norm = SIM_MAX(r_norm, 1e-10);
    tVector3d r_bar = r / r_norm;
    return r_bar / (4 * M_PI * r_norm);
}

void cMonopole::CheckGrad_dPdo()
{
    PushState();
    this->mCenterPos.setRandom();
    tVector3d x_pos = tVector3d::Random();
    // 1. get analytic gradient
    tMatrix3d dPdcenter_ana = EvaulatedPdcenter(x_pos);
    tMatrix3d dPdcenter_num = tMatrix3d::Zero();

    // 2. get numberical gradient
    double eps = 1e-6;
    tVector3d old_p = EvaluatePressure(x_pos);
    for (int i = 0; i < 3; i++)
    {
        mCenterPos[i] += eps;
        tVector3d new_p = EvaluatePressure(x_pos);
        tVector3d dpdo_num = (new_p - old_p) / eps;

        dPdcenter_num.col(i) = dpdo_num;

        mCenterPos[i] -= eps;
    }
    auto diff = dPdcenter_ana - dPdcenter_num;
    double diff_norm = diff.norm();
    std::cout << "diff norm = " << diff_norm << std::endl;
    std::cout << "dPdcenter_num = \n" << dPdcenter_num << std::endl;
    std::cout << "dPdcenter_ana = \n" << dPdcenter_ana << std::endl;
    PopState();
}
void cMonopole::CheckGrad_dPdc()
{
    PushState();
    this->mCenterPos.setRandom();
    tVector3d x_pos = tVector3d::Random();
    // 1. get analytic gradient
    tVector3d dPdcoef_ana = EvaluatedPdcoef(x_pos);

    // 2. get numberical gradient
    double eps = 1e-6;
    tVector3d old_p = EvaluatePressure(x_pos);

    mStrength += eps;
    tVector3d new_p = EvaluatePressure(x_pos);
    tVector3d dPdcoef_num = (new_p - old_p) / eps;

    auto diff = dPdcoef_ana - dPdcoef_num;
    double diff_norm = diff.norm();
    std::cout << "diff norm = " << diff_norm << std::endl;
    std::cout << "dPdcoef_num = \n" << dPdcoef_num << std::endl;
    std::cout << "dPdcoef_ana = \n" << dPdcoef_ana << std::endl;
    PopState();
}

void cMonopole::PushState()
{
    tPreState.mId = mId;
    tPreState.mCenterPos = mCenterPos;
    tPreState.mOmega = mOmega;
    tPreState.mStrength = mStrength;
}
void cMonopole::PopState()
{
    mId = tPreState.mId;
    mCenterPos = tPreState.mCenterPos;
    mOmega = tPreState.mOmega;
    mStrength = tPreState.mStrength;
}