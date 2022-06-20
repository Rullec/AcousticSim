#include "sim/acoustic/Monopole.h"
// #include "utils/DefUtil.h"
#include "utils/LogUtil.h"
#include <iostream>
double gWaveSpeed = 340;
cMonopole::cMonopole(int _id, double omega)
{
    mId = _id;
    mPos.setRandom();
    mA = cMathUtil::RandDouble(0, 1);
    mOmega = omega;
    // CheckCoef();
    // CheckdPdr();
}

// void cMonopole::Init(double strength, const tVector3d &center_pos)
// {
//     mA = strength;
//     mPos = center_pos;
    
// }

/*
c_ij = r.dot(n) * base_item
*/
double cMonopole::CalcCoef(const tVector3d &pos, const tVector3d &normal_)
{
    tVector3d r_vec = pos - mPos;
    tVector3d n = normal_.normalized();

    return CalcBaseItem(pos) * r_vec.dot(n);
}

/*
base_item = -1 / (4 * pi * |r|^3) * (cos k r + kr * sin(kr))
*/
double cMonopole::CalcBaseItem(const tVector3d &pos)
{
    SIM_ASSERT(pos.hasNaN() == false);
    tVector3d r_vec = pos - mPos;
    double r = r_vec.norm();
    r = SIM_MAX(r, 1e-10);

    // k = w / c
    double kr = r * mOmega / gWaveSpeed;

    double cij =
        -1 / (4 * M_PI * r * r * r) * (std::cos(kr) + kr * std::sin(kr));
    return cij;
}

/*
dPdr \in R^3
dPdr = A * \bold r * base_item
*/
tVector3d cMonopole::CalcdPdr(const tVector3d &pos)
{
    tVector3d r_vec = pos - mPos;
    return mA * r_vec * CalcBaseItem(pos);
}

/*
    dp_ij/dn_i = dPdr.dot(\bold n)
*/
double cMonopole::CalcdPdn(const tVector3d &pos, const tVector3d &normal)
{
    tVector3d n = normal.normalized();
    return CalcdPdr(pos).dot(n);
}

void cMonopole::CheckdPdr()
{
    // 1. set up pos
    double time = 0;
    tVector3d test_pos = tVector3d::Random();
    double eps = 1e-6;

    // 2. calculate ana dPdr
    tVector3d dPdr_ana = CalcdPdr(test_pos);
    tVector3d dPdr_num = tVector3d::Zero();

    // 3. calculate old p
    double old_p = EvaluatePressure(test_pos, time);

    // 4. change to new pos, calculate new p
    for (int i = 0; i < 3; i++)
    {
        test_pos[i] += eps;
        double new_p = EvaluatePressure(test_pos, time);
        dPdr_num[i] = (new_p - old_p) / eps;
        test_pos[i] -= eps;
    }
    tVector3d diff = (dPdr_num - dPdr_ana).cwiseAbs();
    std::cout << "dPdr num = " << dPdr_num.transpose() << std::endl;
    std::cout << "dPdr ana = " << dPdr_ana.transpose() << std::endl;
    std::cout << "diff = " << diff.transpose() << std::endl;
    // 5. calculate num dPdr
}
/*
p   = A/(4 * pi * r) * exp(i(w t - kr))
    = A/(4 * pi * r) * cos(wt - kr)
*/
double cMonopole::EvaluatePressure(const tVector3d &pos, double time)
{
    double r = SIM_MAX((pos - mPos).norm(), 1e-10);
    double k = mOmega / gWaveSpeed;
    return mA / (4 * M_PI * r) * std::cos(mOmega * time - k * r);
}

/*
confirm coef:

coef * A = dpdn
*/
void cMonopole::CheckCoef()
{
    tVector3d pos = tVector3d::Random();
    tVector3d n = tVector3d::Random();
    double coef_A = CalcCoef(pos, n) * mA;
    double dpdn = this->CalcdPdn(pos, n);
    double diff = std::fabs(coef_A - dpdn);
    printf("coef A %.3f, dpdn %.3f, diff %.3f\n", coef_A, dpdn, diff);
}
// /*
//  *
//  */
// tVector3d cMonopole::EvaluatePressure(const tVector3d &pos)
// {
//     tVector3d r = pos - this->mCenterPos;
//     double r_norm = r.norm();
//     r_norm = SIM_MAX(r_norm, 1e-10);
//     tVector3d r_bar = r / r_norm;
//     return mStrength / (4 * M_PI * r_norm) * r_bar;
// }

// /*
// - cm / (4 * pi) * (I - 2 r_bar r_bar^T ) / (r^2)
// */
// tMatrix3d cMonopole::EvaulatedPdcenter(const tVector3d &pos)
// {
//     tVector3d r = pos - this->mCenterPos;
//     double r_norm = r.norm();
//     r_norm = SIM_MAX(r_norm, 1e-10);
//     tVector3d r_bar = r / r_norm;
//     tMatrix3d grad = -this->mStrength / (4 * M_PI) *
//                      (tMatrix3d::Identity() - 2 * r_bar * r_bar.transpose())
//                      / (r_norm * r_norm);
//     return grad;
// }

// /*
// r_bar / (4 * pi * r_norm)
// */
// tVector3d cMonopole::EvaluatedPdcoef(const tVector3d &pos)
// {
//     tVector3d r = pos - this->mCenterPos;
//     double r_norm = r.norm();
//     r_norm = SIM_MAX(r_norm, 1e-10);
//     tVector3d r_bar = r / r_norm;
//     return r_bar / (4 * M_PI * r_norm);
// }

// void cMonopole::CheckGrad_dPdo()
// {
//     PushState();
//     this->mCenterPos.setRandom();
//     tVector3d x_pos = tVector3d::Random();
//     // 1. get analytic gradient
//     tMatrix3d dPdcenter_ana = EvaulatedPdcenter(x_pos);
//     tMatrix3d dPdcenter_num = tMatrix3d::Zero();

//     // 2. get numberical gradient
//     double eps = 1e-6;
//     tVector3d old_p = EvaluatePressure(x_pos);
//     for (int i = 0; i < 3; i++)
//     {
//         mCenterPos[i] += eps;
//         tVector3d new_p = EvaluatePressure(x_pos);
//         tVector3d dpdo_num = (new_p - old_p) / eps;

//         dPdcenter_num.col(i) = dpdo_num;

//         mCenterPos[i] -= eps;
//     }
//     auto diff = dPdcenter_ana - dPdcenter_num;
//     double diff_norm = diff.norm();
//     std::cout << "diff norm = " << diff_norm << std::endl;
//     std::cout << "dPdcenter_num = \n" << dPdcenter_num << std::endl;
//     std::cout << "dPdcenter_ana = \n" << dPdcenter_ana << std::endl;
//     PopState();
// }
// void cMonopole::CheckGrad_dPdc()
// {
//     PushState();
//     this->mCenterPos.setRandom();
//     tVector3d x_pos = tVector3d::Random();
//     // 1. get analytic gradient
//     tVector3d dPdcoef_ana = EvaluatedPdcoef(x_pos);

//     // 2. get numberical gradient
//     double eps = 1e-6;
//     tVector3d old_p = EvaluatePressure(x_pos);

//     mStrength += eps;
//     tVector3d new_p = EvaluatePressure(x_pos);
//     tVector3d dPdcoef_num = (new_p - old_p) / eps;

//     auto diff = dPdcoef_ana - dPdcoef_num;
//     double diff_norm = diff.norm();
//     std::cout << "diff norm = " << diff_norm << std::endl;
//     std::cout << "dPdcoef_num = \n" << dPdcoef_num << std::endl;
//     std::cout << "dPdcoef_ana = \n" << dPdcoef_ana << std::endl;
//     PopState();
// }

void cMonopole::PushState()
{
    tPreState.mId = mId;
    tPreState.mPos = mPos;
    tPreState.mOmega = mOmega;
    tPreState.mA = mA;
}
void cMonopole::PopState()
{
    mId = tPreState.mId;
    mPos = tPreState.mPos;
    mOmega = tPreState.mOmega;
    mA = tPreState.mA;
}