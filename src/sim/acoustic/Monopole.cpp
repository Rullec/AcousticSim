#include "sim/acoustic/Monopole.h"
// #include "utils/DefUtil.h"
#include "utils/LogUtil.h"
#include <iostream>
double gWaveSpeed = 340;
cMonopole::cMonopole(int _id, double omega, const tVector3d &pos)
    : cBasePole(ePoleType::eMonopole, _id, omega, pos)
{
    mId = _id;
    mOmega = omega;
}

// void cMonopole::Init(double strength, const tVector3d &center_pos)
// {
//     mA = strength;
//     mPos = center_pos;

// }

/*
base_item = -1 / (4 * pi * |r|^3) * (cos k r + kr * sin(kr))
*/
double cMonopole::CalcBaseItem(const tVector3d &pos)
{
    // SIM_ASSERT(pos.hasNaN() == false);
    // tVector3d r_vec = pos - mPos;
    // double r = r_vec.norm();
    // r = SIM_MAX(r, 1e-10);

    // // k = w / c
    // double kr = r * mOmega / gWaveSpeed;

    // double cij =
    //     -1 / (4 * M_PI * r * r * r) * (std::cos(kr) + kr * std::sin(kr));
    // return cij;
    return 0;
}

/*
dPdr \in R^3
dPdr = - (x - c) / (4 * pi * r^3)
*/
tVector3d cMonopole::CalcdPdx_Re(const tVector3d &pos)
{
    tVector3d r_vec = pos - mPos;
    double r = r_vec.norm();
    r = SIM_MAX(r, 1e-10);

    double kr = r * mOmega / gWaveSpeed;
    // return -r_vec * (std::cos(kr) + kr * std::sin(kr)) / (4 * M_PI * r * r *
    // r);
    return -r_vec / (4 * M_PI * r * r * r);
}

/*
    dp_ij/dn_i = dPdr.dot(\bold n)
*/
double cMonopole::CalcdPdn_Re(const tVector3d &pos, const tVector3d &normal)
{
    tVector3d n = normal.normalized();
    return CalcdPdx_Re(pos).dot(n);
}

double cMonopole::CalcPressureForSoundSynthesis(const tVector3d &pos,
                                                const tVectorXd &weight)
{
    SIM_ASSERT(weight.size() == 1);
    double r = (pos - mPos).norm();
    return weight[0] / (4 * M_PI * r);
}