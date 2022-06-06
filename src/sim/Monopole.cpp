#include "sim/Monopole.h"
#include "utils/DefUtil.h"

cMonopole::cMonopole(int _id)
{
    mId = _id;
    mCenterPos.setZero();
    mStrength = 0;
}

void cMonopole::Init(double strength, const tVector3d &center_pos)
{
    mStrength = strength;
    mCenterPos = center_pos;
}

/*
p(x) = e^{-i * k * r} / (4 * pi * r)
r = |pos - center|

p(x) is a complex number, p(x) = |p(x)| e^{i \phi}

and |p(x)| is sound pressure, e^{i \phi} is phase shift

|p(x)| = 1.0 / (4 * pi * r)
e^{i \phi} = e^{- i * k * r}

p(x, t) = p(x) e^{i \omega t}
*/

double cMonopole::EvaluatePressure(const tVector3d &pos)
{
    double x_minus_r = (pos - mCenterPos).norm();
    return 1.0 / (4 * M_PI * x_minus_r);
}

/*
d|p(x)|/dx = - (x - r) / (4 * pi * |x-r|^3)
*/
tVector3d cMonopole::EvaluatePressureGrad(const tVector3d &pos)
{
    double x_minus_r = (pos - mCenterPos).norm();
    double deno = 4 * M_PI * x_minus_r * x_minus_r * x_minus_r;
    deno = SIM_MAX(1e-10, deno);
    return (mCenterPos - pos) / (deno);
}