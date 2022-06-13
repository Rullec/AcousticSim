#include "ModeVibration.h"

tModeVibration::tModeVibration(const tVector &coef_vec)
    : tModeVibration(coef_vec[0], coef_vec[1], coef_vec[2], coef_vec[3])
{
}

tModeVibration::tModeVibration(double coef, double w, double xi, double wd)
    : mCoef(coef), mW(w), mXi(xi), mWd(wd)
{
}

tVector tModeVibration::GetCoefVec() { return tVector(mCoef, mW, mXi, mWd); }
tVectorXd tModeVibration::GetWave(double duration_sec, double sampling_freq)
{
    double dt = 1.0 / sampling_freq;
    double num_of_steps = duration_sec * sampling_freq;
    // q(t) = coef_j / w_d e^{- \xi * w * t} sin(w_d * t)
    tVectorXd data = tVectorXd::Zero(num_of_steps);
    for (int i = 0; i < num_of_steps; i++)
    {
        double t = i * dt;
        data[i] = mCoef / mWd * std::exp(-mXi * mW * t) * std::sin(mWd * t);
    }
    return data;
}