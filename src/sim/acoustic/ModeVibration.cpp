#include "ModeVibration.h"

tModeVibration::tModeVibration(const tVector &coef_vec)
    : tModeVibration(coef_vec[0], coef_vec[1], coef_vec[2], coef_vec[3])
{
}

tModeVibration::tModeVibration(double coef, double w, double xi, double wd)
    : mCoef(coef), mW(w), mXi(xi), mWd(wd)
{
    mCurDurationSec = 0;
    mCurSamplingFreq = 0;
}

tVector tModeVibration::GetCoefVec() { return tVector(mCoef, mW, mXi, mWd); }
tVectorXd tModeVibration::GetWave(double duration_sec, double sampling_freq)
{
    return mCoef / mWd * GetWaveExpAndSin(duration_sec, sampling_freq);
}

tVectorXd tModeVibration::GetWaveExpAndSin(double duration_sec,
                                           double sampling_freq)
{
    if (std::fabs(mCurDurationSec - duration_sec) > 1e-3 ||
        std::fabs(sampling_freq - mCurSamplingFreq) > 1e-3)
    {

        double dt = 1.0 / sampling_freq;
        double num_of_steps = duration_sec * sampling_freq;
        // q(t) = coef_j / w_d e^{- \xi * w * t} sin(w_d * t)
        mCurDiscreteData = tVectorXd::Zero(num_of_steps);
        for (int i = 0; i < num_of_steps; i++)
        {
            double t = i * dt;
            mCurDiscreteData[i] = std::exp(-mXi * mW * t) * std::sin(mWd * t);
        }
    }
    return mCurDiscreteData;
}