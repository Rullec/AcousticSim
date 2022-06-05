#pragma once
#include "sim/softbody/SoftBodyImplicit.h"

class cAcousticSoftBody : public cSoftBodyImplicit
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    explicit cAcousticSoftBody(int id_);
    virtual void Init(const Json::Value &conf) override;
    virtual void Update(float dt) override;

protected:
    int mAcousticSamplingFreq;
    float mAcousticDuration;
    void SolveForMonopole();
    std::string GetWaveName() const;
    tEigenArr<tMatrixXd> SolveVibration(
                   const tVectorXd &MassDiag, const tSparseMatd &StiffMat,
                   const tVector2f &rayleigh_damping, const tVectorXd &xcur,
                   const tVectorXd &xprev, int sampling_freq, float duration);
};