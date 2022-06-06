#pragma once
#include "sim/softbody/SoftBodyImplicit.h"

SIM_DECLARE_CLASS_AND_PTR(tDiscretedWave);
class cAcousticSoftBody : public cSoftBodyImplicit
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    explicit cAcousticSoftBody(int id_);
    virtual void Init(const Json::Value &conf) override;
    virtual void Update(float dt) override;

protected:
    double mAcousticSamplingFreqHZ;
    float mAcousticDuration;
    double mHammerForce;
    int mNumOfMaxModes;
    tVectorXd mEigenValues;
    tMatrixXd mEigenVecs;
    
    tEigenArr<tVector> mModalVibrationsInfo;
    std::vector<tDiscretedWavePtr> mModalWaves;
    void SolveForMonopole();
    std::string GetWaveName() const;
    std::vector<tVector3d> mSurfaceSoundPressureGrad;

    tEigenArr<tMatrixXd>
    SolveVibration(const tVectorXd &MassDiag, const tSparseMatd &StiffMat,
                   const tVector2f &rayleigh_damping, const tVectorXd &xcur,
                   const tVectorXd &xprev, int sampling_freq, float duration);
    virtual void EigenDecompose(tVectorXd &eigenValues,
                                tMatrixXd &eigenVecs) const;
    virtual void CalculateModalVibration();
    virtual tDiscretedWavePtr CalculateVertexVibration(int v_id);
    int GetNumOfModes() const;
    double GetDt() const;
};