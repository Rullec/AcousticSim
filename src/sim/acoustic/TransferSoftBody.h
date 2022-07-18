#pragma once
#include "sim/kinematic/KinematicBody.h"

SIM_DECLARE_CLASS_AND_PTR(tDiscretedWave);
SIM_DECLARE_CLASS_AND_PTR(tAnalyticWave);
SIM_DECLARE_CLASS_AND_PTR(cArrow);
SIM_DECLARE_CLASS_AND_PTR(cMonopole);
SIM_DECLARE_STRUCT_AND_PTR(tModeVibration);
using cPoleArray = std::vector<cMonopolePtr>;

class cTransferSoftBody : public cKinematicBody
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    explicit cTransferSoftBody(int id_);
    virtual void Init(const Json::Value &conf) override;
    virtual void Update(float dt) override;
    virtual void UpdateImGui() override;
    virtual void SetCameraPos(const tVector3d &pos);

protected:
    tVector3d mCamPos;
    tVector3d mCurSolvedCamPos;
    // double mOmega, mC;
    // tAnalyticWavePtr mAnaWave;
    // cMonopolePtr mPole;

    // modal analysis info
    std::string mModalAnalysisPath;
    std::vector<tModeVibrationPtr>
        mModalVibrationWaveArray; // modal vibration wave for each mode
    tMatrixXd mVertexModesCoef;   // vertex combination coef for each modes
    tVector mAABBMin, mAABBMax;   // object AABB min and max
    int mNumOfPolesPerMode;

    // transfer info
    std::vector<cPoleArray> mPolesArrayArray;
    tEigenArr<tVectorXd> mPolesWeightArray;
    virtual void LoadModalAnalysisResult(const Json::Value &root);

    // init pole position
    virtual void InitPole(cPoleArray &pole, int mode_idx);

    // calculate vertex acceleration on each mode (b)
    virtual tVectorXd CalcVertexAccel(int mode_idx);

    // calculate sound pressure grad (A)
    virtual tMatrixXd CalcSoundPressureGradient(const cPoleArray &pole_array,
                                                int mode_idx);
    // solve pole amp: x = A.inv * b
    virtual tVectorXd SolveForTargetMode(int mode_idx);

    // generate sound from sources
    virtual void
    GenerateSound(const std::vector<cPoleArray> &pole_array,
                       const tEigenArr<tVectorXd> &pole_weight_array);
};