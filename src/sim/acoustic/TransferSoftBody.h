#pragma once
#include "sim/kinematic/KinematicBody.h"

SIM_DECLARE_CLASS_AND_PTR(tDiscretedWave);
SIM_DECLARE_CLASS_AND_PTR(tAnalyticWave);
SIM_DECLARE_CLASS_AND_PTR(cArrow);
SIM_DECLARE_CLASS_AND_PTR(cMonopole);
SIM_DECLARE_STRUCT_AND_PTR(tModeVibration);

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
    double mOmega, mC;
    tAnalyticWavePtr mAnaWave;
    cMonopolePtr mPole;
    std::vector<tModeVibrationPtr> mModalVibrationWaveArray;
    std::string mModalAnalysisResult;
    tMatrixXd mVertexModesCoef;
    std::vector<std::vector<cMonopolePtr>> pole_array_array;
    tEigenArr<tVectorXd> pole_weight_array;
    virtual void ResolveSoundByCamPos();
    virtual void LoadModalAnalysisResult();
    virtual void SetModalAnalysisSound();
    virtual std::vector<cMonopolePtr> SolveForTargetMode(int tar_mode,
                                                         tVectorXd &weight);
    virtual void
    SolveSoundForPoles(const std::vector<std::vector<cMonopolePtr>> &pole_array,
                       const tEigenArr<tVectorXd> &pole_weight_array);
};