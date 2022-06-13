#pragma once
#include "sim/softbody/SoftBodyImplicit.h"

SIM_DECLARE_CLASS_AND_PTR(tDiscretedWave);
SIM_DECLARE_CLASS_AND_PTR(cArrow);
SIM_DECLARE_CLASS_AND_PTR(cMonopole);
SIM_DECLARE_STRUCT_AND_PTR(tModeVibration);

class cTransferSoftBody : public cSoftBodyImplicit
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    explicit cTransferSoftBody(int id_);
    virtual void Init(const Json::Value &conf) override;
    virtual void Update(float dt) override;
    virtual void UpdateImGui() override;
    virtual int GetNumOfDrawTriangles() const override;
    virtual int GetNumOfDrawEdges() const override;
    virtual int GetNumOfDrawVertices() const override;
    virtual void CalcTriangleDrawBuffer(Eigen::Map<tVectorXf> &res,
                                        int &st) const override;
    virtual void CalcEdgeDrawBuffer(Eigen::Map<tVectorXf> &res,
                                    int &st) const override;

protected:
    virtual void SolveMonopole(int mode_idx);
    virtual void LoadModalAnalysisResult();
    virtual void InitSurfaceNormal();
    virtual tVectorXd CalcSoundPressure(int mode_idx);
    virtual void SetModalVibrationSound();
    tEigenArr<tModeVibrationPtr> mModalVibrationArray;
    tMatrixXd mSurfaceVertexDOFCoef;
    std::string mModalAnlysisResultPath;
    bool mEnableDumpedModalSolution; // use undamped & dumped modal analysis
                                     // result
    int mNumOfMonopolesPerFreq;      // number of monopoles
    tEigenArr<tVector> mSurfaceVertexNormalArray;
    struct tPolesPerFreq
    {
        std::vector<cMonopolePtr> mPoles;
        double mOmega;
    };
    std::vector<tPolesPerFreq> mPolesArray;
    virtual void InitPole();
    virtual double GetEnergy(int mode_idx);
    virtual tVectorXd GetGrad(int mode_idx);
    virtual tVectorXd GetX(int mode_idx);
};