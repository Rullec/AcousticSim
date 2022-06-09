#pragma once
#include "sim/softbody/SoftBodyImplicit.h"

SIM_DECLARE_CLASS_AND_PTR(tDiscretedWave);
SIM_DECLARE_CLASS_AND_PTR(cArrow);
class cAcousticSoftBody : public cSoftBodyImplicit
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    explicit cAcousticSoftBody(int id_);
    virtual void Init(const Json::Value &conf) override;
    virtual void Update(float dt) override;
    virtual void UpdateImGui() override;
    virtual int GetNumOfDrawTriangles() const override;
    virtual int GetNumOfDrawEdges() const override;
    virtual int GetNumOfDrawVertices() const override;
    virtual void CalcTriangleDrawBuffer(Eigen::Map<tVectorXf> &res,
                                        int &st) const override;
    virtual void CalcEdgeDrawBuffer(Eigen::Map<tVectorXf> &res, int &st) const override;
protected:
    double mAcousticSamplingFreqHZ;
    float mAcousticDuration;
    double mHammerForce;
    int mNumOfMaxModes;
    tVectorXd mEigenValues;
    tMatrixXd mEigenVecs;
    // material model idx: stvk, neoohookean or linear
    int mAcousticMaterialModelIdx;

    // material parameters idx: wood, steel, plasticity
    std::vector<std::string> mMaterialParamNameLst;
    std::vector<std::string> mMaterialParamPathLst;
    int mMaterialParamIdx;

    // selected vertex
    int mSelectedVertexForHearing;
    tEigenArr<tVector> mModalVibrationsInfo;
    std::vector<tDiscretedWavePtr> mModalWaves;
    void SolveForMonopole();
    std::string GetWaveName() const;
    std::vector<tVector3d> mSurfaceSoundPressureGrad;
    std::vector<int> mIsModeValid;
    tSparseMatd mLinearElasticityStiffnessMatrix;
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
    double CalcTotalMass() const;
    virtual void ChangeMaterial(int old_idx, int new_idx);
    virtual void LoadMaterialParams();
    virtual void DumpResultForTransfer();

    // call stack
    virtual void ResolveModalVibration();
    virtual void GenerateSound();

    std::vector<cArrowPtr> mDrawVertexNormalArray, mDrawTriNormalArray;
    bool mEnableDrawVertexNormal, mEnableDrawTriNormal;
    virtual void InitArrowFromNormal();
};