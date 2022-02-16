#pragma once

#include "SoftBody.h"

class cSoftBodyImplicit : public cSoftBody
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    explicit cSoftBodyImplicit(int id_);
    virtual ~cSoftBodyImplicit();
    virtual void Update(float dt) override;
    virtual void Init(const Json::Value &conf) override;

protected:
    virtual tEigenArr<tMatrix3d> CalcDFDx(int ele_id);
    virtual void CalcDPDF();
    virtual void CheckDFDx();
    virtual void CheckDPDF();
    virtual void UpdateIntForce() override;
    virtual void SolveForNextPos(float dt) override;
    virtual void CalculateStiffnessMatrix();
    tEigenArr<tMatrix3d> mDDsDxTensor_const;    // tensor C_{ijk} for d(Ds)/dx
    virtual void InitDDsDxTensor();
};