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
    virtual void CheckDFDx();
    virtual void CheckDFDF();
    virtual void CheckDEDF();
    // virtual void CheckDPDF();
    // virtual void CheckDPDF_part1();
    // virtual void CheckDPDF_part2();
    virtual void CheckDPDx();
    virtual void CheckElementStiffnessMat();
    virtual void CheckGlobalStiffnessMat();

    
    virtual void SetTetVerticesPos(size_t tet_id, const tVectorXd &tet_vertices_pos);

    virtual void UpdateIntForce() override;
    virtual void SolveForNextPos(float dt) override;
    virtual tMatrixXd CalcElementStiffnessMatrix(int tet_id);
    virtual tMatrixXd CalcGlobalStiffnessMatrix();
    virtual void InitDDsDxTensor();

    tEigenArr<tMatrix3d> mDDsDxTensor_const; // tensor C_{ijk} for d(Ds)/dx
    tMatrixXd mGlobalStiffnessMatrix;
};