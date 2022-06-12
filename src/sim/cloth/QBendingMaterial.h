#pragma once
#include "BaseBendingMaterial.h"
using tMatrix12f = Eigen::Matrix<float, 12, 12>;

class cQBendingMaterial : public cBaseBendingMaterial
{
public:
    explicit cQBendingMaterial();
    virtual void Init(const std::vector<tVertexPtr> &v_array,
                      const std::vector<tEdgePtr> &e_array,
                      const std::vector<tTrianglePtr> &t_array,
                      const tVector3d &bending_stiffness_warpweftbias) override;

    virtual void Update() override;
    double CalcEnergy(const tVectorXd &xcur) override;
    tVectorXd CalcForce(const tVectorXd &xcur) override;
    // std::vector<tMatrix12f> GetEleStiffnessMatrixLst() const;
    // std::vector<tVector4i> GetEdgeConstraintVertex() const;
    virtual void CheckForce() override;
    virtual void CheckStiffnessMatrix() override;

protected:
    // tSparseMatd mStiffnessMat;
    // std::vector<tMatrix12f> mEleKLst;
    // std::vector<tVector4i> mEdgeConstraintVertexLst;
    static tSparseMatd CalcKCoef(const tVector3d &v0, const tVector3d &v1,
                                 const tVector3d &v2, const tVector3d &v3);
};