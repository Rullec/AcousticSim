#pragma once
#include "BaseBendingMaterial.h"
// using tMatrix12f = Eigen::Matrix<float, 12, 12>;
// using tVector12f = Eigen::Matrix<float, 12, 1>;
using tMatrix12d = Eigen::Matrix<double, 12, 12>;
using tVector12d = Eigen::Matrix<double, 12, 1>;

class cDihedralMaterial : public cBaseBendingMaterial
{
public:
    explicit cDihedralMaterial();
    virtual void Init(const std::vector<tVertexPtr> &v_array,
                      const std::vector<tEdgePtr> &e_array,
                      const std::vector<tTrianglePtr> &t_array,
                      const tVector3d &bending_stiffness_warpweftbias);
    virtual void Update() override;
    virtual double CalcEnergy(const tVectorXd &xcur) override;
    virtual tVectorXd CalcForce(const tVectorXd &xcur) override;
    virtual void CheckForce() override;
    virtual void CheckStiffnessMatrix() override;

protected:
    std::vector<float> mRawEdgeLengthArray;
    tEigenArr<tVector2f> mRawHeightArray;
    std::vector<float> mRawTriangleAreaArray;
    tEigenArr<tVector12d> mEleForceArray;
    tVectorXd mGlobalForce;
    void AssignForceAndMatrix();
};