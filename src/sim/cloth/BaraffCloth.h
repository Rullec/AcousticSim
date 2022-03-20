#pragma once
#include "sim/cloth/BaseCloth.h"
#include "utils/SparseUtil.h"
/**
 * \brief           Baraff cloth
 */

SIM_DECLARE_CLASS_AND_PTR(cBaraffMaterial);
SIM_DECLARE_CLASS_AND_PTR(cQBendingMaterial);
typedef Eigen::Matrix<double, 3, 2> tMatrix32d;
class cBaraffCloth : public cBaseCloth
{
public:
    explicit cBaraffCloth(int id_);
    virtual ~cBaraffCloth();
    virtual void Init(const Json::Value &conf) override final;
    virtual void UpdatePos(double dt) override final;
    virtual void ApplyUserPerturbForceOnce(tPerturb *) override;
    virtual void UpdateImGui() override;

protected:
    tMatrixXd mVertexMateralCoords; // cloth's material coordinates
    // tEigenArr<tMatrix32d> mF;       // deformation gradient
    // tEigenArr<tMatrix2d>
    //     mDInv; // the inverse of [B-A, C-A] in each triangle. A, B, C are
    //            // material coordinate for each vertices

    // tVectorXd mJ; // \det(F), the determinant of deformation gradient
    // tEigenArr<tMatrixXd>
    //     mPK1; // first piola kirchhoff tensor, means the gradient of strain
    //           // energy w.r.t the nodal position. the definition of PK1 is
    //           // dependent on the energy definition (consitutional model)
    // tEigenArr<tEigenArr<tMatrixXd>>
    //     mdFdx; // the gradient of deformation gradient
    cBaraffMaterialPtr mMaterial;
    tSparseMatd mStiffnessMatrix;
    float mRayleightA, mRayleightB;
    tVector3f mStretchK, mBendingK; // warp weft bias
    cQBendingMaterialPtr mBendingMaterial;
    float mRayleighA, mRayleighB; // rayleigh damping
    int mDragVertexIdx;
    tVectorXd mMassMatrixDiag;
    // virtual void InitBuffer();
    virtual void InitMass(const Json::Value &conf) override;
    virtual void InitMaterialCoords();
    virtual int GetSingleElementFreedom() const;
    virtual void UpdateCollisionForce(tVectorXd &col_force);
    virtual void CalcIntForce(const tVectorXd &xcur,
                              tVectorXd &int_force) const override final;
    // virtual void CalculateF(); // calculate deformation gradient
    virtual void CalcStiffnessMatrix(const tVectorXd &xcur,
                                     tSparseMatd &K) const;
    virtual void SolveForNextPos(double dt);
    virtual void CheckForce();
    virtual double CalcEnergy(const tVectorXd &xcur);
    virtual void Solve(const tSparseMatd &A, const tVectorXd &b, tVectorXd &x,
                       float threshold, int &iters, float &residual, int fix_vertex = -1);
    // virtual void CheckStiffnessMatrix();
};