#pragma once
#include "sim/cloth/BaseCloth.h"

/**
 * \brief           Baraff cloth
 */

SIM_DECLARE_CLASS_AND_PTR(cBaraffMaterial);
typedef Eigen::Matrix<double, 3, 2> tMatrix32d;
class cBaraffCloth : public cBaseCloth
{
public:
    explicit cBaraffCloth(int id_);
    virtual ~cBaraffCloth();
    virtual void Init(const Json::Value &conf) override final;
    virtual void UpdatePos(double dt) override final;
    virtual void ApplyUserPerturbForceOnce(tPerturb *) override;
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
    tVectorXd mMassMatDiag;
    tSparseMat mStiffnessMatrix;
    double mRayleightA, mRayleightB;

    virtual void InitBuffer();
    virtual void InitMass(const Json::Value &conf) override;
    virtual void InitMaterialCoords();
    virtual int GetSingleElementFreedom() const;
    virtual void UpdateCollisionForce(tVectorXd & col_force);
    virtual void CalcIntForce(const tVectorXd &xcur,
                              tVectorXd &int_force) const override final;
    // virtual void CalculateF(); // calculate deformation gradient
    virtual void CalcStiffnessMatrix(const tVectorXd &xcur,
                                     tSparseMat &K) const;
    virtual void SolveForNextPos(double dt);
};