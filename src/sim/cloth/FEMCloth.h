#pragma once
#include "sim/cloth/BaseCloth.h"

/**
 * \brief           FEM based cloth
 *
 * constitution model (would like to be supported)
 * 1. neo-hookean
 * 2. stvk
 * 3. corotated
 * 4. linear elasticity (only work for small deformation)
 */

typedef Eigen::Matrix<double, 3, 2> tMatrix32d;
class cFEMCloth : public cBaseCloth
{
public:
    explicit cFEMCloth(int id_);
    virtual ~cFEMCloth();
    virtual void Init(const Json::Value &conf) override final;
    virtual void UpdatePos(double dt) override final;

protected:
    tMatrixXd mVertexMateralCoords; // cloth's material coordinates
    tEigenArr<tMatrix32d> mF;       // deformation gradient
    tEigenArr<tMatrix2d>
        mDInv; // the inverse of [B-A, C-A] in each triangle. A, B, C are
               // material coordinate for each vertices

    tVectorXd mJ; // \det(F), the determinant of deformation gradient
    tEigenArr<tMatrixXd>
        mPK1; // first piola kirchhoff tensor, means the gradient of strain
              // energy w.r.t the nodal position. the definition of PK1 is
              // dependent on the energy definition (consitutional model)
    tEigenArr<tEigenArr<tMatrixXd>>
        mdFdx; // the gradient of deformation gradient
    virtual void InitBuffer();
    virtual void InitMaterialCoords();
    virtual int GetSingleElementFreedom() const;
    virtual void CalcIntForce(const tVectorXd &xcur,
                              tVectorXd &int_force) const override final;
    virtual void CalculateF(); // calculate deformation gradient
};