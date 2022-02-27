#pragma once
#include "utils/MathUtil.h"
typedef Eigen::Matrix<double, 3, 2> tMatrix32d;
typedef Eigen::Matrix<double, 9, 1> tVector9d;
typedef Eigen::Matrix<double, 9, 9> tMatrix9d;

class cBaraffMaterial
{
public:
    explicit cBaraffMaterial();
    virtual void SetK(double Kwarp, double Kweft);
    virtual void CheckForce();
    virtual void CheckStiffnessMatrix();
    virtual double CalcEnergy(const tMatrix3d &pos, const tMatrix32d &uv_coords) const;
    virtual tVector9d CalcForce(const tMatrix3d &pos, const tMatrix32d &uv_coords) const;
    virtual tMatrix9d CalcStiffMatrix(const tMatrix3d &pos, const tMatrix32d &rest_texture_coords) const;

protected:
    // virtual tVector2d CalcCondition(const tMatrix32d &F) const;
    virtual void CalcFAndN(const tMatrix3d &pos, const tMatrix32d &uv_coords, tMatrix32d &F, tMatrix32d &N) const;
    virtual tVector2f CalcC(tMatrix32d &F) const;
    virtual tMatrixXd Calcg(const tMatrix32d &N, const tMatrix32d &n) const;
    virtual tMatrix32d Calcn(const tMatrix3d &pos, tMatrix32d &N) const;
    virtual void CalcPi(const tMatrix32d &n, tMatrix3d &P0, tMatrix3d &P1) const;
    tMatrix32d mS; // selection matrix from x to \delta x
    double mBu = 1.0, mBv = 1.0;
    tVector2d mK; // warp and weft stiffness
};