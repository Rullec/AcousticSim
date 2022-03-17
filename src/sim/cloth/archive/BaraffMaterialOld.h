#pragma once
#include "utils/MathUtil.h"
typedef Eigen::Matrix<double, 3, 2> tMatrix32d;
typedef Eigen::Matrix<double, 9, 1> tVector9d;
typedef Eigen::Matrix<double, 9, 9> tMatrix9d;

class cBaraffMaterialUnstable
{
public:
    explicit cBaraffMaterialUnstable();
    // stretch
    virtual void SetStretchK(double Kwarp, double Kweft);
    virtual void SetSheaingK(double Kwarp);
    virtual double CalcEnergy(const tMatrix3d &pos, const tMatrix32d &uv_coords) const;
    virtual tVector9d CalcForce(const tMatrix3d &pos, const tMatrix32d &uv_coords) const;
    virtual tMatrix9d CalcStiffMatrix(const tMatrix3d &pos, const tMatrix32d &rest_texture_coords) const;
    virtual void CheckForce();
    virtual void CheckStiffnessMatrix();

protected:
    // virtual tVector2d CalcCondition(const tMatrix32d &F) const;
    virtual void CalcFAndN(const tMatrix3d &pos, const tMatrix32d &uv_coords, tMatrix32d &F, tMatrix32d &N) const;
    virtual tVector2f CalcC(tMatrix32d &F) const;
    virtual tMatrixXd Calcg(const tMatrix32d &N, const tMatrix32d &n) const;
    virtual tMatrix32d Calcn(const tMatrix3d &pos, tMatrix32d &N) const;
    virtual void CalcPi(const tMatrix32d &n, tMatrix3d &P0, tMatrix3d &P1) const;

    virtual void CalcFAndN_shearing(const tMatrix3d &pos, const tMatrix32d &uv_coords, tMatrix32d &F, tMatrix32d &N) const;
    tMatrix32d mS; // selection matrix from x to \delta x
    double mBu = 1.0, mBv = 1.0, mBshear =  1.0;
    tVector3d mKwarpweftshear; // warp and weft stiffness

    virtual void CheckStretchForce();
    virtual void CheckStretchStiffnessMatrix();
    virtual double CalcStretchEnergy(const tMatrix3d &pos, const tMatrix32d &uv_coords) const;
    virtual tVector9d CalcStretchForce(const tMatrix3d &pos, const tMatrix32d &uv_coords) const;
    virtual tMatrix9d CalcStretchStiffMatrix(const tMatrix3d &pos, const tMatrix32d &rest_texture_coords) const;

    
    // shearing
    virtual void CheckShearingForce();
    virtual void CheckShearingStiffnessMatrix();
    virtual double CalcShearingEnergy(const tMatrix3d &pos, const tMatrix32d &uv_coords) const;
    virtual tVector9d CalcShearingForce(const tMatrix3d &pos, const tMatrix32d &uv_coords) const;
    virtual tMatrix9d CalcShearingStiffMatrix(const tMatrix3d &pos, const tMatrix32d &rest_texture_coords) const;
};