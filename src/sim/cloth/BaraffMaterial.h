#pragma once
#include "utils/DefUtil.h"
#include "utils/MathUtil.h"
#include "utils/SparseUtil.h"

SIM_DECLARE_CLASS_AND_PTR(cBaseObject);
class cBaraffMaterial
{
public:
    cBaraffMaterial();
    virtual void Init(cBaseObjectPtr object, const tVector3d &Kwarpweftbias);
    virtual void SetK(const tVector3d &Kwarpweftbias);
    virtual void Update(); // update

    // virtual double GetEnergy(int tri_id) const;
    // virtual tVector9d GetForce(int tri_id) const;
    // virtual tMatrix9d GetStiffMatrix(int tri_id) const;
    virtual double CalcTotalEnergy() const;
    virtual tVectorXd CalcTotalForce() const;
    virtual tSparseMatd CalcTotalStiffnessMatrix();
    // virtual void CheckForce();
    // virtual void CheckStiffnessMatrix();

protected:
    cBaseObjectPtr mObject;
    tVector3d mKwarpweftbias;
    tEigenArr<tVector3d> mCoefFu_warp_weft, mCoefFv_warp_weft;
    tEigenArr<tVector3d> mCoefFu_diag, mCoefFv_diag;
    
    std::vector<double> mTriangleAreaArray;
    int mNumOfVertices, mNumOfTriangles;
    void UpdateEnergy() const;

    tSparseMatd mHessian;
    std::vector<double> mTriEnergy;
    tEigenArr<tMatrix3d> mTriFint;
    tEigenArr<tMatrix3d> mTriHessian;
};