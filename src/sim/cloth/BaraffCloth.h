#pragma once
#include "sim/cloth/BaseCloth.h"
#include "utils/SparseUtil.h"
/**
 * \brief           Baraff cloth
 */

SIM_DECLARE_CLASS_AND_PTR(cBaraffMaterial);
SIM_DECLARE_CLASS_AND_PTR(cBaseBendingMaterial);
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
    cBaraffMaterialPtr mMaterial;
    tSparseMatd mStiffnessMatrix;
    float mRayleightA, mRayleightB;
    tVector3f mStretchK, mBendingK; // warp weft bias
    cBaseBendingMaterialPtr mBendingMaterial;
    float mRayleighA, mRayleighB; // rayleigh damping
    int mDragVertexIdx;
    tVectorXd mMassMatrixDiag;
    tVectorXd mDx_buf;
    int mMaxCGIters;

    virtual void InitMass(const Json::Value &conf) override;
    virtual void InitMaterialCoords();
    virtual int GetSingleElementFreedom() const;
    virtual void UpdateCollisionForce(tVectorXd &col_force);
    virtual void CalcIntForce(const tVectorXd &xcur,
                              tVectorXd &int_force) const override final;
    virtual void CalcStiffnessMatrix(const tVectorXd &xcur,
                                     tSparseMatd &K) const;
    virtual void SolveForDx(double dt);
    virtual void CheckForce();
    virtual double CalcEnergy(const tVectorXd &xcur);
    virtual void Solve(const tSparseMatd &A, const tVectorXd &b, tVectorXd &x,
                       float threshold, int &iters, float &residual,
                       const std::vector<int> &fix_vertex_array);
    void Repulsion(double dt, tVectorXd &dx) const;
    tSparseMatd PrepareCollisionHessian();

};