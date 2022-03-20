#pragma once
#include "sim/cloth/BaseCloth.h"
#include "sim/gpu_utils/Cuda2DArray.h"
#include "sim/gpu_utils/CudaArray.h"
struct tFixPoint
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    int mVertexId;
    tVector3f mTargetPos;
};

SIM_DECLARE_CLASS_AND_PTR(tFixPoint);
SIM_DECLARE_CLASS_AND_PTR(cQBendingMaterial);
SIM_DECLARE_CLASS_AND_PTR(cPCGSolver);
class cBaraffClothGpu : public cBaseCloth
{
public:
    using tELLLocal2GlobalInfo = tCudaVector32i;

    explicit cBaraffClothGpu(int id_);
    cBaraffClothGpu(const cBaraffClothGpu &) = delete;
    virtual void Init(const Json::Value &conf) override final;
    virtual void UpdatePos(double dt) override final;
    virtual void ApplyUserPerturbForceOnce(tPerturb *) override;
    virtual void UpdateImGui() override;
    virtual void Reset() override;

protected:
    float mRayleighDampingA, mRayleighDampingB;
    tCudaVector3f mStretchK, mBendingK;
    cQBendingMaterialPtr mBendingMaterial;
    // cPCGSolverPtr mSolver;
    std::vector<tFixPointPtr> mConstrainedPt;
    tFixPointPtr mDragPt;

    cCudaArray<tCudaVector3f> mXcurCuda, mXpreCuda;  // node pos on Cuda
    cCudaArray<float> mTriangleInitAreaCuda;         // triangle area
    cCudaArray<tCudaVector3i> mTriangleVertexIdCuda; // triangle vertex id
    cCudaArray<float> mVerticesMassCuda;             // vertex mass
    cCudaArray<int> mConstraintVertexId;

    // =================stretch info===================
    cCudaArray<tCudaVector3f> mCoefFu_warp_weft, mCoefFv_warp_weft;
    cCudaArray<tCudaVector3f> mCoefFu_diag, mCoefFv_diag;

    cCudaArray<tELLLocal2GlobalInfo>
        mELLVidToGlobalVid; // ELL index to global vid
    cCudaArray<tCudaVector3f> mIntForceCuda, mGravityForceCuda;
    cCuda2DArray<tCudaMatrix3f> mStiffnessMatrixCuda;

    cCuda2DArray<tCudaMatrix3f> mSystemMatrixCuda;
    cCudaArray<tCudaVector3f> mSystemRHSCuda;
    cCudaArray<tCudaVector3f> mSolutionCuda;

    virtual void UpdatePosToCuda(const tVectorXd &cpu_pos,
                                 cCudaArray<tCudaVector3f> &cuda_pos) const;
    virtual void UpdatePosToCpu(const cCudaArray<tCudaVector3f> &cuda_pos,
                                tVectorXd &cpu_pos) const;
    virtual void UpdateConstraintVid(cCudaArray<int> &constrained_vid) const;
    virtual void UpdateConstraintVertexPos();
    // virtual void UpdateStiffnessMatrixAndFint();
    // virtual void UpdateLinearSystem();
    // virtual void SolveLinearSystem();
    virtual void InitGeometry(
        const Json::Value &conf); // discretazation from square cloth to
    virtual void InitMass(const Json::Value &conf);
    virtual void InitConstraint(const Json::Value &root);
    virtual void InitELLMatrix();
    virtual void InitStretchCoef();
    virtual void InitGravity();
};