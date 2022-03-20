#include "BaseCloth.h"
#include "sim/gpu_utils/Cuda2DArray.h"
#include "sim/gpu_utils/CudaArray.h"
#include "sim/gpu_utils/CudaMatrix.h"
#include "utils/DefUtil.h"
#include "utils/SparseUtil.h"
// #include "sim/gpu_utils/CudaVector.h"
SIM_DECLARE_CLASS_AND_PTR(cBaraffMaterialUnstable);
SIM_DECLARE_CLASS_AND_PTR(cQBendingMaterial);
typedef tCudaMatrix<float, 9, 2> tCudaMatrix92f;
// typedef tCudaMatrix<float, 9, 9> tCudaMatrix9f;

class cBaraffClothGPU : public cBaseCloth
{
public:
    using tConnectedInfo = tCudaVector32i;
    using tELLLocal2GlobalInfo = tCudaVector32i;
    explicit cBaraffClothGPU(int id_);
    virtual ~cBaraffClothGPU();
    virtual void Init(const Json::Value &conf) override final;
    virtual void UpdatePos(double dt) override final;
    virtual void ApplyUserPerturbForceOnce(tPerturb *) override;
    virtual void UpdateImGui() override;
    virtual void SetPos(const tVectorXd &newpos) override;
    virtual void Reset() override;
    virtual void ClearForce() override;

protected:
    tVector3f mBendingK, mStretchK;
    // -----------CPU data-----------
    // ------ CPU constant data -----
    cBaraffMaterialUnstablePtr mStretchMaterialCPU;
    cQBendingMaterialPtr mQBendingMaterialCPU;
    std::vector<tCudaMatrix32f> mNArrayCpu, mNprimeArrayCpu;
    std::vector<int> mVerticesMassCpu;

    // ------ CPU time-dependent data -----
    std::vector<tCudaVector3f> mXcurCpu, mXpreCpu;
    std::vector<tCudaVector3f> mUserForceCpu;
    // -----------GPU data-----------

    // ------ GPU constant data -----
    cCudaArray<tCudaMatrix32f> mNLstCuda, mNprimeLstCuda;
    cCudaArray<int> mVerticesMassCuda;
    cCudaArray<tCudaVector3i> mTriangleVerticesIdLstCuda;
    cCudaArray<tConnectedInfo>
        mVerticesTriangleIdLstCuda; // vertex-connected triangle id, -1 is
                                    // invalid
    cCudaArray<float> mTriangleAreaCuda;
    // all edges connected with this vertex
    cCudaArray<tConnectedInfo> mVertexConnectedAllEdge,
        mVertexInvolvedInnerEdge; // only inner edges connected with this vertex

    cCudaArray<tCudaMatrix12f> mQBendingElementKLst;
    cCudaArray<tCudaVector4i> mQBendingConstraintVertexLst;
    cCudaArray<tCudaVector3f> mMassMatrixDiagCuda;
    /*
        edge info:
            v0, v1, triangle0, v0_localid, v1_localid,
                    triangle1, v0_localid, v1_localid
    */
    cCudaArray<tCudaVector8i> mEdgeInfo;
    cCudaArray<tELLLocal2GlobalInfo> mELLLocalVidToGlobalVid;
    // ------ GPU time-dependent data ------
    cCudaArray<tCudaVector3f> mXcurCuda, mXpreCuda;
    cCudaArray<tCudaMatrix32f> mFLstCuda,
        mFprimeLstCuda; // F = X * S * D_m^{-1}, N = S * D_m^{-1}
    cCudaArray<tCudaMatrix32f> mnLstCuda,
        mnprimeLstCuda; // Fi / |Fi|, F colwise normalized
    cCudaArray<tCudaVector2f> mCLstCuda, mCprimeLstCuda;  // condition
    cCudaArray<tCudaMatrix92f> mgLstCuda, mgprimeLstCuda; // gi = Ni \otimes ni
    cCudaArray<tCudaMatrix9f> mEleStretchStiffnessMatrixLstCuda;
    cCudaArray<tCudaVector9f> mEleStretchInternalForceCuda;

    cCuda2DArray<tCudaMatrix3f> mGlobalStiffnessMatrix_stretch,
        mGlobalStiffnessMatrix_bending, mGlobalStiffnessMatrix;
    cCudaArray<tCudaVector3f> mGravityCuda, mUserForceCuda, mCollisionForce;
    cCudaArray<tCudaVector3f> mIntForceCuda_stretch, mIntForceCuda_bending,
        mIntForceCuda;
    cCudaArray<tCudaVector3f> mVelCuda;
    cCuda2DArray<tCudaMatrix3f> mSystemMatrixCuda;
    cCudaArray<tCudaVector3f> mSystemRHSCuda;
    cCudaArray<tCudaMatrix3f> mPreconditionerCuda;
    cCudaArray<tCudaVector3f> mSolutionCuda;

    cCudaArray<tCudaVector3f> mPCGResidualCuda;  // r = b - A x
    cCudaArray<tCudaVector3f> mPCGDirectionCuda; // conjugate di
    cCudaArray<tCudaVector3f> mPCGzCuda;         // z = A * di
    cCudaArray<float> mPCG_rMinvr_arrayCuda; // rMinvr, Minv is a preconditioner
    cCudaArray<float> mPCG_dTAd_arrayCuda;   // dT * A * d, A is the system mat

    // 1. CPU fixed point, inited at first
    std::vector<int> mFixedVertexIndicesCPU;
    std::vector<tCudaVector3f> mFixedVertexTargetPosCPU;

    // 2. drag point
    int mDragPointVertexId;            // default -1
    tCudaVector3f mDragPointTargetPos; // drag point pos

    // 2. Update CUDA per frame
    cCudaArray<int> mFixedVertexIndicesCUDA; // fix point indices
    cCudaArray<tCudaVector3f>
        mFixedVertexTargetPosCUDA; // fix point target position
    float mRayleighA, mRayleighB;  // rayleigh damping

    virtual void InitGeometry(const Json::Value &conf) override final;
    virtual void InitMass(const Json::Value &conf) override final;
    virtual void InitGlobalStiffnessMatrix();
    virtual void InitCpuData();
    virtual void UpdateCpuData();
    virtual void InitGpuData();
    virtual void InitQBendingHessian();
    virtual void UpdateGpuData();
    virtual void VerifyData();
    virtual void VerifyLinearSystem(float dt);
    virtual void BuildEdgeInfo();
    tMatrixXf
    FromELLMatrixToEigenMatrix(const cCuda2DArray<tCudaMatrix3f> &ell_mat);
    tVectorXf
    FromCudaVectorToEigenVector(const cCudaArray<tCudaVector3f> &cuda_vec);
    virtual void UpdateStiffnessMatrixAndIntForce();
    virtual void UpdateCollisionForceAndUserForce();
    virtual void UpdateLinearSystem(float dt);
    virtual void SolveLinearSystem();
    virtual void PostSolve();

    virtual void ClearDragPt();
    virtual void UpdateFixedPt();
    // virtual void BuildVertexConnectedEdgeId();
};
