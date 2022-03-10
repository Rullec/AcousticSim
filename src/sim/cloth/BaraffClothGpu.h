#include "BaseCloth.h"
#include "sim/gpu_utils/Cuda2DArray.h"
#include "sim/gpu_utils/CudaArray.h"
#include "sim/gpu_utils/CudaMatrix.h"
#include "utils/DefUtil.h"
#include "utils/SparseUtil.h"
// #include "sim/gpu_utils/CudaVector.h"
SIM_DECLARE_CLASS_AND_PTR(cBaraffMaterial);
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

protected:
    tVector3f mBendingK, mStretchK;
    // -----------CPU data-----------
    // ------ CPU constant data -----
    std::vector<tCudaMatrix32f> mNArrayCpu, mNprimeArrayCpu;
    std::vector<int> mVerticesMassCpu;
    // ------ CPU time-dependent data -----
    std::vector<tCudaVector3f> mXcurCpu, mXpreCpu;

    // -----------GPU data-----------

    // ------ GPU constant data -----
    cCudaArray<tCudaMatrix32f> mNLstCuda, mNprimeLstCuda;
    cCudaArray<int> mVerticesMassCuda;
    cCudaArray<tCudaVector3i> mTriangleVerticesIdLstCuda;
    cCudaArray<tConnectedInfo>
        mVerticesTriangleIdLstCuda; // vertex-connected triangle id, -1 is
                                    // invalid
    // all edges connected with this vertex 
    cCudaArray<tConnectedInfo> mVertexConnectedAllEdge, mVertexInvolvedInnerEdge; // only inner edges connected with this vertex

    cCudaArray<tCudaMatrix12f> mQBendingElementKLst;
    cCudaArray<tCudaVector4i> mQBendingConstraintVertexLst;
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
    cCudaArray<tCudaMatrix9f> mEleStiffnessMatrixLstCuda;
    cCudaArray<tCudaVector3f>
        mEleBendingStiffnessMatrixItemLstCuda; // row, col, value

    cBaraffMaterialPtr mStretchMaterialCPU;
    cQBendingMaterialPtr mQBendingMaterialCPU;
    cCuda2DArray<tCudaMatrix3f> mGlobalStiffnessMatrix;

    virtual void InitGeometry(const Json::Value &conf) override final;
    virtual void InitMass(const Json::Value &conf) override final;
    virtual void InitGlobalStiffnessMatrix();
    virtual void InitCpuData();
    virtual void UpdateCpuData();
    virtual void InitGpuData();
    virtual void InitQBendingHessian();
    virtual void UpdateGpuData();
    virtual void VerifyData();
    virtual void BuildEdgeInfo();
    tMatrixXf GetCudaGlobalStiffnessMatrix();
    // virtual void BuildVertexConnectedEdgeId();
};