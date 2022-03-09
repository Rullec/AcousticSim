#include "BaseCloth.h"
#include "utils/DefUtil.h"
#include "sim/gpu_utils/CudaArray.h"
#include "sim/gpu_utils/CudaELLMatrix.h"
#include "sim/gpu_utils/CudaMatrix.h"
// #include "sim/gpu_utils/CudaVector.h"
SIM_DECLARE_CLASS_AND_PTR(cBaraffMaterial);
SIM_DECLARE_CLASS_AND_PTR(cQBendingMaterial);
typedef tCudaMatrix<float, 9, 2> tCudaMatrix92f;
// typedef tCudaMatrix<float, 9, 9> tCudaMatrix9f;

class cBaraffClothGPU : public cBaseCloth
{
public:
    using tConnectedTriangle = tCudaVector12i;
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
    cCudaArray<tConnectedTriangle>
        mVerticesTriangleIdLstCuda; // the first number is triangle number [1 -
                                    // 11], the latter 11 number are triangle
                                    // id, -1 is invalid
    // ------ GPU time-dependent data ------
    cCudaArray<tCudaVector3f> mXcurCuda, mXpreCuda;
    cCudaArray<tCudaMatrix32f> mFLstCuda,
        mFprimeLstCuda; // F = X * S * D_m^{-1}, N = S * D_m^{-1}
    cCudaArray<tCudaMatrix32f> mnLstCuda,
        mnprimeLstCuda; // Fi / |Fi|, F colwise normalized
    cCudaArray<tCudaVector2f> mCLstCuda, mCprimeLstCuda;   // condition
    cCudaArray<tCudaMatrix92f> mgLstCuda, mgprimeLstCuda; // gi = Ni \otimes ni
    cCudaArray<tCudaMatrix9f> mEleStiffnessMatrixLstCuda;

    cBaraffMaterialPtr mStretchMaterialCPU;
    cQBendingMaterialPtr mQBendingMaterialCPU;
    // tVectorXf mXfpre
    cCudaELLMatrix mGlobalStiffnessMatrix_cpu_buffer;
    cCudaArray<cCudaELLMatrix> mGlobalStiffnessMatrix_cuda;
    // tVectorXf mXfcur, mXfpre, mMassMatrixDiagf;

    virtual void InitGeometry(const Json::Value &conf) override final;
    virtual void InitMass(const Json::Value &conf) override final;
    virtual void InitGlobalStiffnessMatrix();
    virtual void InitCpuData();
    virtual void UpdateCpuData();
    virtual void InitGpuData();
    virtual void UpdateGpuData();
    virtual void VertifyData();
};