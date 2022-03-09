#include "BaraffClothGPU.h"
#include "geometries/Primitives.h"
#include "sim/cloth/BaraffMaterial.h"
#include "sim/cloth/QBendingMaterial.h"
#include "sim/gpu_utils/CudaDef.h"
#include "sim/gpu_utils/CudaMatrix.h"
#include "sim/gpu_utils/CudaMatrixUtil.h"
#include "utils/JsonUtil.h"
#include "utils/RotUtil.h"
#include <iostream>
// strcut Params
// {
//     // std::vector<tMatrix2d>
//     devPtr<CudaMatrix2>  pDsinv;

// }

// __constant__ Params parm;

// __global__ void CalcF(pointer1, pointer2, pointer3, ..)
// {
//     // jisuan
// }

cBaraffClothGPU::cBaraffClothGPU(int id_)
    : cBaseCloth(eClothType::FEM_CLOTH_GPU, id_)
{
}

namespace nBaraffClothGpu
{
extern void Calc(int);
};

void cBaraffClothGPU::Init(const Json::Value &conf)
{
    /*
        1. load parameter:
            bending + stretch stiffness
    */
    mBendingK.noalias() =
        cJsonUtil::ReadVectorJson(
            cJsonUtil::ParseAsValue("cloth_bending_stiffness", conf))
            .cast<float>();
    mStretchK.noalias() =
        cJsonUtil::ReadVectorJson(
            cJsonUtil::ParseAsValue("cloth_stretch_stiffness", conf))
            .segment(0, 3)
            .cast<float>();
    /*
        2. init base class; N; mass
    */
    cBaseCloth::Init(conf);

    // 2.1 init CPU material (for comparision)
    {
        mStretchMaterialCPU = std::make_shared<cBaraffMaterial>();
        mStretchMaterialCPU->Init(shared_from_this(), mStretchK.cast<double>());
        mQBendingMaterialCPU = std::make_shared<cQBendingMaterial>();
        mQBendingMaterialCPU->Init(GetVertexArray(), GetEdgeArray(),
                                   GetTriangleArray(),
                                   mBendingK.cast<double>());
    }

    InitGlobalStiffnessMatrix();
    /*
        3. allocate: CPU buffer and GPU buffer
    */
    InitCpuData();
    UpdateCpuData();

    InitGpuData();
    UpdateGpuData();
    // 5. calculate energy, calculate internal force, calculate hessian
    // check the numerical diff
    // UpdateDataGPU();
    VertifyData();
    // int dims = 5;
    // nBaraffClothGpu::Calc(dims);
    // 1. Ds^{-1}
    // CudaArray<tMatrix2f> mem;
    // mem.upload();

    // param.upload();
    // KERNEL_NAME CUDA_AT(mem.ptr(), )
    exit(1);
}

cBaraffClothGPU::~cBaraffClothGPU() {}
void cBaraffClothGPU::UpdatePos(double dt) {}
void cBaraffClothGPU::ApplyUserPerturbForceOnce(tPerturb *) {}
void cBaraffClothGPU::UpdateImGui() {}

void cBaraffClothGPU::InitMass(const Json::Value &conf)
{
    mClothDensity = cJsonUtil::ParseAsDouble("cloth_density", conf);
    int dof = 3 * GetNumOfVertices();
    mMassMatrixDiag.noalias() = tVectorXd::Zero(dof);
    for (auto &t : mTriangleArrayShared)
    {
        // 1. total area
        auto v0 = mVertexArrayShared[t->mId0];
        auto v1 = mVertexArrayShared[t->mId1];
        auto v2 = mVertexArrayShared[t->mId2];

        double triangle_area =
            cMathUtil::CalcTriangleArea(v0->mPos, v1->mPos, v2->mPos);
        mMassMatrixDiag.segment(3 * t->mId0, 3) +=
            triangle_area / 3 * mClothDensity * tVector3d::Ones();
        mMassMatrixDiag.segment(3 * t->mId1, 3) +=
            triangle_area / 3 * mClothDensity * tVector3d::Ones();
        mMassMatrixDiag.segment(3 * t->mId2, 3) +=
            triangle_area / 3 * mClothDensity * tVector3d::Ones();
    }
    for (int i = 0; i < GetNumOfVertices(); i++)
    {
        mVertexArrayShared[i]->mMass = mMassMatrixDiag[3 * i];
        // std::cout << "v" << i << " mass = " << mVertexArrayShared[i]->mMass
        // << std::endl;
    }
    // mMassMatrixDiagf.noalias() = mMassMatrixDiag.cast<float>();
}

void cBaraffClothGPU::InitGpuData()
{
    mNLstCuda.Upload(mNArrayCpu);
    mNprimeLstCuda.Upload(mNprimeArrayCpu);
    mVerticesMassCuda.Upload(mVerticesMassCpu);

    // triangle's vertices id
    int num_of_tris = GetNumOfTriangles();
    {
        std::vector<tCudaVector3i> triangle_vid_lst(num_of_tris);
        for (size_t i = 0; i < num_of_tris; i++)
        {
            triangle_vid_lst[i] = tCudaVector3i(
                {mTriangleArrayShared[i]->mId0, mTriangleArrayShared[i]->mId1,
                 mTriangleArrayShared[i]->mId2});
        }
        mTriangleVerticesIdLstCuda.Upload(triangle_vid_lst);
    }

    // vertex's included triangle id
    int num_of_v = GetNumOfVertices();
    {
        auto Add = [](tConnectedTriangle &info, int tri_id)
        {
            for (int k = 0; k < tConnectedTriangle::mElements; k++)
            {
                if (info[k] == -1)
                {
                    info[k] = tri_id;
                    return;
                }
            }
            SIM_ERROR("triangle {} vertices is included in more than {} "
                      "triangles, error!",
                      tri_id, tConnectedTriangle::mElements);
        };
        std::vector<tConnectedTriangle> vertex_triangle_id_lst(num_of_v, -1);
        for (size_t tri_id = 0; tri_id < num_of_tris; tri_id++)
        {
            Add(vertex_triangle_id_lst[mTriangleArrayShared[tri_id]->mId0],
                tri_id);
            Add(vertex_triangle_id_lst[mTriangleArrayShared[tri_id]->mId1],
                tri_id);
            Add(vertex_triangle_id_lst[mTriangleArrayShared[tri_id]->mId2],
                tri_id);
        }
        mVerticesTriangleIdLstCuda.Upload(vertex_triangle_id_lst);
    }

    mXcurCuda.Resize(num_of_v);
    mXpreCuda.Resize(num_of_v);
    mFLstCuda.Resize(num_of_tris);
    mFprimeLstCuda.Resize(
        num_of_tris); // F = X * S * D_m^{-1}, N = S * D_m^{-1}
    mnLstCuda.Resize(num_of_tris);
    mnprimeLstCuda.Resize(num_of_tris); // Fi / |Fi|, F colwise normalized
    mCLstCuda.Resize(num_of_tris);
    mCprimeLstCuda.Resize(num_of_tris);
    mgLstCuda.Resize(num_of_tris);
    mgprimeLstCuda.Resize(num_of_tris); // gi = Ni \otimes ni
    mEleStiffnessMatrixLstCuda.Resize(num_of_tris);
    mEleStiffnessMatrixLstCuda.MemsetAsync(tCudaMatrix9f::Zero());
}

void cBaraffClothGPU::InitCpuData()
{
    mNArrayCpu.resize(GetNumOfTriangles());
    mNprimeArrayCpu.resize(GetNumOfTriangles());
    // mNArrayCuda.clear();
    // mNprimeArrayCuda.clear();
    // 1. calculate DmInv for cpu
    {
        tMatrix32f mS;
        mS.setZero();
        mS(0, 0) = -1;
        mS(0, 1) = -1;
        mS(1, 0) = 1;
        mS(2, 1) = 1;

        auto &tri_lst = GetTriangleArray();
        auto &v_lst = GetVertexArray();
        tMatrix2f R_neg45 = cRotUtil::RotMat2D(-M_PI / 4).cast<float>();
        for (int i = 0; i < GetNumOfTriangles(); i++)
        {
            auto t = GetTriangleArray()[i];

            tMatrix2f DmInv = tMatrix2f::Zero();
            DmInv.col(0) = (v_lst[t->mId1]->muv - v_lst[t->mId0]->muv);
            DmInv.col(1) = (v_lst[t->mId2]->muv - v_lst[t->mId0]->muv);
            DmInv = DmInv.inverse().eval();
            DmInv.hasNaN();

            // 2. calculate N' = S * Dminv * R^{-45}
            tMatrix32f result = mS * DmInv;
            mNArrayCpu[i] = cCudaMatrixUtil::EigenMatrixToCudaMatrix(result);
            mNprimeArrayCpu[i] = cCudaMatrixUtil::EigenMatrixToCudaMatrix(
                tMatrix32f(result * R_neg45));
        }
    }
    int num_of_v = GetNumOfVertices();
    mVerticesMassCpu.resize(num_of_v);
    for (int i = 0; i < num_of_v; i++)
    {
        mVerticesMassCpu[i] = GetVertexArray()[i]->mMass;
    }

    // init mXcur
    mXcurCpu.resize(num_of_v);
    mXpreCpu.resize(num_of_v);
    for (int i = 0; i < num_of_v; i++)
    {
        // mXcurCpu[i] = cCudaMatrixUtil::EigenVectorToSimVector<float, 3>(
        // mXcur.segment(3 * i, 3).cast<float>());
        mXcurCpu[i] = cCudaMatrixUtil::EigenMatrixToCudaMatrix(
            tVector3f(mXcur.segment(3 * i, 3).cast<float>()));
        mXpreCpu[i] = mXcurCpu[i];
    }
}

/**
 * \brief
 */
namespace BaraffClothGpu
{

// extern void UpdateF(const cCudaArray<float> &Xcur, const cCudaArray<float>
// &Xpre, const cCudaArray<tMatrix32f> &Nlst, cCudaArray<tMatrix32f> &FLst);
extern void UpdateK(float K0, float K1, cCudaArray<tCudaVector3f> &pos_lst,
                    cCudaArray<tCudaVector3i> &tri_vertices_lst,
                    cCudaArray<tCudaMatrix32f> &N_lst,
                    cCudaArray<tCudaMatrix32f> &F_lst,
                    cCudaArray<tCudaMatrix32f> &n_lst,
                    cCudaArray<tCudaVector2f> &C_lst, // condition
                    cCudaArray<tCudaMatrix92f> &g_lst,
                    cCudaArray<tCudaMatrix9f> &K_lst);
} // namespace BaraffClothGpu
void cBaraffClothGPU::UpdateGpuData()
{
    // 1. update Xcur and Xpre
    {
        mXcurCuda.Upload(mXcurCpu);
        mXpreCuda.Upload(mXpreCpu);
    }

    // 2. update F
    {

        BaraffClothGpu::UpdateK(mStretchK[0], mStretchK[1], mXcurCuda,
                                mTriangleVerticesIdLstCuda, mNLstCuda,
                                mFLstCuda, mnLstCuda, mCLstCuda, mgLstCuda,
                                mEleStiffnessMatrixLstCuda);
        BaraffClothGpu::UpdateK(mStretchK[2], mStretchK[2], mXcurCuda,
                                mTriangleVerticesIdLstCuda, mNprimeLstCuda,
                                mFprimeLstCuda, mnprimeLstCuda, mCprimeLstCuda,
                                mgprimeLstCuda, mEleStiffnessMatrixLstCuda);
    }
}

template <typename T1, typename T2>
bool CheckIsSyncVec(const char *name, const std::vector<T1> &array_cpu,
                    const cCudaArray<T2> &array_gpu)
{
    std::vector<T2> array_cpu_download;
    array_gpu.Download(array_cpu_download);
    std::vector<T1> array_cpu_download_eigen(array_cpu_download.size());

    for (size_t i = 0; i < array_cpu.size(); i++)
    {
        array_cpu_download_eigen[i] =
            cCudaMatrixUtil::CudaMatrixToEigenMatrix(array_cpu_download[i]);
        float diff =
            (array_cpu_download_eigen[i] - array_cpu[i]).cwiseAbs().maxCoeff();
        if (diff > 1e-6 || std::isnan(diff))
        {
            SIM_ERROR("{}: Check sync vec[{}] failed, diff = {}", name, i,
                      diff);
            std::cout << "gpu = " << array_cpu_download_eigen[i] << std::endl;
            std::cout << "cpu = " << array_cpu[i] << std::endl;
            return false;
        }
    }

    SIM_INFO("{}: Check synced", name);
}
template <typename T1, typename T2>
bool CheckIsSyncMat(const char *name, const std::vector<T1> &array_cpu,
                    const cCudaArray<T2> &array_gpu)
{
    std::vector<T2> array_cpu_download;
    array_gpu.Download(array_cpu_download);
    std::vector<T1> array_cpu_download_eigen(array_cpu_download.size());

    for (size_t i = 0; i < array_cpu.size(); i++)
    {
        array_cpu_download_eigen[i] =
            cCudaMatrixUtil::CudaMatrixToEigenMatrix(array_cpu_download[i]);
        float diff =
            (array_cpu_download_eigen[i] - array_cpu[i]).cwiseAbs().maxCoeff();
        if (diff > 1e-6 || std::isnan(diff))
        {
            SIM_ERROR("{}: Check sync vec[{}] failed, diff = {}", name, i,
                      diff);
            std::cout << "cpu = " << array_cpu[i] << std::endl;
            std::cout << "gpu = " << array_cpu_download_eigen[i] << std::endl;
            return false;
        }
    }
    SIM_INFO("{}: Check synced", name);
}

void cBaraffClothGPU::VertifyData()
{
    // 1. check data synchronization
    {
        // 1.1 check xcuda, xpre, massdiag synchronization
        // SIM_ASSERT(CheckIsSync(mXfcur, mXcurCuda));
        // SIM_ASSERT(CheckIsSync(mXfpre, mXpreCuda));
        // SIM_ASSERT(CheckIsSync(mMassMatrixDiagf, mMassDiagCuda));
        // SIM_ASSERT(CheckIsSync(mNArrayCpu, mNLstCuda));
        // SIM_ASSERT(CheckIsSync(mNprimeArrayCpu, mNprimeLstCuda));
    } // 2. check calculated value
    {
        mStretchMaterialCPU->Update();
        // 1. F
        {
            std::vector<tMatrix32f> mFLstCPU, mFprimeLstCPU;
            std::vector<tMatrix32f> nLstCPU;
            std::vector<tMatrix92f> gLstCPU;
            std::vector<tVector2f> CLstCPU;
            std::vector<tMatrix32f> nprimeLstCPU;
            std::vector<tMatrix92f> gprimeLstCPU;
            std::vector<tMatrix9f> KLstCPU;
            std::vector<tVector2f> CprimeLstCPU;
            mStretchMaterialCPU->GetFLst(mFLstCPU);
            mStretchMaterialCPU->GetFprimeLst(mFprimeLstCPU);
            mStretchMaterialCPU->GetnLst(nLstCPU);
            mStretchMaterialCPU->GetgLst(gLstCPU);
            mStretchMaterialCPU->GetCLst(CLstCPU);
            mStretchMaterialCPU->GetnprimeLst(nprimeLstCPU);
            mStretchMaterialCPU->GetgprimeLst(gprimeLstCPU);
            mStretchMaterialCPU->GetCprimeLst(CprimeLstCPU);
            mStretchMaterialCPU->GetEleKLst(KLstCPU);

            CheckIsSyncMat("F", mFLstCPU, mFLstCuda);
            CheckIsSyncMat("Fprime", mFprimeLstCPU, mFprimeLstCuda);

            CheckIsSyncMat("n", nLstCPU, mnLstCuda);
            CheckIsSyncMat("nprime", nprimeLstCPU, mnprimeLstCuda);

            CheckIsSyncMat("g", gLstCPU, mgLstCuda);
            CheckIsSyncMat("gprime", gprimeLstCPU, mgprimeLstCuda);

            CheckIsSyncVec("C", CLstCPU, mCLstCuda);
            CheckIsSyncVec("Cprime", CprimeLstCPU, mCprimeLstCuda);

            CheckIsSyncVec("K", KLstCPU, mEleStiffnessMatrixLstCuda);
            // for (int i = 0; i < GetNumOfTriangles(); i++)
            // {
            //     tMatrix32f cpu = mFLstCPU[i];
            //     tMatrix32f gpu = tCudaMatrixToEigenMatrix(mFLstFromCuda[i]);
            //     tMatrix32f diff = cpu - gpu;
            //     std::cout << diff.norm() << std::endl;
            //     std::cout << cpu << std::endl;
            //     std::cout << gpu << std::endl;
            // }

            // SIM_ASSERT(CheckIsSync(mFLstCPU, mFLstCuda));
            // SIM_ASSERT(CheckIsSync(mFprimeLstCPU, mFprimeLstCuda));
        }
    }
}

void cBaraffClothGPU::InitGeometry(const Json::Value &conf)
{
    cBaseCloth::InitGeometry(conf);
    // mXfcur.noalias() = mXcur.cast<float>();
    // mXfpre.noalias() = mXfcur;
}

void cBaraffClothGPU::SetPos(const tVectorXd &newpos)
{
    cBaseCloth::SetPos(newpos);
}

void cBaraffClothGPU::UpdateCpuData()
{
    for (int i = 0; i < GetNumOfVertices(); i++)
    {
        mXcurCpu[i] = cCudaMatrixUtil::EigenMatrixToCudaMatrix(
            tVector3f(mXcur.segment(3 * i, 3).cast<float>()));
        mXpreCpu[i] = cCudaMatrixUtil::EigenMatrixToCudaMatrix(
            tVector3f(mXpre.segment(3 * i, 3).cast<float>()));
    }
}

void cBaraffClothGPU::InitGlobalStiffnessMatrix()
{
    // 1. build vertice connection
    int num_of_v = GetNumOfVertices();
    std::vector<std::vector<int>> mColumnIdLstPerRow;
    for (int i = 0; i < num_of_v; i++)
        mColumnIdLstPerRow.push_back({});
    int num_of_e = GetNumOfEdges();
    for (int e_id = 0; e_id < num_of_e; e_id++)
    {
        auto cur_e = mEdgeArrayShared[e_id];
        int v0 = cur_e->mId0, v1 = cur_e->mId1;
        mColumnIdLstPerRow[v0].push_back(v1);
        mColumnIdLstPerRow[v1].push_back(v0);
    }
    for (int i = 0; i < num_of_v; i++)
    {
        std::sort(mColumnIdLstPerRow[i].begin(), mColumnIdLstPerRow[i].end());
    }
    mGlobalStiffnessMatrix_cpu_buffer.Init(num_of_v, num_of_v,
                                           mColumnIdLstPerRow);
    mGlobalStiffnessMatrix_cuda.Resize(1);
    mGlobalStiffnessMatrix_cpu_buffer.SetZero();
    mGlobalStiffnessMatrix_cuda.Upload({mGlobalStiffnessMatrix_cpu_buffer});
}