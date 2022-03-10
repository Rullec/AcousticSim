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
    VerifyData();
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

void SortInfo(cBaraffClothGPU::tConnectedInfo &info)
{
    std::vector<int> val = {};
    for (int i = 0; i < info.size(); i++)
    {
        if (info[i] != -1)
        {
            val.push_back(info[i]);
        }
    }
    std::sort(val.begin(), val.end());

    for (int i = 0; i < val.size(); i++)
    {
        info[i] = val[i];
    }
}
void AddInfo(cBaraffClothGPU::tConnectedInfo &info, int add)
{
    for (int k = 0; k < info.size(); k++)
    {
        if (info[k] == -1)
        {
            info[k] = add;
            return;
        }
        else if (info[k] == add)
        {
            return;
        }
    }
    SIM_ERROR("connected info {} try to include more than {} "
              "info, error!",
              info.transpose(), info.size());
    exit(1);
};

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

        std::vector<tConnectedInfo> vertex_triangle_id_lst(num_of_v, -1);
        for (size_t tri_id = 0; tri_id < num_of_tris; tri_id++)
        {
            AddInfo(vertex_triangle_id_lst[mTriangleArrayShared[tri_id]->mId0],
                    tri_id);
            AddInfo(vertex_triangle_id_lst[mTriangleArrayShared[tri_id]->mId1],
                    tri_id);
            AddInfo(vertex_triangle_id_lst[mTriangleArrayShared[tri_id]->mId2],
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

    InitQBendingHessian();
    BuildEdgeInfo();
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
extern void
UpdateK(float K0, float K1, const cCudaArray<tCudaVector3f> &pos_lst,
        const cCudaArray<tCudaVector3i> &tri_vertices_lst,
        cCudaArray<tCudaMatrix32f> &N_lst, cCudaArray<tCudaMatrix32f> &F_lst,
        cCudaArray<tCudaMatrix32f> &n_lst,
        cCudaArray<tCudaVector2f> &C_lst, // condition
        cCudaArray<tCudaMatrix92f> &g_lst, cCudaArray<tCudaMatrix9f> &K_lst);
extern void AssembleStretchStiffnessMatrix(
    const cCudaArray<tCudaMatrix9f> &ele_K_lst_per_triangle,
    const cCudaArray<tCudaVector32i> &vertex_connected_triangle_lst,
    const cCudaArray<tCudaVector3i> &vertex_id_of_triangle_lst,
    const cCudaArray<tCudaVector32i> &ELL_local_vertex_id_to_global_vertex_id,
    cCuda2DArray<tCudaMatrix3f> &global_K);

extern void AssembleBendingStiffnessMatrix(
    const cCudaArray<tCudaVector4i> &bending_vertex_per_edge_lst,
    const cCudaArray<tCudaVector32i> &vertex_connected_edges_lst,
    const cCudaArray<tCudaMatrix12f> &ele_K_lst,
    const cCudaArray<tCudaVector32i> &ELL_local_vertex_id_to_global_vertex_id,
    cCuda2DArray<tCudaMatrix3f> &global_K);
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
        BaraffClothGpu::AssembleStretchStiffnessMatrix(
            mEleStiffnessMatrixLstCuda, mVerticesTriangleIdLstCuda,
            mTriangleVerticesIdLstCuda, mELLLocalVidToGlobalVid,
            this->mGlobalStiffnessMatrix);
        BaraffClothGpu::AssembleBendingStiffnessMatrix(
            mQBendingConstraintVertexLst, mVertexInvolvedInnerEdge,
            mQBendingElementKLst, mELLLocalVidToGlobalVid,
            mGlobalStiffnessMatrix);
        // assemble bending stiffness matrix
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

void cBaraffClothGPU::VerifyData()
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
            // tMatrixXf global_K_cpu =
            //     mStretchMaterialCPU->CalcTotalStiffnessMatrix()
            //         .cast<float>()
            //         .toDense();
            tMatrixXf global_K_cpu =
                (mQBendingMaterialCPU->GetStiffnessMatrix().cast<float>() +
                mStretchMaterialCPU->CalcTotalStiffnessMatrix().cast<float>()).toDense();
            // CheckIsSyncMat("F", mFLstCPU, mFLstCuda);
            // CheckIsSyncMat("Fprime", mFprimeLstCPU, mFprimeLstCuda);

            // CheckIsSyncMat("n", nLstCPU, mnLstCuda);
            // CheckIsSyncMat("nprime", nprimeLstCPU, mnprimeLstCuda);

            // CheckIsSyncMat("g", gLstCPU, mgLstCuda);
            // CheckIsSyncMat("gprime", gprimeLstCPU, mgprimeLstCuda);

            // CheckIsSyncVec("C", CLstCPU, mCLstCuda);
            // CheckIsSyncVec("Cprime", CprimeLstCPU, mCprimeLstCuda);

            // CheckIsSyncVec("K", KLstCPU, mEleStiffnessMatrixLstCuda);

            tMatrixXf global_K_cuda = GetCudaGlobalStiffnessMatrix();
            auto K_diff = (global_K_cuda - global_K_cpu);
            std::cout << "global K diff norm = " << K_diff.norm() << std::endl;
            std::cout << "global K diff abs maxcoef = "
                      << K_diff.cwiseAbs().maxCoeff() << std::endl;
            // std::cout << "global K cuda = \n" << global_K_cuda << std::endl;
            // std::cout << "global K cpu = \n" << global_K_cpu << std::endl;
            // check block
            int num_of_v = GetNumOfVertices();
            for (int i = 0; i < num_of_v; i++)
                for (int j = 0; j < num_of_v; j++)
                {
                    auto block_diff = K_diff.block(3 * i, 3 * j, 3, 3);
                    auto block_cuda = global_K_cuda.block(3 * i, 3 * j, 3, 3);
                    auto block_cpu = global_K_cpu.block(3 * i, 3 * j, 3, 3);
                    float max_coef = block_diff.cwiseAbs().maxCoeff();
                    if (max_coef > 1e-3)
                    {
                        std::cout << "block " << i << " " << j << " diff = \n"
                                  << block_diff << std::endl;
                        std::cout << "cuda = \n" << block_cuda << std::endl;
                        std::cout << "cpu = \n" << block_cpu << std::endl;
                    }
                }
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
#include <set>
void cBaraffClothGPU::InitGlobalStiffnessMatrix()
{
    // 1. build vertice connection
    int num_of_v = GetNumOfVertices();
    std::vector<tELLLocal2GlobalInfo> ELL_local_vid_to_global_vid(
        num_of_v, tELLLocal2GlobalInfo::Ones() * -1);

    for (int i = 0; i < num_of_v; i++)
    {
        ELL_local_vid_to_global_vid[i][0] = i;
    }
    int num_of_e = GetNumOfEdges();

    for (int e_id = 0; e_id < num_of_e; e_id++)
    {
        auto cur_e = mEdgeArrayShared[e_id];
        int v0 = cur_e->mId0, v1 = cur_e->mId1;
        // stretch: one vertex only affect vertices in this triangle
        AddInfo(ELL_local_vid_to_global_vid[v0], v1);
        AddInfo(ELL_local_vid_to_global_vid[v1], v0);
        // bending: one vertex will be affected in all vertices involed in
        // triangles
        if (cur_e->mIsBoundary == false)
        {
            // inner edge, apply inner affect
            std::set<int> cur_set = {
                mTriangleArrayShared[cur_e->mTriangleId0]->mId0,
                mTriangleArrayShared[cur_e->mTriangleId0]->mId1,
                mTriangleArrayShared[cur_e->mTriangleId0]->mId2,
                mTriangleArrayShared[cur_e->mTriangleId1]->mId0,
                mTriangleArrayShared[cur_e->mTriangleId1]->mId1,
                mTriangleArrayShared[cur_e->mTriangleId1]->mId2};
            int k = 0;
            tVector4i involed_vertices = tVector4i::Ones() * -1;
            for (auto &x : cur_set)
            {
                involed_vertices[k] = x;
                assert(k++ <= 4);
            }
            for (int a = 0; a < 4; a++)
                for (int b = 0; b < 4; b++)
                {
                    AddInfo(ELL_local_vid_to_global_vid[involed_vertices[a]],
                            involed_vertices[b]);
                }
        }
    }

    for (int i = 0; i < num_of_v; i++)
    {
        SortInfo(ELL_local_vid_to_global_vid[i]);
        // std::cout << "v" << i << " ELL local to global = "
        //           << ELL_local_vid_to_global_vid[i].transpose() << std::endl;
    }
    mELLLocalVidToGlobalVid.Upload(ELL_local_vid_to_global_vid);
    int max_connected = tELLLocal2GlobalInfo::mElements;
    // for (int i = 0; i < num_of_v; i++)
    // {
    //     std::sort(mColumnIdLstPerRow[i].begin(),
    //     mColumnIdLstPerRow[i].end()); if (mColumnIdLstPerRow[i].size() >
    //     max_connected)
    //     {
    //         printf("max connected vertices %d > %d, please increase the upper
    //         "
    //                "limit\n",
    //                mColumnIdLstPerRow[i].size(), max_connected);
    //     }
    // }
    mGlobalStiffnessMatrix.Resize(num_of_v, max_connected);
    mGlobalStiffnessMatrix.MemsetAsync(tCudaMatrix3f::Zero());
}

void cBaraffClothGPU::InitQBendingHessian()
{

    std::vector<tVector4i> edge_constraint_vertex_lst_eigen =
        mQBendingMaterialCPU->GetEdgeConstraintVertex();
    std::vector<tMatrix12f> eleK_lst_eigen =
        mQBendingMaterialCPU->GetEleStiffnessMatrixLst();
    // std::cout << "edge_constraint_vertex_lst_eigen = "
    //           << edge_constraint_vertex_lst_eigen.size() << std::endl;
    // std::cout << "eleK_lst_eigen = " << eleK_lst_eigen.size() << std::endl;

    std::vector<tCudaVector4i> edge_constraint_vertex_lst_cpu(
        edge_constraint_vertex_lst_eigen.size());
    std::vector<tCudaMatrix12f> eleK_lst_cpu(eleK_lst_eigen.size());
    for (int i = 0; i < eleK_lst_eigen.size(); i++)
    {
        eleK_lst_cpu[i] =
            -1 * cCudaMatrixUtil::EigenMatrixToCudaMatrix(eleK_lst_eigen[i]);

        edge_constraint_vertex_lst_cpu[i] =
            cCudaMatrixUtil::EigenMatrixToCudaMatrix(
                edge_constraint_vertex_lst_eigen[i]);
        // std::cout << "edge constarint " << i << " vertex lst ="
        //           << edge_constraint_vertex_lst_eigen[i].transpose()
        //           << std::endl;
        // std::cout << "ele K " << i << " = \n" << eleK_lst_cpu[i] <<
        // std::endl;
    }

    mQBendingElementKLst.Upload(eleK_lst_cpu);
    mQBendingConstraintVertexLst.Upload(edge_constraint_vertex_lst_cpu);
    // std::cout << "mQBendingElementKLst size = " <<
    // mQBendingElementKLst.Size()
    //           << std::endl;
    // std::cout << "mQBendingConstraintVertexLst size = "
    //           << mQBendingConstraintVertexLst.Size() << std::endl;
    // exit(1);
}

int GetLocalId(tTrianglePtr triangle, int v_id)
{
    if (v_id == triangle->mId0)
        return 0;
    else if (v_id == triangle->mId1)
        return 1;
    else if (v_id == triangle->mId2)
        return 2;
    else
    {
        printf("[error] no v %d in triangle\n", v_id);
        assert(false);
    }
    return -1;
}

void cBaraffClothGPU::BuildEdgeInfo()
{
    int num_of_e = GetNumOfEdges();
    std::vector<tCudaVector8i> edge_info_lst(num_of_e);
    /*
        edge info:
        vec[0]:     vertex0 id,
        vec[1]:     vertex1 id,
        vec[2]:     triangle0 id,
        vec[3]:     v0_local_id in triangle0
        vec[4]:     v1_local_id in triangle0
        vec[5]:     triangle1 id,       (-1 if boundary edge)
        vec[6]:     v0_local_id in triangle1 (-1 if boundary edge)
        vec[7]:     v1_local_id in triangle1 (-1 if boundary edge)
    */
    int num_of_v = this->GetNumOfVertices();
    std::vector<tConnectedInfo> vertex_connected_edges(
        num_of_v, -1 * tConnectedInfo::Ones());
    std::vector<tConnectedInfo> vertex_involved_inner_edges(
        num_of_v, -1 * tConnectedInfo::Ones());

    for (int all_e_id = 0, inner_e_id = 0; all_e_id < num_of_e; all_e_id++)
    {
        tCudaVector8i &cur_edge_info = edge_info_lst[all_e_id];
        cur_edge_info.setValue(-1);
        auto e = mEdgeArrayShared[all_e_id];
        cur_edge_info[0] = e->mId0;
        cur_edge_info[1] = e->mId1;
        cur_edge_info[2] = e->mTriangleId0;
        cur_edge_info[3] =
            GetLocalId(mTriangleArrayShared[e->mTriangleId0], e->mId0);
        cur_edge_info[4] =
            GetLocalId(mTriangleArrayShared[e->mTriangleId0], e->mId1);

        if (e->mTriangleId1 != -1)
        {
            cur_edge_info[5] = e->mTriangleId1;
            cur_edge_info[6] =
                GetLocalId(mTriangleArrayShared[e->mTriangleId1], e->mId0);
            cur_edge_info[7] =
                GetLocalId(mTriangleArrayShared[e->mTriangleId1], e->mId1);
        }

        if (e->mIsBoundary == false)
        {
            // AddInfo(vertex_involved_inner_edges[e->mId0], inner_e_id);
            // AddInfo(vertex_involved_inner_edges[e->mId1], inner_e_id);
            AddInfo(vertex_involved_inner_edges
                        [mTriangleArrayShared[e->mTriangleId0]->mId0],
                    inner_e_id);
            AddInfo(vertex_involved_inner_edges
                        [mTriangleArrayShared[e->mTriangleId0]->mId1],
                    inner_e_id);
            AddInfo(vertex_involved_inner_edges
                        [mTriangleArrayShared[e->mTriangleId0]->mId2],
                    inner_e_id);
            AddInfo(vertex_involved_inner_edges
                        [mTriangleArrayShared[e->mTriangleId1]->mId0],
                    inner_e_id);
            AddInfo(vertex_involved_inner_edges
                        [mTriangleArrayShared[e->mTriangleId1]->mId1],
                    inner_e_id);
            AddInfo(vertex_involved_inner_edges
                        [mTriangleArrayShared[e->mTriangleId1]->mId2],
                    inner_e_id);

            inner_e_id += 1;
        }
        AddInfo(vertex_connected_edges[e->mId0], all_e_id);
        AddInfo(vertex_connected_edges[e->mId1], all_e_id);
    }
    // for (int i = 0; i < num_of_v; i++)
    // {
    //     std::cout << "v" << i << " inner edge = "
    //               << vertex_involved_inner_edges[i].transpose() << std::endl;
    // }
    mEdgeInfo.Upload(edge_info_lst);
    mVertexConnectedAllEdge.Upload(vertex_connected_edges);
    mVertexInvolvedInnerEdge.Upload(vertex_involved_inner_edges);
}

tMatrixXf cBaraffClothGPU::GetCudaGlobalStiffnessMatrix()
{
    int num_of_v = this->GetNumOfVertices();
    int num_of_dof = GetNumOfFreedom();
    // tSparseMatf mat(num_of_dof, num_of_dof);
    tMatrixXf mat(num_of_dof, num_of_dof);
    mat.setZero();
    std::vector<tCudaMatrix3f> cuda_global_K_download;
    mGlobalStiffnessMatrix.Download(cuda_global_K_download);
    std::vector<tELLLocal2GlobalInfo> ell_local_vid_to_global_vid;
    mELLLocalVidToGlobalVid.Download(ell_local_vid_to_global_vid);
    std::vector<tTriplet> content = {};

    for (int i = 0; i < num_of_v; i++)
    {
        for (int j = 0; j < tELLLocal2GlobalInfo::mElements; j++)
        {
            int global_vid = ell_local_vid_to_global_vid[i][j];
            if (global_vid != -1)
            {
                // row major
                mat.block(3 * i, 3 * global_vid, 3, 3) +=
                    cCudaMatrixUtil::CudaMatrixToEigenMatrix(
                        cuda_global_K_download
                            [i * tELLLocal2GlobalInfo::mElements + j]);
                // std::cout << "row " << i << " col " << j << " = \n"
                //           << cuda_global_K_download[i * 12 + j] << std::endl;
            }
        }
    }
    return mat;
}