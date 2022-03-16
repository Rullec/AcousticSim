#include "BaraffClothGPU.h"
#include "geometries/Primitives.h"
#include "sim/cloth/BaraffMaterial.h"
#include "sim/cloth/QBendingMaterial.h"
#include "sim/gpu_utils/CudaDef.h"
#include "sim/gpu_utils/CudaMatrix.h"
#include "sim/gpu_utils/CudaMatrixUtil.h"
#include "utils/EigenUtil.h"
#include "utils/JsonUtil.h"
#include "utils/RotUtil.h"
#include "utils/TimeUtil.hpp"
#include <imgui.h>
#include <iostream>

template <typename int N, int M>
std::vector<tCudaMatrix<float, N, M>>
FetchFromGPU(const cCudaArray<tCudaMatrix<float, N, M>> &cuda_array)
{
    std::vector<tCudaMatrix<float, N, M>> cpu;
    cuda_array.Download(cpu);
    return cpu;
}

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

namespace FRPCGSolver
{
extern void
Solve(const cCudaArray<tCudaMatrix3f> &Winv_preconditioner,
      const cCudaArray<tCudaVector32i> &ELL_local_vertex_id_to_global_vertex_id,
      const cCuda2DArray<tCudaMatrix3f> &A, const cCudaArray<tCudaVector3f> &b,
      cCudaArray<tCudaVector3f> &x, cCudaArray<float> &debugEnergy,
      cCudaArray<float> &pap, cCudaArray<float> &rdotZ,
      cCudaArray<tCudaVector3f> &p, cCudaArray<tCudaVector3f> &Ap,
      cCudaArray<tCudaVector3f> &residual, cCudaArray<tCudaVector3f> &z);
};
namespace PCGSolver
{
extern void TestAtomicAdd(const cCudaArray<float> &data_array,
                          cCudaArray<float> &sum);
extern void
Solve(const cCudaArray<tCudaMatrix3f> &Winv_preconditioner,
      const cCudaArray<tCudaVector32i> &ELL_local_vertex_id_to_global_vertex_id,
      const cCuda2DArray<tCudaMatrix3f> &W, const cCudaArray<tCudaVector3f> &b,
      cCudaArray<tCudaVector3f> &pcg_rbuf, cCudaArray<tCudaVector3f> &pcg_dbuf,
      cCudaArray<tCudaVector3f> &pcg_zbuf, cCudaArray<float> &pcg_rMinvr_array,
      cCudaArray<tCudaVector3f> &x);

extern void CalcJacobiPreconditioner(
    const cCudaArray<tCudaVector32i> &ELL_local_vertex_id_to_global_vertex_id,
    const cCuda2DArray<tCudaMatrix3f> &W,
    cCudaArray<tCudaMatrix3f> &Winv_preconditioner);
extern void CalcBlockJacobiPreconditioner(
    const cCudaArray<tCudaVector32i> &ELL_local_vertex_id_to_global_vertex_id,
    const cCuda2DArray<tCudaMatrix3f> &W,
    cCudaArray<tCudaMatrix3f> &Winv_preconditioner);

extern void CalcNoPreconditioner(
    const cCudaArray<tCudaVector32i> &ELL_local_vertex_id_to_global_vertex_id,
    const cCuda2DArray<tCudaMatrix3f> &W,
    cCudaArray<tCudaMatrix3f> &Winv_preconditioner);

} // namespace PCGSolver

cBaraffClothGPU::cBaraffClothGPU(int id_)
    : cBaseCloth(eClothType::FEM_CLOTH_GPU, id_)
{
    mDragPointVertexId = -1;
    mDragPointTargetPos.setZero();
}
#include "geometries/Triangulator.h"
void cBaraffClothGPU::Init(const Json::Value &conf)
{
    cTimeUtil::Begin("init");
    // {
    //     int test_size = 158;
    //     cCudaArray<float> array_cuda;
    //     std::vector<float> array_cpu(test_size);
    //     float sum = 0;
    //     for (int i = 0; i < test_size; i++)
    //     {
    //         array_cpu[i] = cMathUtil::RandDouble(0, 3);
    //         sum += array_cpu[i];
    //     }
    //     array_cuda.Upload(array_cpu);
    //     cCudaArray<float> sum_cuda;
    //     std::vector<float> sum_cuda_download;
    //     sum_cuda.Upload({0.0});
    //     // calculate it in cuda
    //     PCGSolver::TestAtomicAdd(array_cuda, sum_cuda);
    //     sum_cuda.Download(sum_cuda_download);
    //     std::cout << "sum cuda =  " << sum_cuda_download[0] << std::endl;
    //     std::cout << "sum cpu =  " << sum << std::endl;

    //     exit(1);
    // };

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
        // cTimeUtil::Begin("CPU material");
        mStretchMaterialCPU = std::make_shared<cBaraffMaterial>();
        mStretchMaterialCPU->Init(shared_from_this(), mStretchK.cast<double>());
        mQBendingMaterialCPU = std::make_shared<cQBendingMaterial>();
        mQBendingMaterialCPU->Init(GetVertexArray(), GetEdgeArray(),
                                   GetTriangleArray(),
                                   mBendingK.cast<double>());
        // cTimeUtil::End("CPU material");
    }

    cTimeUtil::Begin("Init global stiffness");
    InitGlobalStiffnessMatrix();
    cTimeUtil::End("Init global stiffness");
    /*
        3. allocate: CPU buffer and GPU buffer
    */
    // cTimeUtil::Begin("Init CPU data");
    InitCpuData();
    UpdateCpuData();
    // cTimeUtil::End("Init CPU data");

    // cTimeUtil::Begin("Init GPU data");
    InitGpuData();
    UpdateGpuData();
    // cTimeUtil::End("Init GPU data");

    // begin to add random noise on vertex data, and then check the result
    {
        // tVectorXd cur_pos = GetPos();
        // cur_pos += tVectorXd::Random(cur_pos.size());
        // SetPos(cur_pos);
        // UpdateStiffnessMatrixAndIntForce();

        // show fint
        // {
        //     auto bending = FetchFromGPU(mIntForceCuda_bending);
        //     auto stretch = FetchFromGPU(mIntForceCuda_stretch);
        //     auto sum = FetchFromGPU(mIntForceCuda);
        //     for (int i = 0; i < bending.size(); i++)
        //     {
        //         std::cout << "---v " << i << "---\n";
        //         std::cout << "bending = " << bending[i].transpose()
        //                   << std::endl;
        //         std::cout << "stretch = " << stretch[i].transpose()
        //                   << std::endl;
        //         std::cout << "sum = " << sum[i].transpose() << std::endl;
        //     }
        // }
    }
    // 5. calculate energy, calculate internal force, calculate hessian
    // check the numerical diff
    // UpdateDataGPU();
    // cTimeUtil::Begin("verify data");
    // VerifyData();
    // cTimeUtil::End("verify data");
    cTimeUtil::End("init");
    // exit(1);
    // int dims = 5;
    // nBaraffClothGpu::Calc(dims);
    // 1. Ds^{-1}
    // CudaArray<tMatrix2f> mem;
    // mem.upload();

    // param.upload();
    // KERNEL_NAME CUDA_AT(mem.ptr(), )
    // exit(1);
}

cBaraffClothGPU::~cBaraffClothGPU() {}

typedef std::pair<std::string, float> tTimeRecms;
static std::vector<tTimeRecms> gProfRec;
void cBaraffClothGPU::UpdatePos(double dt)
{
    // std::cout << "dt = " << dt << std::endl;
    printf("-------new frame---------\n");
    dt = this->mIdealDefaultTimestep;
    UpdateFixedPt();
    // data, CPU->GPU, update pos, pre on CUDA, set dt
    // 1. calculate stiffness matrix
    gProfRec.clear();
    cTimeUtil::Begin("update K and fint");
    UpdateStiffnessMatrixAndIntForce();
    gProfRec.push_back(tTimeRecms("update K and fint",
                                  cTimeUtil::End("update K and fint", true)));
    // 2. system matrix assemble: W x = b
    cTimeUtil::Begin("assemble Ax=b");
    UpdateLinearSystem(dt);
    gProfRec.push_back(
        tTimeRecms("assemble Ax=b", cTimeUtil::End("assemble Ax=b", true)));
    // VerifyLinearSystem(dt);

    // 3. conjugate gradient, solve x

    cTimeUtil::Begin("solve");
    SolveLinearSystem();

    gProfRec.push_back(tTimeRecms("solve", cTimeUtil::End("solve", true)));

    // 4. update CPU rendering data

    cTimeUtil::Begin("post solve");
    PostSolve();
    gProfRec.push_back(
        tTimeRecms("post solve", cTimeUtil::End("post solve", true)));
}
#include "sim/Perturb.h"
void cBaraffClothGPU::ApplyUserPerturbForceOnce(tPerturb *pert)
{
    if (pert == nullptr)
        return;
    int max_id = 0;
    float max_bary = 0;
    auto tri = mTriangleArrayShared[pert->mAffectedTriId];
    std::vector<int> res = {tri->mId0, tri->mId1, tri->mId2};
    for (max_id = 0; max_id < 3; max_id++)
    {
        if (pert->mBarycentricCoords[max_id] > max_bary)
        {
            max_bary = pert->mBarycentricCoords[max_id];
            mDragPointVertexId = res[max_id];
        }
    }
    mDragPointTargetPos = cCudaMatrixUtil::EigenMatrixToCudaMatrix(
        tVector3f(pert->GetGoalPos().segment(0, 3).cast<float>()));
}

void cBaraffClothGPU::UpdateImGui()
{

    for (auto &x : gProfRec)
    {
        ImGui::Text("%s %d ms", x.first.c_str(), int(x.second));
    }
}

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
    mEleStretchStiffnessMatrixLstCuda.Resize(num_of_tris);
    mEleStretchStiffnessMatrixLstCuda.MemsetAsync(tCudaMatrix9f::Zero());

    InitQBendingHessian();
    BuildEdgeInfo();

    // confirm there is no zero
    SIM_ASSERT(mMassMatrixDiag.minCoeff() > 1e-7);

    // update mass matrix by diag
    {
        std::vector<tCudaVector3f> mass_matrix_diag_cpu(num_of_v);
        for (int i = 0; i < num_of_v; i++)
        {
            mass_matrix_diag_cpu[i] = cCudaMatrixUtil::EigenMatrixToCudaMatrix(
                tVector3f(mMassMatrixDiag.segment(3 * i, 3).cast<float>()));
            // std::cout << "mass diag for v " << i << " =  "
            //           << mass_matrix_diag_cpu[i].transpose() << std::endl;
        }
        mMassMatrixDiagCuda.Upload(mass_matrix_diag_cpu);
    }

    // init gravity
    {

        std::vector<tCudaVector3f> gravity_cpu(num_of_v, tCudaVector3f::Zero());
        mGravityForce.resize(GetNumOfFreedom());
        for (int i = 0; i < num_of_v; i++)
        {
            gravity_cpu[i] = tCudaVector3f(
                {0.0f, -9.8f * float(mMassMatrixDiag[3 * i]), 0.0f});
            mGravityForce.segment(3 * i, 3) =
                tVector3d(0.0, -9.8 * mMassMatrixDiag[3 * i], 0.0);
        }
        mGravityCuda.Upload(gravity_cpu);
    }

    // init ele internal force cuda
    {
        mEleStretchInternalForceCuda.Resize(num_of_tris);
        mEleStretchInternalForceCuda.MemsetAsync(tCudaVector9f::Zero());
    }

    // init vertex internal force cuda
    {
        mIntForceCuda_stretch.Resize(num_of_v);
        mIntForceCuda_stretch.MemsetAsync(tCudaVector3f::Zero());

        mIntForceCuda_bending.Resize(num_of_v);
        mIntForceCuda_bending.MemsetAsync(tCudaVector3f::Zero());

        mIntForceCuda.Resize(num_of_v);
        mIntForceCuda.MemsetAsync(tCudaVector3f::Zero());
    }

    // init stretch, bending, total stiffness cuda
    {
        int max_connected = tELLLocal2GlobalInfo::mElements;

        mGlobalStiffnessMatrix_bending.Resize(num_of_v, max_connected);
        mGlobalStiffnessMatrix_bending.MemsetAsync(tCudaMatrix3f::Zero());

        mGlobalStiffnessMatrix_stretch.Resize(num_of_v, max_connected);
        mGlobalStiffnessMatrix_stretch.MemsetAsync(tCudaMatrix3f::Zero());

        mGlobalStiffnessMatrix.Resize(num_of_v, max_connected);
        mGlobalStiffnessMatrix.MemsetAsync(tCudaMatrix3f::Zero());

        mSystemMatrixCuda.Resize(num_of_v, max_connected);
        mSystemMatrixCuda.MemsetAsync(tCudaMatrix3f::Zero());
    }

    // init system matrix and rhs
    {
        int max_connected = tELLLocal2GlobalInfo::mElements;
        mSystemMatrixCuda.Resize(num_of_v, max_connected);
        mSystemRHSCuda.Resize(num_of_v);
    }

    // init velocity on CUDA
    {
        int num_of_v = GetNumOfVertices();
        mVelCuda.Resize(num_of_v);
        mVelCuda.MemsetAsync(tCudaVector3f::Zero());
    }

    // init preconditioner on CUDA
    {
        this->mPreconditionerCuda.Resize(num_of_v);
        mPreconditionerCuda.MemsetAsync(tCudaMatrix3f::Zero());
    }

    // init solution vector on CUDA
    {
        mSolutionCuda.Resize(num_of_v);
        mSolutionCuda.MemsetAsync(tCudaVector3f::Zero());
    }

    // init buffers used in PCG
    {
        mPCGResidualCuda.Resize(num_of_v);  // r = b - A x
        mPCGDirectionCuda.Resize(num_of_v); // conjugate di
        mPCGzCuda.Resize(num_of_v);         // z = A * di
        // iterations: max 1000
        mPCG_rMinvr_arrayCuda.Resize(1000); // rMinvr, Minv is a
        mPCG_dTAd_arrayCuda.Resize(1);      // diT * A * di
    }

    // init the constraint point
    {
        // mFixedVertexIndices.Upload(mConstraint_StaticPointIds);
        mFixedVertexIndicesCPU = {};
        mFixedVertexTargetPosCPU = {};
        for (auto con_vid : mConstraint_StaticPointIds)
        {
            tCudaVector3f pos = cCudaMatrixUtil::EigenMatrixToCudaMatrix(
                tVector3f(GetVertexArray()[con_vid]
                              ->mPos.segment(0, 3)
                              .cast<float>()));
            mFixedVertexIndicesCPU.push_back(con_vid);
            mFixedVertexTargetPosCPU.push_back(pos);
            printf("[debug] fixed point v%d to (%.3f,%.3f,%.3f)\n", con_vid,
                   pos[0], pos[1], pos[2]);
        }
        // mFixedVertexTargetPos.Upload(fixed_pos);
    }
}

void cBaraffClothGPU::ClearDragPt()
{
    mDragPointVertexId = -1;
    mDragPointTargetPos.setZero();
}

void cBaraffClothGPU::UpdateFixedPt()
{
    std::vector<int> fixed_vertex_indices_cpu = mFixedVertexIndicesCPU;
    std::vector<tCudaVector3f> fixed_vertex_target_pos_cpu =
        mFixedVertexTargetPosCPU;
    // if (mDragPointVertexId != -1)
    // {
    //     fixed_vertex_indices_cpu.push_back(mDragPointVertexId);
    //     fixed_vertex_target_pos_cpu.push_back(mDragPointTargetPos);
    // }

    for (int i = 0; i < fixed_vertex_indices_cpu.size(); i++)
    {
        std::cout << "[debug] fix v" << fixed_vertex_indices_cpu[i] << " to "
                  << fixed_vertex_target_pos_cpu[i].transpose() << std::endl;
    }
    mFixedVertexIndicesCUDA.Upload(fixed_vertex_indices_cpu);
    mFixedVertexTargetPosCUDA.Upload(fixed_vertex_target_pos_cpu);
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

    // allocate user force CPU
    {
        mUserForceCpu.resize(num_of_v, tCudaVector3f::Zero());
        mUserForce.resize(3 * num_of_v);
        mUserForce.setZero();
    }
}

/**
 * \brief
 */
namespace BaraffClothGpu
{

extern void UpdateStretch_StiffnessMatrix_Fint(
    float K0, float K1, const cCudaArray<tCudaVector3f> &pos_lst,
    const cCudaArray<tCudaVector3i> &tri_vertices_lst,
    cCudaArray<tCudaMatrix32f> &N_lst, cCudaArray<tCudaMatrix32f> &F_lst,
    cCudaArray<tCudaMatrix32f> &n_lst,
    cCudaArray<tCudaVector2f> &C_lst, // condition
    cCudaArray<tCudaMatrix92f> &g_lst, cCudaArray<tCudaMatrix9f> &K_lst,
    cCudaArray<tCudaVector9f> &fint_lst_in_each_triangle);

extern void AssembleStretch_StiffnessMatrix_Fint(
    const cCudaArray<tCudaMatrix9f> &ele_K_lst_per_triangle,
    const cCudaArray<tCudaVector32i> &vertex_connected_triangle_lst,
    const cCudaArray<tCudaVector3i> &vertex_id_of_triangle_lst,
    const cCudaArray<tCudaVector32i> &ELL_local_vertex_id_to_global_vertex_id,
    const cCudaArray<tCudaVector9f> &ele_fint,
    const cCudaArray<int> &fixed_vertex_indices_lst,
    const cCudaArray<tCudaVector3f> &fixed_vertex_target_pos_lst,
    const cCudaArray<tCudaVector3f> &vertex_pos_lst,
    cCuda2DArray<tCudaMatrix3f> &global_K,
    cCudaArray<tCudaVector3f> &global_fint);

extern void AssembleBendingStiffnessMatrix(
    const cCudaArray<tCudaVector4i> &bending_vertex_per_edge_lst,
    const cCudaArray<tCudaVector32i> &vertex_connected_edges_lst,
    const cCudaArray<tCudaMatrix12f> &ele_K_lst,
    const cCudaArray<tCudaVector32i> &ELL_local_vertex_id_to_global_vertex_id,
    cCuda2DArray<tCudaMatrix3f> &global_K);
extern void UpdateBendingIntForce(
    const cCudaArray<tCudaVector32i> &ELL_local_vertex_id_to_global_vertex_id,
    const cCuda2DArray<tCudaMatrix3f> &bending_hessian,
    const cCudaArray<tCudaVector3f> &vertex_pos_lst,
    cCudaArray<tCudaVector3f> &bending_fint);

extern void AssembleSystemMatrix(
    float dt,
    const cCudaArray<tCudaVector32i> &ELL_local_vertex_id_to_global_vertex_id,
    const cCuda2DArray<tCudaMatrix3f> &K, const cCudaArray<tCudaVector3f> &M,
    cCuda2DArray<tCudaMatrix3f> &W);
extern void CalcCurVelocity(float dt, const cCudaArray<tCudaVector3f> &mXcur,
                            const cCudaArray<tCudaVector3f> &mXpre,
                            cCudaArray<tCudaVector3f> &mVel);

extern void AssembleSystemRHS(
    float dt, const cCudaArray<tCudaVector3f> &Gravity,
    const cCudaArray<tCudaVector3f> &UserForce,
    const cCudaArray<tCudaVector3f> &IntForce,
    const cCudaArray<tCudaVector32i> &ELL_local_vertex_id_to_global_vertex_id,
    const cCudaArray<tCudaVector3f> &M, const cCudaArray<tCudaVector3f> &cur_v,
    cCudaArray<tCudaVector3f> &RHS);

} // namespace BaraffClothGpu
void cBaraffClothGPU::UpdateGpuData()
{
    // 1. update Xcur and Xpre
    {
        mXcurCuda.Upload(mXcurCpu);
        mXpreCuda.Upload(mXpreCpu);
    }

    // 2. update F
    UpdateStiffnessMatrixAndIntForce();
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
            std::vector<tVector9f> EleFintCPU;
            mStretchMaterialCPU->GetFLst(mFLstCPU);
            mStretchMaterialCPU->GetFprimeLst(mFprimeLstCPU);
            mStretchMaterialCPU->GetnLst(nLstCPU);
            mStretchMaterialCPU->GetgLst(gLstCPU);
            mStretchMaterialCPU->GetCLst(CLstCPU);
            mStretchMaterialCPU->GetnprimeLst(nprimeLstCPU);
            mStretchMaterialCPU->GetgprimeLst(gprimeLstCPU);
            mStretchMaterialCPU->GetCprimeLst(CprimeLstCPU);
            mStretchMaterialCPU->GetEleKLst(KLstCPU);
            mStretchMaterialCPU->GetEleFintLst(EleFintCPU);
            // tMatrixXf global_K_cpu =
            //     mStretchMaterialCPU->CalcTotalStiffnessMatrix()
            //         .cast<float>()
            //         .toDense();
            tMatrixXf global_K_cpu =
                (mQBendingMaterialCPU->GetStiffnessMatrix().cast<float>() +
                 mStretchMaterialCPU->CalcTotalStiffnessMatrix().cast<float>())
                    .toDense();
            // CheckIsSyncMat("F", mFLstCPU, mFLstCuda);
            // CheckIsSyncMat("Fprime", mFprimeLstCPU, mFprimeLstCuda);

            // CheckIsSyncMat("n", nLstCPU, mnLstCuda);
            // CheckIsSyncMat("nprime", nprimeLstCPU, mnprimeLstCuda);

            // CheckIsSyncMat("g", gLstCPU, mgLstCuda);
            // CheckIsSyncMat("gprime", gprimeLstCPU, mgprimeLstCuda);

            // CheckIsSyncVec("C", CLstCPU, mCLstCuda);
            // CheckIsSyncVec("Cprime", CprimeLstCPU, mCprimeLstCuda);

            // CheckIsSyncVec("K", KLstCPU, mEleStiffnessMatrixLstCuda);
            {
                tMatrixXf global_K_cuda =
                    FromELLMatrixToEigenMatrix(mGlobalStiffnessMatrix);

                auto K_diff = (global_K_cuda - global_K_cpu);
                std::cout << "global K diff norm = " << K_diff.norm()
                          << std::endl;
                std::cout << "global K diff abs maxcoef = "
                          << K_diff.cwiseAbs().maxCoeff() << std::endl;
                // std::cout << "global K cuda = \n" << global_K_cuda <<
                // std::endl; std::cout << "global K cpu = \n" << global_K_cpu
                // << std::endl; check block
                int num_of_v = GetNumOfVertices();
                for (int i = 0; i < num_of_v; i++)
                    for (int j = 0; j < num_of_v; j++)
                    {
                        auto block_diff = K_diff.block(3 * i, 3 * j, 3, 3);
                        auto block_cuda =
                            global_K_cuda.block(3 * i, 3 * j, 3, 3);
                        auto block_cpu = global_K_cpu.block(3 * i, 3 * j, 3, 3);
                        float max_coef = block_diff.cwiseAbs().maxCoeff();
                        if (max_coef > 1e-3)
                        {
                            std::cout << "block " << i << " " << j
                                      << " diff = \n"
                                      << block_diff << std::endl;
                            std::cout << "cuda = \n" << block_cuda << std::endl;
                            std::cout << "cpu = \n" << block_cpu << std::endl;
                        }
                    }
            }
            std::vector<tVector9f> ele_fint_cuda_eigen(GetNumOfTriangles());
            {
                std::vector<tCudaVector9f> ele_fint_cpu;
                mEleStretchInternalForceCuda.Download(ele_fint_cpu);
                for (int i = 0; i < GetNumOfTriangles(); i++)
                {
                    ele_fint_cuda_eigen[i] =
                        cCudaMatrixUtil::CudaMatrixToEigenMatrix(
                            ele_fint_cpu[i]);
                }
            }

            {
                // tVectorXf fint_cuda =
                //     this->FromCudaVectorToEigenVector(this->mIntForceCuda);
                // tVectorXf fint_cpu =
                //     mStretchMaterialCPU->CalcTotalForce().cast<float>() +
                //     mQBendingMaterialCPU->CalcForce(mXcur).cast<float>();
                // tVectorXf fint_diff = fint_cuda - fint_cpu;
                // float max_coef = fint_diff.cwiseAbs().maxCoeff();
                // std::cout << "fint cpu =  " << fint_cpu.transpose()
                //           << std::endl;
                // std::cout << "fint cuda =  " << fint_cuda.transpose()
                //           << std::endl;
                // std::cout << "fint max diff = " << max_coef << std::endl;
                // for (int i = 0; i < GetNumOfVertices(); i++)
                // {
                //     tVector3f diff = fint_diff.segment(3 * i, 3);
                //     float max_coef = diff.cwiseAbs().maxCoeff();
                //     if (max_coef > 1e-5)
                //     {
                //         std::cout << "v" << i << " fint cpu = "
                //                   << fint_cpu.segment(3 * i, 3).transpose()
                //                   << " cuda = "
                //                   << fint_cuda.segment(3 * i, 3).transpose()
                //                   << std::endl;
                //     }
                // }

                // for (int i = 0; i < GetNumOfTriangles(); i++)
                // {
                //     auto cuda = ele_fint_cuda_eigen[i];
                //     auto cpu = EleFintCPU[i];
                //     auto diff = cuda - cpu;
                //     float max_coef = diff.cwiseAbs().maxCoeff();
                //     if (max_coef > 1e-5)
                //     {
                //         std::cout << "triangle " << i << " cuda = \n"
                //                   << cuda << " cpu = \n"
                //                   << cpu << " max coef diff =  " << max_coef
                //                   << std::endl;
                //     }
                //     else
                //     {
                //         std::cout << "triangle " << i
                //                   << " ele fint verify succ\n";
                //     }
                // }
            }
            // exit(1);
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
    int num_of_v = this->GetNumOfVertices();
    for (int i = 0; i < num_of_v; i++)
    {
        mXcurCpu[i] = cCudaMatrixUtil::EigenMatrixToCudaMatrix(
            tVector3f(newpos.segment(3 * i, 3).cast<float>()));
    }
    mXcurCuda.Upload(mXcurCpu);
}
void cBaraffClothGPU::ClearForce()
{
    cBaseCloth::ClearForce();
    mGravityCuda.MemsetAsync(tCudaVector3f::Zero());
    mUserForceCuda.MemsetAsync(tCudaVector3f::Zero());
    mCollisionForce.MemsetAsync(tCudaVector3f::Zero());
    mIntForceCuda_stretch.MemsetAsync(tCudaVector3f::Zero());
    mIntForceCuda_bending.MemsetAsync(tCudaVector3f::Zero());
    mIntForceCuda.MemsetAsync(tCudaVector3f::Zero());
}
#include <cuda.h>
void cBaraffClothGPU::Reset()
{
    std::cout << "[log] Cloth reset\n";
    SetPos(mClothInitPos);
    mXpre.noalias() = mClothInitPos;
    ClearForce();
    int num_of_v = this->GetNumOfVertices();
    for (int i = 0; i < num_of_v; i++)
    {
        mXpreCpu[i] = cCudaMatrixUtil::EigenMatrixToCudaMatrix(
            tVector3f(mXpre.segment(3 * i, 3).cast<float>()));
    }
    mXpreCuda.Upload(mXpreCpu);
    InitGpuData();
    UpdateGpuData();
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
    // 1. for each edge, vertices included in the connected triangle has effect
    // with each other
    int num_of_v = GetNumOfVertices();
    std::vector<std::set<int>> vertices_included_vertices(num_of_v);
    // 1. add bending constraint
    std::set<int> affected_vertices = {};

    for (int i = 0; i < GetNumOfEdges(); i++)
    {
        affected_vertices.clear();
        auto edge = GetEdgeArray()[i];
        // stretch
        auto t0 = GetTriangleArray()[edge->mTriangleId0];
        affected_vertices.insert(t0->mId0);
        affected_vertices.insert(t0->mId1);
        affected_vertices.insert(t0->mId2);

        // bending
        if (edge->mIsBoundary == false)
        {
            auto t1 = GetTriangleArray()[edge->mTriangleId1];
            affected_vertices.insert(t1->mId0);
            affected_vertices.insert(t1->mId1);
            affected_vertices.insert(t1->mId2);
        }
        // add
        for (auto &v : affected_vertices)
        {
            vertices_included_vertices[v].insert(affected_vertices.begin(),
                                                 affected_vertices.end());
        }
    }

    std::vector<tELLLocal2GlobalInfo> ELL_local_vid_to_global_vid(
        num_of_v, tELLLocal2GlobalInfo::Ones() * -1);
    for (int i = 0; i < num_of_v; i++)
    {
        auto affected_v_lst = vertices_included_vertices[i];
        if (affected_v_lst.size() > tELLLocal2GlobalInfo::mElements)
        {
            printf("[error] vertex %d has more than %d entries in K, exit\n", i,
                   tELLLocal2GlobalInfo::mElements);
            exit(1);
        }
        for (auto &x : affected_v_lst)
        {
            AddInfo(ELL_local_vid_to_global_vid[i], x);
        }
    }
    mELLLocalVidToGlobalVid.Upload(ELL_local_vid_to_global_vid);

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
    std::vector<tConnectedInfo> vertex_connected_edges(num_of_v, -1);
    std::vector<tConnectedInfo> vertex_involved_inner_edges(num_of_v, -1);

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

tMatrixXf cBaraffClothGPU::FromELLMatrixToEigenMatrix(
    const cCuda2DArray<tCudaMatrix3f> &ell_mat)
{
    int num_of_v = this->GetNumOfVertices();
    int num_of_dof = GetNumOfFreedom();
    // tSparseMatf mat(num_of_dof, num_of_dof);
    tMatrixXf mat(num_of_dof, num_of_dof);
    mat.setZero();
    std::vector<tCudaMatrix3f> cuda_global_K_download;
    ell_mat.Download(cuda_global_K_download);
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
namespace GPUMatrixOps
{
extern void ELLMatrixAdd(const cCuda2DArray<tCudaMatrix3f> &mat0,
                         const cCuda2DArray<tCudaMatrix3f> &mat1,
                         cCuda2DArray<tCudaMatrix3f> &mat2);

extern void VectorAdd(const cCudaArray<tCudaVector3f> &vec0,
                      const cCudaArray<tCudaVector3f> &vec1,
                      cCudaArray<tCudaVector3f> &vec2);
} // namespace GPUMatrixOps
/**
 * \brief           update stiffness matrix on GPU
 */
void cBaraffClothGPU::UpdateStiffnessMatrixAndIntForce()
{
    // clear
    {
        // 1. clear stiffness matrix
        mGlobalStiffnessMatrix_stretch.MemsetAsync(tCudaMatrix3f::Zero());
        mGlobalStiffnessMatrix_bending.MemsetAsync(tCudaMatrix3f::Zero());
        mGlobalStiffnessMatrix.MemsetAsync(tCudaMatrix3f::Zero());
        mEleStretchStiffnessMatrixLstCuda.MemsetAsync(tCudaMatrix9f::Zero());

        // 2. clear internal force
        mEleStretchInternalForceCuda.MemsetAsync(tCudaVector9f::Zero());
        mIntForceCuda_stretch.MemsetAsync(tCudaVector3f::Zero());
        mIntForceCuda_bending.MemsetAsync(tCudaVector3f::Zero());
        mIntForceCuda.MemsetAsync(tCudaVector3f::Zero());
    }
    // update stretch stiffness and fint
    BaraffClothGpu::UpdateStretch_StiffnessMatrix_Fint(
        mStretchK[0], mStretchK[1], mXcurCuda, mTriangleVerticesIdLstCuda,
        mNLstCuda, mFLstCuda, mnLstCuda, mCLstCuda, mgLstCuda,
        mEleStretchStiffnessMatrixLstCuda, mEleStretchInternalForceCuda);
    BaraffClothGpu::UpdateStretch_StiffnessMatrix_Fint(
        mStretchK[2], mStretchK[2], mXcurCuda, mTriangleVerticesIdLstCuda,
        mNprimeLstCuda, mFprimeLstCuda, mnprimeLstCuda, mCprimeLstCuda,
        mgprimeLstCuda, mEleStretchStiffnessMatrixLstCuda,
        mEleStretchInternalForceCuda);

    BaraffClothGpu::AssembleStretch_StiffnessMatrix_Fint(
        mEleStretchStiffnessMatrixLstCuda, mVerticesTriangleIdLstCuda,
        mTriangleVerticesIdLstCuda, mELLLocalVidToGlobalVid,
        this->mEleStretchInternalForceCuda, this->mFixedVertexIndicesCUDA,
        this->mFixedVertexTargetPosCUDA, this->mXcurCuda,
        this->mGlobalStiffnessMatrix_stretch, mIntForceCuda_stretch);

    // update bending stiffness and fint
    BaraffClothGpu::AssembleBendingStiffnessMatrix(
        mQBendingConstraintVertexLst, mVertexInvolvedInnerEdge,
        mQBendingElementKLst, mELLLocalVidToGlobalVid,
        mGlobalStiffnessMatrix_bending);

    BaraffClothGpu::UpdateBendingIntForce(
        mELLLocalVidToGlobalVid, mGlobalStiffnessMatrix_bending,
        this->mXcurCuda, mIntForceCuda_bending);

    // sum global stiffness
    GPUMatrixOps::ELLMatrixAdd(mGlobalStiffnessMatrix_bending,
                               mGlobalStiffnessMatrix_stretch,
                               mGlobalStiffnessMatrix);
    // sum global fint
    // mIntForceCuda_stretch.MemsetAsync(tCudaVector3f::Zero());
    // printf("[error] int force is set to zero\n");
    GPUMatrixOps::VectorAdd(mIntForceCuda_bending, mIntForceCuda_stretch,
                            mIntForceCuda);
}

/**
 * \brief           Set up the linear system Wx = b
 *  dt2K = dt * dt * K
 *  W = M - dt2K
 *  b = dt * dt * (mGravityForce + mUserForce + mIntForce)
 *      + dt * (W + dt2K) * V_cur
 *  V_cur = (mXcur - mXpre) / dt
 */
void cBaraffClothGPU::UpdateLinearSystem(float dt)
{
    /* 1. upload all force
            user force
            gravity
            collision force
    */
    UpdateCollisionForceAndUserForce();

    /*
        3. assemble W
    */
    BaraffClothGpu::AssembleSystemMatrix(
        dt, this->mELLLocalVidToGlobalVid, this->mGlobalStiffnessMatrix,
        this->mMassMatrixDiagCuda, this->mSystemMatrixCuda);

    // 4. assemble b
    // 4.1 calculate current velocity
    BaraffClothGpu::CalcCurVelocity(dt, mXcurCuda, mXpreCuda, this->mVelCuda);

    // 4.2 calcualte RHS
    // extern void AssembleSystemRHS(
    // float dt, const cCudaArray<tCudaVector3f> &Gravity,
    // const cCudaArray<tCudaVector3f> &UserForce,
    // const cCudaArray<tCudaVector3f> &IntForce,
    // const cCudaArray<tCudaVector32i>
    // &ELL_local_vertex_id_to_global_vertex_id, const cCudaArray<tCudaVector3f>
    // &W, const cCudaArray<tCudaVector3f> &K, const cCudaArray<tCudaVector3f>
    // &cur_v, cCudaArray<tCudaVector3f> &RHS);
    BaraffClothGpu::AssembleSystemRHS(
        dt, this->mGravityCuda, this->mUserForceCuda, this->mIntForceCuda,
        this->mELLLocalVidToGlobalVid, this->mMassMatrixDiagCuda,
        this->mVelCuda, this->mSystemRHSCuda);
}
// namespace GDSolver
// {
// extern void
// Solve(const cCudaArray<tCudaVector32i>
// &ELL_local_vertex_id_to_global_vertex_id,
//       const cCuda2DArray<tCudaMatrix3f> &A, const cCudaArray<tCudaVector3f>
//       &b, cCudaArray<tCudaVector3f> &x, cCudaArray<tCudaVector3f> &r_buf);
// }
namespace JacobSolver
{
extern void
Solve(const cCudaArray<tCudaVector32i> &ELL_local_vertex_id_to_global_vertex_id,
      const cCuda2DArray<tCudaMatrix3f> &A, const cCudaArray<tCudaVector3f> &b,
      cCudaArray<tCudaVector3f> &x, cCudaArray<tCudaVector3f> &r_buf);
}
void cBaraffClothGPU::SolveLinearSystem()
{
    // JacobSolver::Solve(mELLLocalVidToGlobalVid, mSystemMatrixCuda,
    //                    mSystemRHSCuda, mSolutionCuda, mPCGResidualCuda);
    // GDSolver::Solve(mELLLocalVidToGlobalVid, mSystemMatrixCuda,
    // mSystemRHSCuda,
    //                 mSolutionCuda, mPCGResidualCuda);
    // 1. create the preconditioner
    PCGSolver::CalcBlockJacobiPreconditioner(
        mELLLocalVidToGlobalVid, mSystemMatrixCuda, mPreconditionerCuda);

    // 2. solve
    PCGSolver::Solve(mPreconditionerCuda, mELLLocalVidToGlobalVid,
                     mSystemMatrixCuda, mSystemRHSCuda, this->mPCGResidualCuda,
                     this->mPCGDirectionCuda, this->mPCGzCuda,
                     this->mPCG_rMinvr_arrayCuda, this->mSolutionCuda);
    // {
    // cCudaArray<float> debugEnergy, pap, rdotz;
    // cCudaArray<tCudaVector3f> Ap;
    // Ap.Resize(GetNumOfVertices());

    // FRPCGSolver::Solve(mPreconditionerCuda, mELLLocalVidToGlobalVid,
    //                    mSystemMatrixCuda, mSystemRHSCuda, mSolutionCuda,
    //                    debugEnergy, pap, rdotz, this->mPCGzCuda, Ap,
    //                    this->mPCGResidualCuda, this->mPCGDirectionCuda);
    // }

    // 3. copy back result, verification
    // tMatrixXf W_cpu = this->FromELLMatrixToEigenMatrix(mSystemMatrixCuda);
    // tVectorXf b_cpu = FromCudaVectorToEigenVector(this->mSystemRHSCuda);
    // tVectorXf x_cpu = FromCudaVectorToEigenVector(this->mSolutionCuda);
    // tVectorXf residual = W_cpu * x_cpu - b_cpu;
    // std::cout << "[solve] residual = " << residual.cwiseAbs().maxCoeff()
    //           << std::endl;
    // exit(1);
}

void cBaraffClothGPU::PostSolve()
{
    // 1. copy back the dx data from GPU
    tVectorXf dx = FromCudaVectorToEigenVector(this->mSolutionCuda);
    this->mXpre = mXcur;
    mXcur += dx.cast<double>();

    // 2. update cpu data (mXcur, mXpre)
    this->SetPos(mXcur);

    // 3. update the GPU data
    int num_of_v = this->GetNumOfVertices();
    for (int i = 0; i < num_of_v; i++)
    {
        mXpreCpu[i] = cCudaMatrixUtil::EigenMatrixToCudaMatrix(
            tVector3f(mXpre.segment(3 * i, 3).cast<float>()));
        mXcurCpu[i] = cCudaMatrixUtil::EigenMatrixToCudaMatrix(
            tVector3f(mXcur.segment(3 * i, 3).cast<float>()));
    }
    mXpreCuda.Upload(mXpreCpu);

    ClearDragPt();
}

tVectorXf cBaraffClothGPU::FromCudaVectorToEigenVector(
    const cCudaArray<tCudaVector3f> &cuda_vec)
{
    int num_of_v = this->GetNumOfVertices();
    int num_of_dof = GetNumOfFreedom();
    // tSparseMatf mat(num_of_dof, num_of_dof);
    tVectorXf vec(num_of_dof);
    vec.setZero();

    std::vector<tCudaVector3f> int_force_cpu;
    cuda_vec.Download(int_force_cpu);
    assert(int_force_cpu.size() == num_of_v);
    for (int i = 0; i < num_of_v; i++)
    {
        vec.segment(3 * i, 3).noalias() =
            cCudaMatrixUtil::CudaMatrixToEigenMatrix(int_force_cpu[i]);
    }
    return vec;
}

/**
 *\brief            upload user force and collision force
 */
void cBaraffClothGPU::UpdateCollisionForceAndUserForce()
{
    int num_of_v = GetNumOfVertices();

    mUserForceCuda.Resize(num_of_v);
    mCollisionForce.Resize(num_of_v);

    mUserForceCuda.MemsetAsync(tCudaVector3f::Zero());
    mCollisionForce.MemsetAsync(tCudaVector3f::Zero());

    // begin to add user for
    if (this->mDragPointVertexId != -1)
    {
        this->mUserForceCpu.resize(num_of_v);
        for (int i = 0; i < num_of_v; i++)
            mUserForceCpu[i].setZero();
        float k = 1;
        mUserForceCpu[mDragPointVertexId] +=
            k * (mDragPointTargetPos - mXcurCpu[mDragPointVertexId]);

        std::cout << "[fixed] add force = "
                  << mUserForceCpu[mDragPointVertexId].transpose()
                  << " tar pos = " << mDragPointTargetPos.transpose()
                  << " cur pos = " << mXcurCpu[mDragPointVertexId].transpose()
                  << " cur new = "
                  << mXcur.segment(3 * mDragPointVertexId, 3).transpose()
                  << std::endl;
        mUserForceCuda.Upload(mUserForceCpu);
    }
}

void cBaraffClothGPU::VerifyLinearSystem(float dt)
{
    // 1. get the cpu linear system (barrow from CPU version)
    tMatrixXf W_cpu;
    tVectorXf b_cpu;
    {
        int dof = GetNumOfFreedom();
        tSparseMatf M(dof, dof);
        M.reserve(dof);
        // M.diagonal() = mMassMatDiag;
        for (size_t i = 0; i < dof; i++)
        {
            M.coeffRef(i, i) = mMassMatrixDiag[i];
            // I.coeff(i, i) = 1;
        }
        mStretchMaterialCPU->Update();
        auto global_K_cpu =
            (mQBendingMaterialCPU->GetStiffnessMatrix().cast<float>() +
             mStretchMaterialCPU->CalcTotalStiffnessMatrix().cast<float>())
                .toDense();
        W_cpu = (M - dt * dt * global_K_cpu).cast<float>();

        // set gravity, set user force, set int force
        tVectorXf int_force = (this->mStretchMaterialCPU->CalcTotalForce() +
                               this->mQBendingMaterialCPU->CalcForce(mXcur))
                                  .cast<float>();
        b_cpu = dt * dt *
                    (mGravityForce.cast<float>() + mUserForce.cast<float>() +
                     int_force) +
                dt * (W_cpu + dt * dt * global_K_cpu) *
                    (mXcur - mXpre).cast<float>() / dt;
    }
    // 2. get the cuda version W and residual
    tMatrixXf W_cuda = this->FromELLMatrixToEigenMatrix(mSystemMatrixCuda);
    tMatrixXf b_cuda = FromCudaVectorToEigenVector(mSystemRHSCuda);

    tMatrixXf W_diff = W_cuda - W_cpu;
    tVectorXf b_diff = b_cpu - b_cuda;
    std::cout << "[verify] W diff = " << W_diff.cwiseAbs().maxCoeff()
              << std::endl;
    std::cout << "[verify] b diff = " << b_diff.cwiseAbs().maxCoeff()
              << std::endl;
}
