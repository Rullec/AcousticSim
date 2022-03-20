#include "BaraffClothGpu.h"
#include "geometries/Primitives.h"
#include "sim/Perturb.h"
#include "sim/gpu_solvers/PCGSolver.h"
#include "sim/gpu_utils/CudaMatrixUtil.h"
#include "utils/JsonUtil.h"
#include "utils/RotUtil.h"
#include <set>

cBaraffClothGpu::cBaraffClothGpu(int id_)
    : cBaseCloth(eClothType::FEM_CLOTH_GPU, id_)
{
}

void cBaraffClothGpu::Init(const Json::Value &conf)
{
    // 1. load config
    mBendingK = cCudaMatrixUtil::EigenMatrixToCudaMatrix(
        tVector3f(cJsonUtil::ReadVectorJson(
                      cJsonUtil::ParseAsValue("cloth_bending_stiffness", conf))
                      .cast<float>()));
    mStretchK = cCudaMatrixUtil::EigenMatrixToCudaMatrix(
        tVector3f(cJsonUtil::ReadVectorJson(
                      cJsonUtil::ParseAsValue("cloth_stretch_stiffness", conf))
                      .segment(0, 3)
                      .cast<float>()));

    mRayleighDampingA = cJsonUtil::ParseAsDouble("rayleigh_damping_a", conf);
    mRayleighDampingB = cJsonUtil::ParseAsDouble("rayleigh_damping_b", conf);
    // mSolver = std::make_shared<cPCGSolver>();
    cBaseCloth::Init(conf);

    // 2. init data
    InitELLMatrix();
    InitStretchCoef();
    UpdatePosToCuda(mXcur, mXcurCuda);
    UpdatePosToCuda(mXpre, mXpreCuda);
    InitGravity();
}
namespace BaraffClothGpu
{
extern void UpdateStiffnessMatrixAndFint(
    const cCudaArray<tCudaVector3i> &mTriangleVertexIdCuda,
    const cCudaArray<float> &mTriangleInitAreaCuda,
    const cCudaArray<tCudaVector3f> &mXcurCuda, const tCudaVector3f &mStretchK,
    const cCudaArray<tCudaVector3f> &mCoefFu_warp_weft,
    const cCudaArray<tCudaVector3f> &mCoefFv_warp_weft,
    const cCudaArray<tCudaVector3f> &mCoefFu_diag,
    const cCudaArray<tCudaVector3f> &mCoefFv_diag,
    const cCudaArray<tCudaVector32i> &mELLVidToGlobalVid,
    cCudaArray<tCudaVector3f> &mIntForceCuda,
    cCuda2DArray<tCudaMatrix3f> &mStiffnessMatrixCuda);
extern void AssembleSystemRHS(
    float dt, const cCudaArray<tCudaVector3f> &Gravity,
    const cCudaArray<tCudaVector3f> &IntForce,
    const cCudaArray<tCudaVector32i> &ELL_local_vertex_id_to_global_vertex_id,
    const cCudaArray<float> &vertices_mass,
    const cCudaArray<tCudaVector3f> &x_cur,
    const cCudaArray<tCudaVector3f> &x_pre, cCudaArray<tCudaVector3f> &RHS);
extern void AssembleSystemMatrix(
    float dt, const float rayleigh_damping_a, const float rayleigh_damping_b,
    const cCudaArray<tCudaVector32i> &ELL_local_vertex_id_to_global_vertex_id,
    const cCuda2DArray<tCudaMatrix3f> &K,
    const cCudaArray<float> &vertices_mass, cCuda2DArray<tCudaMatrix3f> &A);
} // namespace BaraffClothGpu

namespace PCGSolver
{
extern void
Solve(const cCudaArray<tCudaVector32i> &ELL_local_vertex_id_to_global_vertex_id,
      const cCuda2DArray<tCudaMatrix3f> &A, const cCudaArray<tCudaVector3f> &b,
      cCudaArray<tCudaVector3f> &x, const cCudaArray<int> &mConstraintVertexId);
};
tVectorXf FromCudaVectorToEigenVector(const cCudaArray<tCudaVector3f> &cuda_vec)
{

    // tSparseMatf mat(num_of_dof, num_of_dof);
    tVectorXf vec(cuda_vec.Size() * 3);
    vec.setZero();

    std::vector<tCudaVector3f> int_force_cpu;
    cuda_vec.Download(int_force_cpu);
    // assert(int_force_cpu.size() == num_of_v);
    for (int i = 0; i < cuda_vec.Size(); i++)
    {
        vec.segment(3 * i, 3).noalias() =
            cCudaMatrixUtil::CudaMatrixToEigenMatrix(int_force_cpu[i]);
    }
    return vec;
}

tMatrixXf FromELLMatrixToEigenMatrix(
    const cCuda2DArray<tCudaMatrix3f> &ell_mat,
    const cCudaArray<tCudaVector32i> &mELLLocalVidToGlobalVid)
{
    int num_of_v = ell_mat.Rows();
    int num_of_dof = 3 * num_of_v;
    // tSparseMatf mat(num_of_dof, num_of_dof);
    tMatrixXf mat(num_of_dof, num_of_dof);
    mat.setZero();
    std::vector<tCudaMatrix3f> cuda_global_K_download;
    ell_mat.Download(cuda_global_K_download);
    std::vector<tCudaVector32i> ell_local_vid_to_global_vid;
    mELLLocalVidToGlobalVid.Download(ell_local_vid_to_global_vid);
    // std::vector<tTriplet> content = {};

    for (int i = 0; i < num_of_v; i++)
    {
        for (int j = 0; j < tCudaVector32i::mElements; j++)
        {
            int global_vid = ell_local_vid_to_global_vid[i][j];
            if (global_vid != -1)
            {
                // row major
                mat.block(3 * i, 3 * global_vid, 3, 3) +=
                    cCudaMatrixUtil::CudaMatrixToEigenMatrix(
                        cuda_global_K_download[i * tCudaVector32i::mElements +
                                               j]);
                // std::cout << "row " << i << " col " << j << " = \n"
                //           << cuda_global_K_download[i * 12 + j] << std::endl;
            }
        }
    }
    return mat;
}
void cBaraffClothGpu::UpdatePos(double dt)
{
    dt = mIdealDefaultTimestep;
    // 1. update fixt pt
    UpdateConstraintVertexPos();
    UpdateConstraintVid(mConstraintVertexId);

    // 2. update K
    // printf("fint size enter = %d\n", mIntForceCuda.Size());
    {
        BaraffClothGpu::UpdateStiffnessMatrixAndFint(
            mTriangleVertexIdCuda, mTriangleInitAreaCuda, mXcurCuda, mStretchK,
            mCoefFu_warp_weft, mCoefFv_warp_weft, mCoefFu_diag, mCoefFv_diag,
            mELLVidToGlobalVid, mIntForceCuda, mStiffnessMatrixCuda);
        // tMatrixXf K = FromELLMatrixToEigenMatrix(mStiffnessMatrixCuda,
        //                                          this->mELLVidToGlobalVid);
        // std::cout << "fint = "
        //           << FromCudaVectorToEigenVector(mIntForceCuda).transpose()
        //           << std::endl;
        // std::cout << " K = \n" << K << std::endl;
    }
    // exit(1);
    // 3. create Ax=b
    {
        BaraffClothGpu::AssembleSystemRHS(
            dt, this->mGravityForceCuda, this->mIntForceCuda,
            this->mELLVidToGlobalVid, this->mVerticesMassCuda, this->mXcurCuda,
            mXpreCuda, this->mSystemRHSCuda);
        // std::cout << "cuda b = "
        //           << FromCudaVectorToEigenVector(mSystemRHSCuda).transpose()
        //           << std::endl;
        BaraffClothGpu::AssembleSystemMatrix(
            dt, this->mRayleighDampingA, this->mRayleighDampingB,
            this->mELLVidToGlobalVid, this->mStiffnessMatrixCuda,
            this->mVerticesMassCuda, this->mSystemMatrixCuda);
        // std::cout << "cuda A = "
        //           << FromELLMatrixToEigenMatrix(mSystemMatrixCuda,
        //                                         this->mELLVidToGlobalVid)
        //           << std::endl;
        // exit(1);
    }
    // 4. solve Ax = b, forward the system, data synchronization
    // mConstraintVertexId.Clear();
    PCGSolver::Solve(mELLVidToGlobalVid, mSystemMatrixCuda, mSystemRHSCuda,
                     this->mSolutionCuda, mConstraintVertexId);

    {
        // solution = dx
        tVectorXf dx = FromCudaVectorToEigenVector(this->mSolutionCuda);
        // std::cout << "dx = " << dx.transpose() << std::endl;
        this->mXpre = mXcur;
        mXcur += dx.cast<double>();
        SetPos(mXcur);
        UpdatePosToCuda(mXcur, mXcurCuda);
        UpdatePosToCuda(mXpre, mXpreCuda);
    }
    // 5. clear drag pt
    mDragPt = nullptr;
    mConstraintVertexId.Resize(0);
}
void cBaraffClothGpu::ApplyUserPerturbForceOnce(tPerturb *pert)
{
    if (pert != nullptr)
    {
        mDragPt = std::make_shared<tFixPoint>();
        mDragPt->mTargetPos =
            pert->GetGoalPos().segment(0, 3).segment(0, 3).cast<float>();

        int ids[3] = {mTriangleArrayShared[pert->mAffectedTriId]->mId0,
                      mTriangleArrayShared[pert->mAffectedTriId]->mId1,
                      mTriangleArrayShared[pert->mAffectedTriId]->mId2};

        int max_id = cMathUtil::Argmax(pert->mBarycentricCoords);
        mDragPt->mVertexId = ids[max_id];
    }
}
void cBaraffClothGpu::UpdateImGui() {}

void cBaraffClothGpu::Reset() { cBaseCloth::Reset(); }

/**
 * \brief           sync a cpu pos data to GPU
 */
void cBaraffClothGpu::UpdatePosToCuda(const tVectorXd &cpu_pos_eigen,
                                      cCudaArray<tCudaVector3f> &cuda_pos) const
{
    int num_of_v = static_cast<int>(cpu_pos_eigen.size() / 3);

    std::vector<tCudaVector3f> cpu_pos(num_of_v);
    for (int i = 0; i < num_of_v; i++)
    {
        cpu_pos[i] = cCudaMatrixUtil::EigenMatrixToCudaMatrix(
            tVector3f(cpu_pos_eigen.segment(3 * i, 3).cast<float>()));
    }
    cuda_pos.Upload(cpu_pos);
}

/**
 * \brief           fetch GPU pos data to cpu
 */
void cBaraffClothGpu::UpdatePosToCpu(const cCudaArray<tCudaVector3f> &cuda_pos,
                                     tVectorXd &cpu_pos_eigen) const
{
    int num_of_v = cuda_pos.Size();

    std::vector<tCudaVector3f> cpu_pos(num_of_v);
    cuda_pos.Download(cpu_pos);
    cpu_pos_eigen.resize(3 * num_of_v);

    for (int i = 0; i < num_of_v; i++)
    {
        cpu_pos_eigen.segment(3 * i, 3) =
            cCudaMatrixUtil::CudaMatrixToEigenMatrix(cpu_pos[i]).cast<double>();
    }
}

void cBaraffClothGpu::InitGeometry(const Json::Value &conf)
{
    cBaseCloth::InitGeometry(conf);
    // update mTriangleAreaCuda
    mTriangleInitAreaCuda.Upload(mTriangleInitArea);

    // update mTriangleVertexIdCuda
    int num_of_tri = GetNumOfTriangles();
    std::vector<tCudaVector3i> tri_vertex_id_cpu(num_of_tri);
    for (int i = 0; i < num_of_tri; i++)
    {
        auto tri = mTriangleArrayShared[i];
        tri_vertex_id_cpu[i][0] = tri->mId0;
        tri_vertex_id_cpu[i][1] = tri->mId1;
        tri_vertex_id_cpu[i][2] = tri->mId2;
        // std::cout << "tri " << i << " v = " << tri_vertex_id_cpu[i].transpose()
        //           << std::endl;
    }
    mTriangleVertexIdCuda.Upload(tri_vertex_id_cpu);

    // update
}

void cBaraffClothGpu::InitMass(const Json::Value &conf)
{
    cBaseCloth::InitMass(conf);
    // upload vertices mass to cuda
    int num_of_v = GetNumOfVertices();
    std::vector<float> vertices_mass(num_of_v);
    for (int i = 0; i < num_of_v; i++)
    {
        vertices_mass[i] = mVertexArrayShared[i]->mMass;
    }
    mVerticesMassCuda.Upload(vertices_mass);
}

void cBaraffClothGpu::InitConstraint(const Json::Value &root)
{
    cBaseCloth::InitConstraint(root);
    mConstrainedPt.clear();
    for (auto &vid : mConstraint_StaticPointIds)
    {
        auto pt = std::make_shared<tFixPoint>();
        pt->mVertexId = vid;
        pt->mTargetPos = this->GetInitPos().segment(3 * vid, 3).cast<float>();
        mConstrainedPt.push_back(pt);
        mDragPt = nullptr;
    }
}
/**
 * \brief           Update constraint vid
 */
void cBaraffClothGpu::UpdateConstraintVid(
    cCudaArray<int> &constrained_vid) const
{

    std::vector<int> cpu_vid = {};
    for (auto &x : mConstrainedPt)
    {
        std::cout << "fix " << x->mVertexId << std::endl;
        cpu_vid.push_back(x->mVertexId);
    }

    if (mDragPt != nullptr)
    {
        std::cout << "fix " << mDragPt->mVertexId << std::endl;
        cpu_vid.push_back(mDragPt->mVertexId);
    }
    constrained_vid.Upload(cpu_vid);
}

/**
 * \brief           update constraint vertex pos to both CPU and GPU
 */
void cBaraffClothGpu::UpdateConstraintVertexPos()
{
    if (mDragPt != nullptr)
        mConstrainedPt.push_back(mDragPt);

    for (auto &x : mConstrainedPt)
    {
        tCudaVector3f pos =
            cCudaMatrixUtil::EigenMatrixToCudaMatrix(x->mTargetPos);
        int vid = x->mVertexId;

        // 1. update cpu data
        mXcur.segment(3 * vid, 3) = x->mTargetPos.cast<double>();
        mXpre.segment(3 * vid, 3) = x->mTargetPos.cast<double>();
        mVertexArrayShared[vid]->mPos.segment(0, 3) =
            x->mTargetPos.cast<double>();

        // 2. update gpu data
        // mXcurCuda.Upload(&pos, vid);
        // mXpreCuda.Upload(&pos, vid);
    }
    this->UpdatePosToCuda(mXcur, mXcurCuda);
    this->UpdatePosToCuda(mXpre, mXpreCuda);
    if (mDragPt != nullptr)
        mConstrainedPt.pop_back();
}

void cBaraffClothGpu::InitELLMatrix()
{
    // 1. for each vertex, get the affect vertex
    int num_of_v = GetNumOfVertices();
    int num_of_e = GetNumOfEdges();
    int num_of_tri = GetNumOfTriangles();
    std::vector<std::set<int>> affected_vertex_array(num_of_v);

    for (int e_id = 0; e_id < num_of_e; e_id++)
    {
        if (mEdgeArrayShared[e_id]->mIsBoundary == true)
            continue;
        tVector4i cur_edge = mEdgeAffectVertexId[e_id];
        for (int i = 0; i < 4; i++)
        {
            int cur_v = cur_edge[i];
            for (int j = 0; j < 4; j++)
                affected_vertex_array[cur_v].insert(cur_edge[j]);
        }
    }

    for (int i = 0; i < num_of_tri; i++)
    {
        auto tri = mTriangleArrayShared[i];
        int v_id[3] = {tri->mId0, tri->mId1, tri->mId2};
        for (int j = 0; j < 3; j++)
        {
            int cur_v = v_id[j];
            affected_vertex_array[cur_v].insert(v_id[j]);
        }
    }

    // 2. set ELL local to global map
    std::vector<tELLLocal2GlobalInfo> ell_map_cpu(
        num_of_v, -1 * tELLLocal2GlobalInfo::Ones());
    OMP_PARALLEL_FOR
    for (int i = 0; i < num_of_v; i++)
    {
        int j = 0;
        for (std::set<int>::iterator it = affected_vertex_array[i].begin();
             it != affected_vertex_array[i].end(); it++, j++)
        {
            ell_map_cpu[i][j] = *it;
        }
        // std::cout << "v " << i << " ELL = " << ell_map_cpu[i].transpose()
        //           << std::endl;
    }
    // exit(1);
    mELLVidToGlobalVid.Upload(ell_map_cpu);

    // 3. allocate matrix
    mStiffnessMatrixCuda.Resize(num_of_v, tELLLocal2GlobalInfo::mElements);
    mSystemMatrixCuda.Resize(num_of_v, tELLLocal2GlobalInfo::mElements);
    mSystemRHSCuda.Resize(num_of_v);
    mIntForceCuda.Resize(num_of_v);
    std::cout << "mIntForceCuda size = " << mIntForceCuda.Size() << std::endl;
    mSolutionCuda.Resize(num_of_v);
}
extern void CalcBaraffStretchCoef(const tMatrix2d &mat, tVector3d &coef_u,
                                  tVector3d &coef_v);

void CalcBaraffStretchCoef(const tMatrix2d &mat, tCudaVector3f &coef_u,
                           tCudaVector3f &coef_v)
{
    tVector3d u, v;
    CalcBaraffStretchCoef(mat, u, v);

    coef_u =
        cCudaMatrixUtil::EigenMatrixToCudaMatrix(tVector3f(u.cast<float>()));

    coef_v =
        cCudaMatrixUtil::EigenMatrixToCudaMatrix(tVector3f(v.cast<float>()));
}

void cBaraffClothGpu::InitStretchCoef()
{
    auto tri_array = GetTriangleArray();
    auto v_array = GetVertexArray();
    int num_of_tri = tri_array.size();
    tMatrix2d Dminv;
    tMatrix2d R45 = cRotUtil::RotMat2D(M_PI / 4);
    std::vector<tCudaVector3f> cpu_coef_Fu_warp_weft(num_of_tri),
        cpu_coef_Fv_warp_weft(num_of_tri);
    std::vector<tCudaVector3f> cpu_coef_Fu_diag(num_of_tri),
        cpu_coef_Fv_diag(num_of_tri);
    for (int i = 0; i < num_of_tri; i++)
    {
        /*
        Dm^{-1} * [R45.T] = [
            a, b
            c, d
        ]
        */
        auto tri = tri_array[i];
        auto v0 = v_array[tri->mId0], v1 = v_array[tri->mId1],
             v2 = v_array[tri->mId2];
        Dminv.col(0) = (v1->muv - v0->muv).cast<double>();
        Dminv.col(1) = (v2->muv - v0->muv).cast<double>();
        Dminv = Dminv.inverse().eval();

        CalcBaraffStretchCoef(Dminv, cpu_coef_Fu_warp_weft[i],
                              cpu_coef_Fv_warp_weft[i]);
        Dminv = Dminv * R45.transpose();
        CalcBaraffStretchCoef(Dminv, cpu_coef_Fu_diag[i], cpu_coef_Fv_diag[i]);
    }
    mCoefFu_warp_weft.Upload(cpu_coef_Fu_warp_weft);
    mCoefFv_warp_weft.Upload(cpu_coef_Fv_warp_weft);

    mCoefFu_diag.Upload(cpu_coef_Fu_diag);
    mCoefFv_diag.Upload(cpu_coef_Fv_diag);
}

void cBaraffClothGpu::InitGravity()
{
    int num_of_v = GetNumOfVertices();
    std::vector<tCudaVector3f> g(num_of_v);
    for (int i = 0; i < num_of_v; i++)
    {
        g[i][1] = -9.8 * mVertexArrayShared[i]->mMass;
        // g[i][1] = 0 * mVertexArrayShared[i]->mMass;
    }
    mGravityForceCuda.Upload(g);
}