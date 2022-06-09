#include "SoftBody.h"
#include "geometries/Primitives.h"
#include "geometries/Tetrahedron.h"
#include "imgui.h"
#include "sim/softbody/BaseMaterial.h"
#include "sim/softbody/MaterialBuilder.h"
#include "utils/JsonUtil.h"
#include "utils/RenderUtil.h"
#include "utils/RotUtil.h"
#include "utils/SparseUtil.h"
#include "utils/TetUtil.h"
#include <iostream>
#include <set>
// extern const tVector gGravity;

cSoftBody::cSoftBody(int id) : cBaseObject(eObjectType::SOFTBODY_TYPE, id) {}

cSoftBody::~cSoftBody() {}

void cSoftBody::InitTetTransform(const Json::Value &conf)
{
    mInitTranslation = cJsonUtil::ReadVectorJson(
                           cJsonUtil::ParseAsValue("init_translation", conf))
                           .segment(0, 3);
    mInitRotation = cJsonUtil::ReadVectorJson(
                        cJsonUtil::ParseAsValue("init_rotation", conf))
                        .segment(0, 3);
    tMatrix transform =
        cRotUtil::AxisAngleToRotmat(cMathUtil::Expand(mInitRotation, 0));
    transform.block(0, 3, 3, 1) = mInitTranslation.segment(0, 3);
    std::cout << "[debug] apply init rotation(AA) = "
              << mInitRotation.transpose() << std::endl;
    std::cout << "[debug] apply init translation = "
              << mInitTranslation.transpose() << std::endl;
    for (auto &v : mVertexArray)
    {
        v->mPos[3] = 1;
        v->mPos = transform * v->mPos;
    }
}
/**
 * \brief           create and allocate the softbody data
 */
void cSoftBody::Init(const Json::Value &conf)
{
    cBaseObject::Init(conf);

    mTetMeshPath = cJsonUtil::ParseAsString("tet_path", conf);

    cTetUtil::LoadTet(mTetMeshPath, mVertexArray, mEdgeArray, mTriangleArray,
                      mTetArrayShared);
    InitTetTransform(conf);

    {
        std::string mat_model_str =
            cJsonUtil::ParseAsString("material_type", conf);
        std::string mat_path = cJsonUtil::ParseAsString("material_path", conf);
        auto mat_type = BuildMaterialTypeFromStr(mat_model_str);
        mMaterial = BuildMaterial(mat_path, mat_type);
    }

    mRho = mMaterial->GetRho();

    mFrictionCoef = cJsonUtil::ParseAsFloat("friction", conf);
    mCollisionK = cJsonUtil::ParseAsFloat("collision_k", conf);

    std::cout << "load from path " << mTetMeshPath << " done\n";
    UpdateTriangleNormal();
    UpdateVertexNormalFromTriangleNormal();
    InitSurface();

    InitInvDm();
    InitPos();
    InitTetVolume();
    InitRawMassMatrix();
    InitDiagLumpedMassMatrix();
    InitForceVector();
    // mXcur[3 * mTetArrayShared[0]->mVertexId[0] + 1] += 0.02;
    mXprev = mXcur;
    mXInit = mXcur;
    SyncPosToVectorArray();
}
void cSoftBody::InitForceVector()
{
    int num_of_v = this->mVertexArray.size();
    mIntForce.noalias() = tVectorXd::Zero(3 * num_of_v);
    mExtForce.noalias() = tVectorXd::Zero(3 * num_of_v);
    mUserForce.noalias() = tVectorXd::Zero(3 * num_of_v);
    mGravityForce.noalias() = tVectorXd::Zero(3 * num_of_v);
    for (int i = 0; i < num_of_v; i++)
    {
        // mGravityForce.segment(3 * i, 3) = gGravity.segment(0,
        // 3).cast<double>() / mInvLumpedMassMatrixDiag[3 * i];
    }
}
void cSoftBody::InitPos()
{
    int num_of_v = this->mVertexArray.size();
    mXcur.resize(3 * num_of_v);
    for (int i = 0; i < num_of_v; i++)
    {
        mXcur.segment(3 * i, 3) =
            mVertexArray[i]->mPos.segment(0, 3).cast<double>();
    }
    mXprev.noalias() = mXcur;
}

void cSoftBody::SetVerticesPos(const tVectorXd &pos)
{
    int num_of_v = this->mVertexArray.size();
    SIM_ASSERT(pos.size() == 3 * num_of_v);
    mXcur.noalias() = pos;
    mXprev.noalias() = pos;
    SyncPosToVectorArray();
}

void cSoftBody::SyncPosToVectorArray()
{
    int num_of_v = this->mVertexArray.size();
    for (int i = 0; i < num_of_v; i++)
    {
        mVertexArray[i]->mPos.segment(0, 3) =
            mXcur.segment(3 * i, 3).cast<double>();
    }
}

void cSoftBody::CalcEdgeDrawBuffer(Eigen::Map<tVectorXf> &res, int &st) const
{
    tVector normal = tVector::Zero();
    tVector black_color = tVector(0, 0, 0, 1);
    for (auto &e : mEdgeArray)
    {
        // 1. get the normal of this edge
        normal = mTriangleArray[e->mTriangleId0]->mNormal;
        if (e->mTriangleId1 != -1)
        {
            normal += mTriangleArray[e->mTriangleId1]->mNormal;
            normal /= 2;
        }

        cRenderUtil::CalcEdgeDrawBufferSingle(mVertexArray[e->mId0],
                                              mVertexArray[e->mId1], normal,
                                              res, st, black_color);
    }
}

/**
 * \brief           Calcualte per triangle normal direction
 */
void cSoftBody::UpdateTriangleNormal()
{
    for (int i = 0; i < this->mTriangleArray.size(); i++)
    {
        auto cur_t = mTriangleArray[i];
        tVector e01 =
            mVertexArray[cur_t->mId1]->mPos - mVertexArray[cur_t->mId0]->mPos;

        tVector e12 =
            mVertexArray[cur_t->mId2]->mPos - mVertexArray[cur_t->mId1]->mPos;
        cur_t->mNormal = e01.cross3(e12).normalized();
    }
}
void cSoftBody::UpdateVertexNormalFromTriangleNormal()
{
    // 1. clear all vertex normal
    // cTimeUtil::Begin("update_v_normal");
    for (auto &x : mVertexArray)
        x->mNormal.setZero();
    // 2. iter each edge
    for (auto &x : mTriangleArray)
    {
        mVertexArray[x->mId0]->mNormal += x->mNormal;
        mVertexArray[x->mId1]->mNormal += x->mNormal;
        mVertexArray[x->mId2]->mNormal += x->mNormal;
    }

    // 3. averge each vertex
    for (int i = 0; i < mVertexArray.size(); i++)
    {
        auto &v = mVertexArray[i];
        v->mNormal.normalize();
    }
    // cTimeUtil::End("update_v_normal");
}

// only get the number of surface triangles
int cSoftBody::GetNumOfTriangles() const
{
    return mSurfaceTriangleIdArray.size();
}

void cSoftBody::CalcTriangleDrawBuffer(Eigen::Map<tVectorXf> &res,
                                       int &st) const
{
    int old_st = st;
    int num_of_surface_tris = mSurfaceTriangleIdArray.size();
    for (int idx = 0; idx < num_of_surface_tris; idx++)
    {
        int tri_id = mSurfaceTriangleIdArray[idx];
        auto tri = mTriangleArray[tri_id];
        cRenderUtil::CalcTriangleDrawBufferSingle(
            mVertexArray[tri->mId0], mVertexArray[tri->mId1],
            mVertexArray[tri->mId2], tri->mColor, res, st);
    }
}
int cSoftBody::GetNumOfEdges() const { return this->mEdgeArray.size(); }
int cSoftBody::GetNumOfVertices() const { return this->mVertexArray.size(); }

/**
 * \brief            update the result
 */
void cSoftBody::Update(float dt)
{
    UpdateDeformationGradient();
    UpdateExtForce();
    UpdateIntForce();
    SolveForNextPos(dt);
    mUserForce.setZero();
}

void cSoftBody::UpdateExtForce()
{
    int num_of_v = this->mVertexArray.size();
    double ground_height = 1e-3;
    double k = mCollisionK;
    float KinectFrictionCoef = mFrictionCoef;
    // mExtForce.fill(5);
    // mExtForce[3 * 1 + 1] = 10;
    for (int i = 0; i < mVertexArray.size(); i++)
    {
        double dist = mVertexArray[i]->mPos[1] - ground_height;

        if (dist < 0)
        {
            mVertexArray[i]->mPos[1] = 0;
            float normal_force_amp = -dist * k;

            tVector3d friction_dir =
                -1 * (mXcur.segment(3 * i, 3) - mXprev.segment(3 * i, 3));
            friction_dir[1] = 0;
            friction_dir.normalize();
            float friction_force_amp =
                KinectFrictionCoef * std::fabs(normal_force_amp);
            tVector3d force = tVector3d::Zero();
            force[1] = normal_force_amp;
            force[0] = friction_dir[0] * friction_force_amp;
            force[2] = friction_dir[2] * friction_force_amp;
            // std::cout << "[col] force " << force.transpose() << " on v" << i
            // << std::endl;
            mExtForce.segment(3 * i, 3) = force;
        }
    }
    // for (int i = 0; i < mXcur.size() / 3; i++)
    // {
    //     if (mXcur[3 * i + 1] < 0)
    //         mXcur[3 * i + 1] = 1e-3;
    // }
    // SyncPosToVectorArray();
}

void cSoftBody::SolveForNextPos(float dt)
{
    // AX = b; A(system matrix), b(residual)
    // std::cout << "mGravityForce = " << mGravityForce.transpose() <<
    // std::endl;
    tVectorXd accel = mInvLumpedMassMatrixDiag.cwiseProduct(
                          mIntForce + mExtForce + mGravityForce + mUserForce) *
                      dt * dt;
    tVectorXd Xnew = accel + 2 * this->mXcur - mXprev;
    tVectorXd Xnew_new = accel + mXcur;

    mXprev = mXcur;
    mXcur = Xnew_new;
    SyncPosToVectorArray();
}

/**
 * \brief           create Dminv
 */
void cSoftBody::InitInvDm()
{
    int num_of_tet = this->mTetArrayShared.size();
    mInvDm.resize(num_of_tet, tMatrix3d::Zero());
    for (int i = 0; i < num_of_tet; i++)
    {
        auto cur_tet = mTetArrayShared[i];
        tVector3d X1 = this->mVertexArray[cur_tet->mVertexId[0]]
                           ->mPos.segment(0, 3)
                           .cast<double>();
        tVector3d X2 = this->mVertexArray[cur_tet->mVertexId[1]]
                           ->mPos.segment(0, 3)
                           .cast<double>();
        tVector3d X3 = this->mVertexArray[cur_tet->mVertexId[2]]
                           ->mPos.segment(0, 3)
                           .cast<double>();
        tVector3d X4 = this->mVertexArray[cur_tet->mVertexId[3]]
                           ->mPos.segment(0, 3)
                           .cast<double>();
        tMatrix3d cur_dm;
        cur_dm.col(0) = X1 - X4;
        cur_dm.col(1) = X2 - X4;
        cur_dm.col(2) = X3 - X4;
        mInvDm[i] = cur_dm.inverse();
        // std::cout << "tet " << i << " Dminv = " << mInvDm[i] << std::endl;
    }
}

void cSoftBody::UpdateDeformationGradient()
{
    int num_of_tet = this->mTetArrayShared.size();

    if (mF.size() != num_of_tet)
    {
        mF.resize(num_of_tet);
    }

    for (int i = 0; i < num_of_tet; i++)
    {
        // 1. calculate Ds (please check the siggraph 2012 note for more
        // details)
        UpdateDeformationGradientForTet(i);
        // std::cout << "calc tet " << i << " F = \n " << mF[i] << std::endl;
    }
}

void cSoftBody::UpdateDeformationGradientForTet(int i)
{
    auto cur_tet = mTetArrayShared[i];
    tVector3d x1 = this->mVertexArray[cur_tet->mVertexId[0]]
                       ->mPos.segment(0, 3)
                       .cast<double>();
    tVector3d x2 = this->mVertexArray[cur_tet->mVertexId[1]]
                       ->mPos.segment(0, 3)
                       .cast<double>();
    tVector3d x3 = this->mVertexArray[cur_tet->mVertexId[2]]
                       ->mPos.segment(0, 3)
                       .cast<double>();
    tVector3d x4 = this->mVertexArray[cur_tet->mVertexId[3]]
                       ->mPos.segment(0, 3)
                       .cast<double>();
    tMatrix3d Ds = tMatrix3d::Zero();
    Ds.col(0) = x1 - x4;
    Ds.col(1) = x2 - x4;
    Ds.col(2) = x3 - x4;
    mF[i].noalias() = Ds * mInvDm[i];
}

void cSoftBody::InitTetVolume()
{
    mInitTetVolume.resize(mTetArrayShared.size());
    for (int i = 0; i < mTetArrayShared.size(); i++)
    {
        auto tet = mTetArrayShared[i];
        mInitTetVolume[i] =
            cTetUtil::CalculateTetVolume(mVertexArray[tet->mVertexId[0]]->mPos,
                                         mVertexArray[tet->mVertexId[1]]->mPos,
                                         mVertexArray[tet->mVertexId[2]]->mPos,
                                         mVertexArray[tet->mVertexId[3]]->mPos);
    }
}

double cSoftBody::CalcEnergy() { return 0; }

/**
 * \brief           the "consistent" mass matrix, as defined in FEM, is a real
 * symmetric matrix shaped as 3Nx3N. Its derivation is mostly based on the
 * principle of virtual work. For more details please check "The Finite Element
 * Process" P165 by Klaus-Jurgen Bathe. (and my own note)
 *
 *      In practice (for a better convergence and simplicity), we use "lumped"
 * mass matrix. The most used lumped scheme is the so-called
 * "row-diagnozation-sum" lump, which reduce all weight in consistent mass
 * matrix to its diagnoal.
 *
 *      Diagonalization Lumped Mass Matrix(DLMM)
 *
 *
 *  M_{total} = rho / 20 * [ \sum_e
 *       V_e * S_e^T *
 *          (2 1 1 1)
 *          (1 2 1 1)
 *          (1 1 2 1)
 *          (1 1 1 2) * S_e
 * ]
 */
void cSoftBody::InitDiagLumpedMassMatrix()
{
    // calcualte sum of volume
    std::cout << "total init volume = " << mInitTetVolume.sum() << std::endl;

    int dof = GetNumOfFreedoms();
    int num_of_v = GetNumOfVertices();
    tVectorXd DLMM_diag = tVectorXd::Zero(num_of_v);

    tMatrix4f ele_mass_template = tMatrix4f::Zero();
    ele_mass_template << 2, 1, 1, 1, 1, 2, 1, 1, 1, 1, 2, 1, 1, 1, 1, 2;
    // std::cout << "ele_mass_template = " << ele_mass_template << std::endl;
    // exit(1);
    for (size_t i = 0; i < GetNumOfTets(); i++)
    {
        auto cur_tet = mTetArrayShared[i];
        // 1. shape the volume
        double tet_volume = mInitTetVolume[i];
        // 2. shape the selection matrix
        double lump_value_row_sum =
            ele_mass_template.row(0).sum(); // lump_value_row_sum = 5
        double lump_value = lump_value_row_sum * tet_volume;
        // 3. add the component
        for (int j = 0; j < 4; j++)
        {
            DLMM_diag[cur_tet->mVertexId[j]] += lump_value;
        }
    }
    // std::cout << "rho = " << mRho << std::endl;
    DLMM_diag = DLMM_diag.eval() / 20 * this->mRho;
    std::cout << "total init mass = " << DLMM_diag.sum() << std::endl;

    // std::cout << "DLMM = " << DLMM_diag.transpose();
    // mInvLumpedMassMatrixDiag = DLMM_diag.cwiseInverse();
    mInvLumpedMassMatrixDiag.noalias() = tVectorXd::Zero(dof);
    for (int i = 0; i < num_of_v; i++)
    {
        double val = 1.0 / DLMM_diag[i];
        mInvLumpedMassMatrixDiag[3 * i + 0] = val;
        mInvLumpedMassMatrixDiag[3 * i + 1] = val;
        mInvLumpedMassMatrixDiag[3 * i + 2] = val;
    }
    std::cout << "new total mass = "
              << mInvLumpedMassMatrixDiag.cwiseInverse().sum() << std::endl;
    // exit(1);
    // std::cout << "mInvLumpedMassMatrixDiag = " <<
    // mInvLumpedMassMatrixDiag.transpose();
    SIM_ASSERT(mInvLumpedMassMatrixDiag.hasNaN() == false);
    // exit(1);
}

void cSoftBody::InitRawMassMatrix()
{
    // calcualte sum of volume

    int dof = GetNumOfFreedoms();
    int num_of_v = GetNumOfVertices();
    tVectorXd DLMM_diag = tVectorXd::Zero(num_of_v);

    tMatrix4f ele_mass_template = tMatrix4f::Zero();
    ele_mass_template << 2, 1, 1, 1, 1, 2, 1, 1, 1, 1, 2, 1, 1, 1, 1, 2;

    ele_mass_template *= mRho / 20;

    /*
     *  M_{total} = rho / 20 * [ \sum_e
     *       V_e * S_e^T *
     *          (2 1 1 1)
     *          (1 2 1 1)
     *          (1 1 2 1)
     *          (1 1 1 2) * S_e
     */
    std::vector<tTriplet> tripletList = {};
    for (size_t i = 0; i < GetNumOfTets(); i++)
    {
        auto cur_tet = mTetArrayShared[i];
        // 1. shape the volume
        double tet_volume = mInitTetVolume[i];

        tMatrix4f tet_ele_mass = tet_volume * ele_mass_template;

        for (int a = 0; a < 4; a++)
            for (int b = 0; b < 4; b++)
            {
                int v0 = cur_tet->mVertexId[a];
                int v1 = cur_tet->mVertexId[b];
                double mass = tet_ele_mass(a, b);
                tripletList.push_back(tTriplet(3 * v0 + 0, 3 * v1 + 0, mass));
                tripletList.push_back(tTriplet(3 * v0 + 1, 3 * v1 + 1, mass));
                tripletList.push_back(tTriplet(3 * v0 + 2, 3 * v1 + 2, mass));
            }
        // 3. add the component
    }

    mRawMassMatrix.setZero();
    mRawMassMatrix.resize(3 * num_of_v, 3 * num_of_v);
    mRawMassMatrix.setFromTriplets(tripletList.begin(), tripletList.end());
    std::cout << "new total mass is " << mRawMassMatrix.sum() << std::endl;
}

int cSoftBody::GetNumOfTets() const { return this->mTetArrayShared.size(); }

int cSoftBody::GetNumOfFreedoms() const { return 3 * GetNumOfVertices(); }
#include "sim/Perturb.h"
/**
 * \brief           apply perturb force for once
 */
void cSoftBody::ApplyUserPerturbForceOnce(tPerturb *pert)
{
    mUserForce.setZero();
    int tri_id = pert->mAffectedTriId;
    tVector3d bary = pert->mBarycentricCoords;
    tVector3d force = pert->GetPerturbForce().cast<double>().segment(0, 3) * 10;
    int v0 = mTriangleArray[tri_id]->mId0;
    int v1 = mTriangleArray[tri_id]->mId1;
    int v2 = mTriangleArray[tri_id]->mId2;
    mUserForce.segment(3 * v0, 3) += force * bary[0];
    mUserForce.segment(3 * v1, 3) += force * bary[1];
    mUserForce.segment(3 * v2, 3) += force * bary[2];
    std::cout << "[user] apply user force on soft body " << force.transpose()
              << std::endl;
}

/**
 * \brief           return material type
 */
cBaseMaterialPtr cSoftBody::GetMaterial() const { return this->mMaterial; }

eMaterialType cSoftBody::GetMaterialType() const
{
    return this->mMaterial->GetType();
}
/**
 * \brief           update imgui
 */
void cSoftBody::UpdateImGui()
{

    static int item_cur_idx = GetMaterial()->GetType();
    auto old_idx = item_cur_idx;
    std::vector<std::string> items_sto = {};
    std::vector<const char *> items = {};
    for (int i = 0; i < eMaterialType::NUM_OF_MATERIAL_MODEL; i++)
    {
        // auto res = BuildMaterialTypeStr(static_cast<eMaterialType>(i));
        // std::cout
        //     << res << std::endl;
        // items_sto.push_back(res);
        items.push_back(gMaterialModelTypeStr[i].c_str());
    }
    ImGui::Combo("Material", &item_cur_idx, items.data(), items.size());
    ImGui::SliderFloat("dampingA", &mMaterial->mRayleighDamplingA, 0, 10);
    ImGui::SliderFloat("dampingB", &mMaterial->mRayleighDamplingB, 0, 1);
    ImGui::SliderFloat("friction", &mFrictionCoef, 0, 1);
    ImGui::SliderFloat("collisionK", &mCollisionK, 0, 2e3);
    // ImGui::SliderFloat("Mu", &gMu, 0.0, 1e5);
    if (old_idx != item_cur_idx)
    {
        // rebuild the material

        mMaterial = BuildMaterial(mMaterial->mMatPath,
                                  static_cast<eMaterialType>(item_cur_idx));
    }
}

// #include "sim/softbody/SoftBody.h"
// #include "geometries/Tetrahedron.h"
// #include "utils/LogUtil.h"

void cSoftBody::UpdateIntForce()
{
    // update internal force
    int num_of_tet = this->mTetArrayShared.size();
    // internal force H = - W P(F) D_m^{-T}, iter on each tet

    mIntForce.setZero();
    tMatrix3d P = tMatrix3d::Zero();
    for (int i = 0; i < num_of_tet; i++)
    {
        auto tet = mTetArrayShared[i];
        // 1.1 get W: tet volume
        double W = mInitTetVolume[i];
        // 1.2 get deformation gradient F
        const tMatrix3d &F = mF[i];
        // 1.3 get P(F) in linear elasticity

        P.noalias() = mMaterial->CalcP(F);
        // 1.4 calculate force on nodes
        tMatrix3d H = -W * P * mInvDm[i].transpose();
        // std::cout << "tet " << i << " H = \n"
        //           << H << std::endl;
        tVector3d f3 = -(H.col(0) + H.col(1) + H.col(2));
        // std::cout << "f0 = " << H.col(0).transpose() << std::endl;
        // std::cout << "f1 = " << H.col(1).transpose() << std::endl;
        // std::cout << "f2 = " << H.col(2).transpose() << std::endl;
        // std::cout << "f3 = " << f3.transpose() << std::endl;
        mIntForce.segment(tet->mVertexId[0] * 3, 3) += H.col(0);
        mIntForce.segment(tet->mVertexId[1] * 3, 3) += H.col(1);
        mIntForce.segment(tet->mVertexId[2] * 3, 3) += H.col(2);
        mIntForce.segment(tet->mVertexId[3] * 3, 3) += f3;
    }
    // std::cout << "fint = " << mIntForce.transpose() << std::endl;
}

tVectorXd cSoftBody::CalcTetIntForce(size_t tet_id)
{
    auto tet = mTetArrayShared[tet_id];
    // 1.1 get W: tet volume
    double W = mInitTetVolume[tet_id];
    // 1.2 get deformation gradient F
    const tMatrix3d &F = mF[tet_id];
    // 1.3 get P(F) in linear elasticity

    tMatrix3d P = mMaterial->CalcP(F);
    // 1.4 calculate force on nodes
    tMatrix3d H = -W * P * mInvDm[tet_id].transpose();
    // std::cout << "tet " << i << " H = \n"
    //           << H << std::endl;
    tVector3d f3 = -(H.col(0) + H.col(1) + H.col(2));
    tVectorXd tet_int_force = tVectorXd::Zero(12);
    tet_int_force.segment(0, 3) = H.col(0);
    tet_int_force.segment(1, 3) = H.col(1);
    tet_int_force.segment(2, 3) = H.col(2);
    tet_int_force.segment(3, 3) = f3;
    return tet_int_force;
}

/**
 * \brief           two constant selection matrices used in stiffness matrix
 *  For more details please check the note "FEM course 第三部分
 * 离散化(刚度矩阵计算)"
 */
void GetSelectionMatrix(tMatrixXd &Sd, tMatrixXd &Sb)
{
    {
        Sb.resize(3, 1);
        Sb << 1, 1, 1;
    }
    {
        tMatrixXd Sa = tMatrixXd::Zero(9, 3), Sc = tMatrixXd::Zero(12, 9);
        Sa.block(0, 0, 3, 1).setOnes();
        Sa.block(3, 1, 3, 1).setOnes();
        Sa.block(6, 2, 3, 1).setOnes();

        Sc.block(0, 0, 3, 3).setIdentity();
        Sc.block(3, 3, 3, 3).setIdentity();
        Sc.block(6, 6, 3, 3).setIdentity();

        tMatrix3d minus_I = -tMatrix3d::Identity();
        Sc.block(9, 0, 3, 3).noalias() = minus_I;
        Sc.block(9, 3, 3, 3).noalias() = minus_I;
        Sc.block(9, 6, 3, 3).noalias() = minus_I;
        Sd.noalias() = Sc * Sa;
    }
}

tVectorXd cSoftBody::CalcTetIntForceBySelectionMatrix(size_t tet_id)
{
    auto tet = mTetArrayShared[tet_id];
    // 1.1 get W: tet volume
    double W = mInitTetVolume[tet_id];
    // 1.2 get deformation gradient F
    const tMatrix3d &F = mF[tet_id];
    // 1.3 get P(F) in linear elasticity

    tMatrix3d P = mMaterial->CalcP(F);
    tMatrixXd Sd, Sb;
    GetSelectionMatrix(Sd, Sb);
    tVectorXd int_force = -W * Sd * P * mInvDm[tet_id].transpose() * Sb;
    return int_force;
}

/**
 * \brief       Given deformation gradient F, calculate P(F)
 */
// tMatrix3d CalcPK1(const tMatrix3d &F)
// {
// tMatrix3d P = tMatrix3d::Zero();
// tMatrix3d I = tMatrix3d::Identity();

// switch (mMaterial)
// {
// // case eMaterialModelType::COROTATED:
// // {
// //     P.setZero();
// //     break;
// // }
// case eMaterialModelType::LINEAR_ELASTICITY:
// {
//     /*
//         P(F) = \mu * (FT + F -2I) + \lambda * tr(F - I) I
//     */
//     P.noalias() = gMu * (F.transpose() + F - 2 * I) + gLambda * (F -
//     I).trace() * I; break;
// }
// // case eMaterialModelType::FIX_COROTATED:
// // {
// //     P.setZero();
// //     break;
// // }
// case eMaterialModelType::STVK:
// {
//     tMatrix3d E = 0.5 * (F.transpose() * F - I);
//     P = F * (2 * gMu * E + gLambda * E.trace() * I);
//     break;
// }
// case eMaterialModelType::NEO_HOOKEAN:
// {
//     /*
//     P(F) = mu * F - mu * F^{-T} + lambda * log(I3) / 2 * F^{-T}
//     */
//     double I3 = 100;
//     tMatrix3d F_inv_T = F.inverse().transpose();
//     P.noalias() = gMu * (F - F_inv_T) + gLambda * std::log(I3) / 2 * F_inv_T;
//     break;
// }
// default:
// {
//     SIM_ERROR("do not support material model {}",
//     BuildMaterialTypeStr(mMaterial)); exit(1);
// }
// break;
// }
// return P;

//     tMatrix3d P = CalcPK1_part1(F) + CalcPK1_part2(F);
//     return P;
// }

// tMatrix3d CalcGreenStrain(const tMatrix3d &F)
// {
//     tMatrix3d I = tMatrix3d::Identity();
//     return 0.5 * (F.transpose() * F - I);
// }
// tMatrix3d CalcPK1_part1(const tMatrix3d &F)
// {
//     tMatrix3d E = CalcGreenStrain(F);
//     return 2 * gMu * F * E;
// }

// tMatrix3d CalcPK1_part2(const tMatrix3d &F)
// {
//     tMatrix3d E = CalcGreenStrain(F);
//     return gLambda * E.trace() * F;
// }

void cSoftBody::Reset()
{
    InitForceVector();
    mXcur = mXInit;
    mXprev = mXInit;
    SyncPosToVectorArray();
}

tVectorXd cSoftBody::GetTetForce(size_t tet_id, const tVectorXd &total_force)
{
    tVectorXd tet_force = tVectorXd::Zero(12);
    for (size_t i = 0; i < 4; i++)
    {
        tet_force.segment(3 * i, 3) =
            total_force.segment(3 * mTetArrayShared[tet_id]->mVertexId[i], 3);
    }
    return tet_force;
}

tVectorXd cSoftBody::GetTetVerticesPos(size_t tet_id,
                                       const tVectorXd &total_pos)
{
    return GetTetForce(tet_id, total_pos);
}

void cSoftBody::InitSurface()
{
    mSurfaceTriangleIdArray.clear();
    mSurfaceVertexIdArray.clear();
    mSurfaceEdgeIdArray.clear();
    int num_of_tris_total = this->mTriangleArray.size();
    std::vector<int> tri_involved_times(num_of_tris_total, 0);
    for (auto &t : this->mTetArrayShared)
    {
        for (int j = 0; j < 4; j++)
            tri_involved_times[t->mTriangleId[j]] += 1;
    }

    std::set<int> involved_vertices = {};
    std::set<int> involved_edges = {};
    for (int id = 0; id < num_of_tris_total; id++)
    {
        if (tri_involved_times[id] == 1)
        {
            auto tri = mTriangleArray[id];
            mSurfaceTriangleIdArray.push_back(id);
            involved_vertices.insert(tri->mId0);
            involved_vertices.insert(tri->mId1);
            involved_vertices.insert(tri->mId2);

            involved_edges.insert(tri->mEId0);
            involved_edges.insert(tri->mEId1);
            involved_edges.insert(tri->mEId2);
        }
    }
    mSurfaceVertexIdArray.resize(involved_vertices.size());
    std::copy(involved_vertices.begin(), involved_vertices.end(),
              mSurfaceVertexIdArray.begin());

    mSurfaceEdgeIdArray.resize(involved_edges.size());
    std::copy(involved_edges.begin(), involved_edges.end(),
              mSurfaceEdgeIdArray.begin());
    printf("[log] surface triangle %d, vertices %d edges %d\n",
           mSurfaceTriangleIdArray.size(), mSurfaceVertexIdArray.size(),
           mSurfaceEdgeIdArray.size());
}