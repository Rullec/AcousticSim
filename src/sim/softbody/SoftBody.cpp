#include "SoftBody.h"
#include "utils/JsonUtil.h"
#include "utils/TetUtil.h"
#include "geometries/Tetrahedron.h"
#include "geometries/Primitives.h"
#include <iostream>
extern const tVector gGravity;

std::string gMaterialModelTypeStr[eMaterialModelType::NUM_OF_MATERIAL_MODEL] =
    {
        "LINEAR_ELASTICITY",
        "COROTATED",
        "FIX_COROTATED",
        "STVK",
        "NEO_HOOKEAN"};
std::string BuildMaterialTypeStr(eMaterialModelType type)
{
    return gMaterialModelTypeStr[type];
}

cSoftBody::cSoftBody(int id) : cBaseObject(eObjectType::SOFTBODY_TYPE, id)
{
}

cSoftBody::~cSoftBody()
{
}

/**
 * \brief           create and allocate the softbody data
*/
void cSoftBody::Init(const Json::Value &conf)
{
    cBaseObject::Init(conf);

    mTetMeshPath = cJsonUtil::ParseAsString("tet_path", conf);
    mRho = cJsonUtil::ParseAsFloat("rho", conf);
    cTetUtil::LoadTet(mTetMeshPath,
                      mVertexArrayShared,
                      mEdgeArrayShared,
                      mTriangleArrayShared,
                      mTetArrayShared);
    mMaterial = eMaterialModelType::STVK;

    std::cout << "load from path " << mTetMeshPath << "done\n";
    UpdateTriangleNormal();
    UpdateVertexNormalFromTriangleNormal();

    InitInvDm();
    int num_of_v = this->mVertexArrayShared.size();
    mIntForce.noalias() = tVectorXd::Zero(3 * num_of_v);
    mExtForce.noalias() = tVectorXd::Zero(3 * num_of_v);
    mUserForce.noalias() = tVectorXd::Zero(3 * num_of_v);
    mGravityForce = tVectorXd::Zero(3 * num_of_v);
    InitPos();
    InitTetVolume();
    InitDiagLumpedMassMatrix();

    // mXcur[3 * 2 + 1] += 0.1;
    mXprev = mXcur;
    SyncPosToVectorArray();
    for (int i = 0; i < num_of_v; i++)
    {
        mGravityForce.segment(3 * i, 3) = gGravity.segment(0, 3).cast<double>() / mInvLumpedMassMatrixDiag[3 * i];
    }

    // exit(1);
}
void cSoftBody::InitPos()
{
    int num_of_v = this->mVertexArrayShared.size();
    mXcur.resize(3 * num_of_v);
    for (int i = 0; i < num_of_v; i++)
    {
        mXcur.segment(3 * i, 3) = mVertexArrayShared[i]->mPos.segment(0, 3).cast<double>();
    }
    mXprev.noalias() = mXcur;
}

void cSoftBody::SetVerticesPos(const tVectorXd &pos)
{
    int num_of_v = this->mVertexArrayShared.size();
    SIM_ASSERT(pos.size() == 3 * num_of_v);
    mXcur.noalias() = pos;
    mXprev.noalias() = pos;
    SyncPosToVectorArray();
}

void cSoftBody::SyncPosToVectorArray()
{
    int num_of_v = this->mVertexArrayShared.size();
    for (int i = 0; i < num_of_v; i++)
    {
        mVertexArrayShared[i]->mPos.segment(0, 3) =
            mXcur.segment(3 * i, 3).cast<double>();
    }
}
extern void CalcTriangleDrawBufferSingle(tVertexPtr v0, tVertexPtr v1, tVertexPtr v2,
                                         Eigen::Map<tVectorXf> &buffer, int &st_pos);
void cSoftBody::CalcTriangleDrawBuffer(Eigen::Map<tVectorXf> &res,
                                       int &st) const
{
    // std::cout << "caclulate triangle draw buffer " << mTriangleArrayShared.size() << std::endl;
    for (auto &x : mTriangleArrayShared)
    {
        CalcTriangleDrawBufferSingle(mVertexArrayShared[x->mId0],
                                     mVertexArrayShared[x->mId1],
                                     mVertexArrayShared[x->mId2], res, st);
    }
}
extern void CalcEdgeDrawBufferSingle(tVertexPtr v0, tVertexPtr v1,
                                     const tVector &edge_normal,
                                     Eigen::Map<tVectorXf> &buffer, int &st_pos,
                                     const tVector &color);
void cSoftBody::CalcEdgeDrawBuffer(Eigen::Map<tVectorXf> &res,
                                   int &st) const
{
    tVector normal = tVector::Zero();
    tVector black_color = tVector(0, 0, 0, 1);
    for (auto &e : mEdgeArrayShared)
    {
        // 1. get the normal of this edge
        normal = mTriangleArrayShared[e->mTriangleId0]->mNormal;
        if (e->mTriangleId1 != -1)
        {
            normal += mTriangleArrayShared[e->mTriangleId1]->mNormal;
            normal /= 2;
        }

        CalcEdgeDrawBufferSingle(mVertexArrayShared[e->mId0], mVertexArrayShared[e->mId1],
                                 normal, res, st, black_color);
    }
}

/**
 * \brief           Calcualte per triangle normal direction
*/
void cSoftBody::UpdateTriangleNormal()
{
    for (int i = 0; i < this->mTriangleArrayShared.size(); i++)
    {
        auto cur_t = mTriangleArrayShared[i];
        tVector e01 = mVertexArrayShared[cur_t->mId1]->mPos - mVertexArrayShared[cur_t->mId0]->mPos;

        tVector e12 = mVertexArrayShared[cur_t->mId2]->mPos - mVertexArrayShared[cur_t->mId1]->mPos;
        cur_t->mNormal = e01.cross3(e12).normalized();
    }
}
void cSoftBody::UpdateVertexNormalFromTriangleNormal()
{
    // 1. clear all vertex normal
    // cTimeUtil::Begin("update_v_normal");
    for (auto &x : mVertexArrayShared)
        x->mNormal.setZero();
    // 2. iter each edge
    for (auto &x : mTriangleArrayShared)
    {
        mVertexArrayShared[x->mId0]->mNormal += x->mNormal;
        mVertexArrayShared[x->mId1]->mNormal += x->mNormal;
        mVertexArrayShared[x->mId2]->mNormal += x->mNormal;
    }

    // 3. averge each vertex
    for (int i = 0; i < mVertexArrayShared.size(); i++)
    {
        auto &v = mVertexArrayShared[i];
        v->mNormal.normalize();
    }
    // cTimeUtil::End("update_v_normal");
}
int cSoftBody::GetNumOfTriangles() const
{
    return this->mTriangleArrayShared.size();
}
int cSoftBody::GetNumOfEdges() const
{
    return this->mEdgeArrayShared.size();
}
int cSoftBody::GetNumOfVertices() const
{
    return this->mVertexArrayShared.size();
}

/**
 * \brief            update the result
*/
void cSoftBody::Update(float dt)
{
    dt = 5e-3;
    UpdateDeformationGradient();
    UpdateIntForce();
    UpdateExtForce();
    SolveForNextPos(dt);
}

void cSoftBody::UpdateExtForce()
{
    int num_of_v = this->mVertexArrayShared.size();
    double ground_height = 1e-2;
    double k = 1e3;
    // mExtForce.fill(5);
    // mExtForce[3 * 1 + 1] = 10;
    for (int i = 0; i < mVertexArrayShared.size(); i++)
    {
        double dist = mVertexArrayShared[i]->mPos[1] - ground_height;
        if (dist < 0)
        {
            mExtForce[3 * i + 1] = -dist * k;
        }
    }
}

void cSoftBody::SolveForNextPos(float dt)
{
    // AX = b; A(system matrix), b(residual)
    // std::cout << "mGravityForce = " << mGravityForce.transpose() << std::endl;
    tVectorXd accel = mInvLumpedMassMatrixDiag.cwiseProduct(mIntForce + mExtForce + mGravityForce + mUserForce) * dt * dt;
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
        tVector3d X1 = this->mVertexArrayShared[cur_tet->mVertexId[0]]->mPos.segment(0, 3).cast<double>();
        tVector3d X2 = this->mVertexArrayShared[cur_tet->mVertexId[1]]->mPos.segment(0, 3).cast<double>();
        tVector3d X3 = this->mVertexArrayShared[cur_tet->mVertexId[2]]->mPos.segment(0, 3).cast<double>();
        tVector3d X4 = this->mVertexArrayShared[cur_tet->mVertexId[3]]->mPos.segment(0, 3).cast<double>();
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
        // 1. calculate Ds (please check the siggraph 2012 note for more details)
        UpdateDeformationGradientForTet(i);
        // std::cout << "calc tet " << i << " F = \n " << mF[i] << std::endl;
    }
}

void cSoftBody::UpdateDeformationGradientForTet(int i)
{
    auto cur_tet = mTetArrayShared[i];
    tVector3d x1 = this->mVertexArrayShared[cur_tet->mVertexId[0]]->mPos.segment(0, 3).cast<double>();
    tVector3d x2 = this->mVertexArrayShared[cur_tet->mVertexId[1]]->mPos.segment(0, 3).cast<double>();
    tVector3d x3 = this->mVertexArrayShared[cur_tet->mVertexId[2]]->mPos.segment(0, 3).cast<double>();
    tVector3d x4 = this->mVertexArrayShared[cur_tet->mVertexId[3]]->mPos.segment(0, 3).cast<double>();
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
        mInitTetVolume[i] = cTetUtil::CalculateTetVolume(
            mVertexArrayShared[tet->mVertexId[0]]->mPos,
            mVertexArrayShared[tet->mVertexId[1]]->mPos,
            mVertexArrayShared[tet->mVertexId[2]]->mPos,
            mVertexArrayShared[tet->mVertexId[3]]->mPos);
    }
}

double cSoftBody::CalcEnergy()
{
    return 0;
}

/**
 * \brief           the "consistent" mass matrix, as defined in FEM, is a real symmetric matrix shaped as 3Nx3N. Its derivation is mostly based on the principle of virtual work. For more details please check "The Finite Element Process" P165 by Klaus-Jurgen Bathe. (and my own note) 
 * 
 *      In practice (for a better convergence and simplicity), we use "lumped" mass matrix. The most used lumped scheme is the so-called "row-diagnozation-sum" lump, which reduce all weight in consistent mass matrix to its diagnoal.
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
        double lump_value_row_sum = ele_mass_template.row(0).sum(); // lump_value_row_sum = 5
        double lump_value = lump_value_row_sum * tet_volume;
        // 3. add the component
        for (int j = 0; j < 4; j++)
        {
            DLMM_diag[cur_tet->mVertexId[j]] += lump_value;
        }
    }
    // std::cout << "rho = " << mRho << std::endl;
    DLMM_diag = DLMM_diag.eval() / 20 * this->mRho;
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
    // std::cout << "mInvLumpedMassMatrixDiag = " << mInvLumpedMassMatrixDiag.transpose();
    SIM_ASSERT(mInvLumpedMassMatrixDiag.hasNaN() == false);
    // exit(1);
}

int cSoftBody::GetNumOfTets() const
{
    return this->mTetArrayShared.size();
}

int cSoftBody::GetNumOfFreedoms() const
{
    return 3 * GetNumOfVertices();
}
#include "sim/Perturb.h"
/**
 * \brief           apply perturb force for once
*/
void cSoftBody::ApplyUserPerturbForceOnce(tPerturb *pert)
{
    mUserForce.setZero();
    int tri_id = pert->mAffectedTriId;
    tVector3d bary = pert->mBarycentricCoords;
    tVector3d force = pert->GetPerturbForce().cast<double>().segment(0, 3) * 100;
    int v0 = mTriangleArrayShared[tri_id]->mId0;
    int v1 = mTriangleArrayShared[tri_id]->mId1;
    int v2 = mTriangleArrayShared[tri_id]->mId2;
    mUserForce.segment(3 * v0, 3) += force * bary[0];
    mUserForce.segment(3 * v1, 3) += force * bary[1];
    mUserForce.segment(3 * v2, 3) += force * bary[2];
    // std::cout << "[user] apply user force on soft body " << force.transpose() << std::endl;
}

/**
 * \brief           return material type
*/
eMaterialModelType cSoftBody::GetMaterial() const
{
    return this->mMaterial;
}

/**
 * \brief           update imgui
*/
#include "imgui.h"
void cSoftBody::UpdateImGUi()
{
    static int item_cur_idx = mMaterial;

    std::vector<const char *> items = {};
    for (int i = 0; i < eMaterialModelType::NUM_OF_MATERIAL_MODEL; i++)
    {
        items.push_back(gMaterialModelTypeStr[i].c_str());
    }

    ImGui::Combo("Material", &item_cur_idx, items.data(), items.size());

    mMaterial = static_cast<eMaterialModelType>(item_cur_idx);
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

        P.noalias() = CalcPK1(F);
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

/**
 * \brief       Given deformation gradient F, calculate P(F)
*/
tMatrix3d cSoftBody::CalcPK1(const tMatrix3d &F)
{
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
    //     P.noalias() = gMu * (F.transpose() + F - 2 * I) + gLambda * (F - I).trace() * I;
    //     break;
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
    //     SIM_ERROR("do not support material model {}", BuildMaterialTypeStr(mMaterial));
    //     exit(1);
    // }
    // break;
    // }
    // return P;

    tMatrix3d P = F * CalcPK1_last(F);
    return P;
}

tMatrix3d cSoftBody::CalcPK1_last(const tMatrix3d &F)
{
    tMatrix3d I = tMatrix3d::Identity();
    tMatrix3d E = 0.5 * (F.transpose() * F - I);
    return 2 * gMu * E + gLambda * E.trace() * I;
}