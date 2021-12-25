#include "SoftBody.h"
#include "utils/JsonUtil.h"
#include "utils/TetUtil.h"
#include "geometries/Tetrahedron.h"
#include "geometries/Primitives.h"
#include <iostream>
extern const tVector gGravity;
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

    cTetUtil::LoadTet(mTetMeshPath,
                      mVertexArrayShared,
                      mEdgeArrayShared,
                      mTriangleArrayShared,
                      mTetArrayShared);

    std::cout << "load from path " << mTetMeshPath << "done\n";
    UpdateTriangleNormal();
    UpdateVertexNormalFromTriangleNormal();

    InitInvDm();
    int num_of_v = this->mVertexArrayShared.size();
    mIntForce.noalias() = tVectorXf::Zero(3 * num_of_v);
    mExtForce.noalias() = tVectorXf::Zero(3 * num_of_v);
    mUserForce.noalias() = tVectorXf::Zero(3 * num_of_v);
    InitPos();
    mInvMassMatrixDiag.noalias() = tVectorXf::Ones(3 * num_of_v) / mVertexArrayShared[0]->mMass;
    mGravityForce = tVectorXf::Zero(3 * num_of_v);
    InitTetVolume();

    mXcur[3 * 2 + 1] += 0.1;
    mXprev = mXcur;
    SyncPosToVectorArray();
    for (int i = 0; i < num_of_v; i++)
    {
        mGravityForce.segment(3 * i, 3) = mVertexArrayShared[i]->mMass * gGravity.segment(0, 3).cast<float>();
    }

    // exit(1);
}
void cSoftBody::InitPos()
{
    int num_of_v = this->mVertexArrayShared.size();
    mXcur.resize(3 * num_of_v);
    for (int i = 0; i < num_of_v; i++)
    {
        mXcur.segment(3 * i, 3) = mVertexArrayShared[i]->mPos.segment(0, 3).cast<float>();
    }
    mXprev.noalias() = mXcur;
}

void cSoftBody::SetVerticesPos(const tVectorXf &pos)
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

float CalculateTetVolume(
    const tVector &pos0,
    const tVector &pos1,
    const tVector &pos2,
    const tVector &pos3)
{
    // 1/6 * (AB X AC) \cdot (AD)
    tVector AB = pos1 - pos0;
    tVector AC = pos2 - pos0;
    tVector AD = pos3 - pos0;
    return (AB.cross3(AC)).dot(AD) / 6.0;
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
    float ground_height = 1e-2;
    float k = 1e3;
    // mExtForce.fill(5);
    mExtForce[3 * 1 + 1] = 10;
    for (int i = 0; i < mVertexArrayShared.size(); i++)
    {
        float dist = mVertexArrayShared[i]->mPos[1] - ground_height;
        if (dist < 0)
        {
            mExtForce[3 * i + 1] = -dist * k;
        }
    }
}

void cSoftBody::SolveForNextPos(float dt)
{
    // AX = b; A(system matrix), b(residual)
    tVectorXf accel = mInvMassMatrixDiag.cwiseProduct(mIntForce + mExtForce + mGravityForce) * dt * dt;
    tVectorXf Xnew = accel + 2 * this->mXcur - mXprev;
    tVectorXf Xnew_new = accel + mXcur;

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
    mInvDm.resize(num_of_tet, tMatrix3f::Zero());
    for (int i = 0; i < num_of_tet; i++)
    {
        auto cur_tet = mTetArrayShared[i];
        tVector3f X1 = this->mVertexArrayShared[cur_tet->mVertexId[0]]->mPos.segment(0, 3).cast<float>();
        tVector3f X2 = this->mVertexArrayShared[cur_tet->mVertexId[1]]->mPos.segment(0, 3).cast<float>();
        tVector3f X3 = this->mVertexArrayShared[cur_tet->mVertexId[2]]->mPos.segment(0, 3).cast<float>();
        tVector3f X4 = this->mVertexArrayShared[cur_tet->mVertexId[3]]->mPos.segment(0, 3).cast<float>();
        tMatrix3f cur_dm;
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
    tMatrix3f Ds = tMatrix3f::Zero();

    for (int i = 0; i < num_of_tet; i++)
    {
        // 1. calculate Ds (please check the siggraph 2012 note for more details)
        auto cur_tet = mTetArrayShared[i];
        tVector3f x1 = this->mVertexArrayShared[cur_tet->mVertexId[0]]->mPos.segment(0, 3).cast<float>();
        tVector3f x2 = this->mVertexArrayShared[cur_tet->mVertexId[1]]->mPos.segment(0, 3).cast<float>();
        tVector3f x3 = this->mVertexArrayShared[cur_tet->mVertexId[2]]->mPos.segment(0, 3).cast<float>();
        tVector3f x4 = this->mVertexArrayShared[cur_tet->mVertexId[3]]->mPos.segment(0, 3).cast<float>();
        Ds.col(0) = x1 - x4;
        Ds.col(1) = x2 - x4;
        Ds.col(2) = x3 - x4;
        mF[i].noalias() = Ds * mInvDm[i];
        // std::cout << "calc tet " << i << " F = \n " << mF[i] << std::endl;
    }
}
void cSoftBody::InitTetVolume()
{
    mInitTetVolume.resize(mTetArrayShared.size());
    for (int i = 0; i < mTetArrayShared.size(); i++)
    {
        auto tet = mTetArrayShared[i];
        mInitTetVolume[i] = CalculateTetVolume(
            mVertexArrayShared[tet->mVertexId[0]]->mPos,
            mVertexArrayShared[tet->mVertexId[1]]->mPos,
            mVertexArrayShared[tet->mVertexId[2]]->mPos,
            mVertexArrayShared[tet->mVertexId[3]]->mPos);
    }
}

void cSoftBody::UpdateIntForce()
{
    // update internal force
    int num_of_tet = this->mTetArrayShared.size();
    // internal force H = - W P(F) D_m^{-T}, iter on each tet
    float mu = 1e5;
    float lambda = 0.5;
    tMatrix3f I = tMatrix3f::Identity();
    mIntForce.setZero();
    for (int i = 0; i < num_of_tet; i++)
    {
        auto tet = mTetArrayShared[i];
        // 1.1 get W: tet volume
        float W = mInitTetVolume[i];
        // 1.2 get deformation gradient F
        const tMatrix3f F = mF[i];
        // 1.3 get P(F) in linear elasticity
        /*
            P(F) = \mu *(FT +F -2I) + \lambda * tr(F - I) I 
        */
        // tMatrix3f P = mu * (F.transpose() + F - 2 * I) + lambda * (F - I).trace() * I;

        // stvk
        // {
        tMatrix3f E = 0.5 * (F.transpose() * F - I);
        tMatrix3f P = F * (2 * mu * E + lambda * E.trace() * I);
        // }
        // 1.4 calculate force on nodes
        tMatrix3f H = -W * P * mInvDm[i].transpose();
        // std::cout << "tet " << i << " H = \n"
        //           << H << std::endl;
        tVector3f f3 = -(H.col(0) + H.col(1) + H.col(2));
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

float cSoftBody::CalcEnergy()
{
    return 0;
}