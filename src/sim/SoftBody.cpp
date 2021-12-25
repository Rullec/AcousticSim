#include "SoftBody.h"
#include "utils/JsonUtil.h"
#include "utils/TetUtil.h"
#include "geometries/Tetrahedron.h"
#include <iostream>
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
    // exit(1);
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