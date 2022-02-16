#include "sim/BaseObject.h"
#include "geometries/Primitives.h"
#include "utils/JsonUtil.h"
#include "utils/LogUtil.h"
#include "utils/MathUtil.h"

#include <string>
std::string gObjectTypeStr[eObjectType::NUM_OBJ_TYPES] = {
    "KinematicBody", "RigidBody", "Cloth", "Fluid",
    "SoftBody", "Acoustic"};

cBaseObject::cBaseObject(eObjectType type, int id_) : mType(type), mObjId(id_)
{
    mObjName = "";
    mEnableDrawBuffer = true;
}

/**
 * \brief           Set object name
 */
void cBaseObject::SetObjName(std::string name) { mObjName = name; }
/**
 * \brief           Get object name
 */
std::string cBaseObject::GetObjName() const { return mObjName; }

cBaseObject::~cBaseObject() {}

eObjectType cBaseObject::BuildObjectType(std::string str)
{
    eObjectType type = eObjectType::INVALID_OBJ_TYPE;
    for (int i = 0; i < eObjectType::NUM_OBJ_TYPES; i++)
    {
        if (gObjectTypeStr[i] == str)
        {
            type = static_cast<eObjectType>(i);
            break;
        }
    }

    SIM_ASSERT(type != eObjectType::INVALID_OBJ_TYPE);
    return type;
}

eObjectType cBaseObject::GetObjectType() const { return this->mType; }

int cBaseObject::GetNumOfTriangles() const { return mTriangleArrayShared.size(); }
int cBaseObject::GetNumOfEdges() const { return mEdgeArrayShared.size(); }
int cBaseObject::GetNumOfVertices() const { return mVertexArrayShared.size(); }
const std::vector<tVertexPtr> &cBaseObject::GetVertexArray() const
{
    return this->mVertexArrayShared;
}
const std::vector<tEdgePtr> &cBaseObject::GetEdgeArray() const
{
    return this->mEdgeArrayShared;
}
const std::vector<tTrianglePtr> &cBaseObject::GetTriangleArray() const
{
    return this->mTriangleArrayShared;
}

std::vector<tVertexPtr> &cBaseObject::GetVertexArrayRef()
{
    return mVertexArrayShared;
}

std::vector<tEdgePtr> &cBaseObject::GetEdgeArrayRef() { return mEdgeArrayShared; }
std::vector<tTrianglePtr> &cBaseObject::GetTriangleArrayRef()
{
    return mTriangleArrayShared;
}

/**
 * \brief           change triangle color
 */
void cBaseObject::ChangeTriangleColor(int tri_id, const tVector3f &color)
{
    tVector new_color = tVector(color[0], color[1], color[2], mColorAlpha);
    mVertexArrayShared[mTriangleArrayShared[tri_id]->mId0]->mColor = new_color;
    mVertexArrayShared[mTriangleArrayShared[tri_id]->mId1]->mColor = new_color;
    mVertexArrayShared[mTriangleArrayShared[tri_id]->mId2]->mColor = new_color;
}

/**
 * \brief           Calculate axis aligned bounding box
 */
#include <cfloat>
void cBaseObject::CalcAABB(tVector &min, tVector &max) const
{
    min = tVector::Ones() * std::numeric_limits<double>::max();
    max = tVector::Ones() * std::numeric_limits<double>::max() * -1;
    for (auto &x : mVertexArrayShared)
    {
        for (int i = 0; i < 3; i++)
        {

            double val = x->mPos[i];
            min[i] = (val < min[i]) ? val : min[i];
            max[i] = (val > max[i]) ? val : max[i];
        }
    }
}

void cBaseObject::Init(const Json::Value &conf)
{
    mObjName = cJsonUtil::ParseAsString(OBJECT_NAME_KEY, conf);
}

/**
 * \brief           Update the normal of triangles
 */
#include "utils/TimeUtil.hpp"
// #include <iostream>
void cBaseObject::UpdateTriangleNormal()
{
    // we assume the rotation axis of v0, v1, v2 is the normal direction here
    // cTimeUtil::Begin("update_normal");
    for (auto &tri : mTriangleArrayShared)
    {
        const tVector &v0 = mVertexArrayShared[tri->mId0]->mPos;
        const tVector &v1 = mVertexArrayShared[tri->mId1]->mPos;
        const tVector &v2 = mVertexArrayShared[tri->mId2]->mPos;
        tri->mNormal = (v1 - v0).cross3(v2 - v1).normalized();
        // std::cout << tri->mNormal.transpose() << std::endl;
    }
    // cTimeUtil::End("update_normal");
}

/**
 * \brief       update the vertex from triangle normal
 */
#include <iostream>
void cBaseObject::UpdateVertexNormalFromTriangleNormal()
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

/**
 * \brief           set the alpha channel for vertex color
 */
void cBaseObject::SetVertexColorAlpha(float val)
{
    mColorAlpha = val;
    for (auto &v : mVertexArrayShared)
    {
        v->mColor[3] = val;
    }
}

/**
 * \brief           get vertex color alpha
 */
float cBaseObject::GetVertexColorAlpha() const { return mColorAlpha; }

/**
 * \brief       calcualte the total area
 */
double cBaseObject::CalcTotalArea() const
{
    float total_area = 0;
    for (auto &t : mTriangleArrayShared)
    {
        total_area += cMathUtil::CalcTriangleArea(mVertexArrayShared[t->mId0]->mPos,
                                                  mVertexArrayShared[t->mId1]->mPos,
                                                  mVertexArrayShared[t->mId2]->mPos);
    }
    return total_area;
}

void cBaseObject::UpdateImGUi()
{
    
}