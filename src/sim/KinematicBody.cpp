#include "KinematicBody.h"
#include "utils/JsonUtil.h"
#include "utils/ObjUtil.h"
#include "utils/DefUtil.h"
#include "utils/RenderUtil.h"
#include "geometries/Primitives.h"
#include <iostream>
std::string gBodyShapeStr[eKinematicBodyShape::NUM_OF_KINEMATIC_SHAPE] = {
    "plane", "cube", "sphere", "capsule", "custom"};
cKinematicBody::cKinematicBody(int id_)
    : cBaseObject(eObjectType::KINEMATICBODY_TYPE, id_)
{
    mIsStatic = true;
    mBodyShape = eKinematicBodyShape::KINEMATIC_INVALID;
    mCustomMeshPath = "";
    mTargetAABBDontUseDirectly = tVector::Zero();
    mScaleDontUseDirectly.setZero();
    mPlaneEquation.setZero();
    mTargetPos.setZero();
    mCurTime = 0;
    mTargetOrientation.setZero();
}

cKinematicBody::~cKinematicBody() {}
void cKinematicBody::Init(const Json::Value &value)
{
    cBaseObject::Init(value);
    std::string type =
        cJsonUtil::ParseAsString(cKinematicBody::TYPE_KEY, value);
    mBodyShape = BuildKinematicBodyShape(type);
    switch (mBodyShape)
    {
    case eKinematicBodyShape::KINEMATIC_CUSTOM:
    {
        mCustomMeshPath =
            cJsonUtil::ParseAsString(cKinematicBody::MESH_PATH_KEY, value);
        cJsonUtil::ReadVectorJson(
            cJsonUtil::ParseAsValue(cKinematicBody::TARGET_AABB_KEY, value),
            mTargetAABBDontUseDirectly);

        cJsonUtil::ReadVectorJson(
            cJsonUtil::ParseAsValue(cKinematicBody::SCALE_KEY, value),
            mScaleDontUseDirectly);

        cJsonUtil::ReadVectorJson(
            cJsonUtil::ParseAsValue(cKinematicBody::TRANSLATION_KEY, value),
            mInitPos);
        cJsonUtil::ReadVectorJson(
            cJsonUtil::ParseAsValue(cKinematicBody::ORIENTATION_KEY, value),
            mInitOrientation);

        if (value.isMember(cKinematicBody::IS_STATIC_KEY) == true)
        {
            mIsStatic =
                cJsonUtil::ParseAsBool(cKinematicBody::IS_STATIC_KEY, value);
            if (mIsStatic == false)
            {
                // parse target translation, orientation
                cJsonUtil::ReadVectorJson(
                    cJsonUtil::ParseAsValue(
                        cKinematicBody::TARGET_ORIENTATION_KEY, value),
                    mTargetOrientation);
                cJsonUtil::ReadVectorJson(
                    cJsonUtil::ParseAsValue(
                        cKinematicBody::TARGET_TRANSLATION_KEY, value),
                    mTargetPos);
                // parse elasped time
                mMovingElaspedTimeSec = cJsonUtil::ParseAsDouble(
                    cKinematicBody::ELASPED_TIME_SEC_KEY, value);
            }
            // std::cout << "target orientation = "
            //           << mTargetOrientation.transpose() << std::endl;
            // std::cout << "target pos = " << mTargetPos.transpose() <<
            // std::endl; std::cout << "elasped time sec = " <<
            // mMovingElaspedTimeSec
            //           << std::endl;
            // exit(1);
        }
        else
        {
            mTargetPos = mInitPos;
            mTargetOrientation = mInitOrientation;
        }
        BuildCustomKinematicBody();
        UpdateCurWorldTransformByTime();
        SetMeshPos();
        break;
    }
    case eKinematicBodyShape::KINEMATIC_PLANE:
    {
        cJsonUtil::ReadVectorJson(
            cJsonUtil::ParseAsValue(cKinematicBody::PLANE_EQUATION_KEY, value),
            mPlaneEquation);
        // std::cout << "plane equation = " << mPlaneEquation.transpose() <<
        // std::endl;
        mPlaneScale =
            cJsonUtil::ParseAsDouble(cKinematicBody::PLANE_SCALE_KEY, value);
        BuildPlane();
        break;
    }
    default:
        SIM_ERROR("Unsupported kinematic shape {}", type);
    }
    tVector min, max;
    CalcAABB(min, max);

    UpdateTriangleNormal();
    UpdateVertexNormalFromTriangleNormal();

    // std::cout << "[debug] obstacle aabb min = " << min.transpose() <<
    // std::endl; std::cout << "[debug] obstacle aabb max = " << max.transpose()
    // << std::endl;
}

eKinematicBodyShape
cKinematicBody::BuildKinematicBodyShape(std::string type_str)
{
    eKinematicBodyShape shape = eKinematicBodyShape::KINEMATIC_INVALID;
    for (int i = 0; i < eKinematicBodyShape::NUM_OF_KINEMATIC_SHAPE; i++)
    {
        if (gBodyShapeStr[i] == type_str)
        {
            shape = static_cast<eKinematicBodyShape>(i);
            break;
        }
    }
    SIM_ASSERT(shape != eKinematicBodyShape::KINEMATIC_INVALID);
    return shape;
}

bool cKinematicBody::IsStatic() const { return mIsStatic; }
#include "geometries/Triangulator.h"

/**
 * \brief           Build plane data strucutre
 */
void cKinematicBody::BuildPlane()
{
    // 1. build legacy XOZ plane, then do a transformation
    // for (int i = 0; i < 4; i++)
    cObjUtil::BuildPlaneGeometryData(mPlaneScale, this->mPlaneEquation,
                                     mVertexArrayShared, mEdgeArrayShared, mTriangleArrayShared);
    for (auto &x : mVertexArrayShared)
        x->mColor = tVector(0.1f, 0.1f, 0.1f, 1);
}

/**
 * \brief               Build custom kinematic body
 *      1. fill the scaled mesh vertices
 *      2. load the obj
 *      3. fill the color setting
 */
void cKinematicBody::BuildCustomKinematicBody()
{
    // std::cout << "mesh path = " << mCustomMeshPath << std::endl;
    cObjUtil::tParams obj_params;
    obj_params.mPath = mCustomMeshPath;
    cObjUtil::LoadObj(obj_params, mVertexArrayShared, mEdgeArrayShared, mTriangleArrayShared);

    // tMatrix trans = GetWorldTransform();
    tVector scale_vec = GetScaleVec();
    mScaledMeshVertices.noalias() = tVectorXd::Zero(mVertexArrayShared.size() * 3);
    for (int i = 0; i < mVertexArrayShared.size(); i++)
    {
        auto &x = mVertexArrayShared[i];
        mScaledMeshVertices.segment(3 * i, 3) =
            (scale_vec.cwiseProduct(x->mPos)).segment(0, 3);
        x->mColor = tVector::Ones() * 0.3;
        x->mColor[3] = 1.0;
    }
    // exit(0);
    cTriangulator::ValidateGeometry(mVertexArrayShared, mEdgeArrayShared, mTriangleArrayShared);
}

/**
 * \brief           init the mesh pos by given init_pos and init_orien
 */
void cKinematicBody::SetMeshPos()
{
    tMatrix trans = GetCurWorldTransform();
    for (int i = 0; i < mVertexArrayShared.size(); i++)
    {
        auto &x = mVertexArrayShared[i];
        x->mPos =
            trans * cMathUtil::Expand(mScaledMeshVertices.segment(3 * i, 3), 1);
    }
}

// int cKinematicBody::GetDrawNumOfTriangles() const
// {
//     return mTriangleArrayShared.size();
// }

// int cKinematicBody::GetDrawNumOfEdges() const { return mEdgeArrayShared.size(); }

void cKinematicBody::CalcTriangleDrawBuffer(Eigen::Map<tVectorXf> &res,
                                            int &st) const
{
    for (auto &x : mTriangleArrayShared)
    {
        cRenderUtil::CalcTriangleDrawBufferSingle(mVertexArrayShared[x->mId0],
                                     mVertexArrayShared[x->mId1],
                                     mVertexArrayShared[x->mId2], res, st);
    }
}

void cKinematicBody::CalcEdgeDrawBuffer(Eigen::Map<tVectorXf> &res,
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

        cRenderUtil::CalcEdgeDrawBufferSingle(mVertexArrayShared[e->mId0], mVertexArrayShared[e->mId1],
                                 normal, res, st, black_color);
    }
}

/**
 * \brief           Get the world transform of this kinematic body
 *          this matrix can convert local pos to world pos in homogeneous coords
 *          world_pos = T * local_pos
 */
tMatrix GetWorldTransform(const tVector &init_pos, const tVector &init_ori)
{
    tMatrix trans = tMatrix::Identity();
    trans.block(0, 3, 3, 1) = init_pos.segment(0, 3);
    trans.block(0, 0, 3, 3) =
        cMathUtil::EulerAnglesToRotMat(init_ori, eRotationOrder::XYZ)
            .topLeftCorner<3, 3>();
    return trans;
}

void cKinematicBody::UpdateCurWorldTransformByTime()
{
    if (IsStatic())
    {
        mCurWorldTransform = GetWorldTransform(mInitPos, mInitOrientation);
    }
    else
    {
        float cur_time = mCurTime;
        if (cur_time > mMovingElaspedTimeSec)
            cur_time = mMovingElaspedTimeSec;
        tVector cur_pos = (mMovingElaspedTimeSec - cur_time) /
                              mMovingElaspedTimeSec * mInitPos +
                          cur_time / mMovingElaspedTimeSec * mTargetPos;
        tVector cur_ori = (mMovingElaspedTimeSec - cur_time) /
                              mMovingElaspedTimeSec * mInitOrientation +
                          cur_time / mMovingElaspedTimeSec * mTargetOrientation;
        mCurWorldTransform = GetWorldTransform(cur_pos, cur_ori);
        // std::cout << "[kin] update cur world transform, cur pos = "
        //           << cur_pos.transpose()
        //           << ", cur ori = " << cur_ori.transpose() << std::endl;
    }
}

tVector cKinematicBody::GetScaleVec() const
{
    float eps = 1e-6;
    if (mTargetAABBDontUseDirectly.norm() > eps &&
        mScaleDontUseDirectly.norm() > eps)
    {
        SIM_ERROR("{} and {} can only be set once for body {}",
                  cKinematicBody::SCALE_KEY, cKinematicBody::TARGET_AABB_KEY,
                  this->mObjName);
        exit(1);
    }

    // use AABB
    tVector scale_vec = tVector::Ones();
    if (mTargetAABBDontUseDirectly.norm() > eps)
    {
        tVector aabb_min, aabb_max;
        CalcAABB(aabb_min, aabb_max);
        tVector aabb = aabb_max - aabb_min;
        // std::cout << "init aabb = " << (aabb_max - aabb_min).transpose() <<
        // std::endl; exit(0);
        for (int i = 0; i < 3; i++)
        {
            if (std::fabs(aabb[i]) < 1e-10)
            {
                continue;
            }
            scale_vec[i] = mTargetAABBDontUseDirectly[i] / aabb[i];
        }
    }

    else
    {
        // use scale_mat
        scale_vec.segment(0, 3) = this->mScaleDontUseDirectly.segment(0, 3);
    }
    return scale_vec;
}
/**
 * \brief           update kinectmatic body
 */
void cKinematicBody::Update(float dt)
{
    mCurTime += dt;
    if (IsStatic() == true)
    {
        return;
    }
    else
    {
        UpdateCurWorldTransformByTime();
        SetMeshPos();
        // std::cout << "update kin body transform, cur = \n"
        //           << mCurWorldTransform << std::endl;
    }
}
tMatrix cKinematicBody::GetCurWorldTransform() const
{
    return this->mCurWorldTransform;
}

void cKinematicBody::Reset() { this->mCurTime = 0; }

tVector cKinematicBody::CalcCOM() const
{
    tVector com = tVector::Zero();
    for (auto &v : this->mVertexArrayShared)
    {
        com += v->mPos;
    }
    com /= mVertexArrayShared.size();
    return com;
}
void cKinematicBody::MoveTranslation(const tVector &shift)
{
    for (auto &v : this->mVertexArrayShared)
    {
        v->mPos += shift;
    }
}

void cKinematicBody::ApplyScale(float scale)
{
    for (auto &v : mVertexArrayShared)
    {
        v->mPos.segment(0, 3) = v->mPos.segment(0, 3) * scale;
    }
}

void cKinematicBody::ApplyUserPerturbForceOnce(tPerturb *)
{
}