#include "geometries/Arrow.h"
#include "geometries/Primitives.h"
#include "utils/ColorUtil.h"
#include "utils/FileUtil.h"
#include "utils/MathUtil.h"
#include "utils/ObjUtil.h"
#include "utils/RenderUtil.h"
#include "utils/RotUtil.h"
#include <iostream>

cArrow::cArrow()
{
    std::string arrow_path = "./data/arrow.obj";
    SIM_ASSERT("arrow resource doesn't exist!" &&
               cFileUtil::ExistsFile(arrow_path));
    cObjUtil::LoadObj(arrow_path, mVertexArray, mEdgeArray, mTriangleArray);
    for (auto &v : mVertexArray)
    {
        v->mPos[3] = 1;
        v->mColor = ColorPurple;
    }
    for (auto &t : mTriangleArray)
    {
        t->mColor = ColorPurple;
    }
    // get length
    tVector aabb_min, aabb_max;
    GetAABB(aabb_min, aabb_max);
    tVector aabb_size = aabb_max - aabb_min;
    mInitLength = aabb_size[1];
    mInitDiameter = aabb_size[0];
    mSt = tVector(0, 0, 0, 1);
    mEd = tVector(0, mInitLength, 0, 1);

    // remember init pos
    mVertexInitPos.clear();
    for (auto &v : mVertexArray)
    {
        mVertexInitPos.push_back(v->mPos);
    }
}

void cArrow::GetAABB(tVector &aabb_min, tVector &aabb_max)
{
    aabb_min.segment(0, 3).array() = std::numeric_limits<double>::max();
    aabb_min[3] = 0;
    aabb_max.segment(0, 3).array() = -std::numeric_limits<double>::max();
    aabb_max[3] = 0;

    for (auto &v : mVertexArray)
    {
        cMathUtil::ContainsAABB(v->mPos, aabb_min, aabb_max);
    }
}

void cArrow::SetStEd(const tVector &st, const tVector &ed)
{

    tMatrix scale = tMatrix::Identity(), shift = tMatrix::Identity(),
            rotate = tMatrix::Identity();
    // 1. first do scale
    {
        double length = (ed - st).norm();
        scale(0, 0) = length;
        scale(1, 1) = length;
        scale(2, 2) = length;
    }

    // 2. do rotation
    {
        tVector ori_dir = tVector(0, 1, 0, 0);
        tVector tar_dir = (ed - st).normalized();
        tVector aa =
            cRotUtil::CalcAxisAngleFromOneVectorToAnother(ori_dir, tar_dir);
        // auto res = cRotUtil::CalcAxisAngleFromOneVectorToAnother(
        //     tVector(0, 1, 0, 0), tVector(0, -1, 0, 0));
        // std::cout << "aa = " << res.transpose() << std::endl;

        rotate = cRotUtil::AxisAngleToRotmat(aa);

        // std::cout << "rotate = " << res.transpose() << std::endl;
        // tVector res_vec = rotate * tVector(0, 1, 0, 1);
        // std::cout << "res = " << res_vec.transpose() << std::endl;
        // exit(1);
    }
    // 3. do shift
    {
        shift.block(0, 3, 3, 1) = st.segment(0, 3);
    }
    tMatrix transform = shift * rotate * scale;

    for (int i = 0; i < mVertexArray.size(); i++)
    {
        auto cur_v = mVertexArray[i];
        cur_v->mPos = transform * mVertexInitPos[i];
    }
}

void cArrow::CalcTriangleDrawBuffer(Eigen::Map<tVectorXf> &res, int &st) const
{
    int old_st = st;
    for (auto &tri : mTriangleArray)
    {
        cRenderUtil::CalcTriangleDrawBufferSingle(
            this->mVertexArray[tri->mId0], this->mVertexArray[tri->mId1],
            this->mVertexArray[tri->mId2], tri->mColor, res, st);
    }
}
void cArrow::CalcEdgeDrawBuffer(Eigen::Map<tVectorXf> &res, int &st) const
{
    tVector normal = tVector::Zero();

    for (int e_id = 0; e_id < mEdgeArray.size(); e_id++)
    {
        auto e = mEdgeArray[e_id];
        tVector color = e->mColor;
        // if (affected_edge.find(e_id) != affected_edge.end())
        // {
        //     color = ColorPurple;
        // }
        // 1. get the normal of this edge
        normal =
            (mVertexArray[e->mId0]->mNormal + mVertexArray[e->mId1]->mNormal) /
            2;

        cRenderUtil::CalcEdgeDrawBufferSingle(mVertexArray[e->mId0],
                                              mVertexArray[e->mId1], normal,
                                              res, st, color);
    }
}

int cArrow::GetNumOfDrawEdges() const { return this->mEdgeArray.size(); }
int cArrow::GetNumOfDrawTriangles() const
{
    return this->mTriangleArray.size();
}