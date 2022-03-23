#include "BVHCollisionDetecter.h"
#include "geometries/ObjectBVH.h"
#include "sim/BaseObject.h"
#include "sim/collision/CollisionInfo.h"
#include "utils/TimeUtil.hpp"
#include <iostream>
#include <queue>
/**
 * \brief       check whether the point & triangle are intersected:
 *      point penetrated (normal product < 0): depth < 0, intersected;
 *      point non penetrated (normal): depth > 0; if depth < eps, intersected;
 * else non intersected
 *
 *
 *      point's projection point is outside of the triangle: not
 * intersected
 */
bool PointTriangleDist(const tVector &pt, const tVector &v0, const tVector &v1,
                       const tVector &v2, float &depth, const float eps = 2e-3)
{
    depth = 0;

    // 1. check inside or outside. if the projection is outside, not intersected
    {
        tVector cross0 = (v1 - v0).cross3(pt - v0);
        tVector cross1 = (v2 - v1).cross3(pt - v1);
        tVector cross2 = (v0 - v2).cross3(pt - v2);
        int sign0 = (cross0.dot(cross1) > 0) ? 1 : -1;
        int sign1 = (cross2.dot(cross1) > 0) ? 1 : -1;
        int sign2 = (cross0.dot(cross2) > 0) ? 1 : -1;
        if (sign0 < 0 || sign1 < 0 || sign2 < 0)
            return false;
    }

    // inside, calculate distance
    // 2. calculate normal
    tVector normal = (v1 - v0).cross3(v2 - v1).normalized();
    normal[3] = 0;
    depth = (pt - v0).dot(normal);
    return std::fabs(depth) <= eps;
}

bool EdgeEdgeDist(const tVector &v0, const tVector &v1, const tVector &v2,
                  const tVector &v3, float &dist, const float eps = 2e-3)
{
    dist = 0;
    tVector v0v1 = v1 - v0, v2v3 = v3 - v2;
    tVector v2v0 = v0 - v2;
    tVector v0v1_dir = v0v1.normalized();
    bool is_parallel = (v0v1_dir.cross3(v2v3.normalized())).norm() < 1e-5;
    if (is_parallel == true)
    {
        // calc dist: project v2v0 to v0v1 direction
        tVector v0v2 = v2 - v0;
        dist = (v0v2 - v0v2.dot(v0v1_dir) * v0v1_dir).norm();
        return dist < eps;
    }
    else
    {

        double w0 = v0v1.squaredNorm(), w1 = -v2v3.dot(v0v1),
               w2 = v2v3.squaredNorm(), w3 = v2v0.dot(v0v1),
               w4 = -v2v0.dot(v2v3);
        double common_coef = -1.0 / (w0 * w2 - w1 * w1);
        double a = common_coef * (w2 * w3 - w1 * w4),
               b = common_coef * (-w1 * w3 + w0 * w4);
        tVector P = v0 + (v1 - v0) * a;
        tVector Q = v2 + (v3 - v2) * b;
        // min dist
        dist = (P - Q).norm();
        return dist < eps;
    }
}

tPointTriangleInfo::tPointTriangleInfo(int obj0, int obj1)
{
    mObj0Id = obj0;
    mObj1Id = obj1;
    mPointTrianglePair.clear();
}

bool IsNodeIntersected(const tBVHNodePtr a, const tBVHNodePtr b)
{
    return a->mAABB.Intersect(b->mAABB);
}

bool DescendA(const tBVHNodePtr a, const tBVHNodePtr b) { return !a->mIsLeaf; }
typedef std::pair<tBVHNodePtr, tBVHNodePtr> tNodePair;
typedef std::pair<int, int> tTrianglePair;
std::vector<tTrianglePair> IntersectTwoBVHs(tBVHNodePtr a, tBVHNodePtr b)
{
    std::vector<tNodePair> stack = {};
    std::vector<tTrianglePair> broadphase_result =
        {}; // possible intersected faces
    while (1)
    {
        // printf("visit a %d b %d\n", a->mId, b->mId);
        if (IsNodeIntersected(a, b))
        {
            if (a->mIsLeaf && b->mIsLeaf)
            {
                // printf("[debug] node leaf intersected!\n");
                // std::cout << "a min = " << a->mAABB.mMin.transpose()
                //           << " max = " << a->mAABB.mMax.transpose() <<
                //           std::endl;
                // std::cout << "b min = " << b->mAABB.mMin.transpose()
                //           << " max = " << b->mAABB.mMax.transpose() <<
                //           std::endl;
                // printf("add %d %d\n", a->mTriangleId, b->mTriangleId);
                // narrow phase
                broadphase_result.push_back(
                    tTrianglePair(a->mTriangleId, b->mTriangleId));
                // print("need narrow phase, trianlge intersection");
            }
            else
            {
                if (DescendA(a, b))
                {
                    stack.push_back(tNodePair(a->mRight, b));
                    a = a->mLeft;
                    continue;
                }
                else
                {
                    stack.push_back(tNodePair(a, b->mRight));
                    b = b->mLeft;
                    continue;
                }
            }
        }
        if (stack.empty() == true)
            break;
        // pop
        a = stack.back().first;
        b = stack.back().second;
        stack.pop_back();
    }
    return broadphase_result;
}

cBVHCollisionDetecter::cBVHCollisionDetecter() {}
cBVHCollisionDetecter::~cBVHCollisionDetecter() {}
void cBVHCollisionDetecter::AddObject(cBaseObjectPtr obj,
                                      bool enable_self_collision)
{
    mColObjs.push_back(obj);
    mEnableSelfCollision.push_back(enable_self_collision);
}
void cBVHCollisionDetecter::Init()
{
    cCollisionDetecter::Init();

    // ! begin to build BVH
    this->mBVHList.clear();
    for (auto &obj : mColObjs)
    {
        auto bvh = std::make_shared<cObjBVH>();
        bvh->Init(obj->GetObjId(), obj->GetVertexArray(), obj->GetEdgeArray(),
                  obj->GetTriangleArray());
        bvh->UpdateAABB();
        mBVHList.push_back(bvh);
    }
}
#include <set>
void CalculatePointTrianglePair(
    const cBaseObjectPtr obj0, const cBaseObjectPtr obj1,
    const std::vector<std::pair<int, int>> &triangle_pair,
    std::vector<std::pair<int, int>> &obj0point_obj1triangle,
    std::vector<std::pair<int, int>> &obj1point_obj0triangle)
{
    const auto t0_array = obj0->GetTriangleArray();
    const auto t1_array = obj1->GetTriangleArray();

    const auto v0_array = obj0->GetVertexArray();
    const auto v1_array = obj1->GetVertexArray();
    std::set<std::pair<int, int>> set_obj0point_obj1_triangle = {};
    std::set<std::pair<int, int>> set_obj1point_obj0_triangle = {};
    for (auto &pair : triangle_pair)
    {
        int t0 = pair.first;
        int t1 = pair.second;
        int t0_v[3] = {t0_array[t0]->mId0, t0_array[t0]->mId1,
                       t0_array[t0]->mId2};
        int t1_v[3] = {t1_array[t1]->mId0, t1_array[t1]->mId1,
                       t1_array[t1]->mId2};

        // obj0 point, obj1 triangle
        for (int i = 0; i < 3; i++)
        {
            set_obj0point_obj1_triangle.insert(
                std::pair<int, int>(t0_v[i], t1));
        }
        for (int i = 0; i < 3; i++)
        {
            set_obj1point_obj0_triangle.insert(
                std::pair<int, int>(t1_v[i], t0));
        }

        // obj1 point, obj0 triangle
    }

    obj0point_obj1triangle = std::vector<std::pair<int, int>>(
        set_obj0point_obj1_triangle.begin(), set_obj0point_obj1_triangle.end());
    obj1point_obj0triangle = std::vector<std::pair<int, int>>(
        set_obj1point_obj0_triangle.begin(), set_obj1point_obj0_triangle.end());
}
void cBVHCollisionDetecter::PerformCD()
{
    mPointTriangleCollisionInfo.clear();
    mObjPointTriangleInfo.clear();
    mObjPointTriangleInfo.resize(this->mColObjs.size());
    // 1. object pair, do collision detection (broadphase)
    cTimeUtil::Begin("broadphase");
    int num_of_obj = mBVHList.size();
    std::vector<tPointTriangleInfo> broadphase_info_lst = {};
    int broadphase_num = 0;
    for (int i = 0; i < num_of_obj; i++)
    {
        auto i_bvh = mBVHList[i];
        for (int j = i + 1; j < num_of_obj; j++)
        {
            auto j_bvh = mBVHList[j];
            auto obj0 = mColObjs[i], obj1 = mColObjs[j];
            // check (i, j)
            tPointTriangleInfo info01(obj0->GetObjId(), obj1->GetObjId()),
                info10(obj1->GetObjId(), obj0->GetObjId());
            auto mTrianglePairInfo =
                IntersectTwoBVHs(i_bvh->GetRootNode(), j_bvh->GetRootNode());
            CalculatePointTrianglePair(obj0, obj1, mTrianglePairInfo,
                                       info01.mPointTrianglePair,
                                       info10.mPointTrianglePair);
            // generate p-t & e-e pair which we need to test

            if (mTrianglePairInfo.size() > 0)
            {
                broadphase_num += info01.mPointTrianglePair.size();
                broadphase_num += info10.mPointTrianglePair.size();
                broadphase_info_lst.push_back(info01);
                broadphase_info_lst.push_back(info10);
            }
            // tBroadphaseInfo info = CheckBetweenTwoBVH(i_bvh, j_bvh);
            // broadphase_info_lst.push_back(info);
        }
    }
    cTimeUtil::End("broadphase");
    std::cout << "broadphase num = " << broadphase_num << std::endl;

    // 2. do narrow phase, detail the information
    // we also want to calculate the
    cTimeUtil::Begin("narrowphase");
    for (auto &info : broadphase_info_lst)
    {
        auto obj0 = mColObjs[info.mObj0Id];
        auto obj1 = mColObjs[info.mObj1Id];
        auto v0_array = obj0->GetVertexArray();
        auto v1_array = obj1->GetVertexArray();
        auto tri1_array = obj1->GetTriangleArray();
        for (auto pt_pair : info.mPointTrianglePair)
        {
            int p_id = pt_pair.first;
            int t_id = pt_pair.second;
            float depth = 0;
            if (true == PointTriangleDist(
                            v0_array[p_id]->mPos,
                            v1_array[tri1_array[t_id]->mId0]->mPos,
                            v1_array[tri1_array[t_id]->mId1]->mPos,
                            v1_array[tri1_array[t_id]->mId2]->mPos, depth))
            {
                tPointTriangleCollisionInfoPtr info =
                    std::make_shared<tPointTriangleCollisionInfo>();
                info->mObj0 = obj0;
                info->mObj1 = obj1;
                info->mVertexId0 = p_id;
                info->mTriangleId1 = t_id;
                mPointTriangleCollisionInfo.push_back(info);
                mObjPointTriangleInfo[obj0->GetObjId()].push_back(info);
                mObjPointTriangleInfo[obj1->GetObjId()].push_back(info);
            }
        }
    }
    std::cout << "narrow phase result = " << mPointTriangleCollisionInfo.size()
              << std::endl;
    for (auto &x : mPointTriangleCollisionInfo)
    {
        printf("[narrow] obj %d v %d collided with obj %d tri %d\n",
               x->mObj0->GetObjId(), x->mVertexId0, x->mObj1->GetObjId(),
               x->mTriangleId1);
    }
    cTimeUtil::End("narrowphase");
}

void cBVHCollisionDetecter::Clear() {}

std::vector<tPointTriangleCollisionInfoPtr>
cBVHCollisionDetecter::GetAllPointTriangleCollisionInfo() const
{
    return this->mPointTriangleCollisionInfo;
}

void cBVHCollisionDetecter::Update() { UpdateBVHAABB(); }

void cBVHCollisionDetecter::UpdateBVHAABB()
{
    // cTimeUtil::Begin("update BVH AABB");
    for (auto &bvh : mBVHList)
    {
        // if the object is static, do not update it
        bvh->UpdateAABB();
        if (bvh->GetNumOfLeaves() < 10)
        {
            bvh->Print();
        }
    }
}

std::vector<tPointTriangleCollisionInfoPtr>
cBVHCollisionDetecter::GetObjPointTriangleCollisionInfo(int obj_id) const
{
    return mObjPointTriangleInfo[obj_id];
}
