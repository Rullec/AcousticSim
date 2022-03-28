#include "BVHCollisionDetecter.h"
#include "geometries/ObjectBVH.h"
#include "sim/BaseObject.h"
#include "sim/collision/CollisionInfo.h"
#include "utils/TimeUtil.hpp"
#include <iostream>
#include <queue>
#include <set>
tEdgeEdgeInfo::tEdgeEdgeInfo(int o0, int o1)
{
    mObj0Id = o0;
    mObj1Id = o1;
    mEdgePair.clear();
    mBaryA.clear();
    mBaryB.clear();
}
tVector3d ClosestPtPointTriangle(const tVector3d &p, const tVector3d &a,
                                 const tVector3d &b, const tVector3d &c,
                                 tVector3f &bary)
{
    // Check if P in vertex region outside A
    tVector3d ab = b - a;
    tVector3d ac = c - a;
    tVector3d ap = p - a;
    float d1 = ab.dot(ap);
    float d2 = ac.dot(ap);
    bary.setZero();
    if (d1 <= 0.0f && d2 <= 0.0f)
    {
        bary[0] = 1;
        return a; // barycentric coordinates (1,0,0)
    }
    // Check if P in vertex region outside B
    tVector3d bp = p - b;
    float d3 = ab.dot(bp);
    float d4 = ac.dot(bp);
    if (d3 >= 0.0f && d4 <= d3)
    {
        bary[1] = 1;
        return b; // barycentric coordinates (0,1,0)
    }
    // Check if P in edge region of AB, if so return projection of P onto AB
    float vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0f && d1 >= 0.0f && d3 <= 0.0f)
    {
        float v = d1 / (d1 - d3);
        bary[0] = 1 - v;
        bary[1] = v;
        return a + v * ab; // barycentric coordinates (1-v,v,0)
    }
    // Check if P in vertex region outside C
    tVector3d cp = p - c;
    float d5 = ab.dot(cp);
    float d6 = ac.dot(cp);
    if (d6 >= 0.0f && d5 <= d6)
    {
        bary[2] = 1;
        return c; // barycentric coordinates (0,0,1)
    }
    // Check if P in edge region of AC, if so return projection of P onto AC
    float vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f)
    {
        float w = d2 / (d2 - d6);
        bary[0] = 1 - w;
        bary[2] = w;
        return a + w * ac; // barycentric coordinates (1-w,0,w)
    }
    // Check if P in edge region of BC, if so return projection of P onto BC
    float va = d3 * d6 - d5 * d4;
    if (va <= 0.0f && (d4 - d3) >= 0.0f && (d5 - d6) >= 0.0f)
    {
        float w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        bary[1] = 1 - w;
        bary[2] = w;
        return b + w * (c - b); // barycentric coordinates (0,1-w,w)
    }
    // P inside face region. Compute Q through its barycentric coordinates
    // (u,v,w)
    float denom = 1.0f / (va + vb + vc);
    float v = vb * denom;
    float w = vc * denom;
    bary[0] = 1 - v - w;
    bary[1] = v;
    bary[2] = w;
    return a + ab * v + ac * w; // = u*a + v*b + w*c, u = va * denom = 1.0f-v-w
}
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

/*
    [p1, q1, p2, q2]
    c1 = p1 + (q1 - p1) * s = (1 - s) * p1 + s * q1
    c2 = p2 + (q2 - p2) * t = (1 - t) * p2 + t * q2
*/
void EdgeEdgeDist(const tVector3d &p1, const tVector3d &q1, const tVector3d &p2,
                  const tVector3d &q2, float &dist, float &s, float &t)
{
    tVector3d d1 = q1 - p1; // Direction vector of segment S1
    tVector3d d2 = q2 - p2; // Direction vector of segment S2
    tVector3d r = p1 - p2;
    float a = d1.dot(d1); // Squared length of segment S1, always nonnegative
    float e = d2.dot(d2); // Squared length of segment S2, always nonnegative
    float f = d2.dot(r);
    // Check if either or both segments degenerate into points
    float EPSILON = 1e-10;
    s = 0, t = 0;
    tVector3d c1, c2;
    if (a <= EPSILON && e <= EPSILON)
    {
        // Both segments degenerate into points
        s = t = 0.0f;
        c1 = p1;
        c2 = p2;
        dist = (c1 - c2).norm();
    }
    if (a <= EPSILON)
    {
        // First segment degenerates into a point
        s = 0.0f;
        t = f / e; // s = 0 => t = (b*s + f) / e = f / e
        t = cMathUtil::Clamp(t, 0.0f, 1.0f);
    }
    else
    {
        float c = d1.dot(r);
        if (e <= EPSILON)
        {
            // Second segment degenerates into a point
            t = 0.0f;
            s = cMathUtil::Clamp(-c / a, 0.0f,
                                 1.0f); // t = 0 => s = (b*t - c) / a = -c / a
        }
        else
        {
            // The general nondegenerate case starts here
            float b = d1.dot(d2);
            float denom = a * e - b * b; // Always nonnegative
            // If segments not parallel, compute closest point on L1 to L2
            // and clamp to segment S1. Else pick arbitrary s (here 0)
            if (denom != 0.0f)
            {
                s = cMathUtil::Clamp((b * f - c * e) / denom, 0.0f, 1.0f);
            }
            else
                s = 0.0f;
            // Compute point on L2 closest to S1(s) using
            // t = Dot((P1 + D1*s) - P2,D2) / Dot(D2,D2) = (b*s + f) / e
            t = (b * s + f) / e;
            // If t in [0,1] done. Else clamp t, recompute s for the new
            // value of t using s = Dot((P2 + D2*t) - P1,D1) / Dot(D1,D1)=
            // (t*b - c) / a and clamp s to [0, 1]
            if (t < 0.0f)
            {
                t = 0.0f;
                s = cMathUtil::Clamp(-c / a, 0.0f, 1.0f);
            }
            else if (t > 1.0f)
            {
                t = 1.0f;
                s = cMathUtil::Clamp((b - c) / a, 0.0f, 1.0f);
            }
        }
    }
    c1 = p1 + d1 * s;
    c2 = p2 + d2 * t;
    dist = (c1 - c2).norm();
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
    std::vector<tTrianglePair> triangle_pair = {}; // possible intersected faces
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
                triangle_pair.push_back(
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
    return triangle_pair;
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
typedef std::pair<int, int> int_pair;
#define BUILD_PAIR(a, b) SIM_MIN(a, b), SIM_MAX(a, b)
std::vector<std::pair<int, int>>
CalculateEdgeEdgePair(const cBaseObjectPtr obj0, const cBaseObjectPtr obj1,
                      const std::vector<std::pair<int, int>> &triangle_pair_lst)
{
    const auto t0_array = obj0->GetTriangleArray();
    const auto t1_array = obj1->GetTriangleArray();
    std::set<int_pair> edge_pair_set = {};
    for (auto &t_pair : triangle_pair_lst)
    {
        int t0_e[3] = {t0_array[t_pair.first]->mEId0,
                       t0_array[t_pair.first]->mEId1,
                       t0_array[t_pair.first]->mEId2};
        int t1_e[3] = {t1_array[t_pair.second]->mEId0,
                       t1_array[t_pair.second]->mEId1,
                       t1_array[t_pair.second]->mEId2};
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                edge_pair_set.insert(int_pair(t0_e[i], t1_e[j]));
            }
    }

    std::vector<std::pair<int, int>> edge_pair(edge_pair_set.begin(),
                                               edge_pair_set.end());
    return edge_pair;
}
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
void DetermineObj0AndObj1(const cBaseObjectPtr &obj0,
                          const cBaseObjectPtr &obj1,
                          bool &obj0_is_cloth_obj1_is_rb)
{
    auto obj0_type = obj0->GetObjectType();
    auto obj1_type = obj1->GetObjectType();
    // obj0: cloth, obj1: kinematic
    if (obj0_type == eObjectType::CLOTH_TYPE &&
        obj1_type == eObjectType::KINEMATICBODY_TYPE)
    {
        obj0_is_cloth_obj1_is_rb = true;
    }
    else if (obj1_type == eObjectType::CLOTH_TYPE &&
             obj0_type == eObjectType::KINEMATICBODY_TYPE)
    {
        obj0_is_cloth_obj1_is_rb = false;
    }
    else
    {
        SIM_ERROR("we cannot handle rigit body contant and cloth-cloth "
                  "intersection\n");
        exit(1);
    }
}
void cBVHCollisionDetecter::PerformCD()
{
    mPointTriangleCollisionInfo.clear();
    mObjPointTriangleInfo.clear();
    mObjPointTriangleInfo.resize(this->mColObjs.size());

    mEdgeEdgeCollisionInfo.clear();
    mObjEdgeEdgeCollisionInfo.clear();
    mObjEdgeEdgeCollisionInfo.resize(this->mColObjs.size());
    // 1. object pair, do collision detection (broadphase)
    cTimeUtil::Begin("broadphase");
    int num_of_obj = mBVHList.size();
    std::vector<tPointTriangleInfo> point_tri_broadphase = {};
    std::vector<tEdgeEdgeInfo> ee_broadphase = {};
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

            tEdgeEdgeInfo ee_info(obj0->GetObjId(), obj1->GetObjId());
            auto mTrianglePairInfo =
                IntersectTwoBVHs(i_bvh->GetRootNode(), j_bvh->GetRootNode());
            CalculatePointTrianglePair(obj0, obj1, mTrianglePairInfo,
                                       info01.mPointTrianglePair,
                                       info10.mPointTrianglePair);
            // for (auto &x : mTrianglePairInfo)
            // {
            //     printf("[borad] tri pair %d %d\n", x.first, x.second);
            // }
            ee_info.mEdgePair =
                CalculateEdgeEdgePair(obj0, obj1, mTrianglePairInfo);
            // generate p-t & e-e pair which we need to test
            // for (auto &e_pair : ee_info.mEdgePair)
            // {
            //     printf("[ee] obj0 %d, obj1 %d, e0 %d e1 %d\n", i, j,
            //            e_pair.first, e_pair.second);
            // }
            // std::cout << "obj0 has " << obj0->GetNumOfEdges() << "
            // edges\n"; std::cout << "obj1 has " << obj1->GetNumOfEdges()
            // << " edges\n";

            if (mTrianglePairInfo.size() > 0)
            {
                broadphase_num += info01.mPointTrianglePair.size();
                broadphase_num += info10.mPointTrianglePair.size();
                broadphase_num += ee_info.mEdgePair.size();
                point_tri_broadphase.push_back(info01);
                point_tri_broadphase.push_back(info10);

                ee_broadphase.push_back(ee_info);
            }

            // tBroadphaseInfo info = CheckBetweenTwoBVH(i_bvh, j_bvh);
            // point_tri_broadphase.push_back(info);
        }
    }
    cTimeUtil::End("broadphase");
    std::cout << "broadphase num = " << broadphase_num << std::endl;

    // 2. narrow phase, point triangle
    cTimeUtil::Begin("narrowphase");
    for (auto &info : point_tri_broadphase)
    {
        auto obj0 = mColObjs[info.mObj0Id];
        auto obj1 = mColObjs[info.mObj1Id];
        bool obj0_is_cloth_obj1_is_rb = false;
        DetermineObj0AndObj1(obj0, obj1, obj0_is_cloth_obj1_is_rb);

        auto v0_array = obj0->GetVertexArray();
        auto v1_array = obj1->GetVertexArray();
        auto tri1_array = obj1->GetTriangleArray();
        for (auto pt_pair : info.mPointTrianglePair)
        {
            int p_id = pt_pair.first;
            int t_id = pt_pair.second;
            float depth = 0;
            tVector3f bary = tVector3f::Zero();

            // calculate the closet point
            {
                const tVector3d pt = v0_array[p_id]->mPos.segment(0, 3);
                tVector3d closest_pt_on_triangle = ClosestPtPointTriangle(
                    pt, v1_array[tri1_array[t_id]->mId0]->mPos.segment(0, 3),
                    v1_array[tri1_array[t_id]->mId1]->mPos.segment(0, 3),
                    v1_array[tri1_array[t_id]->mId2]->mPos.segment(0, 3), bary);
                // judge the dist
                tVector3d outer_normal = tVector3d::Zero();
                double pt_tri_dist = 0;
                tVector3d dist_vec = tVector3d::Zero();
                if (obj0_is_cloth_obj1_is_rb == true)
                {
                    // obj0 is cloth point

                    // obj1 is rigid body triangle, get rigid body normal
                    outer_normal = tri1_array[t_id]->mNormal.segment(0, 3);

                    // so the distance vector: triangle -> point
                    dist_vec = (pt - closest_pt_on_triangle);
                }
                else
                {
                    // obj0 is rigid body point, get rigid body normal
                    tVector3d rb_point_normal =
                        v0_array[p_id]->mNormal.segment(0, 3);

                    // take triangle normal as outer normal direction
                    outer_normal = tri1_array[t_id]->mNormal.segment(0, 3);
                    if (outer_normal.dot(rb_point_normal) < 0)
                        outer_normal *= -1;
                    // obj1 is cloth triangle

                    // so the distance vector: point -> triangle
                    dist_vec = closest_pt_on_triangle - pt;
                }
                pt_tri_dist = dist_vec.norm();
                if (pt_tri_dist < 3e-3)
                {
                    // collided!
                    tPointTriangleCollisionInfoPtr info =
                        std::make_shared<tPointTriangleCollisionInfo>();
                    info->mObj0 = obj0;
                    info->mObj1 = obj1;
                    info->mOuterNormal = outer_normal;
                    info->obj0_is_cloth_obj1_is_rb = obj0_is_cloth_obj1_is_rb;
                    info->mDepth = dist_vec.dot(outer_normal);
                    info->mVertexId0 = p_id;
                    info->mTriangleId1 = t_id;
                    info->mBary = bary.cast<double>();
                    // printf("[pt] p%d and t%d intersect, depth %.3f\n", p_id,
                    //        t_id, info->mDepth);
                    mPointTriangleCollisionInfo.push_back(info);
                    mObjPointTriangleInfo[obj0->GetObjId()].push_back(info);
                    mObjPointTriangleInfo[obj1->GetObjId()].push_back(info);
                }
            }
        }
    }

    std::cout << "narrow phase pt result = "
              << mPointTriangleCollisionInfo.size() << std::endl;
    // for (auto &x : mPointTriangleCollisionInfo)
    // {
    //     printf("[narrow] obj %d v %d collided with obj %d tri %d\n",
    //            x->mObj0->GetObjId(), x->mVertexId0, x->mObj1->GetObjId(),
    //            x->mTriangleId1);
    // }

    // 3. narrowphase, ee
    for (auto &info : ee_broadphase)
    {
        auto obj0 = mColObjs[info.mObj0Id];
        auto obj1 = mColObjs[info.mObj1Id];
        bool obj0_is_cloth_obj1_is_rb = false;
        DetermineObj0AndObj1(obj0, obj1, obj0_is_cloth_obj1_is_rb);
        const auto &e0_array = obj0->GetEdgeArray();
        const auto &e1_array = obj1->GetEdgeArray();
        const auto &v0_array = obj0->GetVertexArray();
        const auto &v1_array = obj1->GetVertexArray();
        const auto &tri0_array = obj0->GetTriangleArray();
        const auto &tri1_array = obj1->GetTriangleArray();
        float dist = 0;
        for (auto &ee : info.mEdgePair)
        {
            int e0 = ee.first;
            int e1 = ee.second;
            dist = 0;
            auto edge0 = e0_array[e0];
            auto edge1 = e1_array[e1];
            // printf("[ee narrow] %d %d %d %d\n", edge0->mId0, edge0->mId1,
            //        edge1->mId0, edge1->mId1);
            const tVector3d &p0 = v0_array[edge0->mId0]->mPos.segment(0, 3);
            const tVector3d &p1 = v0_array[edge0->mId1]->mPos.segment(0, 3);
            const tVector3d &p2 = v1_array[edge1->mId0]->mPos.segment(0, 3);
            const tVector3d &p3 = v1_array[edge1->mId1]->mPos.segment(0, 3);
            float barya = 0, baryb = 0;
            {
                // 1. calculate edge-edge distance

                /*
                    [p1, q1, p2, q2]
                    c1 = p1 + (q1 - p1) * s = (1 - s) * p1 + s * q1
                    c2 = p2 + (q2 - p2) * t = (1 - t) * p2 + t * q2
                */
                EdgeEdgeDist(p0, p1, p2, p3, dist, barya, baryb);
                tVector3d point_on_e0 = (1 - barya) * p0 + barya * p1;
                tVector3d point_on_e1 = (1 - baryb) * p2 + baryb * p3;
                tVector3d dist_vec = tVector3d::Zero();
                tVector3d outer_normal = tVector3d::Zero();
                if (obj0_is_cloth_obj1_is_rb == true)
                {
                    // obj0 is cloth
                    // obj1 is rigid body, take it as the outer normal direction
                    // (the averge of triangle normal)

                    // edge cross
                    outer_normal = (p1 - p0).cross(p3 - p2).normalized();
                    {
                        // get
                        tVector3d obj_outer_normal =
                            tri1_array[edge1->mTriangleId0]->mNormal.segment(0,
                                                                             3);
                        if (edge1->mTriangleId1 != -1)
                        {
                            obj_outer_normal += tri1_array[edge1->mTriangleId1]
                                                    ->mNormal.segment(0, 3);
                        }
                        obj_outer_normal.normalize();
                        if (obj_outer_normal.dot(outer_normal) < 0)
                            outer_normal *= -1;
                    }
                    // tri1_array[edge1->mTriangleId0]->mNormal.segment(0, 3);

                    // outer_normal.normalize();

                    // from rigid body to cloth
                    // from obj1 to obj0
                    dist_vec = point_on_e0 - point_on_e1;
                }
                else
                {
                    // obj0 is rigid body, take it as the outer normal

                    outer_normal = (p1 - p0).cross(p3 - p2).normalized();
                    {
                        // get
                        tVector3d obj_outer_normal =
                            tri0_array[edge0->mTriangleId0]->mNormal.segment(0,
                                                                             3);
                        if (edge0->mTriangleId1 != -1)
                        {
                            obj_outer_normal += tri0_array[edge0->mTriangleId1]
                                                    ->mNormal.segment(0, 3);
                        }
                        obj_outer_normal.normalize();
                        if (obj_outer_normal.dot(outer_normal) < 0)
                            outer_normal *= -1;
                    }

                    // outer_normal =
                    //     tri0_array[edge0->mTriangleId0]->mNormal.segment(0,
                    //     3);
                    // if (edge0->mTriangleId1 != -1)
                    // {
                    //     outer_normal +=
                    //         tri0_array[edge0->mTriangleId1]->mNormal.segment(0,
                    //                                                          3);
                    // }
                    // outer_normal.normalize();

                    // from obj0 to obj1
                    dist_vec = point_on_e1 - point_on_e0;
                }
                double ee_dist = dist_vec.norm();
                double ee_dist_normal_dir = dist_vec.dot(outer_normal);
                // add it!
                // std::cout << "ee_dist = " << ee_dist << std::endl;
                if (ee_dist < 3e-3)
                {
                    auto info = std::make_shared<tEdgeEdgeCollisionInfo>();
                    info->mObj0 = obj0;
                    info->mObj1 = obj1;
                    info->mOuterNormal = outer_normal;
                    info->mDepth = ee_dist_normal_dir;
                    info->obj0_is_cloth_obj1_is_rb = obj0_is_cloth_obj1_is_rb;

                    info->mEdgeId0 = e0;
                    info->mEdgeId1 = e1;
                    info->mBary0 = barya;
                    info->mBary1 = baryb;
                    // printf("[ee] e%d and e%d intersect, depth %.3f\n", e0,
                    // e1,
                    //        ee_dist_normal_dir);
                    mEdgeEdgeCollisionInfo.push_back(info);
                    mObjEdgeEdgeCollisionInfo[obj0->GetObjId()].push_back(info);
                    mObjEdgeEdgeCollisionInfo[obj1->GetObjId()].push_back(info);
                }
            }
        }
    }
    std::cout << "narrow phase ee result = " << mEdgeEdgeCollisionInfo.size()
              << std::endl;
    cTimeUtil::End("narrowphase");
}

void cBVHCollisionDetecter::Clear() {}

void cBVHCollisionDetecter::Update() { UpdateBVHAABB(); }

void cBVHCollisionDetecter::UpdateBVHAABB()
{
    // cTimeUtil::Begin("update BVH AABB");
    for (auto &bvh : mBVHList)
    {
        // if the object is static, do not update it
        bvh->UpdateAABB();
        // if (bvh->GetNumOfLeaves() < 10)
        // {
        //     bvh->Print();
        // }
    }
}
