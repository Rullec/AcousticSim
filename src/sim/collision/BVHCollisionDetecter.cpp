#include "BVHCollisionDetecter.h"
#include "geometries/ObjectBVH.h"
#include "sim/BaseObject.h"
#include "utils/TimeUtil.hpp"
#include <iostream>
#include <queue>
tBroadphaseInfo::tBroadphaseInfo()
{
    mObj0Id = -1;
    mObj1Id = -1;
    mTrianglePairInfo.clear();
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
        if (IsNodeIntersected(a, b))
        {
            // printf("[debug] node intersected!\n");
            // std::cout << "a min = " << a->mAABB.mMin.transpose()
            //           << " max = " << a->mAABB.mMax.transpose() << std::endl;
            // std::cout << "b min = " << b->mAABB.mMin.transpose()
            //           << " max = " << b->mAABB.mMax.transpose() << std::endl;
            if (a->mIsLeaf && b->mIsLeaf)
            {
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
                    a = a->mRight;
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

void cBVHCollisionDetecter::PerformCD()
{
    // 1. object pair, do collision detection (broadphase)
    cTimeUtil::Begin("broadphase");
    int num_of_obj = mBVHList.size();
    std::vector<tBroadphaseInfo> broadphase_info_lst = {};
    for (int i = 0; i < num_of_obj; i++)
    {
        auto i_bvh = mBVHList[i];
        for (int j = i + 1; j < num_of_obj; j++)
        {
            auto j_bvh = mBVHList[j];

            // check (i, j)
            tBroadphaseInfo info;
            info.mObj0Id = this->mColObjs[i]->GetObjId();
            info.mObj1Id = this->mColObjs[j]->GetObjId();
            info.mTrianglePairInfo =
                IntersectTwoBVHs(i_bvh->GetRootNode(), j_bvh->GetRootNode());
            if (info.mTrianglePairInfo.size())
                broadphase_info_lst.push_back(info);
            // tBroadphaseInfo info = CheckBetweenTwoBVH(i_bvh, j_bvh);
            // broadphase_info_lst.push_back(info);
        }
    }
    cTimeUtil::End("broadphase");

    // 2. do narrow phase, detail the information
    // we also want to calculate the
    for (auto &x : broadphase_info_lst)
    {
        std::cout << "narrowphase need to handle " << x.mTrianglePairInfo.size()
                  << std::endl;
    }
}

void cBVHCollisionDetecter::Clear() {}

std::vector<tColPointPtr> cBVHCollisionDetecter::GetContactPoints() const
{
    return {};
}

void cBVHCollisionDetecter::Update() { UpdateBVHAABB(); }

void cBVHCollisionDetecter::UpdateBVHAABB()
{
    // cTimeUtil::Begin("update BVH AABB");
    for (auto &bvh : mBVHList)
    {
        // if the object is static, do not update it
        bvh->UpdateAABB();
    }
    // mBVHList[0]->Print();
}
