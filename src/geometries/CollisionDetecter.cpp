#include "CollisionDetecter.h"
#include "geometries/CollisionInfo.h"
#include "sim/BaseObject.h"
#include "sim/KinematicBody.h"
#include "utils/LogUtil.h"
#include "utils/TimeUtil.hpp"
#include <iostream>
cCollisionDetecter::cCollisionDetecter() {}
cCollisionDetecter::~cCollisionDetecter() {}

void cCollisionDetecter::AddObject(cBaseObjectPtr obj,
                                   bool enable_self_collision /* = false*/)
{
    mColObjs.push_back(obj);
    mEnableSelfCollision.push_back(enable_self_collision);
}

void cCollisionDetecter::Init() { mInited = true; }

void cCollisionDetecter::PerformCD()
{
    if (mInited == false)
    {
        SIM_ERROR("the collision detecter hasn't been inited");
        exit(1);
    }
    cTimeUtil::Begin("DCD");
    // ! 1. clear all buffers
    Clear();

    // ! 2. broadphase collision detection
    BroadphaseCD();

    // ! 3. narrowphase collision detection
    NarrowphaseCD();
    cTimeUtil::End("DCD");
}

void cCollisionDetecter::Clear()
{
    mContactPoints.clear();

    mColCandiadatePairs.clear();
}

std::vector<tColPointPtr> cCollisionDetecter::GetContactPoints() const
{
    return mContactPoints;
}

/**
 * \brief           do broadphase collision (judge AABB intersection)
 */
void cCollisionDetecter::BroadphaseCD()
{
    SIM_ASSERT(this->mColCandiadatePairs.size() == 0);
    // 1. calculate AABB
    int num_of_obj = mColObjs.size();
    tEigenArr<tVector> mAABBmin(num_of_obj), mAABBmax(num_of_obj);

    for (int i = 0; i < num_of_obj; i++)
    {
        mColObjs[i]->CalcAABB(mAABBmin[i], mAABBmax[i]);
    }

    // 2. compare
    for (int i = 0; i < num_of_obj; i++)
    {
        auto obj_type_i = mColObjs[i]->GetObjectType();
        for (int j = i + 1; j < num_of_obj; j++)
        {
            if (true == cMathUtil::IntersectAABB(mAABBmin[i], mAABBmax[i],
                                                 mAABBmin[j], mAABBmax[j]))
            {
                auto obj_type_j = mColObjs[j]->GetObjectType();
                if (obj_type_i <= obj_type_j)
                {

                    mColCandiadatePairs.push_back(tVector2i(i, j));
                }
                else
                {
                    mColCandiadatePairs.push_back(tVector2i(j, i));
                }
                // printf("[debug] broadphse %d and %d collided\n", i, j);
            }
        }
    }
}

/**
 * \brief           do narrowphase collision (no self collision)
 */
void cCollisionDetecter::NarrowphaseCD()
{
    // ! if the broadphase is empty, return
    if (mColCandiadatePairs.size() == 0)
        return;

    // ! begin to check eacy candidate collided obj pair
    for (auto &pair : mColCandiadatePairs)
    {
        int obj0_id = pair[0], obj1_id = pair[1];
        eObjectType type0 = mColObjs[obj0_id]->GetObjectType(),
                    type1 = mColObjs[obj1_id]->GetObjectType();
    }
}
