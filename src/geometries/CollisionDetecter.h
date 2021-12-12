#pragma once
#include "utils/DefUtil.h"
#include "utils/MathUtil.h"

SIM_DECLARE_CLASS_AND_PTR(cBaseObject);
SIM_DECLARE_CLASS_AND_PTR(cKinematicBody);
SIM_DECLARE_CLASS_AND_PTR(tColPoint);

class cCollisionDetecter
{
public:
    cCollisionDetecter();
    virtual ~cCollisionDetecter();
    virtual void AddObject(cBaseObjectPtr obj,
                           bool enable_self_collision = false);
    virtual void PerformCD();
    virtual void Clear();
    virtual std::vector<tColPointPtr> GetContactPoints() const;

protected:
    // permanet info
    std::vector<cBaseObjectPtr> mColObjs;   // collision objects
    std::vector<bool> mEnableSelfCollision; // enable self collision or not

    // buffers (need to be cleared)
    std::vector<tColPointPtr> mContactPoints; // contact points
    tEigenArr<tVector2i> mColCandiadatePairs; // the possible collision obj
                                              // pair, after broad phase
    virtual void BroadphaseCD();
    virtual void NarrowphaseCD();
};