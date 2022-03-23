#pragma once
#include "utils/DefUtil.h"
#include "utils/MathUtil.h"

SIM_DECLARE_CLASS_AND_PTR(cBaseObject);
SIM_DECLARE_CLASS_AND_PTR(cKinematicBody);
SIM_DECLARE_CLASS_AND_PTR(tColPoint);
SIM_DECLARE_STRUCT_AND_PTR(tPointTriangleCollisionInfo);
class cCollisionDetecter
    : public std::enable_shared_from_this<cCollisionDetecter>
{
public:
    cCollisionDetecter();
    virtual ~cCollisionDetecter();
    virtual void AddObject(cBaseObjectPtr obj,
                           bool enable_self_collision = false);
    virtual void Update();
    virtual void Init();
    virtual void PerformCD();
    virtual void Clear();
    virtual std::vector<tPointTriangleCollisionInfoPtr>
    GetAllPointTriangleCollisionInfo() const;
    virtual std::vector<tPointTriangleCollisionInfoPtr>
    GetObjPointTriangleCollisionInfo(int obj_id) const;

protected:
    bool mInited = false;
    // permanet info
    std::vector<cBaseObjectPtr> mColObjs;   // collision objects
    std::vector<bool> mEnableSelfCollision; // enable self collision or not
    std::vector<tPointTriangleCollisionInfoPtr> mPointTriangleCollisionInfo;
    // buffers (need to be cleared)
    // std::vector<tColPointPtr> mContactPoints; // contact points
    // tEigenArr<tVector2i> mColCandiadatePairs; // the possible collision obj
    //                                           // pair, after broad phase
    virtual void BroadphaseCD();
    virtual void NarrowphaseCD();
};