#pragma once
#include "CollisionDetecter.h"

struct tBroadphaseInfo
{
    tBroadphaseInfo();
    int mObj0Id, mObj1Id;
    std::vector<std::pair<int, int>> mTrianglePairInfo;
};

SIM_DECLARE_CLASS_AND_PTR(cObjBVH);
class cBVHCollisionDetecter : public cCollisionDetecter
{
public:
    explicit cBVHCollisionDetecter();
    virtual ~cBVHCollisionDetecter();
    virtual void AddObject(cBaseObjectPtr obj,
                           bool enable_self_collision = false) override;
    virtual void Init() override;
    virtual void Update() override;
    virtual void PerformCD() override;
    virtual void Clear() override;
    virtual std::vector<tColPointPtr> GetContactPoints() const override;

protected:
    std::vector<cObjBVHPtr> mBVHList;
    virtual void UpdateBVHAABB();

    tBroadphaseInfo CheckBetweenTwoBVH(cObjBVHPtr t1, cObjBVHPtr t2) const;
};