#pragma once
#include "CollisionDetecter.h"

SIM_DECLARE_CLASS_AND_PTR(cObjBVH);
class cBVHCollisionDetecter : public cCollisionDetecter
{
public:
    explicit cBVHCollisionDetecter();
    virtual ~cBVHCollisionDetecter();
    virtual void AddObject(cBaseObjectPtr obj,
                           bool enable_self_collision = false) override;
    virtual void Init() override;
    virtual void PerformCD() override;
    virtual void Clear() override;
    virtual std::vector<tColPointPtr> GetContactPoints() const override;

protected:
    std::vector<cObjBVHPtr> mBVHList;
};