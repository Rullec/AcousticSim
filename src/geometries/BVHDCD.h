#pragma once
#include "CollisionDetecter.h"
class cBVHDCD : public cCollisionDetecter
{
public:
    explicit cBVHDCD();
    cBVHDCD(const cBVHDCD &) = delete;
    virtual void AddObject(cBaseObjectPtr obj,
                        bool enable_self_collision = false) override;
    virtual void Init();
    virtual void PerformCD();
    virtual void Clear();
    virtual std::vector<tColPointPtr> GetContactPoints() const override;
}