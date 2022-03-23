#pragma once
#include "CollisionDetecter.h"

struct tPointTriangleInfo
{
    tPointTriangleInfo(int obj0, int obj1);
    int GetNum() const { return mPointTrianglePair.size(); }
    int mObj0Id, mObj1Id;
    std::vector<std::pair<int, int>> mPointTrianglePair; // point id in obj0
};

SIM_DECLARE_STRUCT_AND_PTR(tPointTriangleCollisionInfo);
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
    virtual std::vector<tPointTriangleCollisionInfoPtr>
    GetAllPointTriangleCollisionInfo() const override;
    std::vector<tPointTriangleCollisionInfoPtr>
    GetObjPointTriangleCollisionInfo(int obj_id) const override;

protected:
    std::vector<cObjBVHPtr> mBVHList;
    std::vector<std::vector<tPointTriangleCollisionInfoPtr>>
        mObjPointTriangleInfo;
    virtual void UpdateBVHAABB();

    // tPointTriangleInfo CheckBetweenTwoBVH(cObjBVHPtr t1, cObjBVHPtr t2) const;
};