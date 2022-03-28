#pragma once
#include "CollisionDetecter.h"

struct tPointTriangleInfo
{
    tPointTriangleInfo(int obj0, int obj1);
    int GetNum() const { return mPointTrianglePair.size(); }
    int mObj0Id, mObj1Id;
    std::vector<std::pair<int, int>> mPointTrianglePair; // point id in obj0
    tEigenArr<tVector3f> mBary; // barycentric coords of triangles
};

struct tEdgeEdgeInfo
{
    tEdgeEdgeInfo(int obj0, int obj1);
    int GetNum() const { return mEdgePair.size(); }
    int mObj0Id, mObj1Id;
    std::vector<std::pair<int, int>> mEdgePair;
    std::vector<float> mBaryA,
        mBaryB; // contact point barycentric coords on edgeA and edgeB;
};

SIM_DECLARE_STRUCT_AND_PTR(tPointTriangleCollisionInfo);
SIM_DECLARE_STRUCT_AND_PTR(tEdgeEdgeCollisionInfo);
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

protected:
    std::vector<cObjBVHPtr> mBVHList;
    virtual void UpdateBVHAABB();

    // tPointTriangleInfo CheckBetweenTwoBVH(cObjBVHPtr t1, cObjBVHPtr t2)
    // const;
};