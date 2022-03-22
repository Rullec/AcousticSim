#include "BVHCollisionDetecter.h"
#include "geometries/ObjectBVH.h"
#include "sim/BaseObject.h"
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
        bvh->Print();
        mBVHList.push_back(bvh);
    }
}
void cBVHCollisionDetecter::PerformCD() { cCollisionDetecter::PerformCD(); }
void cBVHCollisionDetecter::Clear() {}
std::vector<tColPointPtr> cBVHCollisionDetecter::GetContactPoints() const {
    return {};
}