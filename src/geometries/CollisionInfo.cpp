#include "geometries/CollisionInfo.h"
#include "sim/BaseObject.h"
#include "sim/KinematicBody.h"
// ! Collision point for each pair of object
tColPoint::tColPoint()
{
    mObjInfo0 = nullptr;
    mObjInfo1 = nullptr;
    mNormal.setZero();
    mDepth = 0;
}

// ! Collision Info for each type of object
tColObjInfo::tColObjInfo(cBaseObjectPtr obj) { mObj = obj; }
tColObjInfo::~tColObjInfo() {}

// ! Collision Info for rigidbody / kinematic body (undeformed body)
tColRigidBodyInfo::tColRigidBodyInfo(cKinematicBodyPtr obj_ptr,
                                     const tVector pos)
    : tColObjInfo(obj_ptr)

{
    mLocalPos = pos;
}
tColRigidBodyInfo::~tColRigidBodyInfo() {}