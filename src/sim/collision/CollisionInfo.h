// #pragma once
// #include "utils/DefUtil.h"
// #include "utils/MathUtil.h"
// #include <memory>
// SIM_DECLARE_CLASS_AND_PTR(cKinematicBody);

// // ! Collision Info for each type of object
// SIM_DECLARE_CLASS_AND_PTR(cBaseObject);
// struct tColObjInfo : public std::enable_shared_from_this<tColObjInfo>
// {
//     tColObjInfo(cBaseObjectPtr obj);
//     virtual ~tColObjInfo() = 0;
//     cBaseObjectPtr mObj;
// };

// SIM_DECLARE_PTR(tColObjInfo);

// // ! Collision Info for rigidbody / kinematic body (undeformed body)
// struct tColRigidBodyInfo : public tColObjInfo
// {
//     tColRigidBodyInfo(cKinematicBodyPtr, const tVector pos);
//     virtual ~tColRigidBodyInfo();
//     tVector mLocalPos;
// };

// // // ! Collision Info for fluid
// // struct tCollisionFluidInfo
// // {
// //     tCollisionFluidInfo();
// //     virtual ~tCollisionFluidInfo();
// // };

// // ! Collision point for each pair of object
// struct tColPoint : std::enable_shared_from_this<tColPoint>
// {
//     tColPoint();
//     EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
//     tColObjInfoPtr mObjInfo0, mObjInfo1;
//     tVector mNormal; // from obj0 to obj1
//     double mDepth;   // minus = penetration
// };

// SIM_DECLARE_PTR(tColPoint);

#include "utils/DefUtil.h"
#include "utils/EigenUtil.h"
#include "utils/LogUtil.h"
SIM_DECLARE_CLASS_AND_PTR(cBaseObject);
struct tClothHumanCollisionInfo
{
    cBaseObjectPtr mObj0, mObj1;

    tVector3d mOuterNormal;
    bool obj0_is_cloth_obj1_is_rb;
    float mDepth;   // penetrated depth; <0 means penetration
};
struct tPointTriangleCollisionInfo : public tClothHumanCollisionInfo
{
    int mVertexId0, mTriangleId1;
    tVector3d mBary;
};
struct tEdgeEdgeCollisionInfo : public tClothHumanCollisionInfo
{
    int mEdgeId0, mEdgeId1;
    float mBary0, mBary1; // edge closest point's barycentric coords
};
SIM_DECLARE_PTR(tPointTriangleCollisionInfo);