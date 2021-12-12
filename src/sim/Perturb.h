#pragma once
#include "utils/DefUtil.h"
#include "utils/MathUtil.h"

/**
 * \brief           user perturb force
 */
struct tTriangle;
struct tVertex;
struct tRay;
SIM_DECLARE_CLASS_AND_PTR(cBaseObject);
struct tPerturb
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    tPerturb();
    void InitTangentRect(const tVector &plane_normal);
    void UpdatePerturbPos(const tVector &cur_camera_pos, const tVector &dir);

    tVector GetPerturbForce() const;
    tVector CalcPerturbPos() const;
    tVector GetGoalPos() const;

    cBaseObjectPtr mObject;
    int mAffectedTriId; // triangle id
    // int mAffectedVerticesId[3];
    // tVertex *mAffectedVertices[3];
    tVector3d mBarycentricCoords; // barycentric coordinates of raw raycast
                                  // point on the affected triangle
protected:
    tVector mPerturbForce;
    tVector mShiftPlaneEquation;
    tVector mGoalPos;
};