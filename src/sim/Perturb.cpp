#include "Perturb.h"
#include "geometries/Primitives.h"
#include "utils/LogUtil.h"
#include <iostream>
tPerturb::tPerturb()
{
    mPerturbForce.setZero();
    mObject = nullptr;
    mAffectedTriId = -1;
    mBarycentricCoords.setZero();
    mShiftPlaneEquation.setZero();
}

void tPerturb::InitTangentRect(const tVector &plane_normal)
{
    // 1. calculate center pos
    tVector center_pos = CalcPerturbPos();
    mShiftPlaneEquation.segment(0, 3) = plane_normal.segment(0, 3);
    mShiftPlaneEquation[3] =
        -plane_normal.segment(0, 3).dot(center_pos.segment(0, 3));
    // std::cout << "[calc] center pos = " << center_pos.transpose() <<
    // std::endl; std::cout << "[calc] plane normal = " <<
    // plane_normal.transpose() << std::endl; std::cout << "[calc] shift plane
    // equation = " << mShiftPlaneEquation.transpose() << std::endl; std::cout
    // << "center pos = " << center_pos.transpose() << std::endl;
    // 2. calculate rectangle vertices
    // each row is a vector
    // tMatrixXd four_vectors = cMathUtil::ExpandFrictionCone(4, plane_normal);
    // // std::cout << "four vecs = \n" << four_vectors << std::endl;
    // double plane_scale = 1e10;
    // four_vectors *= plane_scale;
    // for (int i = 0; i < 4; i++)
    // {
    //     // std::cout << 1 << std::endl;
    //     four_vectors.row(i) += center_pos;
    //     // std::cout << 2 << std::endl;
    //     mRectVertices[i] = new tVertex();
    //     // std::cout << 3 << std::endl;
    //     mRectVertices[i]->mPos = four_vectors.row(i);
    //     // std::cout << 4 << std::endl;
    //     // std::cout << "rect v" << i << " = "
    //     //           << mRectVertices[i]->mPos.transpose() << std::endl;
    // }
}

void tPerturb::UpdatePerturbPos(const tVector &cur_camera_pos, const tVector &dir)
{
    tVector target_goal =

        cMathUtil::RayCastPlane(cur_camera_pos, dir, mShiftPlaneEquation);
    // std::cout << "[perturb] shift plane equation = " <<
    // mShiftPlaneEquation.transpose() << std::endl; std::cout << "[perturb] ray
    // ori = " << cur_camera_pos.transpose() << std::endl; std::cout <<
    // "[perturb] ray dir = " << dir.transpose() << std::endl; std::cout <<
    // "[perturb] inter pos = " << mGoalPos.transpose() << std::endl;
    if (target_goal.hasNaN() == true)
    {
        SIM_WARN("ray has no intersection with the rect range, please "
                  "scale the plane");
        return;
    }
    mGoalPos = target_goal;
    tVector cur_perturb_pos = CalcPerturbPos();
    mPerturbForce = (mGoalPos - cur_perturb_pos) * 10;
    // std::cout << "[update] cur perturb pos = " << cur_perturb_pos.transpose()
    //           << " goal pos = " << goal_pos.transpose() << std::endl;
}
#include "sim/BaseObject.h"
tVector tPerturb::CalcPerturbPos() const
{
    tVector center_pos = tVector::Zero();
    // std::cout << "[calc perurb pos] bary = " <<
    // mBarycentricCoords.transpose() << std::endl;
    auto tri = mObject->GetTriangleArray()[this->mAffectedTriId];
    auto ver_array = mObject->GetVertexArray();
    for (int i = 0; i < 3; i++)
    {
        center_pos += ver_array[tri->mId0]->mPos * mBarycentricCoords[i];
    }
    center_pos[3] = 1;
    return center_pos;
}

tVector tPerturb::GetPerturbForce() const { return mPerturbForce; }

tVector tPerturb::GetGoalPos() const { return this->mGoalPos; }