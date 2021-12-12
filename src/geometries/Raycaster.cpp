#include "Raycaster.h"
#include "utils/JsonUtil.h"
#include "utils/LogUtil.h"
#include <iostream>

cRaycaster::tRaycastResult::tRaycastResult()
{
    mObject = nullptr;
    mLocalTriangleId = -1;
    mIntersectionPoint =
        tVector::Ones() * std::numeric_limits<double>::quiet_NaN();
}

cRaycaster::cRaycaster()
{
    // mEnableOnlyExportCuttedWindow = false;
    mTriangleArray_lst.clear();
    mVertexArray_lst.clear();
    // SIM_INFO("Build raycaster succ");
    // SIM_ASSERT(triangles != nullptr);
}

void cRaycaster::Init(const Json::Value &conf)
{
    
}
// void cRaycaster::AddResources(const std::vector<tTriangle *> triangles,
//                               const std::vector<tVertex *> vertices)
// {
//     mTriangleArray_lst.push_back(triangles);
//     mVertexArray_lst.push_back(vertices);
// }
#include "sim/BaseObject.h"
void cRaycaster::AddResources(cBaseObjectPtr object)
{
    mObjects.push_back(object);

    std::vector<tTriangle *> tri = object->GetTriangleArray();
    std::vector<tVertex *> ver = object->GetVertexArray();

    mTriangleArray_lst.push_back(tri);
    // std::cout << mTriangleArray_lst.size() << std::endl;
    mVertexArray_lst.push_back(ver);
    // std::cout << mVertexArray_lst.size() << std::endl;
}

/**
 * \brief               Do raycast, calculate the intersection
 */
// void cRaycaster::RayCast(const tRay *ray, tTriangle **selected_tri,
//                          int &selected_tri_id, tVector &raycast_point) const
cRaycaster::tRaycastResult cRaycaster::RayCast(const tRay *ray) const
{
    // 1. init
    tRaycastResult res;
    double min_depth = std::numeric_limits<double>::max();
    // raycast_point.noalias() = tVector::Ones() * std::nan("");
    // std::cout << "triangle array lst size = " << mTriangleArray_lst.size()
    //           << std::endl;
    // 2. iterate on each triangle
    for (int obj_id = 0; obj_id < mTriangleArray_lst.size(); obj_id++)
    {
        auto mTriangleArray = mTriangleArray_lst[obj_id];
        auto mVertexArray = mVertexArray_lst[obj_id];
        for (int i = 0; i < mTriangleArray.size(); i++)
        {
            auto &tri = mTriangleArray[i];
            tVector tmp = cMathUtil::RayCastTri(
                ray->mOrigin, ray->mDir, mVertexArray[tri->mId0]->mPos,
                mVertexArray[tri->mId1]->mPos, mVertexArray[tri->mId2]->mPos);

            // if there is an intersection, tmp will have no nan
            if (tmp.hasNaN() == false)
            {
                // std::cout << tmp.transpose() << std::endl;
                double cur_depth = (tmp - ray->mOrigin).segment(0, 3).norm();
                if (cur_depth < min_depth)
                {
                    min_depth = cur_depth;
                    res.mObject = mObjects[obj_id];
                    res.mLocalTriangleId = i;
                    res.mIntersectionPoint = tmp;
                }
            }
        }
    }
    return res;
}

/**
 * \brief           Calculate the depth image
 */
void cRaycaster::CalcDepthMap(const tMatrix2i &cast_range, int height,
                              int width, CameraBasePtr camera, std::string)
{
    SIM_ASSERT("cRaycaster::CalcDepthMap hasn't been finished yet");
}

/**
 * \brief           Calculate the interested window size
 */
void cRaycaster::CalcCastWindowSize(const tMatrix2i &cast_range_window,
                                    int &window_width, int &window_height,
                                    tVector2i &st)
{
    window_width = cast_range_window(0, 1) - cast_range_window(0, 0);
    window_height = cast_range_window(1, 1) - cast_range_window(1, 0);

    st = cast_range_window.col(0).transpose();
}