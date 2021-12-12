#include "CameraBuilder.h"
#include "ArcBallCamera.h"
#include "FPSCamera.h"
#include "OrthoCamera.h"
#include <iostream>

CameraBase *cCameraBuilder::BuildCamera(const Json::Value &camera_json)
{
    CameraBase *cam = nullptr;
    std::string type_str = cJsonUtil::ParseAsString("camera_type", camera_json);
    Json::Value camera_pos_json =
        cJsonUtil::ParseAsValue("camera_pos", camera_json);
    Json::Value camera_focus_json =
        cJsonUtil::ParseAsValue("camera_focus", camera_json);
    SIM_ASSERT(camera_pos_json.size() == 3);
    SIM_ASSERT(camera_focus_json.size() == 3);
    tVector3f mCameraInitFocus = tVector3f(camera_focus_json[0].asFloat(),
                                           camera_focus_json[1].asFloat(),
                                           camera_focus_json[2].asFloat());
    tVector3f mCameraInitPos =
        tVector3f(camera_pos_json[0].asFloat(), camera_pos_json[1].asFloat(),
                  camera_pos_json[2].asFloat());

    float near_plane_dist = cJsonUtil::ParseAsFloat("near", camera_json);
    float far_plane_dist = cJsonUtil::ParseAsFloat("far", camera_json);
    float mCameraInitFov = 0;
    if (camera_json.isMember("fov") == true)
    {
        mCameraInitFov = cJsonUtil::ParseAsFloat("fov", camera_json);
    }

    // cam = new cArcBallCamera(mCameraInitPos, mCameraInitFocus,
    //                          tVector3f(0, 1, 0), mCameraInitFov);
    eCameraType type = CameraBase::BuildCameraTypeFromStr(type_str);
    switch (type)
    {
    case eCameraType::FPS_CAMERA:
    {
        cam =
            new FPSCamera(mCameraInitPos, mCameraInitFocus, tVector3f(0, 1, 0),
                          mCameraInitFov, near_plane_dist, far_plane_dist);

        break;
    }
    case eCameraType::ARCBALL_CAMERA:
    {
        cam = new cArcBallCamera(mCameraInitPos, mCameraInitFocus,
                                 tVector3f(0, 1, 0), mCameraInitFov,
                                 near_plane_dist, far_plane_dist);
        break;
    }
    case eCameraType::ORTHO_CAMERA:
    {

        float init_box =
            cJsonUtil::ParseAsFloat("camera_init_box", camera_json);
        cam = new cOrthoCamera(mCameraInitPos, mCameraInitFocus,
                               tVector3f(0, 1, 0), init_box, near_plane_dist,
                               far_plane_dist);
        break;
    }
    default:
        SIM_ERROR("unsupported camera type {}", type_str);
        exit(1);
        break;
    }
    return cam;
}