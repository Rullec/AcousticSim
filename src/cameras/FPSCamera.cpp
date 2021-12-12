//
// Created by Hanke on 2019-01-31.
//

#include "FPSCamera.h"
#define angle2rad(deg_value) ((deg_value)*gDegreesToRadians)
#define rad2angle(rad_value) ((rad_value)*gRadiansToDegrees)

// FPSCamera::FPSCamera() : CameraBase()
// {
//     mType = eCameraType::FPS_CAMERA;
//     Init();
// }
FPSCamera::~FPSCamera() {}
FPSCamera::FPSCamera(const tVector3f &pos, const tVector3f &center,
                     const tVector3f &up, float fov, float near_plane,
                     float far_plane)
    : CameraBase(pos, center, up, fov, near_plane, far_plane)
{
    mType = eCameraType::FPS_CAMERA;

    pitch = rad2angle(asin(mCamFront(1)));
    float val = asin(mCamFront(0) / cos(angle2rad(pitch)));
    if (val < -1)
        val = -1;
    else if (val > 1)
        val = 1;
    yaw = rad2angle(asin(val));
}

void FPSCamera::MoveForward()
{
    mCamFront.normalize();
    mCamPos += mCamFront * key_acc;
    mCamCenter = mCamPos + mCamFront;
}

void FPSCamera::MoveBackward()
{
    mCamFront.normalize();
    mCamPos -= mCamFront * key_acc;
    mCamCenter = mCamPos + mCamFront;
}

void FPSCamera::MoveLeft()
{
    tVector3f left = mCamUp.cross(mCamFront);
    left.normalize();
    mCamPos += left * key_acc;
    mCamCenter = mCamPos + mCamFront;
}

void FPSCamera::MoveRight()
{
    tVector3f left = -mCamUp.cross(mCamFront);
    left.normalize();
    mCamPos += left * key_acc;
    mCamCenter = mCamPos + mCamFront;
}

void FPSCamera::MoveUp()
{
    mCamPos += mCamUp * key_acc;
    mCamCenter = mCamPos + mCamFront;
}

void FPSCamera::MoveDown()
{
    mCamPos -= mCamUp * key_acc;
    mCamCenter = mCamPos + mCamFront;
}

void FPSCamera::MouseMove(float mouse_x, float mouse_y)
{
    if (first_mouse)
    {
        last_x = mouse_x;
        last_y = mouse_y;
        first_mouse = false;
        // std::cout << "fps first mouse, now mouse = " << mouse_x << " "
        //   << mouse_y << std::endl;
        return;
    }
    // std::cout << "fps move mouse, now mouse = " << mouse_x << " " << mouse_y
    //           << std::endl;
    float x_offset = mouse_x - last_x;
    float y_offset = -mouse_y + last_y;
    // std::cout << "x offset = " << x_offset << " y offset = " << y_offset
    //           << std::endl;
    // // std::cout << "raw yaw = " << yaw << std::endl;
    // // std::cout << "raw picth = " << pitch << std::endl;
    last_x = mouse_x;
    last_y = mouse_y;
    x_offset *= mouse_acc;
    y_offset *= mouse_acc;
    yaw += x_offset;
    pitch += y_offset;

    if (pitch > 89.0f)
        pitch = 89.0f;
    if (pitch < -89.0f)
        pitch = -89.0f;

    // std::cout << "new yaw = " << yaw << std::endl;
    // std::cout << "new picth = " << pitch << std::endl;
    // std::cout << "old front = " << mCamFront.transpose() << std::endl;

    /*

        pitch = rad2angle(asin(mCamFront(1)));
        yaw = rad2angle(asin(asin(mCamFront(0) / cos(angle2rad(pitch)))));

    */
    mCamFront(0) = cos(angle2rad(pitch)) * sin(angle2rad(yaw));
    mCamFront(1) = sin(angle2rad(pitch));
    mCamFront(2) = -cos(angle2rad(pitch)) * cos(angle2rad(yaw));
    mCamFront.normalize();
    // std::cout << "new front = " << mCamFront.transpose() << std::endl;
    mCamCenter = mCamPos + mCamFront;

    // std::cout << "----------------------------\n";
}

void FPSCamera::Init()
{
    mCamPos = tVector3f(0.8, 1, 1.8);
    mCamCenter = tVector3f(0, 0, 0.2);
    mCamUp = tVector3f(0, 1, 0);
    mCamFront = mCamCenter - mCamPos;
    mCamFront.normalize();
    pitch = rad2angle(asin(mCamFront(1)));
    yaw = rad2angle(asin(asin(mCamFront(0) / cos(angle2rad(pitch)))));
}
