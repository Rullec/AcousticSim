#pragma once
//
// Extract to base by Xudong on 2021-01-30
//

// #include "../Utils/EigenUtils.h"
#include "utils/MathUtil.h"
#include <iostream>
#include <memory>
// #include <cmath>

enum eCameraType
{
    FPS_CAMERA = 0, // default wasd FPS camear
    ARCBALL_CAMERA, // arcball camera
    ORTHO_CAMERA,   // orthogonal camera
    NUM_OF_CAMERA_TYPE
};

class CameraBase : std::enable_shared_from_this<CameraBase>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    CameraBase(const tVector3f &pos, const tVector3f &center,
               const tVector3f &up, float fov, float near_plane,
               float far_plane);
    virtual ~CameraBase() = 0;
    eCameraType GetType() const;
    static eCameraType BuildCameraTypeFromStr(std::string name);
    virtual tMatrix4f ViewMatrix();
    virtual tMatrix4f ProjMatrix(float screen_width, float screen_height,
                                 bool is_vulkan = false) const;
    virtual tVector3f GetCameraPos() const;
    virtual tVector3f GetCameraCenter() const;
    virtual tVector3f GetCameraFront() const;
    virtual tVector3f GetCameraUp() const;
    virtual float GetCameraFovDeg() const;
    virtual void MoveForward() = 0;
    virtual void MoveBackward() = 0;
    virtual void MoveLeft() = 0;
    virtual void MoveRight() = 0;
    virtual void MoveUp() = 0;
    virtual void MoveDown() = 0;

    virtual void MouseMove(float mouse_x, float mouse_y) = 0;
    virtual void MiddleKeyMove(float mouse_x, float mouse_y);
    virtual void SetXY(float mouse_x, float mouse_y);
    virtual void ResetFlag();

    virtual void SetKeyAcc(double acc);
    virtual void SetMouseAcc(double acc);
    virtual bool IsFirstMouse() const;

    virtual tVector CalcCursorPointWorldPos(double xpos, double ypos,
                                            int height, int width);
    virtual void RotateAlongYAxis(double angle_deg);
protected:
    eCameraType mType;
    tVector3f mCamPos, mCamCenter, mCamUp, mCamFront;
    float mNearPlane, mFarPlane;
    float mFovDeg;
    float mouse_acc;
    float key_acc;

    float last_x, last_y; // the previous x and y position of mouse
    bool first_mouse;     // whether it's the first mouse event
};