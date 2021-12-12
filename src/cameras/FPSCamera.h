//
// Created by Hanke on 2019-01-31.
//

#ifndef ROBOT_CAMERA_H
#define ROBOT_CAMERA_H
#include "CameraBase.h"

class FPSCamera : public CameraBase
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    FPSCamera(const tVector3f &pos, const tVector3f &center,
              const tVector3f &up, float fov, float near_plane, float far_plane);
    virtual ~FPSCamera();

    virtual void MoveForward() override;
    virtual void MoveBackward() override;
    virtual void MoveLeft() override;
    virtual void MoveRight() override;
    virtual void MoveUp() override;
    virtual void MoveDown() override;
    virtual void MouseMove(float mouse_x, float mouse_y) override;

protected:
    void Init();
    float pitch, yaw;
};

#endif // ROBOT_CAMERA_H
