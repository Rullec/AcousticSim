//
// Created by Xudong on 2021/01/30
//
#pragma once
#include "CameraBase.h"

/**
 * \brief           Arcball camera
 */
class cArcBallCamera : public CameraBase
{
public:
    // cArcBallCamera();
    cArcBallCamera(const tVector3f &pos, const tVector3f &center,
                   const tVector3f &up, float fov, float near_plane,
                   float far_plane);

    virtual ~cArcBallCamera();
    virtual tMatrix4f ViewMatrix();
    virtual void MoveForward() override;
    virtual void MoveBackward() override;
    virtual void MoveLeft() override;
    virtual void MoveRight() override;
    virtual void MoveUp() override;
    virtual void MoveDown() override;
    virtual void MiddleKeyMove(float mouse_x, float mouse_y) override;
    virtual void MouseMove(float mouse_x, float mouse_y) override;

protected:
};