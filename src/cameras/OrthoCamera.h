#pragma once
#include "CameraBase.h"

class cOrthoCamera : public CameraBase
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    cOrthoCamera(const tVector3f &pos, const tVector3f &center,
                 const tVector3f &up, float init_box, float near_plane_dist,
                 float far_plane_dist);
    virtual ~cOrthoCamera();

    virtual void MoveForward() override;
    virtual void MoveBackward() override;
    virtual void MoveLeft() override;
    virtual void MoveRight() override;
    virtual void MoveUp() override;
    virtual void MoveDown() override;

    virtual void MouseMove(float mouse_x, float mouse_y) override;
    virtual tMatrix4f ProjMatrix(float screen_width, float screen_height,
                                 bool is_vulkan = false) const override;
    virtual tVector CalcCursorPointWorldPos(double xpos, double ypos,
                                            int height, int width) override;

protected:
    float mInitBoxSize; // init box width
};