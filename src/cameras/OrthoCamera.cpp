#include "OrthoCamera.h"
#include "utils/LogUtil.h"

cOrthoCamera::cOrthoCamera(const tVector3f &pos_, const tVector3f &center_,
                           const tVector3f &up_, float init_box,
                           float near_plane_dist, float far_plane_dist)
    : CameraBase(pos_, center_, up_, 0, near_plane_dist, far_plane_dist)

{
    mType = eCameraType::ORTHO_CAMERA;
    mInitBoxSize = init_box;
    SIM_ASSERT(mInitBoxSize > 0);
    mouse_acc *= 5e-2;
}

cOrthoCamera::~cOrthoCamera() {}
void cOrthoCamera::MoveForward() { mInitBoxSize *= 0.9; }
void cOrthoCamera::MoveBackward() { mInitBoxSize /= 0.9; }
void cOrthoCamera::MoveLeft() {}
void cOrthoCamera::MoveRight() {}
void cOrthoCamera::MoveUp() {}
void cOrthoCamera::MoveDown() {}
void cOrthoCamera::MouseMove(float mouse_x, float mouse_y)
{
    if (first_mouse)
    {
        last_x = mouse_x;
        last_y = mouse_y;
        first_mouse = false;
        return;
    }
    float x_offset = (last_x - mouse_x) * mInitBoxSize;
    float y_offset = (-mouse_y + last_y) * mInitBoxSize;
    last_x = mouse_x;
    last_y = mouse_y;
    x_offset *= mouse_acc;
    y_offset *= mouse_acc;

    // std::cout << "x offset = " << x_offset << std::endl;
    // std::cout << "y offset = " << y_offset << std::endl;
    tVector3f x_dir = mCamFront.cross(this->mCamUp);
    tVector3f y_dir = -mCamUp;
    x_dir.normalize();
    y_dir.normalize();
    tVector3f shift = x_dir * x_offset + y_dir * y_offset;
    // std::cout << "world shift = " << shift.transpose() << std::endl;
    mCamPos += shift;
    mCamCenter += shift;
}

tMatrix4f cOrthoCamera::ProjMatrix(float screen_width, float screen_height,
                                   bool is_vulkan /*= false*/) const
{
    float gamma = screen_height / screen_width;

    tMatrix4f mat = tMatrix4f::Zero();
    // ortho: refine
    if (is_vulkan == true)
    {
        mat(0, 0) = 1 / (mInitBoxSize);
        mat(0, 3) = 0;
        mat(1, 1) = -1 / (gamma * mInitBoxSize);
        mat(1, 3) = 0;
        mat(2, 2) = -1 / (mFarPlane - mNearPlane);
        mat(2, 3) = -mNearPlane / (mFarPlane - mNearPlane);
        mat(3, 3) = 1;
    }
    else
    {
        SIM_ERROR("invalid");
    }

    return mat;
}

tVector cOrthoCamera::CalcCursorPointWorldPos(double xpos, double ypos,
                                              int height, int width)
{
    tMatrix mat1 = tMatrix::Identity();
    mat1(0, 0) = 1.0 / width;
    mat1(0, 3) = 0.5 / width;
    mat1(1, 1) = 1.0 / height;
    mat1(1, 3) = 0.5 / height;

    tMatrix mat2 = tMatrix::Identity();
    mat2(0, 0) = 2;
    mat2(0, 3) = -1;
    mat2(1, 1) = -2;
    mat2(1, 3) = 1;

    tMatrix mat3 = ProjMatrix(width, height, true).cast<double>().inverse();
    mat3(1, 1) *= -1;
    tMatrix mat4 = ViewMatrix().cast<double>().inverse();
    tVector pixel_coord = mat2 * mat1 * tVector(xpos, ypos, mNearPlane, 1);
    // std::cout << "screen pos = " << xpos << " " << ypos << std::endl;
    // std::cout << "pixel_coord = " << pixel_coord.transpose() << std::endl;
    // tVector pos_in_view_coords = mat3 * pixel_coord;
    // std::cout << "pos_in_view_coords = " << pos_in_view_coords.transpose()
    //           << std::endl;
    tVector world_pos = mat4 * mat3 * pixel_coord;
    // std::cout << "world_pos = " << world_pos.transpose() << std::endl;
    // std::cout << "mat4 = \n" << mat4 << std::endl;
    return world_pos;
}