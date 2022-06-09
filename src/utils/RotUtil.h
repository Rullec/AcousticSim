#include "utils/MathUtil.h"
typedef Eigen::Quaterniond tQuaternion;
typedef Eigen::Affine3d aff3;
typedef Eigen::Affine3f aff3f;

enum eRotationOrder
{
    XYZ = 0, // first X, then Y, then Z. X->Y->Z. R_{total} = Rz * Ry * Rx;
    XZY,
    XYX,
    XZX, // x end
    YXZ,
    YZX,
    YXY,
    YZY, // y end
    ZXY,
    ZYX,
    ZYZ,
    ZXZ, // z end
};

// extern const enum eRotationOrder gRotationOrder;// rotation order. declared
// here and defined in LoboJointV2.cpp
const std::string ROTATION_ORDER_NAME[] = {
    "XYZ", "XZY", "XYX", "XZX", "YXZ", "YZX",
    "YXY", "YZY", "ZXY", "ZYX", "ZYZ", "ZXZ",
};

namespace Eigen
{

/// @brief Returns a perspective transformation matrix like the one from
/// gluPerspective
/// @see http://www.opengl.org/sdk/docs/man2/xhtml/gluPerspective.xml
/// @see glm::perspective
template <typename Scalar>
Eigen::Matrix<Scalar, 4, 4> perspective(Scalar fovy, Scalar aspect,
                                        Scalar zNear, Scalar zFar)
{
    Transform<Scalar, 3, Projective> tr;
    tr.matrix().setZero();
    assert(aspect > 0);
    assert(zFar > zNear);
    assert(zNear > 0);
    Scalar radf = static_cast<Scalar>(M_PI * fovy / 180.0);
    Scalar tan_half_fovy = static_cast<Scalar>(std::tan(radf / 2.0));
    tr(0, 0) = static_cast<Scalar>(1.0 / (aspect * tan_half_fovy));
    tr(1, 1) = static_cast<Scalar>(1.0 / (tan_half_fovy));
    tr(2, 2) = -(zFar + zNear) / (zFar - zNear);
    tr(3, 2) = -1.0;
    tr(2, 3) = static_cast<Scalar>(-(2.0 * zFar * zNear) / (zFar - zNear));
    return tr.matrix();
}

template <typename Scalar>
Eigen::Matrix<Scalar, 4, 4> scale(Scalar x, Scalar y, Scalar z)
{
    Transform<Scalar, 3, Affine> tr;
    tr.matrix().setZero();
    tr(0, 0) = x;
    tr(1, 1) = y;
    tr(2, 2) = z;
    tr(3, 3) = 1;
    return tr.matrix();
}

template <typename Scalar>
Eigen::Matrix<Scalar, 4, 4> translate(Scalar x, Scalar y, Scalar z)
{
    Transform<Scalar, 3, Affine> tr;
    tr.matrix().setIdentity();
    tr(0, 3) = x;
    tr(1, 3) = y;
    tr(2, 3) = z;
    return tr.matrix();
}

/// @brief Returns a view transformation matrix like the one from glu's lookAt
/// @see http://www.opengl.org/sdk/docs/man2/xhtml/gluLookAt.xml
/// @see glm::lookAt
template <typename Derived>
Eigen::Matrix<typename Derived::Scalar, 4, 4>
lookAt(Derived const &eye, Derived const &center, Derived const &up)
{
    typedef Eigen::Matrix<typename Derived::Scalar, 4, 4> Matrix4;
    typedef Eigen::Matrix<typename Derived::Scalar, 3, 1> Vector3;
    Vector3 f = (center - eye).normalized();
    Vector3 u = up.normalized();
    Vector3 s = f.cross(u).normalized();
    u = s.cross(f);
    Matrix4 mat = Matrix4::Zero();
    mat(0, 0) = s.x();
    mat(0, 1) = s.y();
    mat(0, 2) = s.z();
    mat(0, 3) = -s.dot(eye);
    mat(1, 0) = u.x();
    mat(1, 1) = u.y();
    mat(1, 2) = u.z();
    mat(1, 3) = -u.dot(eye);
    mat(2, 0) = -f.x();
    mat(2, 1) = -f.y();
    mat(2, 2) = -f.z();
    mat(2, 3) = f.dot(eye);
    mat.row(3) << 0, 0, 0, 1;
    return mat;
}

/// @see glm::ortho
template <typename Scalar>
Eigen::Matrix<Scalar, 4, 4> ortho(Scalar const &left, Scalar const &right,
                                  Scalar const &bottom, Scalar const &top,
                                  Scalar const &zNear, Scalar const &zFar)
{
    Eigen::Matrix<Scalar, 4, 4> mat = Eigen::Matrix<Scalar, 4, 4>::Identity();
    mat(0, 0) = Scalar(2) / (right - left);
    mat(1, 1) = Scalar(2) / (top - bottom);
    mat(2, 2) = -Scalar(2) / (zFar - zNear);
    mat(3, 0) = -(right + left) / (right - left);
    mat(3, 1) = -(top + bottom) / (top - bottom);
    mat(3, 2) = -(zFar + zNear) / (zFar - zNear);
    return mat;
}

} // namespace Eigen

class cRotUtil
{
public:
    // matrices
    static tMatrix TransformMat(const tVector &translation,
                                const tVector &euler_xyz_orientation);
    static tMatrix TranslateMat(const tVector &trans);
    static tMatrix ScaleMat(double scale);
    static tMatrix ScaleMat(const tVector &scale);
    static tMatrix
    RotateMat(const tVector &euler,
              const eRotationOrder gRotationOrder); // euler angles order rot(Z)
                                                    // * rot(Y) * rot(X)
    static tMatrix RotateMat(const tVector &axis, double theta);
    static tMatrix RotateMat(const tQuaternion &q);
    static tMatrix CrossMat(const tVector &a);
    // inverts a transformation consisting only of rotations and translations
    static tMatrix InvRigidMat(const tMatrix &mat);
    static tVector GetRigidTrans(const tMatrix &mat);
    static tVector InvEuler(const tVector &euler,
                            const eRotationOrder gRotationOrder);
    static void RotMatToAxisAngle(const tMatrix &mat, tVector &out_axis,
                                  double &out_theta);
    static tVector RotMatToEuler(const tMatrix &mat,
                                 const eRotationOrder gRotationOrder);
    static tMatrix2d RotMat2D(double angle);
    static tMatrix AxisAngleToRotmat(const tVector &angvel);
    static tQuaternion RotMatToQuaternion(const tMatrix &mat);
    static tVector EulerangleToAxisAngle(const tVector &euler,
                                         const eRotationOrder gRotationOrder);
    static void EulerToAxisAngle(const tVector &euler, tVector &out_axis,
                                 double &out_theta,
                                 const eRotationOrder gRotationOrder);
    static tVector AxisAngleToEuler(const tVector &axis, double theta);
    static tMatrix DirToRotMat(const tVector &dir, const tVector &up);

    static void DeltaRot(const tVector &axis0, double theta0,
                         const tVector &axis1, double theta1, tVector &out_axis,
                         double &out_theta);
    static tMatrix DeltaRot(const tMatrix &R0, const tMatrix &R1);

    static tQuaternion EulerToQuaternion(const tVector &euler,
                                         const eRotationOrder order);
    static tQuaternion CoefVectorToQuaternion(const tVector &coef);
    static tVector QuaternionToEuler(const tQuaternion &q,
                                     const eRotationOrder gRotationOrder);
    static tQuaternion AxisAngleToQuaternion(const tVector &axis, double theta);
    static tVector QuaternionToAxisAngle(const tQuaternion &q);
    static void QuaternionToAxisAngle(const tQuaternion &q, tVector &out_axis,
                                      double &out_theta);
    static tMatrix BuildQuaternionDiffMat(const tQuaternion &q);
    static tVector CalcQuaternionVel(const tQuaternion &q0,
                                     const tQuaternion &q1, double dt);
    static tVector CalcQuaternionVelRel(const tQuaternion &q0,
                                        const tQuaternion &q1, double dt);
    static tQuaternion VecToQuat(const tVector &v);
    static tVector QuatToVec(const tQuaternion &q);
    static tQuaternion QuatDiff(const tQuaternion &q0, const tQuaternion &q1);
    static double QuatDiffTheta(const tQuaternion &q0, const tQuaternion &q1);
    // static tMatrix Calc_Dq1q0conj_Dq0(const tQuaternion &q0,
    //                                   const tQuaternion &q1);
    // static void TestCalc_Dq1q0conj_Dq0();
    // static tMatrix Calc_DQuaternion_DAxisAngle(const tVector &aa);
    // static tMatrixXd Calc_DQuaterion_DEulerAngles(const tVector
    // &euler_angles,
    //                                               eRotationOrder order);
    // static void TestCalc_DQuaterion_DEulerAngles();
    // static void TestCalc_DQuaterniontDAxisAngle();

    static double QuatTheta(const tQuaternion &dq);
    static tQuaternion VecDiffQuat(const tVector &v0, const tVector &v1);
    static tVector QuatRotVec(const tQuaternion &q, const tVector &dir);
    // static tQuaternion MirrorQuaternion(const tQuaternion &q, eAxis axis);

    // static void QuatSwingTwistDecomposition(const tQuaternion &q,
    //                                         const tVector &dir,
    //                                         tQuaternion &out_swing,
    //                                         tQuaternion &out_twist);
    // static tQuaternion ProjectQuat(const tQuaternion &q, const tVector &dir);

    // static void ButterworthFilter(double dt, double cutoff,
    //                               Eigen::VectorXd &out_x);

    // added by myself
    static tMatrix RotMat(const tQuaternion &quater);
    // static tQuaternion RotMatToQuaternion(const tMatrix &mat);
    static tQuaternion CoefToQuaternion(const tVector &);
    static tQuaternion AxisAngleToQuaternion(const tVector &angvel);
    static tQuaternion EulerAnglesToQuaternion(const tVector &vec,
                                               const eRotationOrder &order);
    static tQuaternion MinusQuaternion(const tQuaternion &quad);
    static tVector QuaternionToCoef(const tQuaternion &quater);
    // static tVector QuaternionToAxisAngle(const tQuaternion &);
    // static tVector CalcAngularVelocity(const tQuaternion &old_rot,
    //                                    const tQuaternion &new_rot,
    //                                    double timestep);
    static tVector CalcAxisAngleFromOneVectorToAnother(const tVector &v0,
                                                       const tVector &v1);

    static tVector CalcAngularVelocityFromAxisAngle(const tQuaternion &old_rot,
                                                    const tQuaternion &new_rot,
                                                    double timestep);
    static tVector QuaternionToEulerAngles(const tQuaternion &,
                                           const eRotationOrder &order);

    static tMatrix EulerAnglesToRotMat(const tVector &euler,
                                       const eRotationOrder &order);
    // static tMatrix EulerAnglesToRotMatDot(const tVector &euler,
    //                                       const eRotationOrder &order);
    // static tVector AngularVelToqdot(const tVector &omega, const tVector
    // &cur_q,
    //                                 const eRotationOrder &order);

protected:
    static tMatrix EulerAngleRotmatX(double x);
    static tMatrix EulerAngleRotmatY(double x);
    static tMatrix EulerAngleRotmatZ(double x);
    static tMatrix EulerAngleRotmatdX(double x);
    static tMatrix EulerAngleRotmatdY(double x);
    static tMatrix EulerAngleRotmatdZ(double x);
};