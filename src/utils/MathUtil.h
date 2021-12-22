#pragma once
#include "Rand.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <random>
#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

const int gInvalidIdx = -1;

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

// for convenience define standard vector for rendering
typedef Eigen::Vector4d tVector;
typedef Eigen::VectorXd tVectorXd;
typedef Eigen::VectorXi tVectorXi;
typedef Eigen::Vector3i tVector3i;
typedef Eigen::Vector4i tVector4i;
typedef Eigen::VectorXf tVectorXf;
typedef Eigen::Vector3d tVector3d;
typedef Eigen::Vector3f tVector3f;
typedef Eigen::Vector4f tVector4f;
typedef Eigen::Vector2f tVector2f;
typedef Eigen::Vector2d tVector2d;
typedef Eigen::Vector2i tVector2i;
typedef Eigen::Matrix2i tMatrix2i;
typedef Eigen::Matrix4d tMatrix;
typedef Eigen::Matrix2f tMatrix2f;
typedef Eigen::Matrix2d tMatrix2d;
typedef Eigen::Matrix3d tMatrix3d;
typedef Eigen::Matrix3f tMatrix3f;
typedef Eigen::MatrixXd tMatrixXd;
typedef Eigen::Matrix4f tMatrix4f;
typedef Eigen::MatrixXf tMatrixXf;
typedef Eigen::MatrixXi tMatrixXi;
typedef Eigen::Quaterniond tQuaternion;
typedef Eigen::Affine3d aff3;
typedef Eigen::Affine3f aff3f;
typedef Eigen::SparseMatrix<double> tSparseMat;
typedef Eigen::Triplet<double> tTriplet;
template <typename T>
using tEigenArr = std::vector<T, Eigen::aligned_allocator<T>>;
typedef tEigenArr<tVector> tVectorArr;

const double gRadiansToDegrees = 57.2957795;
const double gDegreesToRadians = 1.0 / gRadiansToDegrees;
extern const tVector gGravity;
const double gInchesToMeters = 0.0254;
const double gFeetToMeters = 0.3048;

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

class cMathUtil
{
public:
    enum eAxis
    {
        eAxisX,
        eAxisY,
        eAxisZ,
        eAxisMax
    };
    static bool IsPoint(const tVector &vec);
    static tVector VecToPoint(const tVector &vec);
    static int Clamp(int val, int min, int max);
    static void Clamp(const Eigen::VectorXd &min, const Eigen::VectorXd &max,
                      Eigen::VectorXd &out_vec);
    static double Clamp(double val, double min, double max);
    static double Saturate(double val);
    static double Lerp(double t, double val0, double val1);

    static double NormalizeAngle(double theta);

    // rand number
    static double RandDouble();
    static double RandDouble(double min, double max);
    static double RandDoubleNorm(double mean, double stdev);
    static double RandDoubleExp(double lambda);
    static double RandDoubleSeed(double seed);
    static int RandInt();
    static int RandInt(int min, int max);
    static int RandUint();
    static int RandUint(unsigned int min, unsigned int max);
    static int RandIntExclude(int min, int max, int exc);
    static void SeedRand(unsigned long int seed);
    static int RandSign();
    static bool FlipCoin(double p = 0.5);
    static double SmoothStep(double t);

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
    static tMatrix Calc_Dq1q0conj_Dq0(const tQuaternion &q0,
                                      const tQuaternion &q1);
    static void TestCalc_Dq1q0conj_Dq0();
    static tMatrix Calc_DQuaternion_DAxisAngle(const tVector &aa);
    static tMatrixXd Calc_DQuaterion_DEulerAngles(const tVector &euler_angles,
                                                  eRotationOrder order);
    static void TestCalc_DQuaterion_DEulerAngles();
    static void TestCalc_DQuaterniontDAxisAngle();

    static double QuatTheta(const tQuaternion &dq);
    static tQuaternion VecDiffQuat(const tVector &v0, const tVector &v1);
    static tVector QuatRotVec(const tQuaternion &q, const tVector &dir);
    static tQuaternion MirrorQuaternion(const tQuaternion &q, eAxis axis);

    static double Sign(double val);
    static int Sign(int val);

    static double AddAverage(double avg0, int count0, double avg1, int count1);
    static tVector AddAverage(const tVector &avg0, int count0,
                              const tVector &avg1, int count1);
    static void AddAverage(const Eigen::VectorXd &avg0, int count0,
                           const Eigen::VectorXd &avg1, int count1,
                           Eigen::VectorXd &out_result);
    static void CalcSoftmax(const Eigen::VectorXd &vals, double temp,
                            Eigen::VectorXd &out_prob);
    static double EvalGaussian(const Eigen::VectorXd &mean,
                               const Eigen::VectorXd &covar,
                               const Eigen::VectorXd &sample);
    static double EvalGaussian(double mean, double covar, double sample);
    static double CalcGaussianPartition(const Eigen::VectorXd &covar);
    static double EvalGaussianLogp(double mean, double covar, double sample);
    static double EvalGaussianLogp(const Eigen::VectorXd &mean,
                                   const Eigen::VectorXd &covar,
                                   const Eigen::VectorXd &sample);
    static double Sigmoid(double x);
    static double Sigmoid(double x, double gamma, double bias);

    static int SampleDiscreteProb(const Eigen::VectorXd &probs);
    static tVector CalcBarycentric(const tVector &p, const tVector &a,
                                   const tVector &b, const tVector &c);

    static bool ContainsAABB(const tVector &pt, const tVector &aabb_min,
                             const tVector &aabb_max);
    static bool ContainsAABB(const tVector &aabb_min0, const tVector &aabb_max0,
                             const tVector &aabb_min1,
                             const tVector &aabb_max1);
    static bool ContainsAABBXZ(const tVector &pt, const tVector &aabb_min,
                               const tVector &aabb_max);
    static bool ContainsAABBXZ(const tVector &aabb_min0,
                               const tVector &aabb_max0,
                               const tVector &aabb_min1,
                               const tVector &aabb_max1);
    static void CalcAABBIntersection(const tVector &aabb_min0,
                                     const tVector &aabb_max0,
                                     const tVector &aabb_min1,
                                     const tVector &aabb_max1, tVector &out_min,
                                     tVector &out_max);
    static void CalcAABBUnion(const tVector &aabb_min0,
                              const tVector &aabb_max0,
                              const tVector &aabb_min1,
                              const tVector &aabb_max1, tVector &out_min,
                              tVector &out_max);
    static bool IntersectAABB(const tVector &aabb_min0,
                              const tVector &aabb_max0,
                              const tVector &aabb_min1,
                              const tVector &aabb_max1);
    static bool IntersectAABBXZ(const tVector &aabb_min0,
                                const tVector &aabb_max0,
                                const tVector &aabb_min1,
                                const tVector &aabb_max1);

    // check if curr_val and curr_val - delta belong to different intervals
    static bool CheckNextInterval(double delta, double curr_val,
                                  double int_size);

    static tVector SampleRandPt(const tVector &bound_min,
                                const tVector &bound_max);
    // samples a bound within the given bounds with a benter towards the focus
    // pt
    static tVector SampleRandPtBias(const tVector &bound_min,
                                    const tVector &bound_max);
    static tVector SampleRandPtBias(const tVector &bound_min,
                                    const tVector &bound_max,
                                    const tVector &focus);

    static void QuatSwingTwistDecomposition(const tQuaternion &q,
                                            const tVector &dir,
                                            tQuaternion &out_swing,
                                            tQuaternion &out_twist);
    static tQuaternion ProjectQuat(const tQuaternion &q, const tVector &dir);

    static void ButterworthFilter(double dt, double cutoff,
                                  Eigen::VectorXd &out_x);

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
    static tVector CalcAngularVelocity(const tQuaternion &old_rot,
                                       const tQuaternion &new_rot,
                                       double timestep);
    static tVector CalcAngularVelocityFromAxisAngle(const tQuaternion &old_rot,
                                                    const tQuaternion &new_rot,
                                                    double timestep);
    static tVector QuaternionToEulerAngles(const tQuaternion &,
                                           const eRotationOrder &order);

    static tMatrix EulerAnglesToRotMat(const tVector &euler,
                                       const eRotationOrder &order);
    static tMatrix EulerAnglesToRotMatDot(const tVector &euler,
                                          const eRotationOrder &order);
    static tVector AngularVelToqdot(const tVector &omega, const tVector &cur_q,
                                    const eRotationOrder &order);
    static tMatrix VectorToSkewMat(const tVector &);
    static tVector SkewMatToVector(const tMatrix &);
    static bool IsSame(const tVector &v1, const tVector &v2, const double eps);
    static void ThresholdOp(tVectorXd &v, double threshold = 1e-6);
    static tVector CalcAxisAngleFromOneVectorToAnother(const tVector &v0,
                                                       const tVector &v1);
    template <typename T> static const std::string EigenToString(const T &mat)
    {
        std::stringstream ss;
        ss << mat;
        return ss.str();
    }
    static double Truncate(double num, int digits = 5);
    static tMatrixXd ExpandFrictionCone(int num_friction_dirs,
                                        const tVector &normal);
    static tMatrix InverseTransform(const tMatrix &);
    static double CalcConditionNumber(const tMatrixXd &mat);
    static tMatrixXd JacobPreconditioner(const tMatrixXd &mat);
    // static void RoundZero(tMatrixXd &mat, double threshold = 1e-10);

    template <typename T>
    static void RoundZero(T &mat, double threshold = 1e-10)
    {
        mat = (threshold < mat.array().abs()).select(mat, 0.0f);
    }
    template <typename T> static tVector Expand(const T &vec, double n)
    {
        return tVector(vec[0], vec[1], vec[2], n);
    }
    template <typename T> static tMatrix ExpandMat(const T &raw_mat)
    {
        tMatrix mat = tMatrix::Zero();
        mat.block(0, 0, 3, 3) = raw_mat.block(0, 0, 3, 3);
        return mat;
    }
    static tVector RayCastTri(const tVector &ori, const tVector &dir,
                              const tVector &p1, const tVector &p2,
                              const tVector &p3, double eps = 1e-10);
    static tVector RayCastPlane(const tVector &ray_ori, const tVector &ray_dir,
                                const tVector &plane_eqaution,
                                double eps = 1e-10);
    static tMatrixXd
    CartesianProduct(const std::vector<std::vector<double>> &lists);
    static std::vector<std::vector<double>>
    CartesianProductVec(const std::vector<std::vector<double>> &lists);
    static double CalcDistanceFromPointToLine(const tVector3d &point,
                                              const tVector3d &line_origin,
                                              const tVector3d &line_end);
    static tVector CalcNormalFromPlane(const tVector &plane_equation);
    static double EvaluatePlane(const tVector &plane, const tVector &point);
    static double CalcPlanePointDist(const tVector &plane,
                                     const tVector3d &point);
    static tVector SampleFromPlane(const tVector &plane_equation);
    static float CalcTriangleArea(const tVector &p0, const tVector &p1,
                                  const tVector &p2);
    static int RandIntCategorical(const std::vector<double> &prop);

private:
    static cRand gRand;

    template <typename T> static T SignAux(T val)
    {
        return (T(0) < val) - (val < T(0));
    }

    static tMatrix EulerAngleRotmatX(double x);
    static tMatrix EulerAngleRotmatY(double x);
    static tMatrix EulerAngleRotmatZ(double x);
    static tMatrix EulerAngleRotmatdX(double x);
    static tMatrix EulerAngleRotmatdY(double x);
    static tMatrix EulerAngleRotmatdZ(double x);
};
