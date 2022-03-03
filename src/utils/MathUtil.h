#pragma once
#include "Rand.h"
#include <Eigen/Dense>
#include <random>
#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

const int gInvalidIdx = -1;

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
template <typename T>
using tEigenArr = std::vector<T, Eigen::aligned_allocator<T>>;
typedef tEigenArr<tVector> tVectorArr;

const double gRadiansToDegrees = 57.2957795;
const double gDegreesToRadians = 1.0 / gRadiansToDegrees;
const double gInchesToMeters = 0.0254;
const double gFeetToMeters = 0.3048;

class cMathUtil
{
public:
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

    static tMatrix VectorToSkewMat(const tVector &);
    static tVector SkewMatToVector(const tMatrix &);
    static bool IsSame(const tVector &v1, const tVector &v2, const double eps);
    static void ThresholdOp(tVectorXd &v, double threshold = 1e-6);
    static tVector CalcAxisAngleFromOneVectorToAnother(const tVector &v0,
                                                       const tVector &v1);
    template <typename T>
    static const std::string EigenToString(const T &mat)
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
    // static tMatrixXd JacobPreconditioner(const tMatrixXd &mat);
    // static void RoundZero(tMatrixXd &mat, double threshold = 1e-10);

    template <typename T>
    static void RoundZero(T &mat, double threshold = 1e-10)
    {
        mat = (threshold < mat.array().abs()).select(mat, 0.0f);
    }
    template <typename T>
    static tVector Expand(const T &vec, double n)
    {
        return tVector(vec[0], vec[1], vec[2], n);
    }
    template <typename T>
    static tMatrix ExpandMat(const T &raw_mat)
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
    static float CalcTriangleArea3d(const tVector3d &p0, const tVector3d &p1,
                                    const tVector3d &p2);
    static int RandIntCategorical(const std::vector<double> &prop);

private:
    static cRand gRand;

    template <typename T>
    static T SignAux(T val)
    {
        return (T(0) < val) - (val < T(0));
    }

};
