
#ifndef EIGEN_UTIL_H_
#define EIGEN_UTIL_H_
#include <Eigen/Dense>
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
#endif