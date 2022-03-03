#include "LogUtil.h"
#include "MathUtil.h"
#include <iostream>
#include <time.h>

tMatrix cMathUtil::VectorToSkewMat(const tVector &vec)
{
    tMatrix res = tMatrix::Zero();
    double a = vec[0], b = vec[1], c = vec[2];
    res(0, 1) = -c;
    res(0, 2) = b;
    res(1, 0) = c;
    res(1, 2) = -a;
    res(2, 0) = -b;
    res(2, 1) = a;

    return res;
}

tVector cMathUtil::SkewMatToVector(const tMatrix &mat)
{
    // verify mat is a skew matrix
    assert((mat + mat.transpose()).norm() < 1e-10);

    // squeeze a mat to a vector
    tVector res = tVector::Zero();
    res[0] = mat(2, 1);
    res[1] = mat(0, 2);
    res[2] = mat(1, 0);
    return res;
}

bool cMathUtil::IsSame(const tVector &v1, const tVector &v2, const double eps)
{
    for (int i = 0; i < v1.size(); i++)
        if (std::fabs(v1[i] - v2[i]) > eps)
            return false;
    return true;
}

void cMathUtil::ThresholdOp(tVectorXd &v, double threshold)
{
    v = (threshold < v.array().abs()).select(v, 0.0f);
}

tVector cMathUtil::CalcAxisAngleFromOneVectorToAnother(const tVector &v0_,
                                                       const tVector &v1_)
{
    tVector v0 = v0_.normalized(), v1 = v1_.normalized();

    tVector rot_axis = v0.cross3(v1);
    double theta = std::asin(rot_axis.norm()); //[-pi/2, pi/2]

    // if the angle between v0 and v1 > 90
    if (v0.dot(v1) < 0)
    {
        theta = theta > 0 ? (theta + (M_PI / 2 - theta) * 2)
                          : (theta + (-M_PI / 2 - theta) * 2);
    }
    rot_axis = rot_axis.normalized() * std::fabs(theta);
    return rot_axis;
}

double cMathUtil::Truncate(double num, int digits)
{
    return round(num * pow(10, digits)) / pow(10, digits);
}

// // Nx3 friction cone
// // each row is a direction now
// tMatrixXd cMathUtil::ExpandFrictionCone(int num_friction_dirs,
//                                         const tVector &normal_)
// {
//     // 1. check the input
//     tVector normal = normal_;
//     normal[3] = 0;
//     normal.normalize();
//     if (normal.norm() < 1e-6)
//     {
//         std::cout << "[error] ExpandFrictionCone normal = "
//                   << normal_.transpose() << std::endl;
//         exit(0);
//     }

//     // 2. generate a standard friction cone
//     tMatrixXd D = tMatrixXd::Zero(4, num_friction_dirs);
//     double gap = 2 * M_PI / num_friction_dirs;
//     for (int i = 0; i < num_friction_dirs; i++)
//     {
//         D(0, i) = std::cos(gap * i);
//         D(2, i) = std::sin(gap * i);
//     }

//     // 3. rotate the fricition cone
//     tVector Y_normal = tVector(0, 1, 0, 0);
//     tVector axis = Y_normal.cross3(normal).normalized();
//     double theta = std::acos(Y_normal.dot(normal)); // [0, pi]
//     D = cRotUtil::RotMat(cRotUtil::AxisAngleToQuaternion(axis * theta)) * D;
//     D.transposeInPlace();
//     // each row is a direction now
//     return D;
// }

tMatrix cMathUtil::InverseTransform(const tMatrix &raw_trans)
{
    std::cout << "wrong api InverseTransform should not be called\n";
    exit(1);
    tMatrix inv_trans = tMatrix::Identity();
    inv_trans.block(0, 0, 3, 3).transposeInPlace();
    inv_trans.block(0, 3, 3, 1) =
        -inv_trans.block(0, 0, 3, 3) * raw_trans.block(0, 3, 3, 1);
    return inv_trans;
}

double cMathUtil::CalcConditionNumber(const tMatrixXd &mat)
{
    Eigen::EigenSolver<tMatrixXd> solver(mat);
    tVectorXd eigen_values = solver.eigenvalues().real();
    return eigen_values.maxCoeff() / eigen_values.minCoeff();
}

// /**
//  * \brief		Get the jacobian preconditioner P = diag(A)
//  *
//  */
// tMatrixXd cMathUtil::JacobPreconditioner(const tMatrixXd &A)
// {
//     if (A.rows() != A.cols())
//     {
//         std::cout << "cMathUtil::JacobPreconditioner: A is not a square matrix "
//                   << A.rows() << " " << A.cols() << std::endl;
//         exit(1);
//     }
//     tVectorXd diagonal = A.diagonal();
//     if (diagonal.cwiseAbs().minCoeff() < 1e-10)
//     {
//         std::cout
//             << "cMathUtil::JacobPreconditioner: diagnoal is nearly zero for "
//             << diagonal.transpose() << std::endl;
//         exit(1);
//     }

//     return diagonal.cwiseInverse().asDiagonal();
// }

tVector cMathUtil::RayCastTri(const tVector &ori, const tVector &dir,
                              const tVector &p1, const tVector &p2,
                              const tVector &p3, double eps /* = 1e-5*/)
{
    Eigen::Matrix3d mat;
    mat.col(0) = (p1 - p2).segment(0, 3);
    mat.col(1) = (p1 - p3).segment(0, 3);
    mat.col(2) = dir.segment(0, 3);
    Eigen::Vector3d vec;
    vec = (p1 - ori).segment(0, 3);
    Eigen::Vector3d res = mat.inverse() * vec;
    // std::cout << "res = " << res.transpose() << std::endl;
    double beta = res[0], gamma = res[1], t = res[2], alpha = 1 - beta - gamma;
    // std::cout <<"ray cast = " << res.transpose() << std::endl;
    tVector inter =
        tVector(std::nan(""), std::nan(""), std::nan(""), std::nan(""));
    if (0 - eps < alpha && alpha < 1 + eps && 0 - eps < beta &&
        beta < 1 + eps && 0 - eps < gamma && gamma < 1 + eps && t > 0 - eps)
    {
        inter = ori + t * dir;
    }
    return inter;
}

tVector cMathUtil::RayCastPlane(const tVector &ray_ori, const tVector &ray_dir,
                                const tVector &plane_equation,
                                double eps /*= 1e-10*/)
{
    double t = -(plane_equation.segment(0, 3).dot(ray_ori.segment(0, 3)) +
                 plane_equation[3]) /
               (plane_equation.segment(0, 3).dot(ray_dir.segment(0, 3)));
    if (t < eps)
    {
        return tVector::Ones() * std::nan("");
    }
    else
    {
        return ray_ori + ray_dir * t;
    }
}
/**
 * \brief               cartesian product for sets
 */
tMatrixXd
cMathUtil::CartesianProduct(const std::vector<std::vector<double>> &lists)
{
    std::vector<std::vector<double>> result = CartesianProductVec(lists);

    tMatrixXd eigen_res = tMatrixXd::Zero(result.size(), result[0].size());
    for (int i = 0; i < result.size(); i++)
    {
        for (int j = 0; j < result[i].size(); j++)
        {
            eigen_res(i, j) = result[i][j];
        }
    }
    return eigen_res;
}

std::vector<std::vector<double>>
cMathUtil::CartesianProductVec(const std::vector<std::vector<double>> &lists)
{
    std::vector<std::vector<double>> result(0);
    if (std::find_if(std::begin(lists), std::end(lists),
                     [](auto e) -> bool
                     { return e.size() == 0; }) != std::end(lists))
    {
        return result;
    }
    for (auto &e : lists[0])
    {
        result.push_back({e});
    }
    for (size_t i = 1; i < lists.size(); ++i)
    {
        std::vector<std::vector<double>> temp;
        for (auto &e : result)
        {
            for (auto f : lists[i])
            {
                auto e_tmp = e;
                e_tmp.push_back(f);
                temp.push_back(e_tmp);
            }
        }
        result = temp;
    }
    return result;
}

/**
 * \brief           calcualte the distance between a point to a line
 */
double cMathUtil::CalcDistanceFromPointToLine(const tVector3d &point,
                                              const tVector3d &line_origin,
                                              const tVector3d &line_end)
{
    tVector3d origin_2_point = point - line_origin;
    // std::cout << "origin_2_point = " << origin_2_point.transpose() <<
    // std::endl;
    tVector3d origin_2_end = (line_end - line_origin).normalized();
    // std::cout << "origin_2_end = " << origin_2_end.transpose() << std::endl;
    double length = origin_2_point.dot(origin_2_end);

    // std::cout << "length = " << length << std::endl;
    tVector3d origin_2_point_proj = length * origin_2_end;
    // std::cout << "origin_2_point_proj = " << origin_2_point_proj.transpose()
    // << std::endl;
    tVector3d res = origin_2_point - origin_2_point_proj;
    return res.norm();
}

/**
 * \brief           evaulate the plane equation:
 *      ax + by + cd + d = ? given:
 *      a = plane[0]
 *      b = plane[1]
 *      c = plane[2]
 *      d = plane[3]
 */
double cMathUtil::EvaluatePlane(const tVector &plane, const tVector &point)
{
    double sum = plane[3];
    for (int i = 0; i < 3; i++)
        sum += plane[i] * point[i];
    return sum;
}

/**
 * \brief           calcualte the distance between given plane and point
 *
 *      1. select any point in the plane "A"
 *      2. given the point P calculat the vector "PA"
 *      3. given the plane normal "n",  dist = n \cdot PA
 */
double cMathUtil::CalcPlanePointDist(const tVector &plane,
                                     const tVector3d &point_P)
{
    // 1. select random point on A
    if (std::fabs(plane.segment(0, 3).norm()) < 1e-10)
    {
        SIM_WARN("invalid plane {}", plane.transpose());
        return std::nan("");
    }
    tVector3d point_A = SampleFromPlane(plane).segment(0, 3);

    // 2. calcualte the vector PA
    tVector3d vec_PA = point_A - point_P;
    tVector3d vec_n = plane.segment(0, 3).normalized();
    return std::fabs(vec_n.dot(vec_PA));
}

/**
 * \brief           sample from plane
 */
tVector cMathUtil::SampleFromPlane(const tVector &plane_equation)
{
    // 1. find the non-zero index
    int id = -1;
    for (int i = 0; i < 3; i++)
    {
        if (std::fabs(plane_equation[i]) > 1e-10)
        {
            id = i;
            break;
        }
    }

    if (id == 3)
    {
        SIM_ERROR("failed to sample from plane {}", plane_equation.transpose());
        exit(1);
    }

    tVector point = tVector::Random();
    point[3] = 1;
    double residual = -plane_equation[3];
    for (int i = 0; i < 3; i++)
    {
        if (i != id)
        {
            residual += -plane_equation[i] * point[i];
        }
    }
    point[id] = residual / plane_equation[id];
    return point;
};

/**
 * \brief           Calculate normal by given plane equaiton "ax + by + cz + d =
 * 0"
 */
tVector cMathUtil::CalcNormalFromPlane(const tVector &plane_equation)
{
    tVector res = tVector::Zero();
    res.segment(0, 3) = plane_equation.segment(0, 3).normalized();
    return res;

    // tVector
    //     p0 = SampleFromPlane(plane_equation),
    //     p1 = SampleFromPlane(plane_equation),
    //     p2 = SampleFromPlane(plane_equation);

    // tVector v0 = p0 - p1,
    //         v1 = p1 - p2;
    // tVector normal = (v0.cross3(v1)).normalized();
    // if (EvaluatePlane(plane_equation, p0 + normal) < 0)
    //     normal *= -1;
    // // std::cout << "plane = " << plane_equation.transpose() << std::endl;
    // // std::cout << "p0 = " << p0.transpose() << std::endl;
    // // std::cout << "p1 = " << p1.transpose() << std::endl;
    // // std::cout << "p2 = " << p2.transpose() << std::endl;
    // // exit(1);
    // return normal;
}

float cMathUtil::CalcTriangleArea(const tVector &p0, const tVector &p1,
                                  const tVector &p2)
{
    return cMathUtil::CalcTriangleArea3d(
        p0.segment(0, 3),
        p1.segment(0, 3),
        p2.segment(0, 3));
}

float cMathUtil::CalcTriangleArea3d(const tVector3d &p0, const tVector3d &p1,
                                  const tVector3d &p2)
{
    return 0.5 * ((p1 - p0).cross(p2 - p0)).norm();
}