#include "LogUtil.h"
#include "MathUtil.h"
#include <iostream>
#include <time.h>

tVector cMathUtil::CalcBarycentric(const tVector &p, const tVector &a,
                                   const tVector &b, const tVector &c)
{
    tVector v0 = b - a;
    tVector v1 = c - a;
    tVector v2 = p - a;

    double d00 = v0.dot(v0);
    double d01 = v0.dot(v1);
    double d11 = v1.dot(v1);
    double d20 = v2.dot(v0);
    double d21 = v2.dot(v1);
    double denom = d00 * d11 - d01 * d01;
    double v = (d11 * d20 - d01 * d21) / denom;
    double w = (d00 * d21 - d01 * d20) / denom;
    double u = 1.0f - v - w;
    return tVector(u, v, w, 0);
}

bool cMathUtil::ContainsAABB(const tVector &pt, const tVector &aabb_min,
                             const tVector &aabb_max)
{
    bool contains = pt[0] >= aabb_min[0] && pt[1] >= aabb_min[1] &&
                    pt[2] >= aabb_min[2] && pt[0] <= aabb_max[0] &&
                    pt[1] <= aabb_max[1] && pt[2] <= aabb_max[2];
    return contains;
}

bool cMathUtil::ContainsAABB(const tVector &aabb_min0, const tVector &aabb_max0,
                             const tVector &aabb_min1, const tVector &aabb_max1)
{
    return ContainsAABB(aabb_min0, aabb_min1, aabb_max1) &&
           ContainsAABB(aabb_max0, aabb_min1, aabb_max1);
}

bool cMathUtil::ContainsAABBXZ(const tVector &pt, const tVector &aabb_min,
                               const tVector &aabb_max)
{
    bool contains = pt[0] >= aabb_min[0] && pt[2] >= aabb_min[2] &&
                    pt[0] <= aabb_max[0] && pt[2] <= aabb_max[2];
    return contains;
}

bool cMathUtil::ContainsAABBXZ(const tVector &aabb_min0,
                               const tVector &aabb_max0,
                               const tVector &aabb_min1,
                               const tVector &aabb_max1)
{
    return ContainsAABBXZ(aabb_min0, aabb_min1, aabb_max1) &&
           ContainsAABBXZ(aabb_max0, aabb_min1, aabb_max1);
}

void cMathUtil::CalcAABBIntersection(const tVector &aabb_min0,
                                     const tVector &aabb_max0,
                                     const tVector &aabb_min1,
                                     const tVector &aabb_max1, tVector &out_min,
                                     tVector &out_max)
{
    out_min = aabb_min0.cwiseMax(aabb_min1);
    out_max = aabb_max0.cwiseMin(aabb_max1);
    if (out_min[0] > out_max[0])
    {
        out_min[0] = 0;
        out_max[0] = 0;
    }
    if (out_min[1] > out_max[1])
    {
        out_min[1] = 0;
        out_max[1] = 0;
    }
    if (out_min[2] > out_max[2])
    {
        out_min[2] = 0;
        out_max[2] = 0;
    }
}

void cMathUtil::CalcAABBUnion(const tVector &aabb_min0,
                              const tVector &aabb_max0,
                              const tVector &aabb_min1,
                              const tVector &aabb_max1, tVector &out_min,
                              tVector &out_max)
{
    out_min = aabb_min0.cwiseMin(aabb_min1);
    out_max = aabb_max0.cwiseMax(aabb_max1);
}

bool cMathUtil::IntersectAABB(const tVector &aabb_min0,
                              const tVector &aabb_max0,
                              const tVector &aabb_min1,
                              const tVector &aabb_max1)
{
    tVector center0 = 0.5 * (aabb_max0 + aabb_min0);
    tVector center1 = 0.5 * (aabb_max1 + aabb_min1);
    tVector size0 = aabb_max0 - aabb_min0;
    tVector size1 = aabb_max1 - aabb_min1;
    tVector test_len = 0.5 * (size0 + size1);
    tVector delta = center1 - center0;
    bool overlap = (std::abs(delta[0]) <= test_len[0]) &&
                   (std::abs(delta[1]) <= test_len[1]) &&
                   (std::abs(delta[2]) <= test_len[2]);
    return overlap;
}

bool cMathUtil::IntersectAABBXZ(const tVector &aabb_min0,
                                const tVector &aabb_max0,
                                const tVector &aabb_min1,
                                const tVector &aabb_max1)
{
    tVector center0 = 0.5 * (aabb_max0 + aabb_min0);
    tVector center1 = 0.5 * (aabb_max1 + aabb_min1);
    tVector size0 = aabb_max0 - aabb_min0;
    tVector size1 = aabb_max1 - aabb_min1;
    tVector test_len = 0.5 * (size0 + size1);
    tVector delta = center1 - center0;
    bool overlap = (std::abs(delta[0]) <= test_len[0]) &&
                   (std::abs(delta[2]) <= test_len[2]);
    return overlap;
}

bool cMathUtil::CheckNextInterval(double delta, double curr_val,
                                  double int_size)
{
    double pad = 0.001 * delta;
    int curr_count = static_cast<int>(std::floor((curr_val + pad) / int_size));
    int prev_count =
        static_cast<int>(std::floor((curr_val + pad - delta) / int_size));
    bool new_action = (curr_count != prev_count);
    return new_action;
}

tVector cMathUtil::SampleRandPt(const tVector &bound_min,
                                const tVector &bound_max)
{
    tVector pt = tVector(RandDouble(bound_min[0], bound_max[0]),
                         RandDouble(bound_min[1], bound_max[1]),
                         RandDouble(bound_min[2], bound_max[2]), 0);
    return pt;
}

tVector cMathUtil::SampleRandPtBias(const tVector &bound_min,
                                    const tVector &bound_max)
{
    return SampleRandPtBias(bound_min, bound_max,
                            0.5 * (bound_max + bound_min));
}

tVector cMathUtil::SampleRandPtBias(const tVector &bound_min,
                                    const tVector &bound_max,
                                    const tVector &focus)
{
    double t = RandDouble(0, 1);
    tVector size = bound_max - bound_min;
    tVector new_min = focus + (t * 0.5) * size;
    tVector new_max = focus - (t * 0.5) * size;
    tVector offset = (bound_min - new_min).cwiseMax(0);
    offset += (bound_max - new_max).cwiseMin(0);
    new_min += offset;
    new_max += offset;

    return SampleRandPt(new_min, new_max);
}

void cMathUtil::QuatSwingTwistDecomposition(const tQuaternion &q,
                                            const tVector &dir,
                                            tQuaternion &out_swing,
                                            tQuaternion &out_twist)
{
    assert(std::abs(dir.norm() - 1) < 0.000001);
    assert(std::abs(q.norm() - 1) < 0.000001);

    tVector q_axis = tVector(q.x(), q.y(), q.z(), 0);
    double p = q_axis.dot(dir);
    tVector twist_axis = p * dir;
    out_twist = tQuaternion(q.w(), twist_axis[0], twist_axis[1], twist_axis[2]);
    out_twist.normalize();
    out_swing = q * out_twist.conjugate();
}

tQuaternion cMathUtil::ProjectQuat(const tQuaternion &q, const tVector &dir)
{
    assert(std::abs(dir.norm() - 1) < 0.00001);
    tVector ref_axis = tVector::Zero();
    int min_idx = 0;
    dir.cwiseAbs().minCoeff(&min_idx);
    ref_axis[min_idx] = 1;

    tVector rot_dir0 = dir.cross3(ref_axis);
    tVector rot_dir1 = cMathUtil::QuatRotVec(q, rot_dir0);
    rot_dir1 -= rot_dir1.dot(dir) * dir;

    double dir1_norm = rot_dir1.norm();
    tQuaternion p_rot = tQuaternion::Identity();
    if (dir1_norm > 0.0001)
    {
        rot_dir1 /= dir1_norm;
        p_rot = cMathUtil::VecDiffQuat(rot_dir0, rot_dir1);
    }
    return p_rot;
}

void cMathUtil::ButterworthFilter(double dt, double cutoff,
                                  Eigen::VectorXd &out_x)
{
    double sampling_rate = 1 / dt;
    int n = static_cast<int>(out_x.size());

    double wc = std::tan(cutoff * M_PI / sampling_rate);
    double k1 = std::sqrt(2) * wc;
    double k2 = wc * wc;
    double a = k2 / (1 + k1 + k2);
    double b = 2 * a;
    double c = a;
    double k3 = b / k2;
    double d = -2 * a + k3;
    double e = 1 - (2 * a) - k3;

    double xm2 = out_x[0];
    double xm1 = out_x[0];
    double ym2 = out_x[0];
    double ym1 = out_x[0];

    for (int s = 0; s < n; ++s)
    {
        double x = out_x[s];
        double y = a * x + b * xm1 + c * xm2 + d * ym1 + e * ym2;

        out_x[s] = y;
        xm2 = xm1;
        xm1 = x;
        ym2 = ym1;
        ym1 = y;
    }

    double yp2 = out_x[n - 1];
    double yp1 = out_x[n - 1];
    double zp2 = out_x[n - 1];
    double zp1 = out_x[n - 1];

    for (int t = n - 1; t >= 0; --t)
    {
        double y = out_x[t];
        double z = a * y + b * yp1 + c * yp2 + d * zp1 + e * zp2;

        out_x[t] = z;
        yp2 = yp1;
        yp1 = y;
        zp2 = zp1;
        zp1 = z;
    }
}

tMatrix cMathUtil::RotMat(const tQuaternion &quater_)
{
    // https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation#Quaternion-derived_rotation_matrix

    tMatrix res = tMatrix::Zero();
    double w = quater_.w(), x = quater_.x(), y = quater_.y(), z = quater_.z();
    res << 1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w), 0,
        2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w), 0,
        2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y), 0, 0,
        0, 0, 1;
    return res;
}

// tQuaternion cMathUtil::RotMatToQuaternion(const tMatrix &mat)
//{
//	//
// http://www.iri.upc.edu/files/scidoc/2068-Accurate-Computation-of-Quaternions-from-Rotation-Matrices.pdf
//	double eta = 0;
//	double q1, q2, q3, q4;	// = [w, x, y, z]
//
//	// determine q1
//	{
//		double detect_value = mat(0, 0) + mat(1, 1) + mat(2, 2);
//		if (detect_value > eta)
//		{
//			q1 = 0.5 * std::sqrt(1 + detect_value);
//		}
//		else
//		{
//			double numerator = 0;
//			numerator += std::pow(mat(2, 1) - mat(1, 2), 2);
//			numerator += std::pow(mat(0, 2) - mat(2, 0), 2);
//			numerator += std::pow(mat(1, 0) - mat(0, 1), 2);
//			q1 = 0.5 *  std::sqrt(numerator / (3 - detect_value));
//		}
//	}
//
//	// determine q2
//	{
//		double detect_value = mat(0, 0) - mat(1, 1) - mat(2, 2);
//		if (detect_value > eta)
//		{
//			q2 = 0.5 * std::sqrt(1 + detect_value);
//		}
//		else
//		{
//			double numerator = 0;
//			numerator += std::pow(mat(2, 1) - mat(1, 2), 2);
//			numerator += std::pow(mat(0, 1) + mat(1, 0), 2);
//			numerator += std::pow(mat(2, 0) + mat(0, 2), 2);
//			q2 = 0.5 * std::sqrt(numerator / (3 - detect_value));
//		}
//	}
//
//	// determine q3
//	{
//		double detect_value = -mat(0, 0) + mat(1, 1) - mat(2, 2);
//		if (detect_value > eta)
//		{
//			q3 = 0.5 * std::sqrt(1 + detect_value);
//		}
//		else
//		{
//			double numerator = 0;
//			numerator += std::pow(mat(0, 2) - mat(2, 0), 2);
//			numerator += std::pow(mat(0, 1) + mat(1, 0), 2);
//			numerator += std::pow(mat(1, 2) + mat(2, 1), 2);
//			q3 = 0.5 * std::sqrt(numerator / (3 - detect_value));
//		}
//	}
//
//	// determine q4
//	{
//		double detect_value = -mat(0, 0) - mat(1, 1) + mat(2, 2);
//		if (detect_value > eta)
//		{
//			q4 = 0.5 * std::sqrt(1 + detect_value);
//		}
//		else
//		{
//			double numerator = 0;
//			numerator += std::pow(mat(1, 0) - mat(0, 1), 2);
//			numerator += std::pow(mat(2, 0) + mat(0, 2), 2);
//			numerator += std::pow(mat(2, 1) + mat(1, 2), 2);
//			q4 = 0.5 * std::sqrt(numerator / (3 - detect_value));
//		}
//	}
//
//	return tQuaternion(q1, q2, q3, q4);
//}

tVector cMathUtil::QuaternionToCoef(const tQuaternion &quater)
{
    // quaternion -> vec = [x, y, z, w]
    return tVector(quater.x(), quater.y(), quater.z(), quater.w());
}

tQuaternion cMathUtil::CoefToQuaternion(const tVector &vec)
{
    // vec = [x, y, z, w] -> quaternion
    if (vec[3] > 0)
        return tQuaternion(vec[3], vec[0], vec[1], vec[2]);
    else
        return tQuaternion(-vec[3], -vec[0], -vec[1], -vec[2]);
}

tQuaternion cMathUtil::AxisAngleToQuaternion(const tVector &angvel)
{
    double theta = angvel.norm();
    double theta_2 = theta / 2;
    double cos_theta_2 = std::cos(theta_2), sin_theta_2 = std::sin(theta_2);

    tVector norm_angvel = angvel.normalized();
    return tQuaternion(cos_theta_2, norm_angvel[0] * sin_theta_2,
                       norm_angvel[1] * sin_theta_2,
                       norm_angvel[2] * sin_theta_2);
}

// tVector cMathUtil::QuaternionToAxisAngle(const tQuaternion & quater)
//{
//	/* 	quater = [w, x, y, z]
//			w = cos(theta / 2)
//			x = ax * sin(theta/2)
//			y = ay * sin(theta/2)
//			z = az * sin(theta/2)
//		axis angle = theta * [ax, ay, az, 0]
//	*/
//	tVector axis_angle = tVector::Zero();
//
//	double theta = 2 * std::acos(quater.w());
//
//	if (theta < 1e-4) return tVector::Zero();
//
//	//std::cout << theta << " " << std::sin(theta / 2) << std::endl;
//	double ax = quater.x() / std::sin(theta / 2),
//		ay = quater.y() / std::sin(theta / 2),
//		az = quater.z() / std::sin(theta / 2);
//	return theta * tVector(ax, ay, az, 0);
//}

tVector cMathUtil::CalcAngularVelocity(const tQuaternion &old_rot,
                                       const tQuaternion &new_rot,
                                       double timestep)
{
    tQuaternion trans = new_rot * old_rot.conjugate();
    double theta = std::acos(trans.w()) * 2; // std::acos() output range [0, pi]
    if (true == std::isnan(theta))
        return tVector::Zero(); // theta = nan, when w = 1. Omega = 0, 0, 0

    if (theta > 2 * M_PI - theta)
    {
        // theta = theta - 2*pi
        theta = theta - 2 * M_PI; // -pi - pi
        trans.coeffs().segment(0, 3) *= -1;
    }
    else if (std::abs(theta) < 1e-10)
    {
        return tVector::Zero();
    }
    tVector vel = tVector::Zero();
    double coef = theta / (sin(theta / 2) * timestep);
    vel.segment(0, 3) = trans.coeffs().segment(0, 3) * coef;
    return vel;
}

tVector cMathUtil::CalcAngularVelocityFromAxisAngle(const tQuaternion &old_rot,
                                                    const tQuaternion &new_rot,
                                                    double timestep)
{
    std::cout << "cMathUtil::CalcAngularVelocityFromAxisAngle: this func "
                 "hasn't been well-tested, call another one\n";
    exit(1);
    tVector old_aa = cMathUtil::QuaternionToAxisAngle(old_rot),
            new_aa = cMathUtil::QuaternionToAxisAngle(new_rot);
    return (new_aa - old_aa) / timestep;
}

// tVector cMathUtil::QuatRotVec(const tQuaternion & quater, const tVector &
// vec)
//{
//	tVector res = tVector::Zero();
//	res.segment(0, 3) = quater * vec.segment(0, 3);
//	return res;
//}

tVector cMathUtil::QuaternionToEulerAngles(const tQuaternion &q,
                                           const eRotationOrder &order)
{
    tVector res = tVector::Zero();
    double w = q.w(), x = q.x(), y = q.y(), z = q.z();

    // handle the zero quaternion
    if (order == eRotationOrder::XYZ)
    {
        res[0] = std::atan2(2 * (w * x + y * z), 1 - 2 * (x * x + y * y));
        res[1] = std::asin(2 * (w * y - z * x));
        res[2] = std::atan2(2 * (w * z + x * y), 1 - 2 * (y * y + z * z));
        // SIM_INFO("w {} x {} y {} z {}", w, x, y, z);

        // std::cout << "euler angle = " << res.transpose() << std::endl;
    }
    else if (order == eRotationOrder::ZYX)
    {
        res[0] = std::atan2(2 * (w * x - y * z), 1 - 2 * (x * x + y * y));
        res[1] = std::asin(2 * (w * y + z * x));
        res[2] = std::atan2(2 * (w * z - x * y), 1 - 2 * (y * y + z * z));
    }
    else
    {
        std::cout << "[error] tVector cMathUtil::QuaternionToEulerAngles "
                     "Unsupported rotation order = "
                  << order;
        exit(1);
    }
    return res;
}

tQuaternion cMathUtil::EulerAnglesToQuaternion(const tVector &vec,
                                               const eRotationOrder &order)
{
    tQuaternion q[3];
    for (int i = 0; i < 3; i++)
    {
        tVector axis = tVector::Zero();
        axis[i] = 1.0;

        double theta_2 = vec[i] / 2.0;
        axis = axis * std::sin(theta_2);
        axis[3] = std::cos(theta_2);

        q[i] = tQuaternion(axis[3], axis[0], axis[1], axis[2]);
    }

    tQuaternion res;
    if (order == eRotationOrder::XYZ)
    {
        res = q[2] * q[1] * q[0];
    }
    else if (order == eRotationOrder::ZYX)
    {
        res = q[0] * q[1] * q[2];
    }

    res.normalize();
    if (res.w() < 0)
        res = cMathUtil::MinusQuaternion(res);
    return res;
}

tQuaternion cMathUtil::MinusQuaternion(const tQuaternion &quad)
{
    return tQuaternion(-quad.w(), -quad.x(), -quad.y(), -quad.z());
}

tMatrix cMathUtil::EulerAnglesToRotMat(const tVector &euler,
                                       const eRotationOrder &order)
{
    // input euler angles: the rotation theta from parent to local
    // output rot mat: a rot mat that can convert a vector FROM LOCAL FRAME TO
    // PARENT FRAME
    double x = euler[0], y = euler[1], z = euler[2];
    tMatrix mat = tMatrix::Identity();
    if (order == eRotationOrder::XYZ)
    {
        tMatrix x_mat, y_mat, z_mat;
        x_mat = cMathUtil::EulerAngleRotmatX(x);
        y_mat = cMathUtil::EulerAngleRotmatY(y);
        z_mat = cMathUtil::EulerAngleRotmatZ(z);
        mat = z_mat * y_mat * x_mat;
    }
    else if (order == eRotationOrder::ZYX)
    {
        tMatrix x_mat, y_mat, z_mat;
        x_mat = cMathUtil::EulerAngleRotmatX(x);
        y_mat = cMathUtil::EulerAngleRotmatY(y);
        z_mat = cMathUtil::EulerAngleRotmatZ(z);
        mat = x_mat * y_mat * z_mat;
    }
    else
    {
        std::cout << "[error] cMathUtil::EulerAnglesToRotMat(const "
                     "tVector& euler): Unsupported rotation order"
                  << std::endl;
        exit(1);
    }
    return mat;
}

tMatrix cMathUtil::EulerAnglesToRotMatDot(const tVector &euler,
                                          const eRotationOrder &order)
{
    double x = euler[0], y = euler[1], z = euler[2];
    tMatrix mat = tMatrix::Identity();
    if (order == eRotationOrder::XYZ)
    {
        tMatrix Rz = cMathUtil::EulerAngleRotmatZ(z),
                Ry = cMathUtil::EulerAngleRotmatY(y),
                Rx = cMathUtil::EulerAngleRotmatX(x);
        tMatrix Rz_dot = cMathUtil::EulerAngleRotmatdZ(z),
                Ry_dot = cMathUtil::EulerAngleRotmatdY(y),
                Rx_dot = cMathUtil::EulerAngleRotmatdX(x);
        mat = Rz * Ry * Rx_dot + Rz_dot * Ry * Rx + Rz * Ry_dot * Rx;
    }
    else if (order == eRotationOrder::ZYX)
    {
        tMatrix Rz = EulerAngleRotmatZ(z), Ry = EulerAngleRotmatY(y),
                Rx = EulerAngleRotmatX(x);
        tMatrix Rz_dot = EulerAngleRotmatdZ(z), Ry_dot = EulerAngleRotmatdY(y),
                Rx_dot = EulerAngleRotmatdX(x);
        mat = Rx * Ry * Rz_dot + Rx_dot * Ry * Rz + Rx * Ry_dot * Rz;
    }
    else
    {
        std::cout << "[error] cMathUtil::EulerAnglesToRotMatDot(const "
                     "tVector& euler): Unsupported rotation order"
                  << std::endl;
        exit(1);
    }
    return mat;
}

tVector cMathUtil::AngularVelToqdot(const tVector &omega, const tVector &cur_q,
                                    const eRotationOrder &order)
{
    // w = Jw * q'
    // q' = (Jw)^{-1} * omega
    //[w] = R' * R^T

    // step1: get Jw
    // please read P8 formula (30) in C.K Liu's tutorial "A Quick Tutorial on
    // Multibody Dynamics" for more details
    double x = cur_q[0], y = cur_q[1], z = cur_q[2];
    tMatrix Rx = cMathUtil::EulerAngleRotmatX(x),
            Ry = cMathUtil::EulerAngleRotmatY(y),
            Rz = cMathUtil::EulerAngleRotmatZ(z);
    tMatrix Rx_dotx = cMathUtil::EulerAngleRotmatdX(x),
            Ry_doty = cMathUtil::EulerAngleRotmatdY(y),
            Rz_dotz = cMathUtil::EulerAngleRotmatdZ(z);

    if (order == eRotationOrder::XYZ)
    {
        tMatrix R = Rz * Ry * Rx;
        tMatrix dR_dx = Rz * Ry * Rx_dotx, dR_dy = Rz * Ry_doty * Rx,
                dR_dz = Rz_dotz * Ry * Rx;
        tMatrix x_col_mat = dR_dx * R.transpose(),
                y_col_mat = dR_dy * R.transpose(),
                z_col_mat = dR_dz * R.transpose();
        tVector x_col = cMathUtil::SkewMatToVector(x_col_mat);
        tVector y_col = cMathUtil::SkewMatToVector(y_col_mat);
        tVector z_col = cMathUtil::SkewMatToVector(z_col_mat);
        Eigen::Matrix3d Jw = Eigen::Matrix3d::Zero();
        Jw.block(0, 0, 3, 1) = x_col.segment(0, 3);
        Jw.block(0, 1, 3, 1) = y_col.segment(0, 3);
        Jw.block(0, 2, 3, 1) = z_col.segment(0, 3);
        tVector res = tVector::Zero();
        res.segment(0, 3) = Jw.inverse() * omega.segment(0, 3);
        return res;
    }
    else if (order == eRotationOrder::ZYX)
    {
        tMatrix R = Rx * Ry * Rz;
        tMatrix dR_dx = Rx_dotx * Ry * Rz, dR_dy = Rx * Ry_doty * Rz,
                dR_dz = Rx * Ry * Rz_dotz;
        tMatrix x_col_mat = dR_dx * R.transpose(),
                y_col_mat = dR_dy * R.transpose(),
                z_col_mat = dR_dz * R.transpose();
        tVector x_col = cMathUtil::SkewMatToVector(x_col_mat);
        tVector y_col = cMathUtil::SkewMatToVector(y_col_mat);
        tVector z_col = cMathUtil::SkewMatToVector(z_col_mat);
        Eigen::Matrix3d Jw = Eigen::Matrix3d::Zero();
        Jw.block(0, 0, 3, 1) = x_col.segment(0, 3);
        Jw.block(0, 1, 3, 1) = y_col.segment(0, 3);
        Jw.block(0, 2, 3, 1) = z_col.segment(0, 3);
        tVector res = tVector::Zero();
        res.segment(0, 3) = Jw.inverse() * omega.segment(0, 3);
        return res;
    }
    else
    {

        std::cout << "[error] cMathUtil::AngularVelToqdot: Unsupported "
                     "rotation order"
                  << std::endl;
        exit(1);
    }
}

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

// Nx3 friction cone
// each row is a direction now
tMatrixXd cMathUtil::ExpandFrictionCone(int num_friction_dirs,
                                        const tVector &normal_)
{
    // 1. check the input
    tVector normal = normal_;
    normal[3] = 0;
    normal.normalize();
    if (normal.norm() < 1e-6)
    {
        std::cout << "[error] ExpandFrictionCone normal = "
                  << normal_.transpose() << std::endl;
        exit(0);
    }

    // 2. generate a standard friction cone
    tMatrixXd D = tMatrixXd::Zero(4, num_friction_dirs);
    double gap = 2 * M_PI / num_friction_dirs;
    for (int i = 0; i < num_friction_dirs; i++)
    {
        D(0, i) = std::cos(gap * i);
        D(2, i) = std::sin(gap * i);
    }

    // 3. rotate the fricition cone
    tVector Y_normal = tVector(0, 1, 0, 0);
    tVector axis = Y_normal.cross3(normal).normalized();
    double theta = std::acos(Y_normal.dot(normal)); // [0, pi]
    D = cMathUtil::RotMat(cMathUtil::AxisAngleToQuaternion(axis * theta)) * D;
    D.transposeInPlace();
    // each row is a direction now
    return D;
}
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

/**
 * \brief		Get the jacobian preconditioner P = diag(A)
 *
 */
tMatrixXd cMathUtil::JacobPreconditioner(const tMatrixXd &A)
{
    if (A.rows() != A.cols())
    {
        std::cout << "cMathUtil::JacobPreconditioner: A is not a square matrix "
                  << A.rows() << " " << A.cols() << std::endl;
        exit(1);
    }
    tVectorXd diagonal = A.diagonal();
    if (diagonal.cwiseAbs().minCoeff() < 1e-10)
    {
        std::cout
            << "cMathUtil::JacobPreconditioner: diagnoal is nearly zero for "
            << diagonal.transpose() << std::endl;
        exit(1);
    }

    return diagonal.cwiseInverse().asDiagonal();
}

tMatrix cMathUtil::EulerAngleRotmatX(double x)
{
    tMatrix m = tMatrix::Identity();

    double cosx = cos(x);
    double sinx = sin(x);

    m(0, 0) = 1;
    m(1, 1) = cosx;
    m(1, 2) = -sinx;
    m(2, 1) = sinx;
    m(2, 2) = cosx;

    return m;
}
tMatrix cMathUtil::EulerAngleRotmatY(double y)
{
    // return AngleAxisd(y, Vector3d::UnitY()).toRotationMatrix();
    tMatrix m = tMatrix::Identity();

    double cosy = cos(y);
    double siny = sin(y);

    m(1, 1) = 1;
    m(0, 0) = cosy;
    m(0, 2) = siny;
    m(2, 0) = -siny;
    m(2, 2) = cosy;
    return m;
}
tMatrix cMathUtil::EulerAngleRotmatZ(double z)
{
    // return AngleAxisd(z, Vector3d::UnitZ()).toRotationMatrix();
    tMatrix m = tMatrix::Identity();

    double cosz = cos(z);
    double sinz = sin(z);

    m(2, 2) = 1;
    m(0, 0) = cosz;
    m(0, 1) = -sinz;
    m(1, 0) = sinz;
    m(1, 1) = cosz;

    return m;
}
tMatrix cMathUtil::EulerAngleRotmatdX(double x)
{
    tMatrix output = tMatrix::Zero();

    double cosx = cos(x);
    double sinx = sin(x);

    output(1, 1) = -sinx;
    output(1, 2) = -cosx;
    output(2, 1) = cosx;
    output(2, 2) = -sinx;
    return output;
}
tMatrix cMathUtil::EulerAngleRotmatdY(double y)
{
    tMatrix output = tMatrix::Zero();
    double cosy = cos(y);
    double siny = sin(y);

    output(0, 0) = -siny;
    output(0, 2) = cosy;
    output(2, 0) = -cosy;
    output(2, 2) = -siny;
    return output;
}
tMatrix cMathUtil::EulerAngleRotmatdZ(double z)
{
    tMatrix output = tMatrix::Zero();
    double cosz = cos(z);
    double sinz = sin(z);

    output(0, 0) = -sinz;
    output(0, 1) = -cosz;
    output(1, 0) = cosz;
    output(1, 1) = -sinz;
    return output;
}

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
tMatrix cMathUtil::TransformMat(const tVector &translation,
                                const tVector &euler_xyz_orientation)
{
    tMatrix mat = cMathUtil::EulerAnglesToRotMat(euler_xyz_orientation,
                                                 eRotationOrder::XYZ);
    mat.block(0, 3, 3, 1) = translation.segment(0, 3);
    return mat;
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
    return 0.5 * ((p1 - p0).cross3(p2 - p0)).norm();
}