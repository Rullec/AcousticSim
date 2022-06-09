#include "utils/RotUtil.h"
#include "utils/LogUtil.h"
#include <iostream>
tMatrix cRotUtil::TranslateMat(const tVector &trans)
{
    tMatrix mat = tMatrix::Identity();
    mat(0, 3) = trans[0];
    mat(1, 3) = trans[1];
    mat(2, 3) = trans[2];
    return mat;
}

tMatrix cRotUtil::ScaleMat(double scale)
{
    return ScaleMat(tVector::Ones() * scale);
}

tMatrix cRotUtil::ScaleMat(const tVector &scale)
{
    tMatrix mat = tMatrix::Identity();
    mat(0, 0) = scale[0];
    mat(1, 1) = scale[1];
    mat(2, 2) = scale[2];
    return mat;
}

tMatrix cRotUtil::RotateMat(const tVector &euler,
                            const eRotationOrder gRotationOrder)
{
    double x = euler[0];
    double y = euler[1];
    double z = euler[2];

    double sinx = std::sin(x);
    double cosx = std::cos(x);
    double siny = std::sin(y);
    double cosy = std::cos(y);
    double sinz = std::sin(z);
    double cosz = std::cos(z);

    tMatrix mat = tMatrix::Identity();

    if (gRotationOrder == eRotationOrder::XYZ)
    {
        mat(0, 0) = cosy * cosz;
        mat(1, 0) = cosy * sinz;
        mat(2, 0) = -siny;

        mat(0, 1) = sinx * siny * cosz - cosx * sinz;
        mat(1, 1) = sinx * siny * sinz + cosx * cosz;
        mat(2, 1) = sinx * cosy;

        mat(0, 2) = cosx * siny * cosz + sinx * sinz;
        mat(1, 2) = cosx * siny * sinz - sinx * cosz;
        mat(2, 2) = cosx * cosy;
    }
    else
    {
        std::cout << "[error] cRotUtil::RotateMat(const tVector& euler): "
                     "Unsupported rotation order"
                  << std::endl;
        exit(1);
    }
    return mat;
}

tMatrix cRotUtil::RotateMat(const tVector &axis, double theta)
{
    assert(std::abs(axis.squaredNorm() - 1) < 0.0001);
    double c = std::cos(theta);
    double s = std::sin(theta);
    double x = axis[0];
    double y = axis[1];
    double z = axis[2];

    tMatrix mat;
    mat << c + x * x * (1 - c), x * y * (1 - c) - z * s,
        x * z * (1 - c) + y * s, 0, y * x * (1 - c) + z * s,
        c + y * y * (1 - c), y * z * (1 - c) - x * s, 0,
        z * x * (1 - c) - y * s, z * y * (1 - c) + x * s, c + z * z * (1 - c),
        0, 0, 0, 0, 1;

    return mat;
}

tMatrix cRotUtil::RotateMat(const tQuaternion &q)
{
    tMatrix mat = tMatrix::Identity();

    double sqw = q.w() * q.w();
    double sqx = q.x() * q.x();
    double sqy = q.y() * q.y();
    double sqz = q.z() * q.z();
    double invs = 1 / (sqx + sqy + sqz + sqw);

    mat(0, 0) = (sqx - sqy - sqz + sqw) * invs;
    mat(1, 1) = (-sqx + sqy - sqz + sqw) * invs;
    mat(2, 2) = (-sqx - sqy + sqz + sqw) * invs;

    double tmp1 = q.x() * q.y();
    double tmp2 = q.z() * q.w();
    mat(1, 0) = 2.0 * (tmp1 + tmp2) * invs;
    mat(0, 1) = 2.0 * (tmp1 - tmp2) * invs;

    tmp1 = q.x() * q.z();
    tmp2 = q.y() * q.w();
    mat(2, 0) = 2.0 * (tmp1 - tmp2) * invs;
    mat(0, 2) = 2.0 * (tmp1 + tmp2) * invs;

    tmp1 = q.y() * q.z();
    tmp2 = q.x() * q.w();
    mat(2, 1) = 2.0 * (tmp1 + tmp2) * invs;
    mat(1, 2) = 2.0 * (tmp1 - tmp2) * invs;
    return mat;
}

tMatrix cRotUtil::CrossMat(const tVector &a)
{
    tMatrix m;
    m << 0, -a[2], a[1], 0, a[2], 0, -a[0], 0, -a[1], a[0], 0, 0, 0, 0, 0, 1;
    return m;
}

tMatrix cRotUtil::InvRigidMat(const tMatrix &mat)
{
    tMatrix inv_mat = tMatrix::Zero();
    inv_mat.block(0, 0, 3, 3) = mat.block(0, 0, 3, 3).transpose();
    inv_mat.col(3) = -inv_mat * mat.col(3);
    inv_mat(3, 3) = 1;
    return inv_mat;
}

tVector cRotUtil::GetRigidTrans(const tMatrix &mat)
{
    return tVector(mat(0, 3), mat(1, 3), mat(2, 3), 0);
}

tVector cRotUtil::InvEuler(const tVector &euler,
                           const eRotationOrder gRotationOrder)
{
    if (gRotationOrder == eRotationOrder::XYZ)
    {
        tMatrix inv_mat = cRotUtil::RotateMat(tVector(1, 0, 0, 0), -euler[0]) *
                          cRotUtil::RotateMat(tVector(0, 1, 0, 0), -euler[1]) *
                          cRotUtil::RotateMat(tVector(0, 0, 1, 0), -euler[2]);
        tVector inv_euler =
            cRotUtil::RotMatToEuler(inv_mat, eRotationOrder::XYZ);
        return inv_euler;
    }
    else
    {
        std::cout << "[error] cRotUtil::InvEuler: Unsupported rotation order"
                  << std::endl;
        exit(1);
    }
}

void cRotUtil::RotMatToAxisAngle(const tMatrix &mat, tVector &out_axis,
                                 double &out_theta)
{
    double c = (mat(0, 0) + mat(1, 1) + mat(2, 2) - 1) * 0.5;
    c = cMathUtil::Clamp(c, -1.0, 1.0);

    out_theta = std::acos(c);
    if (std::abs(out_theta) < 0.00001)
    {
        out_axis = tVector(0, 0, 1, 0);
    }
    else
    {
        double m21 = mat(2, 1) - mat(1, 2);
        double m02 = mat(0, 2) - mat(2, 0);
        double m10 = mat(1, 0) - mat(0, 1);
        double denom = std::sqrt(m21 * m21 + m02 * m02 + m10 * m10);
        out_axis[0] = m21 / denom;
        out_axis[1] = m02 / denom;
        out_axis[2] = m10 / denom;
        out_axis[3] = 0;
    }
}

tVector cRotUtil::RotMatToEuler(const tMatrix &mat,
                                const eRotationOrder gRotationOrder)
{
    tVector euler;
    if (gRotationOrder == eRotationOrder::XYZ)
    {

        euler[0] = std::atan2(mat(2, 1), mat(2, 2));
        euler[1] = std::atan2(-mat(2, 0), std::sqrt(mat(2, 1) * mat(2, 1) +
                                                    mat(2, 2) * mat(2, 2)));
        euler[2] = std::atan2(mat(1, 0), mat(0, 0));
        euler[3] = 0;
    }
    else
    {
        std::cout << "[error] cRotUtil::RotateMat: Unsupported rotation order"
                  << std::endl;
        exit(1);
    }

    return euler;
}

tMatrix cRotUtil::AxisAngleToRotmat(const tVector &angvel)
{
    return cRotUtil::RotMat(AxisAngleToQuaternion(angvel));
}

tVector cRotUtil::EulerangleToAxisAngle(const tVector &euler,
                                        const eRotationOrder gRotationOrder)
{
    tVector axis = tVector::Zero();
    double angle = 0;
    cRotUtil::EulerToAxisAngle(euler, axis, angle, gRotationOrder);
    return axis * angle;
}
tQuaternion cRotUtil::RotMatToQuaternion(const tMatrix &mat)
{
    double tr = mat(0, 0) + mat(1, 1) + mat(2, 2);
    tQuaternion q;
    if (tr > 0)
    {
        double S = sqrt(tr + 1.0) * 2; // S=4*qw
        q.w() = 0.25 * S;
        q.x() = (mat(2, 1) - mat(1, 2)) / S;
        q.y() = (mat(0, 2) - mat(2, 0)) / S;
        q.z() = (mat(1, 0) - mat(0, 1)) / S;
    }
    else if ((mat(0, 0) > mat(1, 1) && (mat(0, 0) > mat(2, 2))))
    {
        double S = sqrt(1.0 + mat(0, 0) - mat(1, 1) - mat(2, 2)) * 2; // S=4*qx
        q.w() = (mat(2, 1) - mat(1, 2)) / S;
        q.x() = 0.25 * S;
        q.y() = (mat(0, 1) + mat(1, 0)) / S;
        q.z() = (mat(0, 2) + mat(2, 0)) / S;
    }
    else if (mat(1, 1) > mat(2, 2))
    {
        double S = sqrt(1.0 + mat(1, 1) - mat(0, 0) - mat(2, 2)) * 2; // S=4*qy
        q.w() = (mat(0, 2) - mat(2, 0)) / S;
        q.x() = (mat(0, 1) + mat(1, 0)) / S;
        q.y() = 0.25 * S;
        q.z() = (mat(1, 2) + mat(2, 1)) / S;
    }
    else
    {
        double S = sqrt(1.0 + mat(2, 2) - mat(0, 0) - mat(1, 1)) * 2; // S=4*qz
        q.w() = (mat(1, 0) - mat(0, 1)) / S;
        q.x() = (mat(0, 2) + mat(2, 0)) / S;
        q.y() = (mat(1, 2) + mat(2, 1)) / S;
        q.z() = 0.25 * S;
    }

    return q;
}

void cRotUtil::EulerToAxisAngle(const tVector &euler, tVector &out_axis,
                                double &out_theta,
                                const eRotationOrder gRotationOrder)
{

    if (gRotationOrder == eRotationOrder::XYZ)
    {
        double x = euler[0];
        double y = euler[1];
        double z = euler[2];

        double sinx = std::sin(x);
        double cosx = std::cos(x);
        double siny = std::sin(y);
        double cosy = std::cos(y);
        double sinz = std::sin(z);
        double cosz = std::cos(z);

        double c =
            (cosy * cosz + sinx * siny * sinz + cosx * cosz + cosx * cosy - 1) *
            0.5;
        c = cMathUtil::Clamp(c, -1.0, 1.0);

        out_theta = std::acos(c);
        if (std::abs(out_theta) < 0.00001)
        {
            out_axis = tVector(0, 0, 1, 0);
        }
        else
        {
            double m21 = sinx * cosy - cosx * siny * sinz + sinx * cosz;
            double m02 = cosx * siny * cosz + sinx * sinz + siny;
            double m10 = cosy * sinz - sinx * siny * cosz + cosx * sinz;
            double denom = std::sqrt(m21 * m21 + m02 * m02 + m10 * m10);
            out_axis[0] = m21 / denom;
            out_axis[1] = m02 / denom;
            out_axis[2] = m10 / denom;
            out_axis[3] = 0;
        }
    }
    else
    {
        std::cout << "[error] cRotUtil::EulerToAxisAngle: Unsupported "
                     "rotation order"
                  << std::endl;
        exit(1);
    }
}

tVector cRotUtil::AxisAngleToEuler(const tVector &axis, double theta)
{
    tQuaternion q = AxisAngleToQuaternion(axis, theta);
    return QuaternionToEuler(q, eRotationOrder::XYZ);
}

tMatrix cRotUtil::DirToRotMat(const tVector &dir, const tVector &up)
{
    tVector x = up.cross3(dir);
    double x_norm = x.norm();
    if (x_norm == 0)
    {
        x_norm = 1;
        x = (dir.dot(up) >= 0) ? tVector(1, 0, 0, 0) : tVector(-1, 0, 0, 0);
    }
    x /= x_norm;

    tVector y = dir.cross3(x).normalized();
    tVector z = dir;

    tMatrix mat = tMatrix::Identity();

    mat.block(0, 0, 3, 1) = x.segment(0, 3);
    mat.block(0, 1, 3, 1) = y.segment(0, 3);
    mat.block(0, 2, 3, 1) = z.segment(0, 3);
    return mat;
}

void cRotUtil::DeltaRot(const tVector &axis0, double theta0,
                        const tVector &axis1, double theta1, tVector &out_axis,
                        double &out_theta)
{
    tMatrix R0 = RotateMat(axis0, theta0);
    tMatrix R1 = RotateMat(axis1, theta1);
    tMatrix M = DeltaRot(R0, R1);
    RotMatToAxisAngle(M, out_axis, out_theta);
}

tMatrix cRotUtil::DeltaRot(const tMatrix &R0, const tMatrix &R1)
{
    return R1 * R0.transpose();
}

tQuaternion cRotUtil::EulerToQuaternion(const tVector &euler,
                                        const eRotationOrder order)
{
    tVector axis;
    double theta;
    EulerToAxisAngle(euler, axis, theta, order);
    return AxisAngleToQuaternion(axis, theta);
}

tQuaternion cRotUtil::CoefVectorToQuaternion(const tVector &coef)
{
    // coef = [x, y, z, w]
    return tQuaternion(coef[3], coef[0], coef[1], coef[2]);
}

tVector cRotUtil::QuaternionToEuler(const tQuaternion &q,
                                    const eRotationOrder gRotationOrder)
{
    if (gRotationOrder == eRotationOrder::XYZ)
    {
        double sinr = 2.0 * (q.w() * q.x() + q.y() * q.z());
        double cosr = 1.0 - 2.0 * (q.x() * q.x() + q.y() * q.y());
        double x = std::atan2(sinr, cosr);

        double sinp = 2.0 * (q.w() * q.y() - q.z() * q.x());
        double y = 0;
        if (fabs(sinp) >= 1) // north pole and south pole
        {
            y = copysign(M_PI / 2,
                         sinp); // use 90 degrees if out of range
        }
        else
        {
            y = asin(sinp);
        }

        double siny = 2.0 * (q.w() * q.z() + q.x() * q.y());
        double cosy = 1.0 - 2.0 * (q.y() * q.y() + q.z() * q.z());
        double z = std::atan2(siny, cosy);

        return tVector(x, y, z, 0);
    }
    else
    {
        std::cout << "[error] cRotUtil::QuaternionToEuler: Unsupported "
                     "rotation order"
                  << std::endl;
        exit(1);
    }
}

tQuaternion cRotUtil::AxisAngleToQuaternion(const tVector &axis, double theta)
{
    // axis must be normalized
    // std::cout << axis.transpose() << std::endl;
    SIM_ASSERT(std::fabs(axis.norm() - 1) < 1e-10 ||
               std::fabs(axis.norm()) < 1e-10);
    double c = std::cos(theta / 2);
    double s = std::sin(theta / 2);
    tQuaternion q;
    q.w() = c;
    q.x() = s * axis[0];
    q.y() = s * axis[1];
    q.z() = s * axis[2];
    if (q.w() < 0)
        q = cRotUtil::MinusQuaternion(q);
    return q;
}
tVector cRotUtil::QuaternionToAxisAngle(const tQuaternion &q)
{
    tVector out_axis;
    double out_theta;
    QuaternionToAxisAngle(q, out_axis, out_theta);

    out_axis *= out_theta;
    out_axis[3] = 0;
    return out_axis;
}

void cRotUtil::QuaternionToAxisAngle(const tQuaternion &q, tVector &out_axis,
                                     double &out_theta)
{
    out_theta = 0;
    out_axis = tVector(0, 0, 1, 0);

    tQuaternion q1 = q;
    if (q1.w() > 1)
    {
        q1.normalize();
    }

    double sin_theta = std::sqrt(1 - q1.w() * q1.w());
    if (sin_theta > 0.000001)
    {
        out_theta = 2 * std::acos(q1.w());
        out_theta = cMathUtil::NormalizeAngle(out_theta);
        out_axis = tVector(q1.x(), q1.y(), q1.z(), 0) / sin_theta;
    }
}

tMatrix cRotUtil::BuildQuaternionDiffMat(const tQuaternion &q)
{
    // it's right
    tMatrix mat;
    mat << -0.5 * q.x(), -0.5 * q.y(), -0.5 * q.z(), 0, // for w
        0.5 * q.w(), -0.5 * q.z(), 0.5 * q.y(), 0,      // for x
        0.5 * q.z(), 0.5 * q.w(), -0.5 * q.x(), 0,      // for y
        -0.5 * q.y(), 0.5 * q.x(), 0.5 * q.w(), 0;      // for z
    return mat;
}

tVector cRotUtil::CalcQuaternionVel(const tQuaternion &q0,
                                    const tQuaternion &q1, double dt)
{
    tQuaternion q_diff = cRotUtil::QuatDiff(q0, q1);
    tVector axis;
    double theta;
    QuaternionToAxisAngle(q_diff, axis, theta);
    return (theta / dt) * axis;
}

tVector cRotUtil::CalcQuaternionVelRel(const tQuaternion &q0,
                                       const tQuaternion &q1, double dt)
{
    // calculate relative rotational velocity in the coordinate frame of q0
    tQuaternion q_diff = q0.conjugate() * q1;
    tVector axis;
    double theta;
    QuaternionToAxisAngle(q_diff, axis, theta);
    return (theta / dt) * axis;
}

tQuaternion cRotUtil::VecToQuat(const tVector &v)
{
    // v format: [w, x, y, z]
    return tQuaternion(v[0], v[1], v[2], v[3]);
}

tVector cRotUtil::QuatToVec(const tQuaternion &q)
{
    // return value format : [w, x, y, z]
    return tVector(q.w(), q.x(), q.y(), q.z());
}

tQuaternion cRotUtil::QuatDiff(const tQuaternion &q0, const tQuaternion &q1)
{
    return q1 * q0.conjugate();
}

double cRotUtil::QuatDiffTheta(const tQuaternion &q0, const tQuaternion &q1)
{
    tQuaternion dq = QuatDiff(q0, q1);
    return QuatTheta(dq);
}

// given a
double cRotUtil::QuatTheta(const tQuaternion &dq)
{
    double theta = 0;
    tQuaternion q1 = dq;
    if (q1.w() > 1)
    {
        q1.normalize();
    }

    // theta = angle / 2
    double sin_theta = std::sqrt(
        1 -
        q1.w() *
            q1.w()); // sin(theta) which "theta" is the rotation angle/2 in dq
    if (sin_theta > 1e-7)
    {
        theta = 2 * std::acos(q1.w());            // this is angle now
        theta = cMathUtil::NormalizeAngle(theta); // noramlize angle
    }
    return theta;
}

// /**
//  * \brief               Calculate d(q1 * q0.conj) / dq0
//  */
// tMatrix cRotUtil::Calc_Dq1q0conj_Dq0(const tQuaternion &q0,
//                                       const tQuaternion &q1)
// {
//     double a1 = q1.w(), b1 = q1.x(), c1 = q1.y(), d1 = q1.z();
//     tMatrix deriv = tMatrix::Zero();
//     deriv.col(0) = tVector(a1, b1, c1, d1);
//     deriv.col(1) = tVector(b1, -a1, -d1, c1);
//     deriv.col(2) = tVector(c1, d1, -a1, -b1);
//     deriv.col(3) = tVector(d1, -c1, b1, -a1);
//     return deriv;
// }

// /**
//  * \brief           calculate d(Quaternion)/(daxis angle)
//  */
// tMatrix cRotUtil::Calc_DQuaternion_DAxisAngle(const tVector &aa)
// {
//     double theta = aa.norm();
//     tMatrix dQuaterniondAA = tMatrix::Zero();

//     if (std::fabs(theta) < 1e-5)
//     {
//         dQuaterniondAA.row(0) = -1 / 3 * aa.transpose();
//         dQuaterniondAA(1, 0) = 0.5;
//         dQuaterniondAA(2, 1) = 0.5;
//         dQuaterniondAA(3, 2) = 0.5;
//     }
//     else
//     {
//         dQuaterniondAA.row(0) =
//             -0.5 * std::sin(theta / 2) * aa.transpose() / theta;
//         for (int i = 0; i < 3; i++)
//         {
//             tVector daaidaa = tVector::Zero();
//             daaidaa[i] = 1.0;

//             dQuaterniondAA.row(1 + i) =
//                 (daaidaa * theta - aa[i] * aa / theta) / (theta * theta) *
//                     std::sin(theta / 2) +
//                 aa[i] / theta * std::cos(theta / 2) / (2 * theta) * aa;
//         }
//     }

//     // std::cout << "diff mat = \n" << dQuaterniondAA << std::endl;
//     dQuaterniondAA.col(3).setZero();
//     return dQuaterniondAA;
// }

// /**
//  * \brief           calculate d(quaternion)/d(euler_angles)
//  */
// tMatrixXd cRotUtil::Calc_DQuaterion_DEulerAngles(const tVector &euler_angles,
//                                                   eRotationOrder order)
// {
//     tMatrixXd dqdeuler = tMatrixXd::Zero(4, 3);
//     if (order == eRotationOrder ::XYZ)
//     {
//         double e_x = euler_angles[0], e_y = euler_angles[1],
//                e_z = euler_angles[2];
//         double cx = std::cos(e_x / 2), sx = std::sin(e_x / 2);
//         double cy = std::cos(e_y / 2), sy = std::sin(e_y / 2);
//         double cz = std::cos(e_z / 2), sz = std::sin(e_z / 2);
//         dqdeuler.col(0) = 0.5 * tVector(cx * sy * sz - cy * cz * sx,
//                                         sx * sy * sz + cx * cy * cz,
//                                         cx * cy * sz - cz * sx * sy,
//                                         -cx * cz * sy - cy * sx * sz);

//         dqdeuler.col(1) = 0.5 * tVector(cy * sx * sz - cx * cz * sy,
//                                         -cx * cy * sz - cz * sx * sy,
//                                         cx * cy * cz - sx * sy * sz,
//                                         -cy * cz * sx - cx * sy * sz);

//         dqdeuler.col(2) = 0.5 * tVector(cz * sx * sy - cx * cy * sz,
//                                         -cx * cz * sy - cy * sx * sz,
//                                         cy * cz * sx - cx * sy * sz,
//                                         sx * sy * sz + cx * cy * cz);
//     }
//     else
//     {
//         SIM_ERROR("invalid rotation order");
//     }
//     return dqdeuler;
// }

// void cRotUtil::TestCalc_DQuaterion_DEulerAngles()
// {
//     tVector euler_angles = tVector::Random();
//     tQuaternion old_qua =
//         cRotUtil::EulerAnglesToQuaternion(euler_angles, eRotationOrder::XYZ);
//     double eps = 1e-5;
//     tMatrixXd ideal_dqde = cRotUtil::Calc_DQuaterion_DEulerAngles(
//         euler_angles, eRotationOrder::XYZ);
//     // std::cout << "ideal_dqde = \n" << ideal_dqde << std::endl;
//     for (int i = 0; i < 3; i++)
//     {
//         euler_angles[i] += eps;
//         tQuaternion new_qua = cRotUtil::EulerAnglesToQuaternion(
//             euler_angles, eRotationOrder::XYZ);
//         tVector num_dqde =
//             (cRotUtil::QuatToVec(new_qua) - cRotUtil::QuatToVec(old_qua)) /
//             eps;
//         tVector ideal_dqdei = ideal_dqde.col(i);
//         tVector diff = ideal_dqdei - num_dqde;
//         if (diff.norm() > 10 * eps)
//         {
//             std::cout
//                 << "[error] TestCalc_DQuaterion_DEulerAngles fail for col "
//                 << i
//                 << std::endl;
//             std::cout << "ideal = " << ideal_dqdei.transpose() << std::endl;
//             std::cout << "num = " << num_dqde.transpose() << std::endl;
//             std::cout << "diff = " << diff.transpose() << std::endl;

//             exit(0);
//         }
//         euler_angles[i] -= eps;
//     }
//     std::cout << "[log] TestCalc_DQuaterion_DEulerAngles succ\n";
// }
// void cRotUtil::TestCalc_DQuaterniontDAxisAngle()
// {
//     tVector aa = tVector::Random();
//     aa[3] = 0;
//     tQuaternion qua = cRotUtil::AxisAngleToQuaternion(aa);
//     tMatrix dqua_daa = cRotUtil::Calc_DQuaternion_DAxisAngle(aa);
//     double eps = 1e-5;
//     for (int i = 0; i < 3; i++)
//     {
//         aa[i] += eps;
//         tQuaternion new_qua = cRotUtil::AxisAngleToQuaternion(aa);
//         tVector num_deriv_raw = (new_qua.coeffs() - qua.coeffs()) / eps;
//         tVector num_deriv;
//         num_deriv[0] = num_deriv_raw[3];
//         num_deriv.segment(1, 3) = num_deriv_raw.segment(0, 3);
//         tVector ideal_deriv = dqua_daa.col(i);
//         tVector diff = ideal_deriv - num_deriv;
//         if (diff.norm() > 10 * eps)
//         {
//             std::cout << "[error] TestDiffQuaterniontDAxisAngle fail for " <<
//             i
//                       << std::endl;
//             std::cout << i << " diff = " << diff.transpose() << std::endl;
//             std::cout << "ideal = " << ideal_deriv.transpose() << std::endl;
//             std::cout << "num = " << num_deriv.transpose() << std::endl;
//         }
//         aa[i] -= eps;
//     }
//     std::cout << "[log] TestDiffQuaterniontDAxisAngle succ\n";
// }

// void cRotUtil::TestCalc_Dq1q0conj_Dq0()
// {
//     SIM_INFO("Dq1q0conjDq0 begin test!");
//     tQuaternion q1 = tQuaternion::UnitRandom(), q0 =
//     tQuaternion::UnitRandom(); tQuaternion old_q1_q0_conj = q1 *
//     q0.conjugate(); double eps = 1e-5;

//     tMatrix deriv = cRotUtil::Calc_Dq1q0conj_Dq0(q0, q1);
//     for (int i = 0; i < 4; i++)
//     {
//         switch (i)
//         {
//         case 0:
//             q0.w() += eps;
//             break;
//         case 1:
//             q0.x() += eps;
//             break;
//         case 2:
//             q0.y() += eps;
//             break;
//         case 3:
//             q0.z() += eps;
//             break;

//         default:
//             break;
//         }
//         // q0.normalize();
//         tQuaternion new_q1_q0_conj = q1 * q0.conjugate();
//         tVector chaos_order_d =
//             (new_q1_q0_conj.coeffs() - old_q1_q0_conj.coeffs()) / eps;
//         tVector d = tVector(chaos_order_d[3], chaos_order_d[0],
//                             chaos_order_d[1], chaos_order_d[2]);

//         tVector diff = d - deriv.col(i);

//         if (diff.norm() > 10 * eps)
//         {
//             printf("[error] TestDq1q0conjDq0_experimental fail for %d\n", i);
//             std::cout << "d = " << d.transpose() << std::endl;
//             // printf("d= %.5f, %.5f, %.5f, %.5f\n", );
//             std::cout << "ideal d = " << deriv.col(i).transpose() <<
//             std::endl; std::cout << "diff = " << diff.norm() << std::endl;
//             exit(0);
//         }
//         switch (i)
//         {
//         case 0:
//             q0.w() -= eps;
//             break;
//         case 1:
//             q0.x() -= eps;
//             break;
//         case 2:
//             q0.y() -= eps;
//             break;
//         case 3:
//             q0.z() -= eps;
//             break;

//         default:
//             break;
//         }
//     }
//     printf("[log] TestDq1q0conjDq0_experimental succ\n");
// }

tQuaternion cRotUtil::VecDiffQuat(const tVector &v0, const tVector &v1)
{
    return tQuaternion::FromTwoVectors(v0.segment(0, 3), v1.segment(0, 3));
}

tVector cRotUtil::QuatRotVec(const tQuaternion &q, const tVector &dir)
{
    tVector rot_dir = tVector::Zero();
    rot_dir.segment(0, 3) = q * dir.segment(0, 3);
    return rot_dir;
}

tMatrix2d cRotUtil::RotMat2D(double angle)
{
    tMatrix2d rotmat = cRotUtil::EulerAngleRotmatZ(angle).block(0, 0, 2, 2);
    return rotmat;
}

tMatrix cRotUtil::EulerAngleRotmatX(double x)
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
tMatrix cRotUtil::EulerAngleRotmatY(double y)
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
tMatrix cRotUtil::EulerAngleRotmatZ(double z)
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
tMatrix cRotUtil::EulerAngleRotmatdX(double x)
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
tMatrix cRotUtil::EulerAngleRotmatdY(double y)
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
tMatrix cRotUtil::EulerAngleRotmatdZ(double z)
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

tVector cRotUtil::QuaternionToCoef(const tQuaternion &quater)
{
    // quaternion -> vec = [x, y, z, w]
    return tVector(quater.x(), quater.y(), quater.z(), quater.w());
}

tQuaternion cRotUtil::CoefToQuaternion(const tVector &vec)
{
    // vec = [x, y, z, w] -> quaternion
    if (vec[3] > 0)
        return tQuaternion(vec[3], vec[0], vec[1], vec[2]);
    else
        return tQuaternion(-vec[3], -vec[0], -vec[1], -vec[2]);
}

tQuaternion cRotUtil::AxisAngleToQuaternion(const tVector &angvel)
{
    double theta = angvel.norm();
    double theta_2 = theta / 2;
    double cos_theta_2 = std::cos(theta_2), sin_theta_2 = std::sin(theta_2);

    tVector norm_angvel = angvel.normalized();
    return tQuaternion(cos_theta_2, norm_angvel[0] * sin_theta_2,
                       norm_angvel[1] * sin_theta_2,
                       norm_angvel[2] * sin_theta_2);
}

// tVector cRotUtil::QuaternionToAxisAngle(const tQuaternion & quater)
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

// tVector cRotUtil::CalcAngularVelocity(const tQuaternion &old_rot,
//                                        const tQuaternion &new_rot,
//                                        double timestep)
// {
//     tQuaternion trans = new_rot * old_rot.conjugate();
//     double theta = std::acos(trans.w()) * 2; // std::acos() output range [0,
//     pi] if (true == std::isnan(theta))
//         return tVector::Zero(); // theta = nan, when w = 1. Omega = 0, 0, 0

//     if (theta > 2 * M_PI - theta)
//     {
//         // theta = theta - 2*pi
//         theta = theta - 2 * M_PI; // -pi - pi
//         trans.coeffs().segment(0, 3) *= -1;
//     }
//     else if (std::abs(theta) < 1e-10)
//     {
//         return tVector::Zero();
//     }
//     tVector vel = tVector::Zero();
//     double coef = theta / (sin(theta / 2) * timestep);
//     vel.segment(0, 3) = trans.coeffs().segment(0, 3) * coef;
//     return vel;
// }

tVector cRotUtil::CalcAngularVelocityFromAxisAngle(const tQuaternion &old_rot,
                                                   const tQuaternion &new_rot,
                                                   double timestep)
{
    std::cout << "cRotUtil::CalcAngularVelocityFromAxisAngle: this func "
                 "hasn't been well-tested, call another one\n";
    exit(1);
    tVector old_aa = cRotUtil::QuaternionToAxisAngle(old_rot),
            new_aa = cRotUtil::QuaternionToAxisAngle(new_rot);
    return (new_aa - old_aa) / timestep;
}

// tVector cRotUtil::QuatRotVec(const tQuaternion & quater, const tVector &
// vec)
//{
//	tVector res = tVector::Zero();
//	res.segment(0, 3) = quater * vec.segment(0, 3);
//	return res;
//}

tVector cRotUtil::QuaternionToEulerAngles(const tQuaternion &q,
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
        std::cout << "[error] tVector cRotUtil::QuaternionToEulerAngles "
                     "Unsupported rotation order = "
                  << order;
        exit(1);
    }
    return res;
}

tQuaternion cRotUtil::EulerAnglesToQuaternion(const tVector &vec,
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
        res = cRotUtil::MinusQuaternion(res);
    return res;
}

tQuaternion cRotUtil::MinusQuaternion(const tQuaternion &quad)
{
    return tQuaternion(-quad.w(), -quad.x(), -quad.y(), -quad.z());
}

tMatrix cRotUtil::EulerAnglesToRotMat(const tVector &euler,
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
        x_mat = cRotUtil::EulerAngleRotmatX(x);
        y_mat = cRotUtil::EulerAngleRotmatY(y);
        z_mat = cRotUtil::EulerAngleRotmatZ(z);
        mat = z_mat * y_mat * x_mat;
    }
    else if (order == eRotationOrder::ZYX)
    {
        tMatrix x_mat, y_mat, z_mat;
        x_mat = cRotUtil::EulerAngleRotmatX(x);
        y_mat = cRotUtil::EulerAngleRotmatY(y);
        z_mat = cRotUtil::EulerAngleRotmatZ(z);
        mat = x_mat * y_mat * z_mat;
    }
    else
    {
        std::cout << "[error] cRotUtil::EulerAnglesToRotMat(const "
                     "tVector& euler): Unsupported rotation order"
                  << std::endl;
        exit(1);
    }
    return mat;
}

// tMatrix cRotUtil::EulerAnglesToRotMatDot(const tVector &euler,
//                                           const eRotationOrder &order)
// {
//     double x = euler[0], y = euler[1], z = euler[2];
//     tMatrix mat = tMatrix::Identity();
//     if (order == eRotationOrder::XYZ)
//     {
//         tMatrix Rz = cRotUtil::EulerAngleRotmatZ(z),
//                 Ry = cRotUtil::EulerAngleRotmatY(y),
//                 Rx = cRotUtil::EulerAngleRotmatX(x);
//         tMatrix Rz_dot = cRotUtil::EulerAngleRotmatdZ(z),
//                 Ry_dot = cRotUtil::EulerAngleRotmatdY(y),
//                 Rx_dot = cRotUtil::EulerAngleRotmatdX(x);
//         mat = Rz * Ry * Rx_dot + Rz_dot * Ry * Rx + Rz * Ry_dot * Rx;
//     }
//     else if (order == eRotationOrder::ZYX)
//     {
//         tMatrix Rz = EulerAngleRotmatZ(z), Ry = EulerAngleRotmatY(y),
//                 Rx = EulerAngleRotmatX(x);
//         tMatrix Rz_dot = EulerAngleRotmatdZ(z), Ry_dot =
//         EulerAngleRotmatdY(y),
//                 Rx_dot = EulerAngleRotmatdX(x);
//         mat = Rx * Ry * Rz_dot + Rx_dot * Ry * Rz + Rx * Ry_dot * Rz;
//     }
//     else
//     {
//         std::cout << "[error] cRotUtil::EulerAnglesToRotMatDot(const "
//                      "tVector& euler): Unsupported rotation order"
//                   << std::endl;
//         exit(1);
//     }
//     return mat;
// }

// tVector cRotUtil::AngularVelToqdot(const tVector &omega, const tVector
// &cur_q,
//                                     const eRotationOrder &order)
// {
//     // w = Jw * q'
//     // q' = (Jw)^{-1} * omega
//     //[w] = R' * R^T

//     // step1: get Jw
//     // please read P8 formula (30) in C.K Liu's tutorial "A Quick Tutorial on
//     // Multibody Dynamics" for more details
//     double x = cur_q[0], y = cur_q[1], z = cur_q[2];
//     tMatrix Rx = cRotUtil::EulerAngleRotmatX(x),
//             Ry = cRotUtil::EulerAngleRotmatY(y),
//             Rz = cRotUtil::EulerAngleRotmatZ(z);
//     tMatrix Rx_dotx = cRotUtil::EulerAngleRotmatdX(x),
//             Ry_doty = cRotUtil::EulerAngleRotmatdY(y),
//             Rz_dotz = cRotUtil::EulerAngleRotmatdZ(z);

//     if (order == eRotationOrder::XYZ)
//     {
//         tMatrix R = Rz * Ry * Rx;
//         tMatrix dR_dx = Rz * Ry * Rx_dotx, dR_dy = Rz * Ry_doty * Rx,
//                 dR_dz = Rz_dotz * Ry * Rx;
//         tMatrix x_col_mat = dR_dx * R.transpose(),
//                 y_col_mat = dR_dy * R.transpose(),
//                 z_col_mat = dR_dz * R.transpose();
//         tVector x_col = cMathUtil::SkewMatToVector(x_col_mat);
//         tVector y_col = cMathUtil::SkewMatToVector(y_col_mat);
//         tVector z_col = cMathUtil::SkewMatToVector(z_col_mat);
//         Eigen::Matrix3d Jw = Eigen::Matrix3d::Zero();
//         Jw.block(0, 0, 3, 1) = x_col.segment(0, 3);
//         Jw.block(0, 1, 3, 1) = y_col.segment(0, 3);
//         Jw.block(0, 2, 3, 1) = z_col.segment(0, 3);
//         tVector res = tVector::Zero();
//         res.segment(0, 3) = Jw.inverse() * omega.segment(0, 3);
//         return res;
//     }
//     else if (order == eRotationOrder::ZYX)
//     {
//         tMatrix R = Rx * Ry * Rz;
//         tMatrix dR_dx = Rx_dotx * Ry * Rz, dR_dy = Rx * Ry_doty * Rz,
//                 dR_dz = Rx * Ry * Rz_dotz;
//         tMatrix x_col_mat = dR_dx * R.transpose(),
//                 y_col_mat = dR_dy * R.transpose(),
//                 z_col_mat = dR_dz * R.transpose();
//         tVector x_col = cMathUtil::SkewMatToVector(x_col_mat);
//         tVector y_col = cMathUtil::SkewMatToVector(y_col_mat);
//         tVector z_col = cMathUtil::SkewMatToVector(z_col_mat);
//         Eigen::Matrix3d Jw = Eigen::Matrix3d::Zero();
//         Jw.block(0, 0, 3, 1) = x_col.segment(0, 3);
//         Jw.block(0, 1, 3, 1) = y_col.segment(0, 3);
//         Jw.block(0, 2, 3, 1) = z_col.segment(0, 3);
//         tVector res = tVector::Zero();
//         res.segment(0, 3) = Jw.inverse() * omega.segment(0, 3);
//         return res;
//     }
//     else
//     {

//         std::cout << "[error] cRotUtil::AngularVelToqdot: Unsupported "
//                      "rotation order"
//                   << std::endl;
//         exit(1);
//     }
// }

// void cRotUtil::QuatSwingTwistDecomposition(const tQuaternion &q,
//                                             const tVector &dir,
//                                             tQuaternion &out_swing,
//                                             tQuaternion &out_twist)
// {
//     assert(std::abs(dir.norm() - 1) < 0.000001);
//     assert(std::abs(q.norm() - 1) < 0.000001);

//     tVector q_axis = tVector(q.x(), q.y(), q.z(), 0);
//     double p = q_axis.dot(dir);
//     tVector twist_axis = p * dir;
//     out_twist = tQuaternion(q.w(), twist_axis[0], twist_axis[1],
//     twist_axis[2]); out_twist.normalize(); out_swing = q *
//     out_twist.conjugate();
// }

// tQuaternion cRotUtil::ProjectQuat(const tQuaternion &q, const tVector &dir)
// {
//     assert(std::abs(dir.norm() - 1) < 0.00001);
//     tVector ref_axis = tVector::Zero();
//     int min_idx = 0;
//     dir.cwiseAbs().minCoeff(&min_idx);
//     ref_axis[min_idx] = 1;

//     tVector rot_dir0 = dir.cross3(ref_axis);
//     tVector rot_dir1 = cRotUtil::QuatRotVec(q, rot_dir0);
//     rot_dir1 -= rot_dir1.dot(dir) * dir;

//     double dir1_norm = rot_dir1.norm();
//     tQuaternion p_rot = tQuaternion::Identity();
//     if (dir1_norm > 0.0001)
//     {
//         rot_dir1 /= dir1_norm;
//         p_rot = cRotUtil::VecDiffQuat(rot_dir0, rot_dir1);
//     }
//     return p_rot;
// }

// void cRotUtil::ButterworthFilter(double dt, double cutoff,
//                                   Eigen::VectorXd &out_x)
// {
//     double sampling_rate = 1 / dt;
//     int n = static_cast<int>(out_x.size());

//     double wc = std::tan(cutoff * M_PI / sampling_rate);
//     double k1 = std::sqrt(2) * wc;
//     double k2 = wc * wc;
//     double a = k2 / (1 + k1 + k2);
//     double b = 2 * a;
//     double c = a;
//     double k3 = b / k2;
//     double d = -2 * a + k3;
//     double e = 1 - (2 * a) - k3;

//     double xm2 = out_x[0];
//     double xm1 = out_x[0];
//     double ym2 = out_x[0];
//     double ym1 = out_x[0];

//     for (int s = 0; s < n; ++s)
//     {
//         double x = out_x[s];
//         double y = a * x + b * xm1 + c * xm2 + d * ym1 + e * ym2;

//         out_x[s] = y;
//         xm2 = xm1;
//         xm1 = x;
//         ym2 = ym1;
//         ym1 = y;
//     }

//     double yp2 = out_x[n - 1];
//     double yp1 = out_x[n - 1];
//     double zp2 = out_x[n - 1];
//     double zp1 = out_x[n - 1];

//     for (int t = n - 1; t >= 0; --t)
//     {
//         double y = out_x[t];
//         double z = a * y + b * yp1 + c * yp2 + d * zp1 + e * zp2;

//         out_x[t] = z;
//         yp2 = yp1;
//         yp1 = y;
//         zp2 = zp1;
//         zp1 = z;
//     }
// }

tMatrix cRotUtil::RotMat(const tQuaternion &quater_)
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

tMatrix cRotUtil::TransformMat(const tVector &translation,
                               const tVector &euler_xyz_orientation)
{
    tMatrix mat = cRotUtil::EulerAnglesToRotMat(euler_xyz_orientation,
                                                eRotationOrder::XYZ);
    mat.block(0, 3, 3, 1) = translation.segment(0, 3);
    return mat;
}

tVector cRotUtil::CalcAxisAngleFromOneVectorToAnother(const tVector &v0_,
                                                      const tVector &v1_)
{

    tVector v0 = v0_.normalized(), v1 = v1_.normalized();
    if(v0.norm() < 1e-10 || v1.norm() < 1e-10)
    {
        return tVector::Zero();
    }
    tVector rot_axis = v0.cross3(v1);
    if ((v0 + v1).norm() < 1e-10)
    {
        // if v0 and v1 are oppo
        // 1. find a vector in their tangent space
        rot_axis = tVector::Random();
        rot_axis[3] = 0;
        double sum = v0.segment(0, 3).sum();
        for (int i = 0; i < 3; i++)
        {
            // find nonzero
            if (std::fabs(v0[i]) > 0.56) // greater than sqrt(3) / 3
            {
                //
                rot_axis[i] = -(rot_axis.dot(v0) - v0[i] * rot_axis[i]) / v0[i];
                break;
            }
        }
        rot_axis.normalize();
        rot_axis *= M_PI;
    }
    else
    {

        double theta = std::asin(rot_axis.norm()); //[-pi/2, pi/2]

        // if the angle between v0 and v1 > 90
        if (v0.dot(v1) < 0)
        {
            theta = theta > 0 ? (theta + (M_PI / 2 - theta) * 2)
                              : (theta + (-M_PI / 2 - theta) * 2);
        }
        rot_axis = rot_axis.normalized() * std::fabs(theta);
    }
    return rot_axis;
}