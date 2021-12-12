#include "MathUtil.h"
#include "LogUtil.h"
#include <iostream>
#include <time.h>
// const enum eRotationOrder gRotationOrder = eRotationOrder::XYZ;
const tVector gGravity = tVector(0, -9.8, 0, 0);
// const tVector gGravity = tVector(0, 0, 0, 0);
cRand cMathUtil::gRand = cRand();

bool cMathUtil::IsPoint(const tVector &vec)
{
    return std::fabs(vec[3] - 1.0) < 1e-10;
}
tVector cMathUtil::VecToPoint(const tVector &vec)
{
    tVector new_vec = vec;
    new_vec[3] = 1;
    return new_vec;
}
int cMathUtil::Clamp(int val, int min, int max)
{
    return std::max(min, std::min(val, max));
}

void cMathUtil::Clamp(const Eigen::VectorXd &min, const Eigen::VectorXd &max,
                      Eigen::VectorXd &out_vec)
{
    out_vec = out_vec.cwiseMin(max).cwiseMax(min);
}

double cMathUtil::Clamp(double val, double min, double max)
{
    return std::max(min, std::min(val, max));
}

double cMathUtil::Saturate(double val) { return Clamp(val, 0.0, 1.0); }

double cMathUtil::Lerp(double t, double val0, double val1)
{
    return (1 - t) * val0 + t * val1;
}

double cMathUtil::NormalizeAngle(double theta)
{
    // normalizes theta to be between [-pi, pi]
    double norm_theta = fmod(theta, 2 * M_PI);
    if (norm_theta > M_PI)
    {
        norm_theta = -2 * M_PI + norm_theta;
    }
    else if (norm_theta < -M_PI)
    {
        norm_theta = 2 * M_PI + norm_theta;
    }
    return norm_theta;
}

double cMathUtil::RandDouble() { return RandDouble(0, 1); }

double cMathUtil::RandDouble(double min, double max)
{
    return gRand.RandDouble(min, max);
}

double cMathUtil::RandDoubleNorm(double mean, double stdev)
{
    return gRand.RandDoubleNorm(mean, stdev);
}

double cMathUtil::RandDoubleExp(double lambda)
{
    return gRand.RandDoubleExp(lambda);
}

double cMathUtil::RandDoubleSeed(double seed)
{
    unsigned int int_seed = *reinterpret_cast<unsigned int *>(&seed);
    std::default_random_engine rand_gen(int_seed);
    std::uniform_real_distribution<double> dist;
    return dist(rand_gen);
}

int cMathUtil::RandInt() { return gRand.RandInt(); }

int cMathUtil::RandInt(int min, int max) { return gRand.RandInt(min, max); }

int cMathUtil::RandUint() { return gRand.RandUint(); }

int cMathUtil::RandUint(unsigned int min, unsigned int max)
{
    return gRand.RandUint(min, max);
}

int cMathUtil::RandIntExclude(int min, int max, int exc)
{
    return gRand.RandIntExclude(min, max, exc);
}

void cMathUtil::SeedRand(unsigned long int seed)
{
    gRand.Seed(seed);
    srand(gRand.RandInt());
}

int cMathUtil::RandSign() { return gRand.RandSign(); }

double cMathUtil::SmoothStep(double t)
{
    double val = t * t * t * (t * (t * 6 - 15) + 10);
    return val;
}

bool cMathUtil::FlipCoin(double p) { return gRand.FlipCoin(p); }

tMatrix cMathUtil::TranslateMat(const tVector &trans)
{
    tMatrix mat = tMatrix::Identity();
    mat(0, 3) = trans[0];
    mat(1, 3) = trans[1];
    mat(2, 3) = trans[2];
    return mat;
}

tMatrix cMathUtil::ScaleMat(double scale)
{
    return ScaleMat(tVector::Ones() * scale);
}

tMatrix cMathUtil::ScaleMat(const tVector &scale)
{
    tMatrix mat = tMatrix::Identity();
    mat(0, 0) = scale[0];
    mat(1, 1) = scale[1];
    mat(2, 2) = scale[2];
    return mat;
}

tMatrix cMathUtil::RotateMat(const tVector &euler,
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
        std::cout << "[error] cMathUtil::RotateMat(const tVector& euler): "
                     "Unsupported rotation order"
                  << std::endl;
        exit(1);
    }
    return mat;
}

tMatrix cMathUtil::RotateMat(const tVector &axis, double theta)
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

tMatrix cMathUtil::RotateMat(const tQuaternion &q)
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

tMatrix cMathUtil::CrossMat(const tVector &a)
{
    tMatrix m;
    m << 0, -a[2], a[1], 0, a[2], 0, -a[0], 0, -a[1], a[0], 0, 0, 0, 0, 0, 1;
    return m;
}

tMatrix cMathUtil::InvRigidMat(const tMatrix &mat)
{
    tMatrix inv_mat = tMatrix::Zero();
    inv_mat.block(0, 0, 3, 3) = mat.block(0, 0, 3, 3).transpose();
    inv_mat.col(3) = -inv_mat * mat.col(3);
    inv_mat(3, 3) = 1;
    return inv_mat;
}

tVector cMathUtil::GetRigidTrans(const tMatrix &mat)
{
    return tVector(mat(0, 3), mat(1, 3), mat(2, 3), 0);
}

tVector cMathUtil::InvEuler(const tVector &euler,
                            const eRotationOrder gRotationOrder)
{
    if (gRotationOrder == eRotationOrder::XYZ)
    {
        tMatrix inv_mat = cMathUtil::RotateMat(tVector(1, 0, 0, 0), -euler[0]) *
                          cMathUtil::RotateMat(tVector(0, 1, 0, 0), -euler[1]) *
                          cMathUtil::RotateMat(tVector(0, 0, 1, 0), -euler[2]);
        tVector inv_euler =
            cMathUtil::RotMatToEuler(inv_mat, eRotationOrder::XYZ);
        return inv_euler;
    }
    else
    {
        std::cout << "[error] cMathUtil::InvEuler: Unsupported rotation order"
                  << std::endl;
        exit(1);
    }
}

void cMathUtil::RotMatToAxisAngle(const tMatrix &mat, tVector &out_axis,
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

tVector cMathUtil::RotMatToEuler(const tMatrix &mat,
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
        std::cout << "[error] cMathUtil::RotateMat: Unsupported rotation order"
                  << std::endl;
        exit(1);
    }

    return euler;
}

tMatrix cMathUtil::AxisAngleToRotmat(const tVector &angvel)
{
    return cMathUtil::RotMat(AxisAngleToQuaternion(angvel));
}

tVector cMathUtil::EulerangleToAxisAngle(const tVector &euler,
                                         const eRotationOrder gRotationOrder)
{
    tVector axis = tVector::Zero();
    double angle = 0;
    cMathUtil::EulerToAxisAngle(euler, axis, angle, gRotationOrder);
    return axis * angle;
}
tQuaternion cMathUtil::RotMatToQuaternion(const tMatrix &mat)
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

void cMathUtil::EulerToAxisAngle(const tVector &euler, tVector &out_axis,
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
        c = Clamp(c, -1.0, 1.0);

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
        std::cout << "[error] cMathUtil::EulerToAxisAngle: Unsupported "
                     "rotation order"
                  << std::endl;
        exit(1);
    }
}

tVector cMathUtil::AxisAngleToEuler(const tVector &axis, double theta)
{
    tQuaternion q = AxisAngleToQuaternion(axis, theta);
    return QuaternionToEuler(q, eRotationOrder::XYZ);
}

tMatrix cMathUtil::DirToRotMat(const tVector &dir, const tVector &up)
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

void cMathUtil::DeltaRot(const tVector &axis0, double theta0,
                         const tVector &axis1, double theta1, tVector &out_axis,
                         double &out_theta)
{
    tMatrix R0 = RotateMat(axis0, theta0);
    tMatrix R1 = RotateMat(axis1, theta1);
    tMatrix M = DeltaRot(R0, R1);
    RotMatToAxisAngle(M, out_axis, out_theta);
}

tMatrix cMathUtil::DeltaRot(const tMatrix &R0, const tMatrix &R1)
{
    return R1 * R0.transpose();
}

tQuaternion cMathUtil::EulerToQuaternion(const tVector &euler,
                                         const eRotationOrder order)
{
    tVector axis;
    double theta;
    EulerToAxisAngle(euler, axis, theta, order);
    return AxisAngleToQuaternion(axis, theta);
}

tQuaternion cMathUtil::CoefVectorToQuaternion(const tVector &coef)
{
    // coef = [x, y, z, w]
    return tQuaternion(coef[3], coef[0], coef[1], coef[2]);
}

tVector cMathUtil::QuaternionToEuler(const tQuaternion &q,
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
        std::cout << "[error] cMathUtil::QuaternionToEuler: Unsupported "
                     "rotation order"
                  << std::endl;
        exit(1);
    }
}

tQuaternion cMathUtil::AxisAngleToQuaternion(const tVector &axis, double theta)
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
        q = cMathUtil::MinusQuaternion(q);
    return q;
}
tVector cMathUtil::QuaternionToAxisAngle(const tQuaternion &q)
{
    tVector out_axis;
    double out_theta;
    QuaternionToAxisAngle(q, out_axis, out_theta);

    out_axis *= out_theta;
    out_axis[3] = 0;
    return out_axis;
}

void cMathUtil::QuaternionToAxisAngle(const tQuaternion &q, tVector &out_axis,
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

tMatrix cMathUtil::BuildQuaternionDiffMat(const tQuaternion &q)
{
    // it's right
    tMatrix mat;
    mat << -0.5 * q.x(), -0.5 * q.y(), -0.5 * q.z(), 0, // for w
        0.5 * q.w(), -0.5 * q.z(), 0.5 * q.y(), 0,      // for x
        0.5 * q.z(), 0.5 * q.w(), -0.5 * q.x(), 0,      // for y
        -0.5 * q.y(), 0.5 * q.x(), 0.5 * q.w(), 0;      // for z
    return mat;
}

tVector cMathUtil::CalcQuaternionVel(const tQuaternion &q0,
                                     const tQuaternion &q1, double dt)
{
    tQuaternion q_diff = cMathUtil::QuatDiff(q0, q1);
    tVector axis;
    double theta;
    QuaternionToAxisAngle(q_diff, axis, theta);
    return (theta / dt) * axis;
}

tVector cMathUtil::CalcQuaternionVelRel(const tQuaternion &q0,
                                        const tQuaternion &q1, double dt)
{
    // calculate relative rotational velocity in the coordinate frame of q0
    tQuaternion q_diff = q0.conjugate() * q1;
    tVector axis;
    double theta;
    QuaternionToAxisAngle(q_diff, axis, theta);
    return (theta / dt) * axis;
}

tQuaternion cMathUtil::VecToQuat(const tVector &v)
{
    // v format: [w, x, y, z]
    return tQuaternion(v[0], v[1], v[2], v[3]);
}

tVector cMathUtil::QuatToVec(const tQuaternion &q)
{
    // return value format : [w, x, y, z]
    return tVector(q.w(), q.x(), q.y(), q.z());
}

tQuaternion cMathUtil::QuatDiff(const tQuaternion &q0, const tQuaternion &q1)
{
    return q1 * q0.conjugate();
}

double cMathUtil::QuatDiffTheta(const tQuaternion &q0, const tQuaternion &q1)
{
    tQuaternion dq = QuatDiff(q0, q1);
    return QuatTheta(dq);
}

// given a
double cMathUtil::QuatTheta(const tQuaternion &dq)
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

/**
 * \brief               Calculate d(q1 * q0.conj) / dq0
 */
tMatrix cMathUtil::Calc_Dq1q0conj_Dq0(const tQuaternion &q0,
                                      const tQuaternion &q1)
{
    double a1 = q1.w(), b1 = q1.x(), c1 = q1.y(), d1 = q1.z();
    tMatrix deriv = tMatrix::Zero();
    deriv.col(0) = tVector(a1, b1, c1, d1);
    deriv.col(1) = tVector(b1, -a1, -d1, c1);
    deriv.col(2) = tVector(c1, d1, -a1, -b1);
    deriv.col(3) = tVector(d1, -c1, b1, -a1);
    return deriv;
}

/**
 * \brief           calculate d(Quaternion)/(daxis angle)
 */
tMatrix cMathUtil::Calc_DQuaternion_DAxisAngle(const tVector &aa)
{
    double theta = aa.norm();
    tMatrix dQuaterniondAA = tMatrix::Zero();

    if (std::fabs(theta) < 1e-5)
    {
        dQuaterniondAA.row(0) = -1 / 3 * aa.transpose();
        dQuaterniondAA(1, 0) = 0.5;
        dQuaterniondAA(2, 1) = 0.5;
        dQuaterniondAA(3, 2) = 0.5;
    }
    else
    {
        dQuaterniondAA.row(0) =
            -0.5 * std::sin(theta / 2) * aa.transpose() / theta;
        for (int i = 0; i < 3; i++)
        {
            tVector daaidaa = tVector::Zero();
            daaidaa[i] = 1.0;

            dQuaterniondAA.row(1 + i) =
                (daaidaa * theta - aa[i] * aa / theta) / (theta * theta) *
                    std::sin(theta / 2) +
                aa[i] / theta * std::cos(theta / 2) / (2 * theta) * aa;
        }
    }

    // std::cout << "diff mat = \n" << dQuaterniondAA << std::endl;
    dQuaterniondAA.col(3).setZero();
    return dQuaterniondAA;
}

/**
 * \brief           calculate d(quaternion)/d(euler_angles)
 */
tMatrixXd cMathUtil::Calc_DQuaterion_DEulerAngles(const tVector &euler_angles,
                                                  eRotationOrder order)
{
    tMatrixXd dqdeuler = tMatrixXd::Zero(4, 3);
    if (order == eRotationOrder ::XYZ)
    {
        double e_x = euler_angles[0], e_y = euler_angles[1],
               e_z = euler_angles[2];
        double cx = std::cos(e_x / 2), sx = std::sin(e_x / 2);
        double cy = std::cos(e_y / 2), sy = std::sin(e_y / 2);
        double cz = std::cos(e_z / 2), sz = std::sin(e_z / 2);
        dqdeuler.col(0) = 0.5 * tVector(cx * sy * sz - cy * cz * sx,
                                        sx * sy * sz + cx * cy * cz,
                                        cx * cy * sz - cz * sx * sy,
                                        -cx * cz * sy - cy * sx * sz);

        dqdeuler.col(1) = 0.5 * tVector(cy * sx * sz - cx * cz * sy,
                                        -cx * cy * sz - cz * sx * sy,
                                        cx * cy * cz - sx * sy * sz,
                                        -cy * cz * sx - cx * sy * sz);

        dqdeuler.col(2) = 0.5 * tVector(cz * sx * sy - cx * cy * sz,
                                        -cx * cz * sy - cy * sx * sz,
                                        cy * cz * sx - cx * sy * sz,
                                        sx * sy * sz + cx * cy * cz);
    }
    else
    {
        SIM_ERROR("invalid rotation order");
    }
    return dqdeuler;
}

void cMathUtil::TestCalc_DQuaterion_DEulerAngles()
{
    tVector euler_angles = tVector::Random();
    tQuaternion old_qua =
        cMathUtil::EulerAnglesToQuaternion(euler_angles, eRotationOrder::XYZ);
    double eps = 1e-5;
    tMatrixXd ideal_dqde = cMathUtil::Calc_DQuaterion_DEulerAngles(
        euler_angles, eRotationOrder::XYZ);
    // std::cout << "ideal_dqde = \n" << ideal_dqde << std::endl;
    for (int i = 0; i < 3; i++)
    {
        euler_angles[i] += eps;
        tQuaternion new_qua = cMathUtil::EulerAnglesToQuaternion(
            euler_angles, eRotationOrder::XYZ);
        tVector num_dqde =
            (cMathUtil::QuatToVec(new_qua) - cMathUtil::QuatToVec(old_qua)) /
            eps;
        tVector ideal_dqdei = ideal_dqde.col(i);
        tVector diff = ideal_dqdei - num_dqde;
        if (diff.norm() > 10 * eps)
        {
            std::cout
                << "[error] TestCalc_DQuaterion_DEulerAngles fail for col " << i
                << std::endl;
            std::cout << "ideal = " << ideal_dqdei.transpose() << std::endl;
            std::cout << "num = " << num_dqde.transpose() << std::endl;
            std::cout << "diff = " << diff.transpose() << std::endl;

            exit(0);
        }
        euler_angles[i] -= eps;
    }
    std::cout << "[log] TestCalc_DQuaterion_DEulerAngles succ\n";
}
void cMathUtil::TestCalc_DQuaterniontDAxisAngle()
{
    tVector aa = tVector::Random();
    aa[3] = 0;
    tQuaternion qua = cMathUtil::AxisAngleToQuaternion(aa);
    tMatrix dqua_daa = cMathUtil::Calc_DQuaternion_DAxisAngle(aa);
    double eps = 1e-5;
    for (int i = 0; i < 3; i++)
    {
        aa[i] += eps;
        tQuaternion new_qua = cMathUtil::AxisAngleToQuaternion(aa);
        tVector num_deriv_raw = (new_qua.coeffs() - qua.coeffs()) / eps;
        tVector num_deriv;
        num_deriv[0] = num_deriv_raw[3];
        num_deriv.segment(1, 3) = num_deriv_raw.segment(0, 3);
        tVector ideal_deriv = dqua_daa.col(i);
        tVector diff = ideal_deriv - num_deriv;
        if (diff.norm() > 10 * eps)
        {
            std::cout << "[error] TestDiffQuaterniontDAxisAngle fail for " << i
                      << std::endl;
            std::cout << i << " diff = " << diff.transpose() << std::endl;
            std::cout << "ideal = " << ideal_deriv.transpose() << std::endl;
            std::cout << "num = " << num_deriv.transpose() << std::endl;
        }
        aa[i] -= eps;
    }
    std::cout << "[log] TestDiffQuaterniontDAxisAngle succ\n";
}

void cMathUtil::TestCalc_Dq1q0conj_Dq0()
{
    SIM_INFO("Dq1q0conjDq0 begin test!");
    tQuaternion q1 = tQuaternion::UnitRandom(), q0 = tQuaternion::UnitRandom();
    tQuaternion old_q1_q0_conj = q1 * q0.conjugate();
    double eps = 1e-5;

    tMatrix deriv = cMathUtil::Calc_Dq1q0conj_Dq0(q0, q1);
    for (int i = 0; i < 4; i++)
    {
        switch (i)
        {
        case 0:
            q0.w() += eps;
            break;
        case 1:
            q0.x() += eps;
            break;
        case 2:
            q0.y() += eps;
            break;
        case 3:
            q0.z() += eps;
            break;

        default:
            break;
        }
        // q0.normalize();
        tQuaternion new_q1_q0_conj = q1 * q0.conjugate();
        tVector chaos_order_d =
            (new_q1_q0_conj.coeffs() - old_q1_q0_conj.coeffs()) / eps;
        tVector d = tVector(chaos_order_d[3], chaos_order_d[0],
                            chaos_order_d[1], chaos_order_d[2]);

        tVector diff = d - deriv.col(i);

        if (diff.norm() > 10 * eps)
        {
            printf("[error] TestDq1q0conjDq0_experimental fail for %d\n", i);
            std::cout << "d = " << d.transpose() << std::endl;
            // printf("d= %.5f, %.5f, %.5f, %.5f\n", );
            std::cout << "ideal d = " << deriv.col(i).transpose() << std::endl;
            std::cout << "diff = " << diff.norm() << std::endl;
            exit(0);
        }
        switch (i)
        {
        case 0:
            q0.w() -= eps;
            break;
        case 1:
            q0.x() -= eps;
            break;
        case 2:
            q0.y() -= eps;
            break;
        case 3:
            q0.z() -= eps;
            break;

        default:
            break;
        }
    }
    printf("[log] TestDq1q0conjDq0_experimental succ\n");
}

tQuaternion cMathUtil::VecDiffQuat(const tVector &v0, const tVector &v1)
{
    return tQuaternion::FromTwoVectors(v0.segment(0, 3), v1.segment(0, 3));
}

tVector cMathUtil::QuatRotVec(const tQuaternion &q, const tVector &dir)
{
    tVector rot_dir = tVector::Zero();
    rot_dir.segment(0, 3) = q * dir.segment(0, 3);
    return rot_dir;
}

tQuaternion cMathUtil::MirrorQuaternion(const tQuaternion &q, eAxis axis)
{
    tQuaternion mirror_q;
    mirror_q.w() = q.w();
    mirror_q.x() = (axis == eAxisX) ? q.x() : -q.x();
    mirror_q.y() = (axis == eAxisY) ? q.y() : -q.y();
    mirror_q.z() = (axis == eAxisZ) ? q.z() : -q.z();
    return mirror_q;
}

double cMathUtil::Sign(double val) { return SignAux<double>(val); }

int cMathUtil::Sign(int val) { return SignAux<int>(val); }

double cMathUtil::AddAverage(double avg0, int count0, double avg1, int count1)
{
    double total = count0 + count1;
    return (count0 / total) * avg0 + (count1 / total) * avg1;
}

tVector cMathUtil::AddAverage(const tVector &avg0, int count0,
                              const tVector &avg1, int count1)
{
    double total = count0 + count1;
    return (count0 / total) * avg0 + (count1 / total) * avg1;
}

void cMathUtil::AddAverage(const Eigen::VectorXd &avg0, int count0,
                           const Eigen::VectorXd &avg1, int count1,
                           Eigen::VectorXd &out_result)
{
    double total = count0 + count1;
    out_result = (count0 / total) * avg0 + (count1 / total) * avg1;
}

void cMathUtil::CalcSoftmax(const Eigen::VectorXd &vals, double temp,
                            Eigen::VectorXd &out_prob)
{
    assert(out_prob.size() == vals.size());
    int num_vals = static_cast<int>(vals.size());
    double sum = 0;
    double max_val = vals.maxCoeff();
    for (int i = 0; i < num_vals; ++i)
    {
        double val = vals[i];
        val = std::exp((val - max_val) / temp);
        out_prob[i] = val;
        sum += val;
    }

    out_prob /= sum;
}

double cMathUtil::EvalGaussian(const Eigen::VectorXd &mean,
                               const Eigen::VectorXd &covar,
                               const Eigen::VectorXd &sample)
{
    assert(mean.size() == covar.size());
    assert(sample.size() == covar.size());

    Eigen::VectorXd diff = sample - mean;
    double exp_val = diff.dot(diff.cwiseQuotient(covar));
    double likelihood = std::exp(-0.5 * exp_val);

    double partition = CalcGaussianPartition(covar);
    likelihood /= partition;
    return likelihood;
}

double cMathUtil::EvalGaussian(double mean, double covar, double sample)
{
    double diff = sample - mean;
    double exp_val = diff * diff / covar;
    double norm = 1 / std::sqrt(2 * M_PI * covar);
    double likelihood = norm * std::exp(-0.5 * exp_val);
    return likelihood;
}

double cMathUtil::CalcGaussianPartition(const Eigen::VectorXd &covar)
{
    int data_size = static_cast<int>(covar.size());
    double det = covar.prod();
    double partition = std::sqrt(std::pow(2 * M_PI, data_size) * det);
    return partition;
}

double cMathUtil::EvalGaussianLogp(const Eigen::VectorXd &mean,
                                   const Eigen::VectorXd &covar,
                                   const Eigen::VectorXd &sample)
{
    int data_size = static_cast<int>(covar.size());

    Eigen::VectorXd diff = sample - mean;
    double logp = -0.5 * diff.dot(diff.cwiseQuotient(covar));
    double det = covar.prod();
    logp += -0.5 * (data_size * std::log(2 * M_PI) + std::log(det));

    return logp;
}

double cMathUtil::EvalGaussianLogp(double mean, double covar, double sample)
{
    double diff = sample - mean;
    double logp = -0.5 * diff * diff / covar;
    logp += -0.5 * (std::log(2 * M_PI) + std::log(covar));
    return logp;
}

double cMathUtil::Sigmoid(double x) { return Sigmoid(x, 1, 0); }

double cMathUtil::Sigmoid(double x, double gamma, double bias)
{
    double exp = -gamma * (x + bias);
    double val = 1 / (1 + std::exp(exp));
    return val;
}

int cMathUtil::SampleDiscreteProb(const Eigen::VectorXd &probs)
{
    assert(std::abs(probs.sum() - 1) < 0.00001);
    double rand = RandDouble();

    int rand_idx = gInvalidIdx;
    int num_probs = static_cast<int>(probs.size());
    for (int i = 0; i < num_probs; ++i)
    {
        double curr_prob = probs[i];
        rand -= curr_prob;

        if (rand <= 0)
        {
            rand_idx = i;
            break;
        }
    }
    return rand_idx;
}

tMatrix2d cMathUtil::RotMat2D(double angle)
{
    tMatrix2d rotmat = cMathUtil::EulerAngleRotmatZ(angle).block(0, 0, 2, 2);
    return rotmat;
}

/**
 * \briewf          categorical random distribution
 */
int cMathUtil::RandIntCategorical(const std::vector<double> &prop)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::discrete_distribution<> d(prop.begin(), prop.end());
    int num = d(gen);
    return num;
}