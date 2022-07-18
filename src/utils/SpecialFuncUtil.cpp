#include "utils/SpecialFuncUtil.h"
#include "utils/LogUtil.h"
#include "utils/MathUtil.h"
#include <cmath>
#include <iostream>

const doublec gI(0.0, 1.0);

int factor(int m)
{
    SIM_ASSERT(m >= 0);
    int res = 1;
    while (m > 1)
    {
        res *= m;
        m--;
    }
    return res;
}
/*
P_l^m

1. m \in [-l, l]
2. if l is negative, P_{-l}^m = P_{l-1}^m
3. if m is negative, P_l^{-m} = (-1)^m * (l-m)! / (l + m)! P_l^m
*/
double cSpecFunc::AssociatedLegendrePoly(int l, int m, double x)
{
    if (l < 0)
    {
        return cSpecFunc::AssociatedLegendrePoly(-1 * l - 1, m, x);
    }
    else if (m < 0)
    {
        SIM_ASSERT(m <= l && m >= -l);
        int negative_m = m;
        int real_m = -m;
        double coef = (real_m % 2 == 0 ? 1 : -1) * factor(l - real_m) * 1.0 /
                      factor(l + real_m);
        return coef * cSpecFunc::AssociatedLegendrePoly(l, real_m, x);
    }
    else if (x < 0)
    {
        return ((l + m) % 2 == 0 ? 1 : -1) *
               cSpecFunc::AssociatedLegendrePoly(l, m, -x);
    }
    else
    {
        // proper result
        SIM_ASSERT(m <= l && m >= -l);
        // l = 0
        double coef = std::pow(2, l) * std::pow(1 - x * x, m / 2.0);
        double val = 0;
        for (int k = m; k <= l; k++)
        {
            double coef0 = cSpecFunc::GeneralizedBinomialCoef(l, k);
            double coef1 =
                cSpecFunc::GeneralizedBinomialCoef((l + k - 1) / 2.0, l);
            double x_pow = (k == m) ? 1 : std::pow(x, k - m);
            val += factor(k) * 1.0 / factor(k - m) * x_pow * coef0 * coef1;
        }
        return coef * val;
    }
}
doublec cSpecFunc::SphericalHarmonics(int n, int m, double theta, double phi)
{
    int m_abs =  m > 0 ? m : -m;
    SIM_ASSERT(n >= 0);
    SIM_ASSERT(m_abs <= n);
    double coef =
        (2 * n + 1) * factor(n - m_abs) * 1.0 / factor(n + m_abs) / (4 * M_PI);
    double asso_legendre =
        cSpecFunc::AssociatedLegendrePoly(n, m, std::cos(theta));
    doublec e_imphi(0, m_abs * phi);
    e_imphi = std::exp(e_imphi);
    return coef * asso_legendre * e_imphi;
}

doublec cSpecFunc::SecondHankel(int n, double val)
{
    SIM_ASSERT(n >= 0);
    if (n == 0)
    {
        return gI / val * std::exp(-gI * val);
    }
    else if (n == 1)
    {
        return -1.0 * (val - gI) / (val * val) * std::exp(-gI * val);
    }
    else
    {
        return (2 * n + 1) / (val)*cSpecFunc::SecondHankel(n - 1, val) -
               cSpecFunc::SecondHankel(n - 2, val);
    }
}

doublec cSpecFunc::MultipoleCartesian(int n, int m, double x, double y,
                                      double z)
{
    return 0;
}

doublec cSpecFunc::MultipoleSpherical(int n, int m, double k, double r,
                                      double theta, double phi)
{
    return cSpecFunc::SecondHankel(n, k * r) *
           cSpecFunc::SphericalHarmonics(n, m, theta, phi);
}

void cSpecFunc::Test()
{
    // SIM_ASSERT(cSpecFunc::TestAssociatedLegendre());
    // SIM_ASSERT(cSpecFunc::TestCartersian2Spherical());
    // cSpecFunc::TestAssoLegendreGrad();
    cSpecFunc::TestHankelGrad();
}

tVector3d cSpecFunc::ConvertCartesianToSpherical(const tVector3d &cartesian_pos)
{
    /*
    all theta and phi are [-pi, pi]
    */

    double r = cartesian_pos.norm();
    if (r == 0)
    {
        return tVector3d(0, 0, 0);
    }
    else
    {
        // phi = arctan(y / x)
        double x = cartesian_pos[0];
        double y = cartesian_pos[1];
        double phi = std::atan2(y, x);
        double sin_theta = 0;
        double cos_phi = std::cos(phi);
        double sin_phi = std::sin(phi);
        if (cos_phi != 0)
        {
            sin_theta = x / (r * cos_phi);
        }
        else
        {
            sin_theta = y / (r * sin_phi);
        }
        double cos_theta = cartesian_pos[2] / r;
        double theta = std::atan2(sin_theta, cos_theta);
        // sin(theta) = x / (r * cos(phi))
        // cos(theta) = z / r
        return tVector3d(r, theta, phi);
    }
}

tVector3d cSpecFunc::ConvertSphericalToCartesian(const tVector3d &sph_pos)
{
    double r = sph_pos[0];
    double theta = sph_pos[1];
    double phi = sph_pos[2];
    return tVector3d(r * std::sin(theta) * std::cos(phi),
                     r * std::sin(theta) * std::sin(phi), r * std::cos(theta));
}

double cSpecFunc::GeneralizedBinomialCoef(double alpha, int k)
{
    if (k == 0)
    {
        return 1;
    }
    else if (k > 0)
    {
        double fenzi = alpha;
        double fenmu = 1.0;
        for (int i = 2; i <= k; i++)
        {
            fenzi *= (alpha - (i - 1));
            fenmu *= i;
        }
        return fenzi / fenmu;
    }
    else
    {
        printf("[error] k %d, illegal\n", k);
        exit(1);
    }
}

bool cSpecFunc::TestAssociatedLegendre()
{
#ifndef __APPLE__
    for (int n = 0; n < 5; n++)
    {
        for (int m = -n; m <= n; m++)
        {
            double x = cMathUtil::RandDouble(-1, 1);
            // double x = 0.5;
            double std_val = std::assoc_legendref(n, m, x);
            double self_val = cSpecFunc::AssociatedLegendrePoly(n, m, x);
            double diff = std::fabs(std_val - self_val);
            if (diff > 1e-5)
            {
                printf(
                    "[error] n %d m %d x %.3f std %.3f self %.3f, diff %.1e\n",
                    n, m, x, std_val, self_val, diff);
                return false;
            }
        }
    }
#endif
    printf("[log] test associated legendre poly succ\n");
    return true;
}

bool cSpecFunc::TestCartersian2Spherical()
{
    for (int i = 0; i < 10; i++)
    {
        // tVector3d pos = tVector3d::Random();
        tVector3d pos(1, 1, 1);
        tVector3d sph = cSpecFunc::ConvertCartesianToSpherical(pos);
        tVector3d pos_new = cSpecFunc::ConvertSphericalToCartesian(sph);
        tVector3d diff = (pos_new - pos).cwiseAbs();
        if (diff.norm() > 1e-6)
        {
            std::cout << "[error] test cartesian to spherical!\n";
            std::cout << "raw pos " << pos.transpose() << std::endl;
            std::cout << "sph " << sph.transpose() << std::endl;
            std::cout << "new pos " << pos_new.transpose() << std::endl;
            std::cout << "diff " << diff.transpose() << std::endl;
            return false;
        }
    }
    printf("[log] test convertion between cartesian and spherical succ\n");
    return true;
}

/*
dP_n^m / du =
    (m * u * P_n^m - (n + m) * (n - m + 1) * sqrt(1 -u * u) P_{n}^{m-1})
    /
    (1 - u * u)
*/
double cSpecFunc::AssoLegendreGrad(int n, int m, double u)
{
    if (m == 0 && n == 0)
        return 0;
    else if (m < 0)
    {
        // if m < 0, *= -1
        m *= -1;
        double coef =
            (m % 2 == 0 ? 1 : -1) * factor(n - m) * 1.0 / factor(n + m);
        return coef * cSpecFunc::AssoLegendreGrad(n, m, u);
    }
    else
    {
        SIM_ASSERT(std::fabs(m) <= std::fabs(n));
        SIM_ASSERT(std::fabs(m - 1) <= std::fabs(n));
        double part1 = m * u * cSpecFunc::AssociatedLegendrePoly(n, m, u);
        double part2 = -(n + m) * (n - m + 1) * std::sqrt(1 - u * u) *
                       cSpecFunc::AssociatedLegendrePoly(n, m - 1, u);

        return (part1 + part2) / (1 - u * u);
    }
}

tVector2cd cSpecFunc::SphericalHarmonicsGrad() { return tVector2cd::Zero(); }

/*

dh_n(z) = 0.5 * (
        h_{n-1}(z)
        - (h_n(z) + z * h_{n+1}(z) ) / z
        )
*/
doublec cSpecFunc::SecondHankelGrad(int n, double val)
{

    // 1. val == 0
    if (val == 0)
    {
        // change val to 1e-6 for stability
        val = 1e-6;
    }

    // 2. if n == 0 ?
    if (n == 0)
    {
        return (val - gI) / (val * val) * std::exp(-val * gI);
    }
    else
    {
        return 0.5 * (cSpecFunc::SecondHankel(n - 1, val) -
                      (cSpecFunc::SecondHankel(n, val) +
                       val * cSpecFunc::SecondHankel(n + 1, val)) /
                          val);
    }
}

bool cSpecFunc::TestAssoLegendreGrad()
{
    for (int n = 0; n < 4; n++)
    {
        for (int m = -n; m <= n; m++)
        {
            double x = cMathUtil::RandDouble(-1, 1);
            double grad_ana = cSpecFunc::AssoLegendreGrad(n, m, x);
            double eps = 1e-8;
            double grad_num =
                (cSpecFunc::AssociatedLegendrePoly(n, m, x + eps) -
                 cSpecFunc::AssociatedLegendrePoly(n, m, x)) /
                eps;
            double grad_diff = std::fabs(grad_ana - grad_num);
            if (grad_diff > 1e-4)
            {
                printf("[error] n %d m %d x %.3f asso legendre num grad %.3f, "
                       "ana grad "
                       "%.3f, diff %.3e\n",
                       n, m, x, grad_num, grad_ana, grad_diff);
                return false;
            }
        }
    }
    printf("[log] test asso legendre grad succ\n");
    return true;
}

bool cSpecFunc::TestHankelGrad()
{
    double eps = 1e-6;
    for (int n = 0; n < 5; n++)
    {
        double x = cMathUtil::RandDouble(-10, 10);
        doublec grad_ana = cSpecFunc::SecondHankelGrad(n, x);
        doublec old_val = cSpecFunc::SecondHankel(n, x);
        doublec grad_num =
            (cSpecFunc::SecondHankel(n, x + eps) - old_val) / eps;
        doublec grad_diff = grad_num - grad_ana;
        std::cout << "n = " << n << " grad ana = " << grad_ana
                  << " grad num = " << grad_num << " grad diff = " << grad_diff
                  << std::endl;
    }
    return true;
}
bool cSpecFunc::TestSphericalHarmonicsGrad() { return true; }
