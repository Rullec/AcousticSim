#pragma once
#include "utils/EigenUtil.h"
#include <complex>
typedef std::complex<double> doublec;
/*
special functions support
Associated legendre poly
Spherical Harmonic
Spherical Hankel
Basis
*/

/*
a + bi = e^{ln(sqrt(a^2 + b^2)) + i * tan(b, a)}
*/
class cSpecialFunctionUtil
{
    static tVector3d ConvertSphericalToCartesian(double r, double theta,
                                                 double phi);
    static tVector3d ConvertCartesianToSpherical(double x, double y, double z);
    static doublec Multipole(int n, int m, double k, double r, double theta,
                             double phi);
    static doublec SphericalHankel(int n, double r);
    static doublec SphericalHarmonic(int n, int m, double theta, double phi);
    static void AssociatedLegendrePoly(int n, int m, double theta);

    static bool CheckSphericalHankel(int n, int m);
    static bool CheckSphericalHarmonic(int n, int m);
    static bool CheckAssociatedLegendrePoly(int n, int m);
};