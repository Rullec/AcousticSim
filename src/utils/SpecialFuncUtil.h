#pragma once
#include "utils/ComplexUtil.h"
#include "utils/EigenUtil.h"

class cSpecFunc
{
public:
    // 0 order evaluation
    static double AssociatedLegendrePoly(int n, int m, double theta);
    static doublec SphericalHarmonics(int n, int m, double theta, double phi);
    static doublec SecondHankel(int n, double val);
    static doublec MultipoleCartesian(int n, int m, double x, double y,
                                      double z);
    static doublec MultipoleSpherical(int n, int m, double k, double r,
                                      double theta, double phi);

    // 1 order evaluation
    static tVector3cd MultipoleCartesianGrad();
    static tVector3cd MultipoleSphericalGrad();
    static double AssoLegendreGrad(int n, int m, double theta);
    static tVector2cd SphericalHarmonicsGrad();
    static doublec SecondHankelGrad(int n, double val);

    // cartesian - spherical convertion
    static tVector3d
    ConvertCartesianToSpherical(const tVector3d &cartesian_pos);
    static tVector3d ConvertSphericalToCartesian(const tVector3d &sph_pos);

    // test
    static void Test();
    static bool TestAssociatedLegendre();
    static bool TestCartersian2Spherical();
    static bool TestMultipoleGrad();
    static bool TestAssoLegendreGrad();
    static bool TestHankelGrad();
    static bool TestSphericalHarmonicsGrad();

    // helper
    static double GeneralizedBinomialCoef(double alpha, int k);
};