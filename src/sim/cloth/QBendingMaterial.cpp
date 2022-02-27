#include "QBendingMaterial.h"
#include <iostream>
#include "utils/LogUtil.h"
tVector CalculateCotangentCoeff(const tVector3d &x0, const tVector3d &x1, const tVector3d &x2,
                                const tVector3d &x3)
{
    const tVector3d &e0 = x0 - x1, &e1 = x0 - x2, &e2 = x0 - x3, &e3 = x1 - x2,
                    &e4 = x1 - x3;

    const double &e0_norm = e0.norm(), e1_norm = e1.norm(), e2_norm = e2.norm(),
                 e3_norm = e3.norm(), e4_norm = e4.norm();

    const double &t01 = std::acos(std::fabs(e0.dot(e1)) / (e0_norm * e1_norm)),
                 &t02 = std::acos(std::fabs(e0.dot(e2)) / (e0_norm * e2_norm)),
                 &t03 = std::acos(std::fabs(e0.dot(e3)) / (e0_norm * e3_norm)),
                 &t04 = std::acos(std::fabs(e0.dot(e4)) / (e0_norm * e4_norm));

    const double &c01 = 1.0 / std::tan(t01), &c02 = 1.0 / std::tan(t02),
                 &c03 = 1.0 / std::tan(t03), &c04 = 1.0 / std::tan(t04);
    return tVector(c03 + c04, c01 + c02, -c01 - c03, -c02 - c04);
}

tVectorXd Concatenate(const tVector3d &v0,
                      const tVector3d &v1,
                      const tVector3d &v2,
                      const tVector3d &v3)
{
    tVectorXd rest = tVectorXd::Zero(12);
    rest.segment(0, 3) = v0.segment(0, 3);
    rest.segment(3, 3) = v1.segment(0, 3);
    rest.segment(6, 3) = v2.segment(0, 3);
    rest.segment(9, 3) = v3.segment(0, 3);
    return rest;
}

/**
 * \brief       calcualt energy
 *      E = 0.5 * x^T * Q * x
 *      Q = 3.0 * k / (A0 + A1) * K^T * K
 * 
*/
double cQBendingMaterial::CalcEnergy(
    const tVector3d &v0,
    const tVector3d &v1,
    const tVector3d &v2,
    const tVector3d &v3, double stiffness)
{
    tVectorXd x = Concatenate(v0, v1, v2, v3);
    return 0.5 * x.transpose() * CalcStiffnessMatrix(v0, v1, v2, v3, stiffness) * x;
}

/**
 * \brief           f = -dE/dx = - Q x
*/
tVectorXd cQBendingMaterial::CalcForce(const tVector3d &v0,
                                       const tVector3d &v1,
                                       const tVector3d &v2,
                                       const tVector3d &v3, double K)
{
    return -CalcStiffnessMatrix(v0, v1, v2, v3, K) * Concatenate(v0, v1, v2, v3);
}

/**
 * \brief           stiffness matrix:
 *      Q = 3 / (A0 + A1) * K^T * K
 *      the first v0, v1 is located in the shared edge
*/
tMatrixXd cQBendingMaterial::CalcStiffnessMatrix(const tVector3d &v0,
                                                 const tVector3d &v1,
                                                 const tVector3d &v2,
                                                 const tVector3d &v3, double stiff)
{
    float area_total = cMathUtil::CalcTriangleArea3d(v0, v1, v2) +
                       cMathUtil::CalcTriangleArea3d(v0, v1, v3);
    tMatrixXd K = CalcKCoef(v0, v1, v2, v3);
    return 3 * stiff * K.transpose() * K / area_total;
}

tMatrixXd cQBendingMaterial::CalcKCoef(const tVector3d &v0,
                                       const tVector3d &v1,
                                       const tVector3d &v2,
                                       const tVector3d &v3)
{
    tVector coef = CalculateCotangentCoeff(v0, v1, v2, v3);
    tMatrixXd K = tMatrixXd::Zero(3, 12);
    for (int i = 0; i < 4; i++)
    {
        K.block(0, 3 * i, 3, 3) = tMatrix3d::Identity() * coef[i];
    }

    return K;
}

void GetIndividualVertexPos(const tVectorXd &v,
                            tVector3d &v0,
                            tVector3d &v1,
                            tVector3d &v2,
                            tVector3d &v3)
{
    SIM_ASSERT(v.size() == 12);
    v0 = v.segment(0, 3);
    v1 = v.segment(3, 3);
    v2 = v.segment(6, 3);
    v3 = v.segment(9, 3);
}

// void cQBendingMaterial::CheckForce()
// {
//     double K = 2;
//     tVectorXd v = tVectorXd::Random(12);
//     tVector3d v0, v1, v2, v3;
//     GetIndividualVertexPos(v, v0, v1, v2, v3);
//     double E_old = cQBendingMaterial::CalcEnergy(v0, v1, v2, v3, K);
//     tVectorXd f_ana = cQBendingMaterial::CalcForce(v0, v1, v2, v3, K);
//     tVectorXd f_num = f_ana;
//     double eps = 1e-5;
//     for (size_t i = 0; i < v.size(); i++)
//     {
//         v[i] += eps;
//         GetIndividualVertexPos(v, v0, v1, v2, v3);
//         double E_new = cQBendingMaterial::CalcEnergy(v0, v1, v2, v3, K);
//         f_num[i] = (E_new - E_old) / eps;
//         v[i] -= eps;
//     }
//     auto diff = f_num - f_ana;
//     double diff_norm = diff.norm();
//     std::cout << "f ana = " << f_ana.transpose() << std::endl;
//     std::cout << "f num = " << f_num.transpose() << std::endl;
//     std::cout << "f diff = " << diff.transpose() << std::endl;
//     std::cout << "f diff norm = " << diff_norm << std::endl;
// }
// void cQBendingMaterial::CheckStiffnessMatrix()
// {
// }