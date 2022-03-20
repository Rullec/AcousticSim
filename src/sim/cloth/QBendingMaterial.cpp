#include "QBendingMaterial.h"
#include "geometries/Primitives.h"
#include "utils/LogUtil.h"
#include <iostream>
#include <set>
// #define _CRTDBG_MAP_ALLOC
// #include <stdlib.h>
// #include <crtdbg.h>

// AB and AC
double CalcCotangent(const tVector3d &A, const tVector3d &B, const tVector3d &C)
{
    const tVector3d AB = B - A;
    const tVector3d AC = C - A;
    float ABAC = AB.dot(AC);
    float ABAB = AB.dot(AB);
    float ACAC = AC.dot(AC);
    return ABAC / std::sqrt(ABAB * ACAC - ABAC * ABAC);
}
tVector CalcCotangentCoeff(const tVector3d &x0, const tVector3d &x1,
                           const tVector3d &x2, const tVector3d &x3)
{
    double c01 = CalcCotangent(x0, x1, x2);
    double c02 = CalcCotangent(x0, x1, x3);
    double c03 = CalcCotangent(x1, x0, x2);
    double c04 = CalcCotangent(x1, x0, x3);

    // clamp cotangent value to 1 degree, otherwise the value will be very big
    float threshold = 1.0 / std::sin(1.0 / 180 * M_PI);
    cMathUtil::Clamp(c01, -threshold, threshold);
    cMathUtil::Clamp(c02, -threshold, threshold);
    cMathUtil::Clamp(c03, -threshold, threshold);
    cMathUtil::Clamp(c04, -threshold, threshold);

    tVector res = tVector(c03 + c04, c01 + c02, -c01 - c03, -c02 - c04);
    if (res.hasNaN() == true)
    {
        SIM_ASSERT("CalcCotangentCoeff has Nan");
        exit(1);
    }
    // std::cout << "succ to calc cotangent coef\n";
    return res;
}

// tVectorXd Concatenate(const tVector3d &v0,
//                       const tVector3d &v1,
//                       const tVector3d &v2,
//                       const tVector3d &v3)
// {
//     tVectorXd rest = tVectorXd::Zero(12);
//     rest.segment(0, 3) = v0.segment(0, 3);
//     rest.segment(3, 3) = v1.segment(0, 3);
//     rest.segment(6, 3) = v2.segment(0, 3);
//     rest.segment(9, 3) = v3.segment(0, 3);
//     return rest;
// }

// /**
//  * \brief       calcualt energy
//  *      E = 0.5 * x^T * Q * x
//  *      Q = 3.0 * k / (A0 + A1) * K^T * K
//  *
// */
// double cQBendingMaterial::CalcEnergy(
//     const tVector3d &v0,
//     const tVector3d &v1,
//     const tVector3d &v2,
//     const tVector3d &v3, double stiffness)
// {
//     tVectorXd x = Concatenate(v0, v1, v2, v3);
//     return 0.5 * x.transpose() * CalcStiffnessMatrix(v0, v1, v2, v3,
//     stiffness) * x;
// }

// /**
//  * \brief           f = -dE/dx = - Q x
// */
// tVectorXd cQBendingMaterial::CalcForce(const tVector3d &v0,
//                                        const tVector3d &v1,
//                                        const tVector3d &v2,
//                                        const tVector3d &v3, double K)
// {
//     return -CalcStiffnessMatrix(v0, v1, v2, v3, K) * Concatenate(v0, v1, v2,
//     v3);
// }

// /**
//  * \brief           stiffness matrix:
//  *      Q = 3 / (A0 + A1) * K^T * K
//  *      the first v0, v1 is located in the shared edge
// */
// tMatrixXd cQBendingMaterial::CalcStiffnessMatrix(const tVector3d &v0,
//                                                  const tVector3d &v1,
//                                                  const tVector3d &v2,
//                                                  const tVector3d &v3, double
//                                                  stiff)
// {
//     float area_total = cMathUtil::CalcTriangleArea3d(v0, v1, v2) +
//                        cMathUtil::CalcTriangleArea3d(v0, v1, v3);
//     tMatrixXd K = CalcKCoef(v0, v1, v2, v3);
//     return 3 * stiff * K.transpose() * K / area_total;
// }

// void GetIndividualVertexPos(const tVectorXd &v,
//                             tVector3d &v0,
//                             tVector3d &v1,
//                             tVector3d &v2,
//                             tVector3d &v3)
// {
//     SIM_ASSERT(v.size() == 12);
//     v0 = v.segment(0, 3);
//     v1 = v.segment(3, 3);
//     v2 = v.segment(6, 3);
//     v3 = v.segment(9, 3);
// }

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

cQBendingMaterial::cQBendingMaterial() {}

int SelectAnotherVerteix(tTrianglePtr tri, int v0, int v1)
{
    SIM_ASSERT(tri != nullptr);
    std::set<int> vid_set = {tri->mId0, tri->mId1, tri->mId2};
    // printf("[debug] select another vertex in triangle 3 vertices (%d, %d, %d)
    // besides %d %d\n", tri->mId0, tri->mId1, tri->mId2, v0, v1);
    vid_set.erase(vid_set.find(v0));
    vid_set.erase(vid_set.find(v1));
    return *vid_set.begin();
};
void cQBendingMaterial::Init(const std::vector<tVertexPtr> &v_array,
                             const std::vector<tEdgePtr> &e_array,
                             const std::vector<tTrianglePtr> &t_array,
                             const tVector3d &bending_stiffness_warpweftbias)
{
    double K = bending_stiffness_warpweftbias.mean();
    int num_of_v = v_array.size();
    int num_of_dof = 3 * num_of_v;
    // create the stiffness matrix
    mStiffnessMat.resize(num_of_dof, num_of_dof);
    std::vector<tTriplet> triplet_lst = {};
    // 1. calculate triangle area
    std::vector<double> triangle_area_lst = {};
    for (auto &t : t_array)
    {
        triangle_area_lst.push_back(cMathUtil::CalcTriangleArea3d(
            v_array[t->mId0]->mPos.segment(0, 3),
            v_array[t->mId1]->mPos.segment(0, 3),
            v_array[t->mId2]->mPos.segment(0, 3)));
    }
    // 2. calculate coef
    // mEleKLst.resize(e_array.size());
    // tMatrix12f ele_K;
    mEdgeConstraintVertexLst.clear();
    mEleKLst.clear();
    for (auto &e : e_array)
    {
        if (e->mIsBoundary == true)
            continue;
        int v0 = e->mId0, v1 = e->mId1;
        int v2 = SelectAnotherVerteix(t_array[e->mTriangleId0], v0, v1),
            v3 = SelectAnotherVerteix(t_array[e->mTriangleId1], v0, v1);
        int v_lst[4] = {v0, v1, v2, v3};
        mEdgeConstraintVertexLst.push_back(tVector4i(v0, v1, v2, v3));

        // calculate K
        tVector coef = CalcCotangentCoeff(
            v_array[v0]->mPos.segment(0, 3), v_array[v1]->mPos.segment(0, 3),
            v_array[v2]->mPos.segment(0, 3), v_array[v3]->mPos.segment(0, 3));

        double prefix_coef = K * 3.0 /
                             (triangle_area_lst[e->mTriangleId0] +
                              triangle_area_lst[e->mTriangleId1]);
        tMatrix12f eleK = tMatrix12f::Zero();
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
            {
                // handle (i, j) block 3 * 3
                double cur_coef = coef[i] * coef[j] * prefix_coef;
                // printf("[bending] cpu, coef ij %.1f, %.1f, prefix 0.1f");
                for (int k = 0; k < 3; k++)
                {
                    triplet_lst.push_back(
                        tTriplet(3 * v_lst[i] + k, 3 * v_lst[j] + k, cur_coef));
                    eleK(3 * i + k, 3 * j + k) += cur_coef;
                }
                // triplet_lst.push_back(tTriplet(3 * v_lst[i] + 1, 3 * v_lst[j]
                // + 1, cur_coef)); triplet_lst.push_back(tTriplet(3 * v_lst[i]
                // + 2, 3 * v_lst[j] + 2, cur_coef));
            }
        this->mEleKLst.push_back(eleK);
        // calcualte stiffness
        // get result
    }
    mStiffnessMat.setFromTriplets(triplet_lst.begin(), triplet_lst.end());
    // std::cout << "mStiffnessMat = \n" << mStiffnessMat << std::endl;
}

double cQBendingMaterial::CalcEnergy(const tVectorXd &xcur)
{
    return 0.5 * xcur.transpose() * mStiffnessMat * xcur;
}
tVectorXd cQBendingMaterial::CalcForce(const tVectorXd &xcur)
{
    tVectorXd total_force = -mStiffnessMat * xcur;
    // std::cout << "[mat] bending total force max = "
    //           << total_force.cwiseAbs().maxCoeff() << std::endl;
    return total_force;
}
tSparseMatd cQBendingMaterial::GetStiffnessMatrix() const
{
    return -1 * mStiffnessMat;
}

std::vector<tMatrix12f> cQBendingMaterial::GetEleStiffnessMatrixLst() const
{
    return this->mEleKLst;
}
std::vector<tVector4i> cQBendingMaterial::GetEdgeConstraintVertex() const
{
    return mEdgeConstraintVertexLst;
}