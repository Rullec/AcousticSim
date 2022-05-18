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

cQBendingMaterial::cQBendingMaterial() {}

void cQBendingMaterial::Init(const std::vector<tVertexPtr> &v_array,
                             const std::vector<tEdgePtr> &e_array,
                             const std::vector<tTrianglePtr> &t_array,
                             const tVector3d &bending_stiffness_warpweftbias)
{
    cBaseMaterial::Init(v_array, e_array, t_array,
                        bending_stiffness_warpweftbias);

    int num_of_v = v_array.size();
    int num_of_dof = 3 * num_of_v;
    // create the stiffness matrix
    
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
    
    // tMatrix12f ele_K;

    mEleKLst.clear();
    for (int i = 0; i < e_array.size(); i++)
    {
        auto e = e_array[i];
        if (e->mIsBoundary == true)
            continue;
        // int v0 = e->mId0, v1 = e->mId1;
        // int v2 = SelectAnotherVertex(t_array[e->mTriangleId0], v0, v1),
        //     v3 = SelectAnotherVertex(t_array[e->mTriangleId1], v0, v1);
        // int v_lst[4] = {v0, v1, v2, v3};
        // mEdgeConstraintVertexLst.push_back(tVector4i(v0, v1, v2, v3));
        tVector4i v_lst = mEdgeConstraintVertexLst[i];

        // calculate K
        tVector coef =
            CalcCotangentCoeff(v_array[v_lst[0]]->mPos.segment(0, 3),
                               v_array[v_lst[1]]->mPos.segment(0, 3),
                               v_array[v_lst[2]]->mPos.segment(0, 3),
                               v_array[v_lst[3]]->mPos.segment(0, 3));

        double prefix_coef = mEdgeK * 3.0 /
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

// std::vector<tMatrix12f> cQBendingMaterial::GetEleStiffnessMatrixLst() const
// {
//     return this->mEleKLst;
// }
// std::vector<tVector4i> cQBendingMaterial::GetEdgeConstraintVertex() const
// {
//     return mEdgeConstraintVertexLst;
// }

void cQBendingMaterial::CheckForce() {}

void cQBendingMaterial::CheckStiffnessMatrix() {}

void cQBendingMaterial::Update() {}