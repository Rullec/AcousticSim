#include "SemiCloth.h"
#include "utils/JsonUtil.h"
#include "geometries/Primitives.h"
#include <set>

cSemiCloth::cSemiCloth(int id_) : cBaseCloth(eClothType::SEMI_IMPLICIT_CLOTH, id_) {}
cSemiCloth::~cSemiCloth() {}
void cSemiCloth::Init(const Json::Value &conf)
{
    mEnableQBending = cJsonUtil::ParseAsBool(ENABLEQBENDING_KEY, conf);
    mBendingStiffness = cJsonUtil::ParseAsBool(BENDINGSTIFFNESS_KEY, conf);
    cBaseCloth::Init(conf);

    if (mEnableQBending)
        InitBendingHessian();
    CalcNodePositionVector(mXpre);
    mXcur.noalias() = mXpre;
}

#include <iostream>
void cSemiCloth::UpdatePos(double dt)
{
    tVectorXd mXnext = tVectorXd::Zero(GetNumOfFreedom());

    // 2. calculate force

    CalcIntForce(mXcur, mIntForce);
    CalcExtForce(mExtForce);
    CalcDampingForce((mXcur - mXpre) / mIdealDefaultTimestep, mDampingForce);

    // std::cout << "before x = " << mXcur.transpose() << std::endl;
    // std::cout << "fint = " << mIntForce.transpose() << std::endl;
    // std::cout << "fext = " << mExtForce.transpose() << std::endl;
    // std::cout << "fdamp = " << mDampingForce.transpose() << std::endl;
    mXnext = CalcNextPositionSemiImplicit();

    // std::cout << "mXnext = " << mXnext.transpose() << std::endl;
    mXpre.noalias() = mXcur;
    mXcur.noalias() = mXnext;
    SetPos(mXcur);
}

/**
 * \brief       Given total force, calcualte the next vertices' position
 */
tVectorXd cSemiCloth::CalcNextPositionSemiImplicit() const
{
    /*
        semi implicit
        X_next = dt2 * Minv * Ftotal + 2 * Xcur - Xpre
    */

    double dt2 = mIdealDefaultTimestep * mIdealDefaultTimestep;
    tVectorXd next_pos = dt2 * mInvMassMatrixDiag.cwiseProduct(
                                   mIntForce + mExtForce + mDampingForce) +
                         2 * mXcur - mXpre;

    return next_pos;
}

void cSemiCloth::InitGeometry(const Json::Value &conf)
{
    cBaseCloth::InitGeometry(conf);
    // int gap = mSubdivision + 1;

    // set up the vertex pos data
    // in XOY plane

    mStiffness = cJsonUtil::ParseAsDouble("stiffness", conf);
    for (auto &x : mEdgeArray)
        x->mK_spring = mStiffness;
}

/**
 * \brief           Calculate the bending hessian (constant in inextensible
 * assumption)
 *
 *
 *          E_{bending} = 0.5 xT * Q x
 *
 * Q \in [node_dof, node_dof]
 *
 * According to the discrete mean curvature normal operator, for a edge "ei"
 * with two adjoint triangles, the bending energy hessian matrix can be defined
 * as:
 *
 * "ei" connect two points: x0 and x1.
 *      triangle1: x0 x1 x2
 *      triangle2: x0 x1 x3
 *
 * Five edges in these two triangles
 *          e0 = x0 - x1
 *          e1 = x0 - x2
 *          e2 = x0 - x3
 *          e3 = x1 - x2
 *          e4 = x1 - x3
 * Four angles between adjoint edges
 *          t01 = theta<e0, e1>
 *          t02 = theta<e0, e2>
 *          t03 = theta<e0, e3>
 *          t04 = theta<e0, e4>
 *
 * cotangent values:
 *          c01 = cot(t01)
 *          c02 = cot(t02)
 *          c03 = cot(t03)
 *          c04 = cot(t04)
 *
 * Ki = [(c03​+c04​)I3, ​​(c01​+c02​)I3​​,
 * −(c01​+c03​)I3​​, −(c02​+c04​)I3​​]  \in R^{3 \times 12}
 *    = [Ki0, Ki1, Ki2, Ki3]
 * E = 0.5 x^T KiT * Ki x
 * Qi = KiT * Ki \in R^{12 \times 12}
 *    = [
 *  Ki0 * K.T
 *  Ki1 * K.T
 *  Ki2 * K.T
 *  Ki3 * K.T
 * ]  is assigned to
 * [
 *  (x0_id, x0_id) & (x0_id, x1_id) & (x0_id, x2_id) & (x0_id, x3_id) \\
 *  (x1_id, x0_id) & (x1_id, x1_id) & (x1_id, x2_id) & (x1_id, x3_id) \\
 *  (x2_id, x0_id) & (x2_id, x1_id) & (x2_id, x2_id) & (x2_id, x3_id) \\
 *  (x3_id, x0_id) & (x3_id, x1_id) & (x3_id, x2_id) & (x3_id, x3_id) \\
 * ]
 *
 * It's a symmetric matrix
 *
 */
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

tVector CalculateCotangentCoeff(const tVector &x0, tVector &x1, tVector &x2,
                                tVector &x3)
{
    const tVector &e0 = x0 - x1, &e1 = x0 - x2, &e2 = x0 - x3, &e3 = x1 - x2,
                  &e4 = x1 - x3;
    // std::cout << "e0 = " << e0.transpose() << std::endl;
    // std::cout << "e1 = " << e1.transpose() << std::endl;
    // std::cout << "e2 = " << e2.transpose() << std::endl;
    // std::cout << "e3 = " << e3.transpose() << std::endl;
    // std::cout << "e4 = " << e4.transpose() << std::endl;
    const double &e0_norm = e0.norm(), e1_norm = e1.norm(), e2_norm = e2.norm(),
                 e3_norm = e3.norm(), e4_norm = e4.norm();
    // printf("[debug] norm: e0 = %.3f, e1 = %.3f, e2 = %.3f, e3 = %.3f, e4 =
    // %.3f\n",
    //        e0_norm,
    //        e1_norm,
    //        e2_norm,
    //        e3_norm,
    //        e4_norm);
    const double &t01 = std::acos(std::fabs(e0.dot(e1)) / (e0_norm * e1_norm)),
                 &t02 = std::acos(std::fabs(e0.dot(e2)) / (e0_norm * e2_norm)),
                 &t03 = std::acos(std::fabs(e0.dot(e3)) / (e0_norm * e3_norm)),
                 &t04 = std::acos(std::fabs(e0.dot(e4)) / (e0_norm * e4_norm));
    // const double &t01 = std::acos(e0.dot(e1) / (e0_norm * e1_norm)),
    //              &t02 = std::acos(e0.dot(e2) / (e0_norm * e2_norm)),
    //              &t03 = std::acos(e0.dot(e3) / (e0_norm * e3_norm)),
    //              &t04 = std::acos(e0.dot(e4) / (e0_norm * e4_norm));
    // printf("[debug] theta: t01 = %.3f, t02 = %.3f, t03 = %.3f, t04 = %.3f\n",
    //        t01, t02, t03, t04);
    const double &c01 = 1.0 / std::tan(t01), &c02 = 1.0 / std::tan(t02),
                 &c03 = 1.0 / std::tan(t03), &c04 = 1.0 / std::tan(t04);
    return tVector(c03 + c04, c01 + c02, -c01 - c03, -c02 - c04);
}
void cSemiCloth::InitBendingHessian()
{
    // std::cout << "---------\n";
    int dof = GetNumOfFreedom();
    mBendingHessianQ.resize(dof, dof);

    std::vector<tTriplet> Q_trilet_array(0);
    for (int i = 0; i < mEdgeArray.size(); i++)
    {
        const auto &e = this->mEdgeArray[i];
        if (e->mIsBoundary == false)
        {
            // printf("[debug] bending, tri %d and tri %d, shared edge: %d\n",
            //        e->mTriangleId0, e->mTriangleId1, i);
            int vid[4] = {e->mId0, e->mId1,
                          SelectAnotherVerteix(mTriangleArray[e->mTriangleId0],
                                               e->mId0, e->mId1),
                          SelectAnotherVerteix(mTriangleArray[e->mTriangleId1],
                                               e->mId0, e->mId1)};
            // printf("[debug] bending, tri %d and tri %d, shared edge: %d,
            // total vertices: %d %d %d %d\n",
            //        e->mTriangleId0, e->mTriangleId1, i, vid[0], vid[1],
            //        vid[2], vid[3]);
            tVector cot_vec = CalculateCotangentCoeff(
                mVertexArray[vid[0]]->mPos, mVertexArray[vid[1]]->mPos,
                mVertexArray[vid[2]]->mPos, mVertexArray[vid[3]]->mPos);
            // std::cout << "cot vec = " << cot_vec.transpose() << std::endl;
            for (int row_id = 0; row_id < 4; row_id++)
                for (int col_id = 0; col_id < 4; col_id++)
                {
                    double square = 1.0 / mTriangleArray.size() * 2;
                    double value = mBendingStiffness * cot_vec[row_id] *
                                   cot_vec[col_id] * 3 / square;
                    for (int j = 0; j < 3; j++)
                    {
                        Q_trilet_array.push_back(tTriplet(
                            3 * vid[row_id] + j, 3 * vid[col_id] + j, value));
                    }
                }
        }
    }
    mBendingHessianQ.setFromTriplets(Q_trilet_array.begin(),
                                     Q_trilet_array.end());

    // std::cout << "bending Q constant done = \n"
    //           << mBendingHessianQ << std::endl;
    printf("bending Q constant done %d/%d, bending stiffnes = %.3f\n",
           Q_trilet_array.size(), mBendingHessianQ.size(), mBendingStiffness);
    // exit(0);
}