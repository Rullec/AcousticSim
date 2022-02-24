#include "ImplicitCloth.h"
#include "utils/JsonUtil.h"
#include <iostream>
#include "geometries/Primitives.h"

cImplicitCloth::cImplicitCloth(int id_)
    : cBaseCloth(eClothType::IMPLICIT_CLOTH, id_)
{
}
cImplicitCloth::~cImplicitCloth() {}

void cImplicitCloth::Init(const Json::Value &conf)
{
    mMaxNewtonIters = cJsonUtil::ParseAsInt("max_newton_iters", conf);
    mStiffness = cJsonUtil::ParseAsInt("stiffness", conf);
    cBaseCloth::Init(conf);
    CalcNodePositionVector(mXpre);
    mXcur.noalias() = mXpre;

    for (auto &x : mEdgeArrayShared)
        x->mK_spring = mStiffness;
}

void cImplicitCloth::InitGeometry(const Json::Value &conf)
{
    cBaseCloth::InitGeometry(conf);
    for (auto &x : mEdgeArrayShared)
        x->mK_spring = mStiffness;
}
void cImplicitCloth::UpdatePos(double dt)
{
    // std::cout << "mInt force = " << mIntForce.transpose() << std::endl;
    // std::cout << "mExt force = " << mExtForce.transpose() << std::endl;
    // 3. forward simulation
    tVectorXd mXnext = tVectorXd::Zero(GetNumOfFreedom());

    // 2. calculate force
    CalcIntForce(mXcur, mIntForce);
    CalcExtForce(mExtForce);
    CalcDampingForce((mXcur - mXpre) / mIdealDefaultTimestep, mDampingForce);

    // std::cout << "before x = " << mXcur.transpose() << std::endl;
    // std::cout << "fint = " << mIntForce.transpose() << std::endl;
    // std::cout << "fext = " << mExtForce.transpose() << std::endl;
    // std::cout << "fdamp = " << mDampingForce.transpose() << std::endl;
    mXnext = CalcNextPositionImplicit();

    // std::cout << "mXnext = " << mXnext.transpose() << std::endl;
    mXpre.noalias() = mXcur;
    mXcur.noalias() = mXnext;
    SetPos(mXcur);
}

tVectorXd cImplicitCloth::CalcNextPositionImplicit()
{
    /*
        Let X = Xnext is what we want have, the equation is:

        G(X) = dt2 * F_int(X) - M * X + M (2 * Xcur - Xpre) + dt2 * Fext = 0

        Let x0 = Xcur as an init solution (we can certainly have some other init
       x0) Let i = 0
        1. Calc res = G(Xi), if i > max_iters or |res| < thre, to step5; else
       step 2
        2. Xi <- Xi - (dGdX).inv() * res, i+=1, to step1
        5. return res
    */

    int max_iters = mMaxNewtonIters;
    double cnvg_thre = 1e-5; // convergence threshold
    tVectorXd x0 = mXcur;
    int dof = GetNumOfFreedom();
    tVectorXd res = tVectorXd::Zero(dof);
    tSparseMat dGdx(dof, dof);
    int cur_iter = 0;

    Eigen::SparseLU<tSparseMat> solver;

    while (true)
    {
        mIntForce.setZero();
        mDampingForce.setZero();
        // step 1
        CalcGxImplicit(x0, res, mIntForce, mExtForce, mDampingForce);
        double res_norm = res.norm();
        if (cur_iter > max_iters || res_norm < cnvg_thre)
            break;

        // step 2
        // CalcdGxdxImplicit(x0, dGdx);
        CalcdGxdxImplicitSparse(x0, dGdx);
        // TestdGxdxImplicit(x0, dGdx);

        // if (cur_iter == 3)
        // {
        //     tSparseMat mat(dof, dof);
        //     CalcdGxdxImplicitSparse(x0, mat);
        //     std::cout << "diff norm = " << (mat - dGdx).norm() << std::endl;
        //     std::cout << "dense = \n"
        //               << dGdx << std::endl;
        //     std::cout << "sparse = \n"
        //               << mat << std::endl;
        //     exit(0);
        // }
        if (cur_iter == 0)
            solver.analyzePattern(dGdx);

        solver.factorize(dGdx);
        // Use the factors to solve the linear system
        // x = ;
        x0 = x0 - solver.solve(res);
        cur_iter++;
        // printf("iter %d diff %.5f\n", cur_iter, res_norm);
        if (x0.hasNaN())
        {
            std::cout << "x0 has Nan = " << x0.transpose() << std::endl;
            std::cout << "x0 init = " << mXcur.transpose() << std::endl;
            // std::cout << "dGdx = \n"
            //           << dGdx << std::endl;
            // std::cout << "dGdx inv = \n"
            //           << dGdx.inverse() << std::endl;

            exit(0);
        }
        // std::cout << "iter " << cur_iter << " res = " << res.norm()
        //           << std::endl;
        // exit(0);
    }

    // printf("[debug] newton solver done, iters %d/%d, res norm %.5f\n",
    // cur_iter,
    //        max_iters, res.norm());
    // exit(0);
    return x0;
}
void cImplicitCloth::CalcGxImplicit(const tVectorXd &xcur, tVectorXd &Gx,
                                    tVectorXd &fint_buf, tVectorXd &fext_buf,
                                    tVectorXd &fdamp_buffer) const
{
    CalcIntForce(xcur, fint_buf);
    // CalcExtForce(fext_buf);
    CalcDampingForce((mXcur - mXpre) / mIdealDefaultTimestep, fdamp_buffer);
    Gx.noalias() =
        2 * mXcur - mXpre - xcur +
        mIdealDefaultTimestep * mIdealDefaultTimestep *
            mInvMassMatrixDiag.cwiseProduct(fext_buf + fint_buf + fdamp_buffer);
}

void cImplicitCloth::CalcdGxdxImplicit(const tVectorXd &x,
                                       tMatrixXd &dGdx) const
{

    dGdx.noalias() = tMatrixXd::Zero(GetNumOfFreedom(), GetNumOfFreedom());
    int id0, id1;
    double dist;
    tVector3d pos0, pos1;
    for (auto &spr : this->mEdgeArrayShared)
    {
        id0 = spr->mId0;
        id1 = spr->mId1;

        /*
            For a "l"th spring, the node index of its two ends is "i" and "j"
            let \bar{f}_i is the force contribution to "i"th node
            let \bar{f}_j is the force contribution to "j"th node
            l_{raw} is the raw length of this spring

            Calc:
                d(\bar{f}_i) / dxi, d(\bar{f}_i) / dxj

            dfi/dxi = k * I3
                    - k * l_{raw} *
                        (
                            I_3 * |xi - xj| - (xi - xj) * (xi - xj)^T / |xi -
           xj|
                        )
                        /
                        (xi - xj)^2

            dfi/dxj = - dfi/dxi
            dfj/dxi = - dfi/dxi
            dfj/dxj = - dfi/dxj = dfidxi

            In summary:
            dFdX should be a real symmetric matrix
        */
        pos0 = x.segment(3 * id0, 3);
        pos1 = x.segment(3 * id1, 3);
        double dist = (pos0 - pos1).norm();
        const tMatrix3d &I3 = tMatrix3d::Identity(3, 3);
        tMatrix3d dfidxi =
            spr->mK_spring * I3 -
            spr->mK_spring * spr->mRawLength *
                (I3 * dist - (pos0 - pos1) * (pos0 - pos1).transpose() / dist) /
                (dist * dist);

        if (dfidxi.hasNaN() == true)
        {
            printf("df%ddx%d has Nan, id0 = %d, id1 = %d, dist = %.5f\n", id0,
                   id0, id0, id1, dist);
            std::cout << "dfidxi = \n" << dfidxi << std::endl;
            std::cout << "pos0 = " << pos0.transpose() << std::endl;
            std::cout << "pos1 = " << pos1.transpose() << std::endl;
            std::cout << "x = " << x.transpose() << std::endl;
            exit(0);
        }
        dGdx.block(3 * id0, 3 * id0, 3, 3) += dfidxi;
        dGdx.block(3 * id1, 3 * id1, 3, 3) += dfidxi;
        dGdx.block(3 * id0, 3 * id1, 3, 3) += -dfidxi;
        dGdx.block(3 * id1, 3 * id0, 3, 3) += -dfidxi;
    }

    // std::cout << "dFdx = \n" << dGdx << std::endl;
    // dGdx = dt2 * Minv * dFdX - I
    dGdx = mIdealDefaultTimestep * mIdealDefaultTimestep *
               mInvMassMatrixDiag.asDiagonal().toDenseMatrix() * dGdx -
           tMatrixXd::Identity(GetNumOfFreedom(), GetNumOfFreedom());
}
void cImplicitCloth::CalcdGxdxImplicitSparse(const tVectorXd &x,
                                             tSparseMat &dGdx) const
{
    int dof = GetNumOfFreedom();
    Eigen::SparseMatrix<double> spMat(dof, dof);
    int id0, id1;
    double dist;
    tVector3d pos0, pos1;
    double dt2 = mIdealDefaultTimestep * mIdealDefaultTimestep;
    tEigenArr<tTriplet> tri_lst(0);
    for (auto &spr : this->mEdgeArrayShared)
    {
        id0 = spr->mId0;
        id1 = spr->mId1;

        pos0 = x.segment(3 * id0, 3);
        pos1 = x.segment(3 * id1, 3);
        double dist = (pos0 - pos1).norm();
        const tMatrix3d &I3 = tMatrix3d::Identity(3, 3);
        tMatrix3d dfidxi =
            spr->mK_spring * I3 -
            spr->mK_spring * spr->mRawLength *
                (I3 * dist - (pos0 - pos1) * (pos0 - pos1).transpose() / dist) /
                (dist * dist);

        if (dfidxi.hasNaN() == true)
        {
            printf("df%ddx%d has Nan, id0 = %d, id1 = %d, dist = %.5f\n", id0,
                   id0, id0, id1, dist);
            std::cout << "dfidxi = \n" << dfidxi << std::endl;
            std::cout << "pos0 = " << pos0.transpose() << std::endl;
            std::cout << "pos1 = " << pos1.transpose() << std::endl;
            std::cout << "x = " << x.transpose() << std::endl;
            exit(0);
        }
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                double val = dfidxi(i, j);
                {
                    tri_lst.push_back(tTriplet(3 * id0 + i, 3 * id0 + j, val));
                    tri_lst.push_back(tTriplet(3 * id1 + i, 3 * id1 + j, val));
                    tri_lst.push_back(tTriplet(3 * id0 + i, 3 * id1 + j, -val));
                    tri_lst.push_back(tTriplet(3 * id1 + i, 3 * id0 + j, -val));
                }
            }
    }
    dGdx.setFromTriplets(tri_lst.begin(), tri_lst.end());
    dGdx = dt2 * mInvMassMatrixDiag.asDiagonal() * dGdx;

    for (int i = 0; i < GetNumOfFreedom(); i++)
    {
        dGdx.coeffRef(i, i) -= 1;
    }
    //     tri_lst.push_back(tTriplet(i, i, -1));
    // dGdx.resize(GetNumOfFreedom(), GetNumOfFreedom());
    // dGdx.setFromTriplets(tri_lst.begin(), tri_lst.end());
    /*
        dt2 * inv
    */
    // std::cout << "dFdx = \n" << dGdx << std::endl;
    // dGdx = dt2 * Minv * dFdX - I
    // dGdx = mIdealDefaultTimestep * mIdealDefaultTimestep *
    // mInvMassMatrixDiag.asDiagonal().toDenseMatrix()
    // *
    //            dGdx -
    //        tMatrixXd::Identity(GetNumOfFreedom(), GetNumOfFreedom());
}
void cImplicitCloth::TestdGxdxImplicit(const tVectorXd &x0,
                                       const tMatrixXd &Gx_ana)
{
    double eps = 1e-6;
    tVectorXd Gx_old, Gx_new;
    tVectorXd x_now = x0;
    CalcGxImplicit(x_now, Gx_old, mIntForce, mExtForce, mDampingForce);

    for (int i = 0; i < x0.size(); i++)
    {
        x_now[i] += eps;
        CalcGxImplicit(x_now, Gx_new, mIntForce, mExtForce, mDampingForce);
        tVectorXd dGdxi_num = (Gx_new - Gx_old) / eps;
        tVectorXd dGdxi_ana = Gx_ana.col(i);
        tVectorXd diff = dGdxi_ana - dGdxi_num;
        double diff_norm = diff.norm();
        if (diff_norm > 10 * eps)
        {
            printf("[error] error test dGdx/dx%d failed, diff norm %.5f\n", i,
                   diff_norm);
            exit(0);
        }

        x_now[i] -= eps;
    }
    SIM_INFO("test dG/dx succ");
}