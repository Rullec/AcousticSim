#include "PBDCloth.h"
#include "utils/JsonUtil.h"
#include "geometries/Primitives.h"
cPBDCloth::cPBDCloth(int id_) : cBaseCloth(eClothType::PBD_CLOTH, id_)
{
    mVcur.resize(0);
    mBendingMatrixKArray.clear();
}

cPBDCloth::~cPBDCloth() {}
void cPBDCloth::Init(const Json::Value &conf)
{
    // load config
    Json::Value pbd_config = cJsonUtil::ParseAsValue("pbd_config", conf);
    mItersPBD = cJsonUtil::ParseAsInt("max_pbd_iters", pbd_config);
    mStiffnessPBD = cJsonUtil::ParseAsDouble("stiffness_pbd", pbd_config);
    mEnableParallelPBD =
        cJsonUtil::ParseAsBool("enable_parallel_pbd", pbd_config);
    mEnableBendingPBD =
        cJsonUtil::ParseAsBool("enable_bending_pbd", pbd_config);
    mBendingStiffnessPBD =
        cJsonUtil::ParseAsDouble("bending_stiffness_pbd", pbd_config);
    printf("[pbd] enable bending %d, bendign stiffness %.4f\n",
           mEnableBendingPBD, mBendingStiffnessPBD);

    //
    cBaseCloth::Init(conf);
    mVcur.noalias() = tVectorXd::Zero(GetNumOfFreedom());

    if (mEnableBendingPBD == true)
    {
        InitBendingMatrixPBD();
    }

    CalcNodePositionVector(mXpre);
    mXcur.noalias() = mXpre;
}

/**
 * \brief           Calculat ethe bending constant matrix "K" for each edge
 *
 * if one edge is located in the boundary, its value will be "Nan"
 */

extern int SelectAnotherVerteix(tTrianglePtr tri, int v0, int v1);

extern tVector CalculateCotangentCoeff(const tVector &x0, tVector &x1,
                                       tVector &x2, tVector &x3);
void cPBDCloth::InitBendingMatrixPBD()
{
    mBendingMatrixKArray.resize(mEdgeArrayShared.size(),
                                tVector::Ones() * std::nan(""));
    for (int i = 0; i < mEdgeArrayShared.size(); i++)
    {
        const auto &e = mEdgeArrayShared[i];
        if (e->mIsBoundary == false)
        {
            // it's an interior boundary, begin to calculate the K
            int vid[4] = {e->mId0, e->mId1,
                          SelectAnotherVerteix(mTriangleArrayShared[e->mTriangleId0],
                                               e->mId0, e->mId1),
                          SelectAnotherVerteix(mTriangleArrayShared[e->mTriangleId1],
                                               e->mId0, e->mId1)};
            // printf("[debug] bending, tri %d and tri %d, shared edge: %d,
            // total vertices: %d %d %d %d\n",
            //        e->mTriangleId0, e->mTriangleId1, i, vid[0], vid[1],
            //        vid[2], vid[3]);
            mBendingMatrixKArray[i].noalias() = CalculateCotangentCoeff(
                mVertexArrayShared[vid[0]]->mPos, mVertexArrayShared[vid[1]]->mPos,
                mVertexArrayShared[vid[2]]->mPos, mVertexArrayShared[vid[3]]->mPos);
            // std::cout << "[debug] pbd K for edge " << i << " = " <<
            // mBendingMatrixKArray[i].transpose() << std::endl;
        }
    }
    // exit(0);
}

void cPBDCloth::UpdatePos(double dt)
{
    CalcExtForce(mExtForce);
    UpdateVelAndPosUnconstrained(mExtForce);
    /*
        3. collision detect
            build constraint
    */
    // ConstraintSetupPBD();

    // 4. solve constraint
    ConstraintProcessPBD();

    // 5. post process vel
    PostProcessPBD();
}

/**
 * \brief           Update the unconstrained vel and pos
 */
void cPBDCloth::UpdateVelAndPosUnconstrained(const tVectorXd &fext)
{
    // std::cout << "fext = " << fext.transpose() << std::endl;
    mVcur += mInvMassMatrixDiag.cwiseProduct(fext) * mIdealDefaultTimestep;
    // std::cout << "mVcur = " << mVcur.transpose() << std::endl;
    mXcur += mVcur * mIdealDefaultTimestep;
    // std::cout << "mXcur = " << mXcur.transpose() << std::endl;
}

/**
 * \brief           solve the stretch (spring) constraint for PBD
 * \param final_k   the determined constraint stiffness
 */
void cPBDCloth::StretchConstraintProcessPBD(double final_k)
{

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (int e_id = 0; e_id < mEdgeArrayShared.size(); e_id++)
    {
        auto e = mEdgeArrayShared[e_id];
        int id0 = e->mId0, id1 = e->mId1;
        const tVector3d &p1 = mXcur.segment(3 * id0, 3),
                        &p2 = mXcur.segment(3 * id1, 3);
        double raw = e->mRawLength;
        // std::cout << "raw = " << raw << std::endl;
        double dist = (p1 - p2).norm();
        double w1 = mInvMassMatrixDiag[3 * e->mId0],
               w2 = mInvMassMatrixDiag[3 * e->mId1];
        double w_sum = w1 + w2;
        double coef1 = -w1 / w_sum * final_k, coef2 = w2 / w_sum * final_k;
        if (w_sum == 0)
        {
            continue;
        }
        const tVector3d &delta_p1 = coef1 * (dist - raw) * (p1 - p2) / dist,
                        &delta_p2 = coef2 * (dist - raw) * (p1 - p2) / dist;
        {
            mXcur[3 * id0 + 0] += delta_p1[0];
            mXcur[3 * id0 + 1] += delta_p1[1];
            mXcur[3 * id0 + 2] += delta_p1[2];

            mXcur[3 * id1 + 0] += delta_p2[0];
            mXcur[3 * id1 + 1] += delta_p2[1];
            mXcur[3 * id1 + 2] += delta_p2[2];
        }
    }
}

/**
 * \brief           given raw vertex vector p, solve the constraint and get the
 * new p
 */
void cPBDCloth::ConstraintProcessPBD()
{
    const int iters = mItersPBD;
    double stretch_k = 1 - std::pow((1 - mStiffnessPBD), 1.0 / iters);
    double bending_k = 1 - std::pow((1 - mBendingStiffnessPBD), 1.0 / iters);
    // std::cout << "X cur = " << mXcur.transpose() << std::endl;
    // std::cout << "mEdgeArray size = " << mEdgeArrayShared.size() << std::endl;
    for (int i = 0; i < iters; i++)
    {
        StretchConstraintProcessPBD(stretch_k);

        if (mEnableBendingPBD == true)
        {
            BendingConstraintProcessPBD(bending_k);
        }
        // std::cout << "X cur = " << mXcur.transpose() << std::endl;
        // exit(0);
    }
}

/**
 * \brief           process the bendign constraint
 * \param final_k   determined PBD stiffness
 */
void cPBDCloth::BendingConstraintProcessPBD(double final_k)
{
    tMatrixXd K = tMatrixXd::Zero(3, 12);
    // for (auto &e : mEdgeArray)
    for (int _i = 0; _i < mEdgeArrayShared.size(); _i++)
    {
        const auto &e = mEdgeArrayShared[_i];
        if (e->mIsBoundary == false)
        {
            // for interior edges, calculate the nodal displacement "delta"
            /*
                for a given interior edge (i.e. a pair of adjoint triangles),
                the bending constraint C(x) = 0 is\in R^3
                C(x) = [C_0, C_1, C_2] = [C_i], i = {0, 1, 2}
                Here we solve C_i(x) \in R individually.
            */
            int vid[4] = {e->mId0, e->mId1,
                          SelectAnotherVerteix(mTriangleArrayShared[e->mTriangleId0],
                                               e->mId0, e->mId1),
                          SelectAnotherVerteix(mTriangleArrayShared[e->mTriangleId1],
                                               e->mId0, e->mId1)};
            tVectorXd x_total = tVectorXd::Zero(12);

            const tVector &K_vec = mBendingMatrixKArray[_i];
            for (int j = 0; j < 4; j++)
            {
                K.block(0, 3 * j, 3, 3).noalias() =
                    tMatrix3d::Identity() * K_vec[j];
                x_total.segment(3 * j, 3).noalias() =
                    mXcur.segment(vid[j] * 3, 3);
            }
            // std::cout << "K = \n"
            //           << K << std::endl;
            // std::cout << "x_total = " << x_total.transpose() << std::endl;
            tVector inv_mass_vec = tVector(
                mInvMassMatrixDiag[3 * vid[0]], mInvMassMatrixDiag[3 * vid[1]],
                mInvMassMatrixDiag[3 * vid[2]], mInvMassMatrixDiag[3 * vid[3]]);
            // std::cout << "inv mass vec = " << inv_mass_vec.transpose() <<
            // std::endl;
            for (int j = 0; j < 3; j++) // j is the constraint id in [0,2]
            {
                /*
                begin to Solve C_j
               */
                // SIM_ERROR("hasn't been impled");
                // 1. calculate C_j(x)
                double C_j = K.row(j).dot(x_total);
                // printf("[debug] C_%d = %.3f\n", j, C_j);
                // 1. calculate \nabla_{x_k} C_j(x), k=[0, 1, 2, 3]
                tVector3d DCj_Dxk[4] = {
                    K.row(j).segment(3 * 0, 3), K.row(j).segment(3 * 1, 3),
                    K.row(j).segment(3 * 2, 3), K.row(j).segment(3 * 3, 3)};

                // 2. calculate s_j
                double s_j = 0;
                for (int k = 0; k < 4; k++)
                {
                    s_j += (inv_mass_vec[k] * DCj_Dxk[k].norm());
                }
                s_j = C_j / s_j;
                // printf("[debug] s%d = %.3f\n", j, s_j);
                // 3. calculate \Delta x_i^j
                for (int i = 0; i < 4; i++)
                {
                    int id = vid[i];
                    tVector3d delta =
                        -s_j * inv_mass_vec[i] * DCj_Dxk[i] * final_k;
                    // std::cout << "DC" << j << "Dx" << i << " = " <<
                    // delta.transpose() << std::endl;
                    mXcur.segment(id * 3, 3) += delta;
                }
            }
        }
    }
}
/**
 * \brief
 */
void cPBDCloth::PostProcessPBD()
{
    // SIM_WARN("PostProcessPBD hasn't been impled");
    mVcur = (mXcur - mXpre) / mIdealDefaultTimestep;
    SetPos(mXcur);
    mXpre = mXcur;
}

void cPBDCloth::CalcExtForce(tVectorXd &ext_force) const
{
    cBaseCloth::CalcExtForce(ext_force);
    ext_force += -mDamping * this->mVcur;
    // std::cout << "damping = " << (-mDamping * this->mVcur).transpose() <<
    // std::endl;
}
