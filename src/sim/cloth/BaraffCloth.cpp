#include "BaraffCloth.h"
#include "BaraffMaterial.h"
#include "geometries/Primitives.h"
#include "imgui.h"
#include "sim/Perturb.h"
#include "sim/cloth/DihedralBending.h"
#include "sim/cloth/QBendingMaterial.h"
#include "utils/JsonUtil.h"
#include "utils/TimeUtil.hpp"
#include <iostream>
float gTargetV = 0;
cBaraffCloth::cBaraffCloth(int id_) : cBaseCloth(eClothType::FEM_CLOTH, id_)
{
    // mF.clear();
    // mJ.resize(0);
    // mPK1.clear();
    // mdFdx.clear();
}

cBaraffCloth::~cBaraffCloth() {}
cBaseBendingMaterialPtr mDiMaterial = std::make_shared<cDihedralMaterial>();
void cBaraffCloth::Init(const Json::Value &conf)
{
    mRayleightA = cJsonUtil::ParseAsDouble("rayleigh_damping_a", conf);
    mRayleightB = cJsonUtil::ParseAsDouble("rayleigh_damping_b", conf);
    bool mEnableGravity = cJsonUtil::ParseAsDouble("enable_gravity", conf);
    mMaxCGIters = cJsonUtil::ParseAsInt("max_iters", conf);
    cBaseCloth::Init(conf);
    mXcur.noalias() = mClothInitPos;
    mXpre.noalias() = mClothInitPos;
    // InitBuffer();
    InitMaterialCoords();

    // init material
    // mMaterial = std::make_shared<cBaraffMaterialUnstable>();
    // mBendingMaterial = std::make_shared<cQBendingMaterial>();
    mBendingMaterial = std::make_shared<cDihedralMaterial>();
    mBendingK = tVector3f::Ones() *
                cJsonUtil::ReadVectorJson(
                    cJsonUtil::ParseAsValue("cloth_bending_stiffness", conf))
                    .cast<float>()[0];
    mBendingMaterial->Init(GetVertexArray(), GetEdgeArray(), GetTriangleArray(),
                           mBendingK.cast<double>());
    mDiMaterial->Init(GetVertexArray(), GetEdgeArray(), GetTriangleArray(),
                      mBendingK.cast<double>());
    mStretchK = cJsonUtil::ReadVectorJson(
                    cJsonUtil::ParseAsValue("cloth_stretch_stiffness", conf))
                    .segment(0, 3)
                    .cast<float>();
    // mMaterial->SetStretchK(mStretchK[0], mStretchK[1]);
    // mMaterial->SetSheaingK(mStretchK[2]);
    mMaterial = std::make_shared<cBaraffMaterial>();
    mMaterial->Init(shared_from_this(), mStretchK.cast<double>());
    mDragVertexIdx = -1;
    mYDragForce = 0;
    // tVectorXd random_move = tVectorXd::Random(GetNumOfFreedom());
    // mXcur += random_move;
    // SetPos(mXcur);
    // mMaterial->CheckForce();
    // mMaterial->CheckStiffnessMatrix();
    // mXcur -= random_move;
    // SetPos(mXcur);
    // std::cout << "bending energy = " << mBendingMaterial->CalcEnergy(mXcur);
    // mXcur[0] += 0.1;
    // std::cout << "bending energy = " << mBendingMaterial->CalcEnergy(mXcur);
    // mXcur[2] += 0.1;
    // std::cout << "bending energy = " << mBendingMaterial->CalcEnergy(mXcur);
    // exit(1);
    // allocate forces
    {
        this->ClearForce();
        this->mGravityForce.noalias() = tVectorXd::Zero(GetNumOfFreedom());
        for (int i = 0; i < GetNumOfVertices(); i++)
        {
            mGravityForce[3 * i + 1] =
                -mVertexArray[i]->mMass * 9.8 * mEnableGravity;
        }
        std::cout << "mEnableGravity = " << mEnableGravity << std::endl;
        // std::cout << mGravityForce.transpose();
        // exit(1);
        // mGravityForce[0] = 1;
    }
}

/**
 * \brief       Update the nodal position
 */
typedef std::pair<std::string, float> tTimeRecms;
static std::vector<tTimeRecms> gProfRec;
void cBaraffCloth::UpdatePos(double dt)
{
    gProfRec.clear();
    dt = mIdealDefaultTimestep;
    // CheckForce();
    // 1. calculate stiffness matrix
    cTimeUtil::Begin("update material");
    {
        // mMaterial->Update(true, true, true);
        // cTimeUtil::Begin("update ")
        mMaterial->Update();
        mBendingMaterial->Update();
        mDiMaterial->Update();
        // cTimeUtil::End("update")
        // exit(1);
    }

    gProfRec.push_back(
        tTimeRecms("update material", cTimeUtil::End("update material", true)));
    // std::cout << "-------\n";
    // UpdateCollisionForce(mCollisionForce);
    cTimeUtil::Begin("calc_K");
    CalcStiffnessMatrix(mXcur, mStiffnessMatrix);

    {
        // mStiffnessMatrix = mMaterial->CalcTotalStiffnessMatrix() +
        // this->mBendingMaterial->GetStiffnessMatrix(); auto diff_norm = (new_K
        // - mStiffnessMatrix).cwiseAbs().norm(); std::cout << "new_K = " <<
        // new_K.toDense() << std::endl; std::cout << "old_K = " <<
        // mStiffnessMatrix.toDense() << std::endl; std::cout << "[K] diff_norm
        // = " << diff_norm << std::endl;
    }

    gProfRec.push_back(tTimeRecms("calc_K", cTimeUtil::End("calc_K", true)));

    cTimeUtil::Begin("calc_fint");

    CalcIntForce(mXcur, mIntForce);

    ApplyYDragForce(mUserForce);
    // std::cout << "[fint] diff norm = " << (new_fint - mIntForce).norm() <<
    // std::endl;

    gProfRec.push_back(
        tTimeRecms("calc_fint", cTimeUtil::End("calc_fint", true)));

    // for (int i = 0; i < mVertexArray.size(); i++)
    // {
    //     std::cout << "fint" << i << " = " << mIntForce.segment(3 * i,
    //     3).transpose() << std::endl;
    // std::cout << "x" << i << " = " << mXcur.segment(3 * i, 3).transpose() <<
    // std::endl;
    // }
    // std::cout << "user force = " << mUserForce.transpose() << std::endl;
    cTimeUtil::Begin("solve");
    SolveForDx(dt);
    // mDx_buf.setZero();
    Repulsion(dt, this->mDx_buf);
    gProfRec.push_back(tTimeRecms("solve", cTimeUtil::End("solve", true)));

    mXpre.noalias() = mXcur;
    mXcur.noalias() = mXpre + mDx_buf;

    SetPos(mXcur);
    this->mUserForce.setZero();
    mDragVertexIdx = -1;
    this->mPointTriangleCollisionInfo.clear();
    this->mEdgeEdgeCollisionInfo.clear();
    // check
    // this->mMaterial->CheckForce();
    // this->mMaterial->CheckStiffnessMatrix();
    // if (frame > 5)
    //     exit(1);
}
static int frame = 0;
void cBaraffCloth::ApplyUserPerturbForceOnce(tPerturb *pert)
{
    if (pert)
    {
        int v0 = mTriangleArray[pert->mAffectedTriId]->mId0;
        int v1 = mTriangleArray[pert->mAffectedTriId]->mId1;
        int v2 = mTriangleArray[pert->mAffectedTriId]->mId2;

        mXcur.segment(3 * v0, 3) = pert->GetGoalPos().segment(0, 3);
        mXpre.segment(3 * v0, 3) = pert->GetGoalPos().segment(0, 3);
        this->SetPos(mXcur);
        mDragVertexIdx = v0;
    }
}

/**
 * \brief           Init matrix coords
 */
void cBaraffCloth::InitMaterialCoords()
{
    // calculate material coordinates
    mVertexMateralCoords.noalias() = tMatrixXd::Zero(GetNumOfVertices(), 2);
    for (int i = 0; i < mVertexArray.size(); i++)
    {
        mVertexMateralCoords.row(i).noalias() =
            mVertexArray[i]->muv.cast<double>();
    }
}
/**
 * \brief       the freedom of a triangle, 3 nodes in cartesian space, the dof =
 * 3*3 = 9
 */
int cBaraffCloth::GetSingleElementFreedom() const { return 9; }

/**
 * \brief           calculate internal force in
 */
void cBaraffCloth::CalcIntForce(const tVectorXd &xcur,
                                tVectorXd &int_force) const
{
    tVectorXd bending_force = mBendingMaterial->CalcForce(mXcur);
    int_force = mMaterial->CalcTotalForce() + bending_force;
}

void cBaraffCloth::InitMass(const Json::Value &conf)
{
    mClothDensity = cJsonUtil::ParseAsDouble("cloth_density", conf);
    int dof = 3 * GetNumOfVertices();
    mMassMatrixDiag.noalias() = tVectorXd::Zero(dof);
    for (auto &t : mTriangleArray)
    {
        // 1. total area
        auto v0 = mVertexArray[t->mId0];
        auto v1 = mVertexArray[t->mId1];
        auto v2 = mVertexArray[t->mId2];

        double triangle_area =
            cMathUtil::CalcTriangleArea(v0->mPos, v1->mPos, v2->mPos);
        mMassMatrixDiag.segment(3 * t->mId0, 3) +=
            triangle_area / 3 * mClothDensity * tVector3d::Ones();
        mMassMatrixDiag.segment(3 * t->mId1, 3) +=
            triangle_area / 3 * mClothDensity * tVector3d::Ones();
        mMassMatrixDiag.segment(3 * t->mId2, 3) +=
            triangle_area / 3 * mClothDensity * tVector3d::Ones();
    }
    for (int i = 0; i < GetNumOfVertices(); i++)
    {
        mVertexArray[i]->mMass = mMassMatrixDiag[3 * i];
        // std::cout << "v" << i << " mass = " << mVertexArray[i]->mMass
        // << std::endl;
    }
}

void cBaraffCloth::CalcStiffnessMatrix(const tVectorXd &xcur,
                                       tSparseMatd &K_global) const
{
    tSparseMatd bending_H = this->mBendingMaterial->GetStiffnessMatrix();
    // tSparseMatd di_H = mDiMaterial->GetStiffnessMatrix();
    K_global = mMaterial->CalcTotalStiffnessMatrix();
    K_global += bending_H;
    // std::cout << "bending_H = \n" << bending_H.toDense() << std::endl;
    // std::cout << "qbendning H = \n" << qbending_h << std::endl;
    // std::cout << "di H = \n" << di_H << std::endl;
    // std::cout << "K = \n" << K_global << std::endl;
    // exit(1);
    // +
    //            this->mBendingMaterial->GetStiffnessMatrix();
}

/**
 * \brief           solve for next pos by baraff scheme
 *  let A = (I + dt * damping_a) * I + dt * (damping_b - dt) * MInv * K
 *  let b = dt * Minv (fext + fint + (dt - damping_b) * K * vt) - damping_a * vt
 *  we will have A delta_v = b
 */
#include <Eigen/SparseCholesky>
#include <fstream>
void cBaraffCloth::SolveForDx(double dt)
{
    int dof = GetNumOfFreedom();
    tSparseMatd M(dof, dof);
    M.reserve(dof);
    // M.diagonal() = mMassMatDiag;
    for (size_t i = 0; i < dof; i++)
    {
        M.coeffRef(i, i) = mMassMatrixDiag[i];
        // I.coeff(i, i) = 1;
    }

    tSparseMatd W =
        (1 + dt * mRayleightA) * M + dt * (mRayleightB - dt) * mStiffnessMatrix;
    // add effect for constraint points
    // for (auto &i : this->mConstraint_StaticPointIds)
    // {
    //     for (int j = 0; j < 3; j++)
    //         W.coeffRef(3 * i + j, 3 * i + j) += 1e12;
    // }
    tVectorXd b = dt * dt * (mGravityForce + mUserForce + mIntForce) +
                  dt * M * (mXcur - mXpre) / dt;
    mDx_buf.noalias() = mXcur - mXpre;
    float threshold = 1e-12, residual = 0;
    int iters = 0;
    // std::cout << "cpu A = \n" << W.toDense() << std::endl;
    // std::cout << "cpu b = " << b.transpose() << std::endl;
    // exit(1);
    std::vector<int> fix_array = GetConstraintStaticVertices();
    if (mDragVertexIdx != -1)
        fix_array.push_back(mDragVertexIdx);
    W += PrepareCollisionHessian();

    // std::cout << "W = \n" << W.toDense() << std::endl;
    // std::cout << "b = \n" << b.transpose() << std::endl;
    Solve(W, b, mDx_buf, threshold, iters, residual, fix_array);
}

void cBaraffCloth::UpdateCollisionForce(tVectorXd &col_force)
{
    int num_of_v = this->mVertexArray.size();
    double ground_height = 1e-3;
    double k = 1e2;
    float KinectFrictionCoef = 0.5;
    // mExtForce.fill(5);
    // mExtForce[3 * 1 + 1] = 10;
    col_force.setZero();
    for (int i = 0; i < mVertexArray.size(); i++)
    {
        double dist = mVertexArray[i]->mPos[1] - ground_height;

        if (dist < 0)
        {
            mVertexArray[i]->mPos[1] = 0;
            float normal_force_amp = -dist * k;

            tVector3d friction_dir =
                -1 * (mXcur.segment(3 * i, 3) - mXpre.segment(3 * i, 3));
            friction_dir[1] = 0;
            friction_dir.normalize();
            float friction_force_amp =
                KinectFrictionCoef * std::fabs(normal_force_amp);
            tVector3d force = tVector3d::Zero();
            force[1] = normal_force_amp;
            force[0] = friction_dir[0] * friction_force_amp;
            force[2] = friction_dir[2] * friction_force_amp;
            // std::cout << "[col] force " << force.transpose() << " on v" << i
            // << std::endl;
            col_force.segment(3 * i, 3) = force;
        }
    }
}

#include "InnerForceCalculator.h"
void cBaraffCloth::UpdateImGui()
{
    ImGui::SliderFloat("dampingA", &mRayleightA, 0, 1);
    ImGui::SliderFloat("dampingB", &mRayleightB, -1e-2, 0);

    tVector3f stretch = mStretchK.cast<float>();
    tVector3f bending = mBendingK.cast<float>();
    // float stretch_warp = mStretchK[0],
    //       stretch_weft = mStretchK[1],
    //       stretch_bias = mStretchK[2];

    ImGui::SliderFloat("Kwarp", &stretch[0], 1e-2, 1e4);
    ImGui::SliderFloat("Kweft", &stretch[1], 1e-2, 1e4);
    ImGui::SliderFloat("Kbias", &stretch[2], 0, 1e4);

    if ((stretch - mStretchK).norm() > 1e-6)
    {
        mStretchK = stretch;
        mMaterial->SetK(mStretchK.cast<double>());
        std::cout << "change bending K = " << mStretchK.transpose()
                  << std::endl;
    }

    ImGui::SliderFloat("Bending warp", &bending[0], 1e-5, 1e-2, "%.5f");
    ImGui::SliderFloat("Bending Kweft", &bending[1], 1e-5, 1e-2, "%.5f");
    ImGui::SliderFloat("Bending Kbias", &bending[2], 1e-5, 1e-2, "%.5f");
    if ((bending - mBendingK).norm() > 1e-6)
    {
        mBendingK = bending;
        mBendingMaterial->Init(GetVertexArray(), GetEdgeArray(),
                               GetTriangleArray(), mBendingK.cast<double>());
        std::cout << "change bending K = " << mBendingK.transpose()
                  << std::endl;
    }

    float cloth_v = GetClothShape()[0];
    float old_tar = gTargetV;
    ImGui::SliderFloat("inner force target v", &gTargetV, -cloth_v / 2,
                       cloth_v / 2, "%.3f");
    if (ImGui::Button("calc inner force"))
    {
        // begin to calculate
        auto res = cInnerForceCalculator::FindInvolvedVertices(
            mVertexArray, mEdgeArray, mTriangleArray, this->mClothUVmin,
            this->mClothUVmax, gTargetV);
        std::cout << "involved num of v = " << res.size() << std::endl;
    }

    ImGui::SliderFloat("drag force Y", &mYDragForce, 0, 1, "%.2f");

    tVector aabb_min, aabb_max;
    CalcAABB(aabb_min, aabb_max);
    tVector aabb_size = aabb_max - aabb_min;

    ImGui::Text("aabb size %.3f %.3f %.3f", aabb_size[0], aabb_size[1],
                aabb_size[2]);

    for (auto &x : gProfRec)
    {
        ImGui::Text("%s %d ms", x.first.c_str(), int(x.second));
    }
}

void cBaraffCloth::CheckForce()
{
    // old energy
    tVectorXd cur_x = mXcur;
    double E_old = CalcEnergy(cur_x);

    // ana force
    tVectorXd ana_fint;

    CalcIntForce(cur_x, ana_fint);
    int dof = GetNumOfFreedom();
    tVectorXd num_fint = tVectorXd::Zero(dof);
    // new energy
    double eps = 1e-5;
    for (int i = 0; i < dof; i++)
    {
        cur_x[i] += eps;
        double E_new = CalcEnergy(cur_x);
        num_fint[i] = -(E_new - E_old) / eps;
        cur_x[i] -= eps;
    }
    tVectorXd f_diff = num_fint - ana_fint;
    // std::cout << "[check force] num_fint = " << num_fint.transpose() <<
    // std::endl; std::cout << "[check force] ana_fint = " <<
    // ana_fint.transpose() << std::endl; std::cout << "[check force] f diff = "
    // << f_diff.transpose() << std::endl;
    std::cout << "[check force] f diff norm = " << f_diff.norm() << std::endl;
}

void GetTrianglePosMatAndUV(tVector3d v0, tVector3d v1, tVector3d v2,
                            tVector2f u0, tVector2f u1, tVector2f u2,
                            tMatrix3d &pos_mat, tMatrix32d &uv_mat)
{
    pos_mat.col(0) = v0;
    pos_mat.col(1) = v1;
    pos_mat.col(2) = v2;
    uv_mat.row(0) = u0.cast<double>();
    uv_mat.row(1) = u1.cast<double>();
    uv_mat.row(2) = u2.cast<double>();
}

double cBaraffCloth::CalcEnergy(const tVectorXd &xcur)
{
    double e_total = 0;
    for (auto &t : this->mTriangleArray)
    {
        tMatrix3d pos_mat;
        tMatrix32d uv_mat;
        GetTrianglePosMatAndUV(
            xcur.segment(3 * t->mId0, 3), xcur.segment(3 * t->mId1, 3),
            xcur.segment(3 * t->mId2, 3), mVertexArray[t->mId0]->muv,
            mVertexArray[t->mId1]->muv, mVertexArray[t->mId2]->muv, pos_mat,
            uv_mat);
        e_total += mMaterial->CalcTotalEnergy();

        // pos_mat.col(0) =->mPos.segment(0, 3);
        // pos_mat.col(1) =->mPos.segment(0, 3);
        // pos_mat.col(2) =->mPos.segment(0, 3);
    }
    return e_total;
}

void Filter(const std::vector<int> &constraint_pts, tVectorXd &vec)
{
    for (auto &i : constraint_pts)
    {
        vec.segment(3 * i, 3).setZero();
    }
}

tVectorXd ApplyMatmul(const tSparseMatd &A, const tVectorXd &b,
                      const std::vector<int> &fix_vertex_array)
{
    auto should_remove =
        [](const std::vector<int> &fix_vertex_array, int cur_idx)
    {
        for (auto &drag_pt_idx : fix_vertex_array)
        {
            bool is_return = (cur_idx == 3 * drag_pt_idx + 0) ||
                             (cur_idx == 3 * drag_pt_idx + 1) ||
                             (cur_idx == 3 * drag_pt_idx + 2);
            if (is_return)
                return true;
        }
        return false;
    };
    tVectorXd res = tVectorXd::Zero(b.size());
    OMP_PARALLEL_FOR
    for (int k = 0; k < A.outerSize(); ++k)
    {
        if (should_remove(fix_vertex_array, k))
            continue;
        for (tSparseMatd::InnerIterator it(A, k); it; ++it)
        {
            // std::cout << it.row() << "\t";
            // std::cout << it.col() << "\t";
            // std::cout << it.value() << std::endl;
            if (should_remove(fix_vertex_array, it.col()))
            {
                // std::cout << "drag pt " << drag_pt << " ignore " << it.col()
                //           << std::endl;
                continue;
            }
            res[k] += b[it.col()] * it.value();
        }
    }
    return res;
}
void cBaraffCloth::Solve(const tSparseMatd &A, const tVectorXd &b, tVectorXd &x,
                         float threshold, int &iters, float &residual,
                         const std::vector<int> &fix_vertex_array /* = -1*/)
{
    // std::cout << "Eigen::nbThreads() = " << Eigen::nbThreads() << std::endl;
    int recalc_gap = 20;
    int max_iters = mMaxCGIters;
    // tVectorXd r = b - A * x;
    tVectorXd r = b - ApplyMatmul(A, x, fix_vertex_array);
    for (auto fix_vertex : fix_vertex_array)
    {
        r.segment(3 * fix_vertex, 3).setZero();
        // std::cout << "r[init] = " << r.segment(3 * fix_vertex, 3).transpose()
        //           << std::endl;
        x.segment(3 * fix_vertex, 3).setZero();
    }
    tVectorXd PInv = A.diagonal().cwiseInverse();
    tVectorXd d = PInv.cwiseProduct(r);
    iters = 1;
    tVectorXd Adi;
    tVectorXd rnext;
    tVectorXd PInvr, PInvrnext;
    while (iters < max_iters)
    {
        // printf("-----iter %d-----\n", iters);
        Adi.noalias() = ApplyMatmul(A, d, fix_vertex_array);
        // if (Adi.hasNaN())
        // {
        //     std::cout << "Adi has Nan\n";
        //     exit(1);
        // }
        // if (fix_vertex != -1)
        // {
        //     std::cout << "Adi[fixed] = "
        //               << Adi.segment(3 * fix_vertex, 3).transpose()
        //               << std::endl;
        // }
        PInvr.noalias() = PInv.cwiseProduct(r);
        // Filter(mConstraint_StaticPointIds, q);
        double deno = d.dot(Adi);
        if (deno == 0)
        {
            std::cout << "[warn] deno = 0, break\n";

            break;
        }
        double alpha = r.dot(PInvr) / deno;
        // if (std::isnan(alpha) || std::isinf(alpha))
        // {
        //     std::cout << "alpha is nan or inf " << alpha << std::endl;
        //     std::cout << "deno = " << deno << std::endl;
        //     exit(1);
        // }
        x += alpha * d;

        // if (iters % recalc_gap == 0)
        // {
        //     r = b - A * x;
        //     Filter(this->mConstraint_StaticPointIds, r);
        // }
        // else
        rnext.noalias() = r - alpha * Adi;
        // if (rnext.hasNaN())
        // {

        //     std::cout << "rnext has Nan\n";
        //     std::cout << "r = " << r.transpose() << std::endl;
        //     std::cout << "Adi = " << Adi.transpose() << std::endl;
        //     std::cout << "alpha = " << alpha << std::endl;
        //     exit(1);
        // }
        PInvrnext.noalias() = PInv.cwiseProduct(rnext);
        // if (PInvrnext.hasNaN())
        // {
        //     std::cout << "PInvrnext has Nan\n";
        //     exit(1);
        // }
        double beta = rnext.dot(PInvrnext) / (r.dot(PInvr));
        if (beta == 0)
        {
            std::cout << "[warn] beta = 0, break\n";
            break;
            // break;
        }
        d = PInvrnext + beta * d;
        // if (d.hasNaN())
        // {
        //     std::cout << "d has Nan\n";
        //     exit(1);
        // }
        r.noalias() = rnext;
        // if (r.hasNaN())
        // {
        //     std::cout << "r has Nan\n";
        //     exit(1);
        // }
        iters += 1;
    }
    std::cout << "iters " << iters << " residual = " << r.norm() << std::endl;
    if (x.hasNaN())
    {
        std::cout << "x has Nan\n";
        exit(1);
    }
}
#include "sim/collision/CollisionInfo.h"
/**
 * \brief       Apply repulsion force by Bridson 2002
 *
 * inelastic collision: I = - normal * |Ic|
 * non-penetration: I = - normal * |Ir|
 *
 * let v_thre = 0.5
 * |Ic| = 0, vn > v_thre
 *        m * |vn| / 2, vn < v_thre
 *
 * |Ir| = dt * k * d
 */
#include <set>
void cBaraffCloth::Repulsion(double dt, tVectorXd &dx) const
{
    int iters = 1;
    float scale = 0.1;
    for (int i = 0; i < iters; i++)
    {
        double v_thre = 0.01;
        double K = mStretchK.mean() * 1;
        // double K = mMassMatrixDiag.sum() / 3 / (GetNumOfVertices() * dt *
        // dt);
        tVectorXd ddx = tVectorXd::Zero(dx.size());
        tVectorXd ddx_times = tVectorXd::Zero(dx.size() / 3);
        std::set<int> changed_v = {};
        for (int _idx = 0; _idx < mEdgeEdgeCollisionInfo.size(); _idx++)
        {
            // printf("----------ee %d---------\n", _idx);
            auto info = mEdgeEdgeCollisionInfo[_idx];
            /*
            1. calculate |Ic|
                1.1 calculate vn
                1.2 calculate Ic
            */
            tVector3d outer_normal = info->mOuterNormal;
            tVector3d vrel = tVector3d::Zero();
            double mass = 0;
            int edge_id = -1;
            if (info->obj0_is_cloth_obj1_is_rb == true)
            {
                // obj0 is cloth, get vertex vel
                edge_id = info->mEdgeId0;
                int v0 = mEdgeArray[edge_id]->mId0;
                int v1 = mEdgeArray[edge_id]->mId1;
                vrel = ((1 - info->mBary0) * dx.segment(3 * v0, 3) +
                        info->mBary0 * dx.segment(3 * v1, 3)) /
                       dt;
                mass = mMassMatrixDiag[3 * v0] + mMassMatrixDiag[3 * v1];
            }
            else
            {
                // obj1 is cloth, get vertex vel
                edge_id = info->mEdgeId1;
                int v0 = mEdgeArray[edge_id]->mId0;
                int v1 = mEdgeArray[edge_id]->mId1;
                vrel = ((1 - info->mBary1) * dx.segment(3 * v0, 3) +
                        info->mBary1 * dx.segment(3 * v1, 3)) /
                       dt;
                mass = mMassMatrixDiag[3 * v0] + mMassMatrixDiag[3 * v1];
            }
            // std::cout << "mass = " << mass << std::endl;
            // std::cout << "vrel = " << vrel.transpose() << std::endl;
            // std::cout << "normal = " << outer_normal.transpose() <<
            // std::endl;
            double vn = vrel.dot(outer_normal);
            double Ic = vn > v_thre ? 0 : mass * std::fabs(vn) / 2;
            // std::cout << "Ic = " << Ic << std::endl;
            /*
            2. calculate |Ir|
            */
            double d = (3e-3 - info->mDepth);

            double Ir = 0;

            // if vn is too big (vn < - 0.1 * d /dt)
            float inhibition_coef = 0.1;
            // if the normal velocity can resolve 10% overlap, we don't add more
            // impulse. otherwise it must be bounce up
            if (vn > inhibition_coef * d / dt)
            {
                // printf("vn %.3f is bigger (< %.3f), Ir = 0\n", vn,
                //        -inhibition_coef * d / dt);
                Ir = 0;
            }
            else
            {
                // printf("------vn %.3f is small or negative, we begin to add "
                //        "spring "
                //        "force------\n",
                //        vn);
                // if vn is not so big, we apply Ir
                // Ir = SIM_MIN(dt * K * d, mass * (inhibition_coef * d / dt -
                // vn));
                double Ir0 = dt * K * d;
                double Ir1 = mass * (inhibition_coef * d / dt -
                                     vn); // if we apply this impulse, our next
                                          // vel = inhibition_coef * d / dt
                // printf("spring impulse %.5f, inhibition impulse %.5f, we take
                // "
                //        "%s\n",
                //        Ir0, Ir1, Ir0 < Ir1 ? "spring" : "inhibition");
                Ir = SIM_MIN(Ir0, Ir1);
                // if(Ir0 < Ir1)
                // {

                // }
                // std::cout << "[Ir] ir = " << Ir << " Ic = " << Ic <<
                // std::endl; Ir = dt * K * d;
            }
            // std::cout << "Ir = " << Ir << std::endl;
            // double Ir = 0;
            // std::cout << "Ir = " << Ir << std::endl;
            double I = Ic + Ir;
            // std::cout << "I = " << I << std::endl;
            double I_tilde = 2 * I /
                             (1 + info->mBary0 * info->mBary0 +
                              (1 - info->mBary0) * (1 - info->mBary0) +
                              info->mBary1 * info->mBary1 +
                              (1 - info->mBary1) * (1 - info->mBary1));

            // Delta v0 = (1-a) * I_tilde / m * normal
            // Delta v1 = a * I_tilde / m * normal
            // Delta v2 = (1-b) * I_tilde / m * normal
            // Delta v3 = b * I_tilde / m * normal
            // \Delta (dx/dt) = \Delta V
            // dx = dt * \Delta v
            if (info->obj0_is_cloth_obj1_is_rb == true)
            {
                // obj0 is cloth
                int v0 = mEdgeArray[info->mEdgeId0]->mId0;
                int v1 = mEdgeArray[info->mEdgeId0]->mId1;
                dx.segment(3 * v0, 3) += scale * dt * (1 - info->mBary0) *
                                         I_tilde / mass * outer_normal;
                ddx_times[v0] += 1;
                dx.segment(3 * v1, 3) +=
                    scale * dt * info->mBary0 * I_tilde / mass * outer_normal;
                // std::cout << "ddx for v" << v0 << " = "
                //           << dt * (1 - info->mBary0) * I_tilde / mass *
                //                  outer_normal.transpose()
                //           << std::endl;
                // std::cout << "ddx for v" << v1 << " = "
                //           << dt * info->mBary0 * I_tilde / mass *
                //                  outer_normal.transpose()
                //           << std::endl;
                ddx_times[v1] += 1;
                changed_v.insert(v0);
                changed_v.insert(v1);
            }
            else
            {
                // obj1 is cloth
                int v0 = mEdgeArray[info->mEdgeId1]->mId0;
                int v1 = mEdgeArray[info->mEdgeId1]->mId1;
                dx.segment(3 * v0, 3) += scale * dt * (1 - info->mBary1) *
                                         I_tilde / mass * outer_normal;
                ddx_times[v0] += 1;
                dx.segment(3 * v1, 3) +=
                    scale * dt * info->mBary1 * I_tilde / mass * outer_normal;

                // std::cout << "ddx for v" << v0 << " = "
                //           << dt * (1 - info->mBary1) * I_tilde / mass *
                //                  outer_normal.transpose()
                //           << std::endl;
                // std::cout << "ddx for v" << v1 << " = "
                //           << dt * info->mBary1 * I_tilde / mass *
                //                  outer_normal.transpose()
                //           << std::endl;

                ddx_times[v1] += 1;
                changed_v.insert(v0);
                changed_v.insert(v1);
            }
        }
        for (int _idx = 0; _idx < mPointTriangleCollisionInfo.size(); _idx++)
        {
            /*
            1. calculate |Ic|
                1.1 calculate vn
                1.2 calculate Ic
            */
            // printf("----------pt %d---------\n", _idx);
            auto info = mPointTriangleCollisionInfo[_idx];
            tVector3d vrel =
                tVector3d::Zero(); // vrel = vcloth - vbody = vcloth
            double mass = 0;
            if (info->obj0_is_cloth_obj1_is_rb == true)
            {
                // obj0 is cloth, get vertex vel
                vrel = dx.segment(3 * info->mVertexId0, 3) / dt;
                mass = this->mMassMatrixDiag[3 * info->mVertexId0];
            }
            else
            {
                auto tri = mTriangleArray[info->mTriangleId1];
                // obj1 is cloth get triangle vel and weight
                vrel = (dx.segment(3 * tri->mId0, 3) * info->mBary[0] +
                        dx.segment(3 * tri->mId1, 3) * info->mBary[1] +
                        dx.segment(3 * tri->mId2, 3) * info->mBary[2]) /
                       dt;
                mass = mMassMatrixDiag[3 * tri->mId0] +
                       mMassMatrixDiag[3 * tri->mId1] +
                       mMassMatrixDiag[3 * tri->mId2];
            }
            // std::cout << "vrel = " << vrel.transpose() << std::endl;
            // std::cout << "p(mv) = " << mass * vrel.transpose() << std::endl;
            double vn = vrel.dot(info->mOuterNormal);
            // double Ic = vn > v_thre ? 0 : mass * std::fabs(vn) / 2;
            double Ic = vn > v_thre ? 0 : mass * std::fabs(vn) / 2;

            // std::cout << "normal = " << info->mOuterNormal.transpose() <<
            // std::endl; std::cout << "mass = " << mass << std::endl; std::cout
            // << "vel = " << vrel.transpose() << std::endl; std::cout << "Ic =
            // " << Ic
            // << std::endl;
            /*
            2. calculate |Ir|
            */
            double d = (3e-3 - info->mDepth);
            SIM_ASSERT(d > 0);
            double Ir = 0;
            // if the velocity will bounce, we apply no
            if (vn > 0.1 * d / dt)
            {
                // printf("vn %.3f is bigger (< %.3f), Ir = 0\n", vn,
                //        -0.1 * d / dt);
                Ir = 0;
            }
            else
            {
                // Ir = SIM_MIN(dt * K * d, mass * (0.1 * d / dt - vn));
                // printf("------vn %.3f is small or negative, we begin to add "
                //        "spring "
                //        "force------\n",
                //        vn);
                // if vn is not so big, we apply Ir
                // Ir = SIM_MIN(dt * K * d, mass * (inhibition_coef * d / dt -
                // vn));
                double Ir0 = dt * K * d;
                double Ir1 = mass * (0.1 * d / dt -
                                     vn); // if we apply this impulse, our next
                                          // vel = inhibition_coef * d / dt
                // printf("spring impulse %.5f, inhibition impulse %.5f, we take
                // "
                //        "%s\n",
                //        Ir0, Ir1, Ir0 < Ir1 ? "spring" : "inhibition");
                Ir = SIM_MIN(Ir0, Ir1);
            }
            SIM_ASSERT(Ir >= 0);
            SIM_ASSERT(Ic >= 0);
            // std::cout << "Ir = " << Ir << std::endl;
            double I = Ic + Ir;
            // // std::cout << "I = " << I << std::endl;
            double I_tilde = 2 * I / (1 + info->mBary.squaredNorm());
            // // std::cout << "Itilde = " << I_tilde << std::endl;
            if (info->obj0_is_cloth_obj1_is_rb == true)
            {
                // obj0 is cloth, is vertex, apply this impulse to vertex
                // dx = dt * dv = - dt * I_tilde / m * normal
                tVector3d incre = dt * I_tilde / mass * info->mOuterNormal;
                dx.segment(3 * info->mVertexId0, 3) += scale * incre;
                ddx_times[info->mVertexId0] += 1;
                changed_v.insert(info->mVertexId0);
            }
            else
            {
                // obj1 is cloth, is triangle, apply this impulse to triangle
                // dx = dt * dv = dt * wi * I / m * N
                tVector3d incre0, incre1, incre2;
                auto tri = mTriangleArray[info->mTriangleId1];
                {
                    incre0 = dt * info->mBary[0] * I_tilde / mass *
                             info->mOuterNormal;
                    incre1 = dt * info->mBary[1] * I_tilde / mass *
                             info->mOuterNormal;
                    incre2 = dt * info->mBary[2] * I_tilde / mass *
                             info->mOuterNormal;
                }

                dx.segment(3 * tri->mId0, 3) += scale * incre0;
                dx.segment(3 * tri->mId1, 3) += scale * incre1;
                dx.segment(3 * tri->mId2, 3) += scale * incre2;

                ddx_times[tri->mId0] += 1;
                ddx_times[tri->mId1] += 1;
                ddx_times[tri->mId2] += 1;
                changed_v.insert(tri->mId0);
                changed_v.insert(tri->mId1);
                changed_v.insert(tri->mId2);
                // printf("[repulsion] add dx to v%d: %.3f %.3f %.3f\n",
                // tri->mId2,
                //        incre2[0], incre2[1], incre2[2]);
            }
        }
        // printf("-------------end-------------\n");
        for (auto v : changed_v)
        {
            // printf("-----apply v%d-----\n", v);
            // int times = 2;

            // tVector3d incre = ddx.segment(3 * v, 3) / times;
            // std::cout << "old dx = " << dx.segment(3 * v, 3).transpose()
            //           << " times = " << times << " incre = " <<
            //           incre.transpose()
            //           << std::endl;
            // std::cout << "v" << v
            //           << " cur dx = " << dx.segment(3 * v, 3).transpose()
            //           << std::endl;
            // dx.segment(3 * v, 3) += incre;
            // std::cout << "after = " << dx.segment(3 * v, 3).transpose()
            //           << std::endl;
        }
    }
}

tSparseMatd cBaraffCloth::PrepareCollisionHessian()
{
    // v-t
    int num_of_dof = GetNumOfFreedom();
    tSparseMatd Hessian(num_of_dof, num_of_dof);
    std::vector<tTriplet> triplets = {};
    for (int _idx = 0; _idx < mPointTriangleCollisionInfo.size(); _idx++)
    {
        /*
        1. calculate |Ic|
            1.1 calculate vn
            1.2 calculate Ic
        */
        // printf("----------pt %d---------\n", _idx);
        auto info = mPointTriangleCollisionInfo[_idx];
        tVector3d vrel = tVector3d::Zero(); // vrel = vcloth - vbody = vcloth
        double mass = 0;
        tVectorXd dCdx = tVectorXd::Zero(12);
        std::vector<int> local_id_2_global_id(4, -1);
        if (info->obj0_is_cloth_obj1_is_rb == true)
        {
            // obj0 is cloth, get vertex vel
            // C = N^T * (info->mVertexId0 - b)
            // \nabla_x C = dCdx = N
            dCdx.segment(0, 3) = info->mOuterNormal;
            local_id_2_global_id[0] = info->mVertexId0;
        }
        else
        {
            auto tri = mTriangleArray[info->mTriangleId1];

            dCdx.segment(3, 3) = -info->mBary[0] * info->mOuterNormal;
            dCdx.segment(6, 3) = -info->mBary[1] * info->mOuterNormal;
            dCdx.segment(9, 3) = -info->mBary[2] * info->mOuterNormal;
            local_id_2_global_id[1] = tri->mId0;
            local_id_2_global_id[2] = tri->mId1;
            local_id_2_global_id[3] = tri->mId2;
        }
        tMatrixXd dCTdC = dCdx * dCdx.transpose();
        for (int i = 0; i < 4; i++)
        {
            int global_i = local_id_2_global_id[i];
            if (global_i == -1)
                continue;
            for (int j = 0; j < 4; j++)
            {
                int global_j = local_id_2_global_id[j];
                if (global_j == -1)
                    continue;

                tMatrix3d block = 1e-2 * dCTdC.block(3 * i, 3 * j, 3, 3);
                for (int a = 0; a < 3; a++)
                    for (int b = 0; b < 3; b++)
                    {
                        triplets.push_back(tTriplet(
                            3 * global_i + a, 3 * global_j + b, block(a, b)));
                    }
            }
        }
    }
    Hessian.setFromTriplets(triplets.begin(), triplets.end());
    return Hessian;
}
