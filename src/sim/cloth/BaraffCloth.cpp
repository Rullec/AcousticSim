#include <iostream>
#include "BaraffCloth.h"
#include "utils/TimeUtil.hpp"
#include "geometries/Primitives.h"
#include "BaraffMaterial.h"
#include "utils/JsonUtil.h"
#include "sim/Perturb.h"
int frame = 0;
cBaraffCloth::cBaraffCloth(int id_) : cBaseCloth(eClothType::FEM_CLOTH, id_)
{
    // mF.clear();
    // mJ.resize(0);
    // mPK1.clear();
    // mdFdx.clear();
}

cBaraffCloth::~cBaraffCloth() {}

void cBaraffCloth::Init(const Json::Value &conf)
{
    mRayleightA = cJsonUtil::ParseAsDouble("rayleigh_damping_a", conf);
    mRayleightB = cJsonUtil::ParseAsDouble("rayleigh_damping_b", conf);
    cBaseCloth::Init(conf);
    mXcur.noalias() = mClothInitPos;
    mXpre.noalias() = mClothInitPos;
    InitBuffer();

    // init material
    mMaterial = std::make_shared<cBaraffMaterial>();
    mMaterial->SetK(1e5, 1e5);
    mMaterial->CheckStiffnessMatrix();
    mMaterial->CheckForce();

    // allocate forces
    {
        this->ClearForce();
        this->mGravityForce.noalias() = tVectorXd::Zero(GetNumOfFreedom());
        for (int i = 0; i < GetNumOfVertices(); i++)
        {
            // mGravityForce[3 * i + 1] = -9.8;
        }
        // mGravityForce[0] = 1;
    }
}

/**
 * \brief       Update the nodal position
 */

void cBaraffCloth::UpdatePos(double dt)
{

    std::cout << "--------------frame " << frame++ << "-----------\n";
    dt = 5e-3;
    // 1. calculate stiffness matrix
    UpdateCollisionForce(mCollisionForce);
    CalcStiffnessMatrix(mXcur, mStiffnessMatrix);

    CalcIntForce(mXcur, mIntForce);
    for (int i = 0; i < mVertexArrayShared.size(); i++)
    {
        std::cout << "fint" << i << " = " << mIntForce.segment(3 * i, 3).transpose() << std::endl;
        // std::cout << "x" << i << " = " << mXcur.segment(3 * i, 3).transpose() << std::endl;
    }
    // std::cout << "user force = " << mUserForce.transpose() << std::endl;
    SolveForNextPos(dt);
    SetPos(mXcur);
    this->mUserForce.setZero();
    // if (frame > 5)
    //     exit(1);
}

void cBaraffCloth::ApplyUserPerturbForceOnce(tPerturb *pert)
{
    this->mUserForce.setZero();
    if (pert)
    {
        auto t = this->mTriangleArrayShared[pert->mAffectedTriId];
        tVector3d force = pert->GetPerturbForce().segment(0, 3);
        mUserForce.segment(3 * t->mId0, 3) += force * pert->mBarycentricCoords[0];
        mUserForce.segment(3 * t->mId1, 3) += force * pert->mBarycentricCoords[1];
        mUserForce.segment(3 * t->mId2, 3) += force * pert->mBarycentricCoords[2];
    }
}

/**
 * \brief       Init FEM buffer,
 */
void cBaraffCloth::InitBuffer()
{
    // int element_size = this->GetSingleElementFreedom();
    // int num_of_triangles = GetNumOfTriangles();
    // mF.resize(num_of_triangles, tMatrix32d::Zero());
    // mJ.noalias() = tVectorXd::Zero(num_of_triangles);
    // mPK1.resize(num_of_triangles, tMatrixXd::Zero(element_size, element_size));
    // mdFdx.resize(
    //     num_of_triangles,
    //     tEigenArr<tMatrixXd>(element_size,
    //                          tMatrixXd::Zero(element_size, element_size)));

    InitMaterialCoords();
}

/**
 * \brief           Init matrix coords
 */
void cBaraffCloth::InitMaterialCoords()
{
    // calculate material coordinates
    mVertexMateralCoords.noalias() = tMatrixXd::Zero(GetNumOfVertices(), 2);
    for (int i = 0; i < mVertexArrayShared.size(); i++)
    {
        mVertexMateralCoords.row(i).noalias() =
            mVertexArrayShared[i]->muv.cast<double>();
    }

    // // calculate the DInv (used in material point)
    // mDInv.resize(mTriangleArrayShared.size(), tMatrix2d::Zero());

    // tMatrix2d mat1;
    // for (int i = 0; i < this->mTriangleArrayShared.size(); i++)
    // {
    //     auto tri = mTriangleArrayShared[i];
    //     mat1.col(0).noalias() = (mVertexMateralCoords.row(tri->mId1) -
    //                              mVertexMateralCoords.row(tri->mId0))
    //                                 .transpose();
    //     mat1.col(1).noalias() = (mVertexMateralCoords.row(tri->mId2) -
    //                              mVertexMateralCoords.row(tri->mId0))
    //                                 .transpose();

    //     /*
    //         mat0 = [b-a; c-a;] \in R 3 \times 2
    //         mat1 = [B-A; C-A;] \in R 2 \times 2
    //         mat0 = F * mat1, F \in R 3 \times 2
    //     */
    //     mDInv[i] = mat1.inverse();
    // }
}
/**
 * \brief       the freedom of a triangle, 3 nodes in cartesian space, the dof =
 * 3*3 = 9
 */
int cBaraffCloth::GetSingleElementFreedom() const { return 9; }

/**
 * \brief           calculate internal force in
 */
void cBaraffCloth::CalcIntForce(const tVectorXd &xcur, tVectorXd &int_force) const
{
    auto get_pos_mat = [](const tVectorXd &xcur,
                          int id0, int id1, int id2) -> tMatrix3d
    {
        tMatrix3d val = tMatrix3d::Zero();
        val.col(0) = xcur.segment(3 * id0, 3);
        val.col(1) = xcur.segment(3 * id1, 3);
        val.col(2) = xcur.segment(3 * id2, 3);
        return val;
    };
    // auto get_uv_mat;

    int dof = 3 * GetNumOfVertices();
    int_force.noalias() = tVectorXd::Zero(dof);
    for (auto &t : this->mTriangleArrayShared)
    {
        // 1. calc ele stiffness
        // 1.1 assemble current pos,
        // 1.2 give the 3x2d texture coords
        int v_id[3] = {t->mId0,
                       t->mId1,
                       t->mId2};
        // int id0 = t->mId0,
        //     id1 = t->mId1,
        //     id2 = t->mId2;
        tMatrix32d uv_coords;
        uv_coords.row(0) = mVertexArrayShared[v_id[0]]->muv.cast<double>();
        uv_coords.row(1) = mVertexArrayShared[v_id[1]]->muv.cast<double>();
        uv_coords.row(2) = mVertexArrayShared[v_id[2]]->muv.cast<double>();

        tVectorXd force = mMaterial->CalcForce(get_pos_mat(xcur, v_id[0], v_id[1], v_id[2]), uv_coords);

        int_force.segment(3 * v_id[0], 3) += force.segment(0, 3);
        int_force.segment(3 * v_id[1], 3) += force.segment(3, 3);
        int_force.segment(3 * v_id[2], 3) += force.segment(6, 3);
    }
}

/**
 * \brief           calculate deformation gradient
 */
// void cBaraffCloth::CalculateF()
// {
// cTimeUtil::Begin("CalculateF");
// tMatrix32d mat0 = tMatrix32d::Zero();
// tMatrix2d mat1 = tMatrix2d::Zero();
// for (int i = 0; i < this->mTriangleArrayShared.size(); i++)
// {
//     auto tri = mTriangleArrayShared[i];
//     const tVector3d &a = mXcur.segment(3 * tri->mId0, 3);
//     const tVector3d &b = mXcur.segment(3 * tri->mId1, 3);
//     const tVector3d &c = mXcur.segment(3 * tri->mId2, 3);

//     mat0.col(0).noalias() = b - a;
//     mat0.col(1).noalias() = c - a;

//     /*
//         mat0 = [b-a; c-a;] \in R 3 \times 2
//         mat1 = [B-A; C-A;] \in R 2 \times 2
//         mat0 = F * mat1, F \in R 3 \times 2
//     */
//     mF[i] = mat0 * mDInv[i];
// }
// cTimeUtil::End("CalculateF");
// }

void cBaraffCloth::InitMass(const Json::Value &conf)
{
    mClothDensity = cJsonUtil::ParseAsDouble("cloth_density", conf);
    int dof = 3 * GetNumOfVertices();
    mInvMassMatrixDiag.noalias() = tVectorXd::Zero(dof);
    for (auto &t : mTriangleArrayShared)
    {
        // 1. total area
        auto v0 = mVertexArrayShared[t->mId0];
        auto v1 = mVertexArrayShared[t->mId1];
        auto v2 = mVertexArrayShared[t->mId2];

        double triangle_area = cMathUtil::CalcTriangleArea(v0->mPos,
                                                           v1->mPos,
                                                           v2->mPos);
        mInvMassMatrixDiag.segment(3 * t->mId0, 3) = triangle_area / 3 * mClothDensity * tVector3d::Ones();
        mInvMassMatrixDiag.segment(3 * t->mId1, 3) = triangle_area / 3 * mClothDensity * tVector3d::Ones();
        mInvMassMatrixDiag.segment(3 * t->mId2, 3) = triangle_area / 3 * mClothDensity * tVector3d::Ones();
    }
    mInvMassMatrixDiag = mInvMassMatrixDiag.cwiseInverse();
}

void cBaraffCloth::CalcStiffnessMatrix(const tVectorXd &xcur,
                                       tSparseMat &K_global) const
{
    auto get_pos_mat = [](const tVectorXd &xcur,
                          int id0, int id1, int id2) -> tMatrix3d
    {
        tMatrix3d val = tMatrix3d::Zero();
        val.col(0) = xcur.segment(3 * id0, 3);
        val.col(1) = xcur.segment(3 * id1, 3);
        val.col(2) = xcur.segment(3 * id2, 3);
        return val;
    };
    // auto get_uv_mat;

    tEigenArr<tTriplet>
        sub_triples;
    for (auto &t : this->mTriangleArrayShared)
    {
        // 1. calc ele stiffness
        // 1.1 assemble current pos,
        // 1.2 give the 3x2d texture coords
        int v_id[3] = {t->mId0,
                       t->mId1,
                       t->mId2};
        // int id0 = t->mId0,
        //     id1 = t->mId1,
        //     id2 = t->mId2;
        tMatrix32d uv_coords;
        uv_coords.row(0) = mVertexArrayShared[v_id[0]]->muv.cast<double>();
        uv_coords.row(1) = mVertexArrayShared[v_id[1]]->muv.cast<double>();
        uv_coords.row(2) = mVertexArrayShared[v_id[2]]->muv.cast<double>();

        auto ele_K = mMaterial->CalcStiffMatrix(get_pos_mat(xcur, v_id[0], v_id[1], v_id[2]), uv_coords);

        // 2. assemble
        for (size_t i = 0; i < 3; i++)
        {
            size_t global_vi_idx = v_id[i];
            for (size_t j = 0; j < 3; j++)
            {
                size_t global_vj_idx = v_id[j];

                const tMatrix3d &ele_K_part = ele_K.block(3 * i, 3 * j, 3, 3);
                size_t i3 = 3 * global_vi_idx, j3 = 3 * global_vj_idx;
                sub_triples.push_back(tTriplet(i3, j3, ele_K_part(0, 0)));
                sub_triples.push_back(tTriplet(i3, j3 + 1, ele_K_part(0, 1)));
                sub_triples.push_back(tTriplet(i3, j3 + 2, ele_K_part(0, 2)));

                sub_triples.push_back(tTriplet(i3 + 1, j3, ele_K_part(1, 0)));
                sub_triples.push_back(tTriplet(i3 + 1, j3 + 1, ele_K_part(1, 1)));
                sub_triples.push_back(tTriplet(i3 + 1, j3 + 2, ele_K_part(1, 2)));

                sub_triples.push_back(tTriplet(i3 + 2, j3, ele_K_part(2, 0)));
                sub_triples.push_back(tTriplet(i3 + 2, j3 + 1, ele_K_part(2, 1)));
                sub_triples.push_back(tTriplet(i3 + 2, j3 + 2, ele_K_part(2, 2)));
            }
        }
    }
    int dof = 3 * GetNumOfVertices();
    K_global.resize(dof, dof);
    K_global.setFromTriplets(sub_triples.begin(), sub_triples.end());
}

/**
 * \brief           solve for next pos by baraff scheme
 *  let A = (I + dt * damping_a) * I + dt * (damping_b - dt) * MInv * K
 *  let b = dt * Minv (fext + fint + (dt - damping_b) * K * vt) - damping_a * vt
 *  we will have A delta_v = b
*/
void cBaraffCloth::SolveForNextPos(double dt)
{
    // {
    //     int dof = GetNumOfFreedom();
    //     tVectorXd vt = (mXcur - mXpre) / dt;
    //     tSparseMat A = dt * (mRayleightB - dt) * mStiffnessMatrix;
    //     {
    //         for (size_t row = 0; row < dof; row++)
    //             A.row(row) *= mInvMassMatrixDiag[row];
    //         A.diagonal() += tVectorXd::Ones(dof) * (1 + dt * this->mRayleightA);
    //     }

    //     tVectorXd b = dt * mInvMassMatrixDiag.cwiseProduct(mGravityForce + mUserForce + mIntForce + mCollisionForce + (dt - mRayleightB) * mStiffnessMatrix * vt) - mRayleightA * vt;

    //     Eigen::ConjugateGradient<tSparseMat, Eigen::Upper> solver;
    //     tVectorXd delta_vel = solver.compute(A).solve(b);
    //     // tVectorXd delta_vel = solver.compute(A).solve(b);
    //     // tVectorXd delta_vel = A.toDense().inverse() * b;
    //     mXpre.noalias() = mXcur;
    //     mXcur = mXcur + (vt + delta_vel) * dt;
    // }
    {
        int dof = GetNumOfFreedom();
        tVectorXd vt = (mXcur - mXpre) / dt;
        tSparseMat A = -dt * dt * mStiffnessMatrix;
        {
            for (size_t row = 0; row < dof; row++)
                A.row(row) *= mInvMassMatrixDiag[row];
            A.diagonal() += tVectorXd::Ones(dof);
        }

        tVectorXd b = dt * mInvMassMatrixDiag.cwiseProduct(
            mGravityForce + mUserForce + mIntForce + mCollisionForce + dt * mStiffnessMatrix * vt
            );

        // Eigen::ConjugateGradient<tSparseMat, Eigen::Upper> solver;
        // tVectorXd delta_vel = solver.compute(A).solve(b);
        // tVectorXd delta_vel = solver.compute(A).solve(b);
        tVectorXd delta_vel = A.toDense().inverse() * b;
        tVectorXd v_next = vt + delta_vel;
        mXpre.noalias() = mXcur;
        mXcur = mXcur + v_next * dt;
    }

    // {
    //     /*
    //         explicit:
    //         M dv = dt * (f_total)
    //     */
    //     tVectorXd vt = (mXcur - mXpre) / dt;
    //     // tVectorXd dv = dt * this->mInvMassMatrixDiag.cwiseProduct((mGravityForce + mUserForce + mIntForce + mCollisionForce));
    //     tVectorXd dv = dt * this->mInvMassMatrixDiag.cwiseProduct(mIntForce + mUserForce);

    //     mXpre.noalias() = mXcur;
    //     mXcur = mXcur + (vt + dv) * dt;
    // }
}

void cBaraffCloth::UpdateCollisionForce(tVectorXd &col_force)
{
    int num_of_v = this->mVertexArrayShared.size();
    double ground_height = 1e-3;
    double k = 1e2;
    float KinectFrictionCoef = 0.5;
    // mExtForce.fill(5);
    // mExtForce[3 * 1 + 1] = 10;
    col_force.setZero();
    for (int i = 0; i < mVertexArrayShared.size(); i++)
    {
        double dist = mVertexArrayShared[i]->mPos[1] - ground_height;

        if (dist < 0)
        {
            mVertexArrayShared[i]->mPos[1] = 0;
            float normal_force_amp = -dist * k;

            tVector3d friction_dir = -1 * (mXcur.segment(3 * i, 3) - mXpre.segment(3 * i, 3));
            friction_dir[1] = 0;
            friction_dir.normalize();
            float friction_force_amp = KinectFrictionCoef * std::fabs(normal_force_amp);
            tVector3d force = tVector3d::Zero();
            force[1] = normal_force_amp;
            force[0] = friction_dir[0] * friction_force_amp;
            force[2] = friction_dir[2] * friction_force_amp;
            // std::cout << "[col] force " << force.transpose() << " on v" << i << std::endl;
            col_force.segment(3 * i, 3) = force;
        }
    }
}