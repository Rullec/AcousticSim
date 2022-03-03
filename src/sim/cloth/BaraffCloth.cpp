#include <iostream>
#include "BaraffCloth.h"
#include "sim/cloth/QBendingMaterial.h"
#include "utils/TimeUtil.hpp"
#include "geometries/Primitives.h"
#include "BaraffMaterial.h"
#include "utils/JsonUtil.h"
#include "sim/Perturb.h"
#include "imgui.h"
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
    // mMaterial = std::make_shared<cBaraffMaterial>();
    mBendingMaterial = std::make_shared<cQBendingMaterial>();
    mBendingK = cJsonUtil::ReadVectorJson(cJsonUtil::ParseAsValue("cloth_bending_stiffness", conf)).cast<float>();
    mBendingMaterial->Init(
        GetVertexArray(),
        GetEdgeArray(),
        GetTriangleArray(),
        mBendingK.cast<double>());
    mStretchK = cJsonUtil::ReadVectorJson(cJsonUtil::ParseAsValue("cloth_stretch_stiffness", conf)).segment(0, 3).cast<float>();
    // mMaterial->SetStretchK(mStretchK[0], mStretchK[1]);
    // mMaterial->SetSheaingK(mStretchK[2]);
    mMaterial = std::make_shared<cBaraffMaterial>();
    mMaterial->Init(shared_from_this(), mStretchK.cast<double>());
    // mMaterial->CheckStretchStiffnessMatrix();
    // mMaterial->CheckStretchForce();

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
            // mGravityForce[3 * i + 1] = -mVertexArrayShared[i]->mMass * 9.8;
        }
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
    std::cout << "--------------frame " << frame++ << "-----------\n";
    gProfRec.clear();
    dt = mIdealDefaultTimestep;
    // CheckForce();
    // 1. calculate stiffness matrix
    cTimeUtil::Begin("update material");
    {
        mMaterial->Update(true, true, true);
    }

    gProfRec.push_back(tTimeRecms("update material", cTimeUtil::End("update material", true)));
    // std::cout << "-------\n";
    // UpdateCollisionForce(mCollisionForce);
    cTimeUtil::Begin("calc_K");
    CalcStiffnessMatrix(mXcur, mStiffnessMatrix);

    {
        // mStiffnessMatrix = mMaterial->CalcTotalStiffnessMatrix() + this->mBendingMaterial->GetStiffnessMatrix();
        // auto diff_norm = (new_K - mStiffnessMatrix).cwiseAbs().norm();
        // std::cout << "new_K = " << new_K.toDense() << std::endl;
        // std::cout << "old_K = " << mStiffnessMatrix.toDense() << std::endl;
        // std::cout << "[K] diff_norm = " << diff_norm << std::endl;
    }

    gProfRec.push_back(tTimeRecms("calc_K", cTimeUtil::End("calc_K", true)));

    cTimeUtil::Begin("calc_fint");

    CalcIntForce(mXcur, mIntForce);
    
    // std::cout << "[fint] diff norm = " << (new_fint - mIntForce).norm() << std::endl;

    gProfRec.push_back(tTimeRecms("calc_fint", cTimeUtil::End("calc_fint", true)));
    // for (int i = 0; i < mVertexArrayShared.size(); i++)
    // {
    //     std::cout << "fint" << i << " = " << mIntForce.segment(3 * i, 3).transpose() << std::endl;
    // std::cout << "x" << i << " = " << mXcur.segment(3 * i, 3).transpose() << std::endl;
    // }
    // std::cout << "user force = " << mUserForce.transpose() << std::endl;
    cTimeUtil::Begin("solve");
    SolveForNextPos(dt);
    gProfRec.push_back(tTimeRecms("solve", cTimeUtil::End("solve", true)));
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
    int_force = mMaterial->CalcTotalForce() + mBendingMaterial->CalcForce(mXcur);
    // auto get_pos_mat = [](const tVectorXd &xcur,
    //                       int id0, int id1, int id2) -> tMatrix3d
    // {
    //     tMatrix3d val = tMatrix3d::Zero();
    //     val.col(0) = xcur.segment(3 * id0, 3);
    //     val.col(1) = xcur.segment(3 * id1, 3);
    //     val.col(2) = xcur.segment(3 * id2, 3);
    //     return val;
    // };
    // // auto get_uv_mat;

    // int dof = 3 * GetNumOfVertices();
    // int_force.noalias() = tVectorXd::Zero(dof);
    // for (int i = 0; i < GetNumOfTriangles(); i++)
    // {
    //     auto &t = this->mTriangleArrayShared[i];
    //     // 1. calc ele stiffness
    //     // 1.1 assemble current pos,
    //     // 1.2 give the 3x2d texture coords
    //     int v_id[3] = {t->mId0,
    //                    t->mId1,
    //                    t->mId2};
    //     // int id0 = t->mId0,
    //     //     id1 = t->mId1,
    //     //     id2 = t->mId2;
    //     tMatrix32d uv_coords;
    //     uv_coords.row(0) = mVertexArrayShared[v_id[0]]->muv.cast<double>();
    //     uv_coords.row(1) = mVertexArrayShared[v_id[1]]->muv.cast<double>();
    //     uv_coords.row(2) = mVertexArrayShared[v_id[2]]->muv.cast<double>();

    //     tVectorXd force = mMaterial->CalcForce(get_pos_mat(xcur, v_id[0], v_id[1], v_id[2]), uv_coords);
    //     // std::cout << "[old] triangle " << i << " force = " << force.transpose() << std::endl;
    //     int_force.segment(3 * v_id[0], 3) += force.segment(0, 3);
    //     int_force.segment(3 * v_id[1], 3) += force.segment(3, 3);
    //     int_force.segment(3 * v_id[2], 3) += force.segment(6, 3);
    // }

    // int_force += mBendingMaterial->CalcForce(mXcur);
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
    mMassMatrixDiag.noalias() = tVectorXd::Zero(dof);
    for (auto &t : mTriangleArrayShared)
    {
        // 1. total area
        auto v0 = mVertexArrayShared[t->mId0];
        auto v1 = mVertexArrayShared[t->mId1];
        auto v2 = mVertexArrayShared[t->mId2];

        double triangle_area = cMathUtil::CalcTriangleArea(v0->mPos,
                                                           v1->mPos,
                                                           v2->mPos);
        mMassMatrixDiag.segment(3 * t->mId0, 3) += triangle_area / 3 * mClothDensity * tVector3d::Ones();
        mMassMatrixDiag.segment(3 * t->mId1, 3) += triangle_area / 3 * mClothDensity * tVector3d::Ones();
        mMassMatrixDiag.segment(3 * t->mId2, 3) += triangle_area / 3 * mClothDensity * tVector3d::Ones();
    }
    for (int i = 0; i < GetNumOfVertices(); i++)
    {
        mVertexArrayShared[i]->mMass = mMassMatrixDiag[3 * i];
        // std::cout << "v" << i << " mass = " << mVertexArrayShared[i]->mMass << std::endl;
    }
}

void cBaraffCloth::CalcStiffnessMatrix(const tVectorXd &xcur,
                                       tSparseMat &K_global) const
{
    K_global = mMaterial->CalcTotalStiffnessMatrix() + this->mBendingMaterial->GetStiffnessMatrix();
//     auto get_pos_mat = [](const tVectorXd &xcur,
//                           int id0, int id1, int id2) -> tMatrix3d
//     {
//         tMatrix3d val = tMatrix3d::Zero();
//         val.col(0) = xcur.segment(3 * id0, 3);
//         val.col(1) = xcur.segment(3 * id1, 3);
//         val.col(2) = xcur.segment(3 * id2, 3);
//         return val;
//     };
//     // auto get_uv_mat;

//     tEigenArr<tTriplet>
//         total_triplets;
//     total_triplets.reserve(mTriangleArrayShared.size() * 3 * 3 * 9);

//     // __pragma(omp parallel for num_threads(3))
//     int th_id, nthreads;
// #pragma omp parallel for num_threads(8)
//     for (int idx = 0; idx < mTriangleArrayShared.size(); idx++)
//     {
//         // std::cout << "omp num thread = " << omp_get_num_threads() << std::endl;
//         // printf("Hello World from thread %d\n", omp_get_thread_num());
//         auto t = this->mTriangleArrayShared[idx];
//         std::vector<tTriplet> sub_triples = {};
//         sub_triples.reserve(3 * 3 * 9);
//         // 1. calc ele stiffness
//         // 1.1 assemble current pos,
//         // 1.2 give the 3x2d texture coords
//         int v_id[3] = {t->mId0,
//                        t->mId1,
//                        t->mId2};
//         // int id0 = t->mId0,
//         //     id1 = t->mId1,
//         //     id2 = t->mId2;
//         tMatrix32d uv_coords;
//         uv_coords.row(0) = mVertexArrayShared[v_id[0]]->muv.cast<double>();
//         uv_coords.row(1) = mVertexArrayShared[v_id[1]]->muv.cast<double>();
//         uv_coords.row(2) = mVertexArrayShared[v_id[2]]->muv.cast<double>();

//         auto ele_K = mMaterial->CalcStiffMatrix(get_pos_mat(xcur, v_id[0], v_id[1], v_id[2]), uv_coords);
//         // std::cout << "[old] triangle " << idx << " K = \n"
//         //   << ele_K << std::endl;
//         // 2. assemble
//         for (size_t i = 0; i < 3; i++)
//         {
//             size_t global_vi_idx = v_id[i];
//             for (size_t j = 0; j < 3; j++)
//             {
//                 size_t global_vj_idx = v_id[j];

//                 const tMatrix3d &ele_K_part = ele_K.block(3 * i, 3 * j, 3, 3);
//                 size_t i3 = 3 * global_vi_idx, j3 = 3 * global_vj_idx;
//                 sub_triples.push_back(tTriplet(i3, j3, ele_K_part(0, 0)));
//                 sub_triples.push_back(tTriplet(i3, j3 + 1, ele_K_part(0, 1)));
//                 sub_triples.push_back(tTriplet(i3, j3 + 2, ele_K_part(0, 2)));

//                 sub_triples.push_back(tTriplet(i3 + 1, j3, ele_K_part(1, 0)));
//                 sub_triples.push_back(tTriplet(i3 + 1, j3 + 1, ele_K_part(1, 1)));
//                 sub_triples.push_back(tTriplet(i3 + 1, j3 + 2, ele_K_part(1, 2)));

//                 sub_triples.push_back(tTriplet(i3 + 2, j3, ele_K_part(2, 0)));
//                 sub_triples.push_back(tTriplet(i3 + 2, j3 + 1, ele_K_part(2, 1)));
//                 sub_triples.push_back(tTriplet(i3 + 2, j3 + 2, ele_K_part(2, 2)));
//             }
//         }

// #pragma omp critical
//         total_triplets.insert(total_triplets.end(), sub_triples.begin(), sub_triples.end());
//     }
//     int dof = 3 * GetNumOfVertices();
//     K_global.resize(dof, dof);

//     K_global.setFromTriplets(total_triplets.begin(), total_triplets.end());
//     K_global += mBendingMaterial->GetStiffnessMatrix();
}

/**
 * \brief           solve for next pos by baraff scheme
 *  let A = (I + dt * damping_a) * I + dt * (damping_b - dt) * MInv * K
 *  let b = dt * Minv (fext + fint + (dt - damping_b) * K * vt) - damping_a * vt
 *  we will have A delta_v = b
*/
void cBaraffCloth::SolveForNextPos(double dt)
{
    int dof = GetNumOfFreedom();
    tSparseMat M(dof, dof);
    M.reserve(dof);
    // M.diagonal() = mMassMatDiag;
    for (size_t i = 0; i < dof; i++)
    {
        M.coeffRef(i, i) = mMassMatrixDiag[i];
        // I.coeff(i, i) = 1;
    }
    // std::cout << "M = " << M.rows() << " " << M.cols() << std::endl;
    // std::cout << "mStiffnessMatrix = " << mStiffnessMatrix.rows() << " " << mStiffnessMatrix.cols() << std::endl;
    tSparseMat W = M - dt * dt * mStiffnessMatrix;
    tVectorXd b = dt * dt * (mGravityForce + mUserForce + mIntForce) + dt * (W + dt * dt * mStiffnessMatrix) * (mXcur - mXpre) / dt;
    tVectorXd dx = mXcur - mXpre;
    float threshold = 1e-12, residual = 0;
    int iters = 0;
    Solve(W, b, dx, threshold, iters, residual);

    mXpre = mXcur;
    mXcur = mXcur + dx;
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

void cBaraffCloth::UpdateImGui()
{
    ImGui::SliderFloat("dampingA", &mRayleightA, 0, 10);
    ImGui::SliderFloat("dampingB", &mRayleightB, -1, 1);

    tVector3f stretch = mStretchK.cast<float>();
    tVector3f bending = mBendingK.cast<float>();
    // float stretch_warp = mStretchK[0],
    //       stretch_weft = mStretchK[1],
    //       stretch_bias = mStretchK[2];

    ImGui::SliderFloat("Kwarp", &stretch[0], 1e1, 1e4);
    ImGui::SliderFloat("Kweft", &stretch[1], 1e1, 1e4);
    ImGui::SliderFloat("Kbias", &stretch[2], 1e1, 1e4);

    if ((stretch - mStretchK).norm() > 1e-6)
    {
        mStretchK = stretch;
        mMaterial->SetK(mStretchK.cast<double>());
        std::cout << "change bending K = " << mStretchK.transpose() << std::endl;
    }

    ImGui::SliderFloat("Bending warp", &bending[0], 1e-5, 1e-2, "%.5f");
    ImGui::SliderFloat("Bending Kweft", &bending[1], 1e-5, 1e-2, "%.5f");
    ImGui::SliderFloat("Bending Kbias", &bending[2], 1e-5, 1e-2, "%.5f");
    if ((bending - mBendingK).norm() > 1e-6)
    {
        mBendingK = bending;
        mBendingMaterial->Init(
            GetVertexArray(),
            GetEdgeArray(),
            GetTriangleArray(),
            mBendingK.cast<double>());
        std::cout << "change bending K = " << mBendingK.transpose() << std::endl;
    }

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
    // std::cout << "[check force] num_fint = " << num_fint.transpose() << std::endl;
    // std::cout << "[check force] ana_fint = " << ana_fint.transpose() << std::endl;
    // std::cout << "[check force] f diff = " << f_diff.transpose() << std::endl;
    std::cout << "[check force] f diff norm = " << f_diff.norm() << std::endl;
}

void GetTrianglePosMatAndUV(
    tVector3d v0, tVector3d v1, tVector3d v2,
    tVector2f u0, tVector2f u1, tVector2f u2,
    tMatrix3d &pos_mat,
    tMatrix32d &uv_mat)
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
    for (auto &t : this->mTriangleArrayShared)
    {
        tMatrix3d pos_mat;
        tMatrix32d uv_mat;
        GetTrianglePosMatAndUV(xcur.segment(3 * t->mId0, 3),
                               xcur.segment(3 * t->mId1, 3),
                               xcur.segment(3 * t->mId2, 3),
                               mVertexArrayShared[t->mId0]->muv,
                               mVertexArrayShared[t->mId1]->muv,
                               mVertexArrayShared[t->mId2]->muv,
                               pos_mat, uv_mat);
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
tSparseMat mBlockJacobianPreconditioner;
void cBaraffCloth::Solve(const tSparseMat &A, const tVectorXd &b, tVectorXd &x, float threshold, int &iters, float &residual)
{
    // std::cout << "Eigen::nbThreads() = " << Eigen::nbThreads() << std::endl;
    int recalc_gap = 20;
    int max_iters = 100;
    tVectorXd r = b - A * x;

    tVectorXd PInv = A.diagonal().cwiseInverse();

    {
        int dof = A.rows();
        mBlockJacobianPreconditioner.resize(dof, dof);
        for (int i = 0; i < dof / 3; i++)
        {
            tMatrix3d Ainv = A.block(3 * i, 3 * i, 3, 3).toDense().inverse();
        }
    }

    // tVectorXd PInv = tVectorXd::Ones(b.size());
    // Filter(this->mConstraint_StaticPointIds, r);
    tVectorXd d = PInv.cwiseProduct(r);

    // double delta_new = r.dot(c);
    iters = 1;
    // double residual0 = residual;
    // tVectorXd q = tVectorXd::Zero(r.size());
    // tVectorXd s;
    // while (delta_new > threshold && iters < max_iters)
    tVectorXd Adi;
    tVectorXd rnext;
    tVectorXd PInvr, PInvrnext;
    while (r.norm() > 1e-4 && iters < max_iters)
    {
        Adi.noalias() = A * d;
        PInvr.noalias() = PInv.cwiseProduct(r);
        // Filter(mConstraint_StaticPointIds, q);
        double alpha = r.dot(PInvr) / (d.dot(Adi));
        x += alpha * d;

        // if (iters % recalc_gap == 0)
        // {
        //     r = b - A * x;
        //     Filter(this->mConstraint_StaticPointIds, r);
        // }
        // else
        rnext.noalias() = r - alpha * Adi;
        PInvrnext.noalias() = PInv.cwiseProduct(rnext);
        double beta = rnext.dot(PInvrnext) / (r.dot(PInvr));
        d = PInvrnext + beta * d;

        r.noalias() = rnext;
        if (iters % (max_iters / 10) == 0 || iters == 1)
        {
            std::cout << "iters " << iters << " residual = " << r.norm() << std::endl;
        }
        iters += 1;
    }
    std::cout << "iters " << iters << " residual = " << r.norm() << std::endl;
}