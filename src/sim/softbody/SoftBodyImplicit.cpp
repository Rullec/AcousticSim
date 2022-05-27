#include "SoftBodyImplicit.h"
#include "utils/TimeUtil.hpp"
#include <iostream>
#include "utils/LogUtil.h"
#include "geometries/Tetrahedron.h"
#include "sim/softbody/ThreeOrderTensor.h"
#include "sim/softbody/FourOrderTensor.h"
cFourOrderTensor CalcDEDF(const tMatrix3d &F);
// cFourOrderTensor CalcDPDF(const tMatrix3d &F);
cFourOrderTensor CalcDPDF_part1(const tMatrix3d &F);
cFourOrderTensor CalcDPDF_part2(const tMatrix3d &F);
cFourOrderTensor CalcDTrEIDF(const tMatrix3d &F);
cFourOrderTensor CalcDFDF(const tMatrix3d &F);

typedef std::function<tMatrixXd(const tMatrixXd &)> tForwardCalcFunc;
extern cFourOrderTensor CalcNumericalDerivaitve_MatrixWRTMatrix(tForwardCalcFunc func, const tMatrixXd &X_raw, const tVector2i output_dims, float eps = 1e-5);

cSoftBodyImplicit::cSoftBodyImplicit(int id_) : cSoftBody(id_)
{
}

cSoftBodyImplicit::~cSoftBodyImplicit()
{
}

typedef std::pair<std::string, float> tTimeRecms;
static std::vector<tTimeRecms> gProfRec;
void cSoftBodyImplicit::Update(float dt)
{
    gProfRec.clear();
    dt = 5e-3;
    cTimeUtil::Begin("update_F");
    UpdateDeformationGradient();
    gProfRec.push_back(tTimeRecms("update_F", cTimeUtil::End("update_F", true)));

    UpdateExtForce();

    // cTimeUtil::Begin("update_K_dense");
    // mGlobalStiffnessMatrix.noalias() = CalcGlobalStiffnessMatrix();
    // cTimeUtil::End("update_K_dense");

    cTimeUtil::Begin("update_K_sparse");
    mGlobalStiffnessSparseMatrix = CalcGlobalStiffnessSparseMatrix();
    // std::cout << "K diff = " << (mGlobalStiffnessMatrix - sparse_mat).norm() << std::endl;
    gProfRec.push_back(tTimeRecms("update_K_sparse", cTimeUtil::End("update_K_sparse", true)));

    cTimeUtil::Begin("update_fint");
    UpdateIntForce();
    gProfRec.push_back(tTimeRecms("update_fint", cTimeUtil::End("update_fint", true)));

    // cTimeUtil::Begin("solve_pos_dense");
    // SolveForNextPos(dt);
    // cTimeUtil::End("solve_pos_dense");

    cTimeUtil::Begin("solve_pos_sparse");
    SolveForNextPosSparse(dt);
    gProfRec.push_back(tTimeRecms("solve_pos_sparse", cTimeUtil::End("solve_pos_sparse", true)));

    SyncPosToVectorArray();
    mUserForce.setZero();
}
void cSoftBodyImplicit::Init(const Json::Value &conf)
{
    InitDDsDxTensor();
    cSoftBody::Init(conf);
    // CheckDPDx();
    // exit(1);
    // CheckDFDx();
    // CheckDFDF();
    // CheckDEDF();
    // CheckDPDF_part1();
    // CheckDPDF_part2();
    // CheckSelectionMat();
    // CheckElementStiffnessMat();
    // CheckGlobalStiffnessMat();
    // CheckDPDF();
    // CheckDFDF();
}

tEigenArr<tMatrix3d> cSoftBodyImplicit::CalcDFDx(int ele_id)
{
    tEigenArr<tMatrix3d> dfdx = this->mDDsDxTensor_const;
    for (int i = 0; i < dfdx.size(); i++)
    {
        dfdx[i] = (dfdx[i] * mInvDm[ele_id]).eval();
    }
    return dfdx;
}

// cFourOrderTensor CalcDPDF(const tMatrix3d &F)
// {

//     return CalcDPDF_part1(F) + CalcDPDF_part2(F);
// }

void cSoftBodyImplicit::UpdateIntForce()
{
    // SIM_ERROR("for implicit scheme, we do not calculate int force explicitly.");
    cSoftBody::UpdateIntForce();
    // exit(1);
    // CalculateStiffnessMatrix();
}

void cSoftBodyImplicit::SolveForNextPos(float dt)
{
    // tVectorXd MassDiag = mInvLumpedMassMatrixDiag.cwiseInverse();

    // tMatrixXd A = dt * (mRayleighDamplingB - dt) * this->mGlobalStiffnessMatrix;
    // A.diagonal() += (1 + mRayleighDamplingA * dt) * MassDiag;

    // tVectorXd b = MassDiag.cwiseProduct((mXcur - mXprev) / dt) + dt * (mExtForce + mIntForce + mGravityForce + mUserForce);
    // tVectorXd vel = A.inverse() * b;
    // mXprev = mXcur;
    // mXcur = dt * vel + mXcur;
    // if (mXcur.hasNaN())
    // {
    //     SIM_ERROR("Xcur hasNan, exit");
    //     exit(1);
    // }
}
#include <Eigen/IterativeLinearSolvers>
void cSoftBodyImplicit::SolveForNextPosSparse(float dt)
{
    tVectorXd MassDiag = mInvLumpedMassMatrixDiag.cwiseInverse();

    tSparseMatd A = dt * (mRayleighDamplingB - dt) * this->mGlobalStiffnessSparseMatrix;
    // std::cout << "A.diagonal() = " << A.diagonal().size() << std::endl;
    // std::cout << "MassDiag = " << MassDiag.size() << std::endl;
    A.diagonal() += (1 + mRayleighDamplingA * dt) * MassDiag;

    tVectorXd b = MassDiag.cwiseProduct((mXcur - mXprev) / dt) + dt * (mExtForce + mIntForce + mGravityForce + mUserForce);
    // tVectorXd vel = A.inverse() * b;

    Eigen::ConjugateGradient<tSparseMatd, Eigen::Upper> solver;
    tVectorXd vel = solver.compute(A).solve(b);

    mXprev = mXcur;
    mXcur = dt * vel + mXcur;
    if (mXcur.hasNaN())
    {
        SIM_ERROR("Xcur hasNan, exit");
        exit(1);
    }
    SyncPosToVectorArray();
}

/**
 * \brief           calculate dPdx:
 *  {dPdx}_{ijm} = {dP/dF}_{ijkl} * {dF/dx}_{klm}
*/
cThreeOrderTensor CalcDPDx(const cFourOrderTensor &dPdF,
                           const cThreeOrderTensor &dFdx)
{
    size_t dimI = dPdF.GetShape()[0];
    size_t dimJ = dPdF.GetShape()[1];
    size_t dimM = dFdx.GetShape()[2];

    tEigenArr<tMatrixXd> dpdx_lst(0);
    // we will get {dPdx}_{ijm}
    for (size_t m = 0; m < dimM; m++)
    {
        tMatrixXd dPdx_comp = tMatrixXd::Zero(dimI, dimJ);
        const tMatrixXd &dFdx_comp = dFdx(m);
        for (size_t i = 0; i < dimI; i++)
            for (size_t j = 0; j < dimJ; j++)
            {
                dPdx_comp(i, j) = dFdx_comp.cwiseProduct(dPdF(i, j)).sum();
            }
        dpdx_lst.push_back(dPdx_comp);
    }
    return cThreeOrderTensor(dpdx_lst);
}
/**
 * \brief           calculate sitffness matrix K
 * 
 *  1. calculate dPdF
 *  2. calculate dFdx
 *  3. calculate some selection matrix (for each tet)
 *  4. finish the element stiffness mat
 *  5. assembly
*/
// void cSoftBodyImplicit::CalculateStiffnessMatrix()
// {
//     // auto & tet : this->mTetArrayShared
//     for (size_t tet_id = 0; tet_id < mTetArrayShared.size(); tet_id++)
//     {
//         auto &tet = mTetArrayShared[tet_id];
//         auto &F = mF[tet_id];
//         const cFourOrderTensor &dPdF = CalcDPDF(F);
//         const cThreeOrderTensor &dFdx = cThreeOrderTensor(CalcDFDx(tet_id));
//         const cThreeOrderTensor &dPdx = CalcDPDx(dPdF, dFdx);
//     }
// }

void cSoftBodyImplicit::CheckDFDx()
{
    UpdateDeformationGradient();
    auto F_old_lst = this->mF;

    for (int ele_id = 0; ele_id < GetNumOfTets(); ele_id++)
    {
        // calculate dfdx analytically
        tEigenArr<tMatrix3d> dfdx_ana_lst = CalcDFDx(ele_id);
        tMatrix3d F_old = mF[ele_id];

        // calculate dfdx numerically
        auto cur_tet = mTetArrayShared[ele_id];
        tVector4i v_id = cur_tet->mVertexId;
        float eps = 1e-5;
        // calculate for a point
        for (int j = 0; j < 4; j++)
        {
            int cur_v = v_id[j];
            for (int comp = 0; comp < 3; comp += 1)
            {
                // calculate new deformation gradient
                mXcur[3 * cur_v + comp] += eps;
                SyncPosToVectorArray();
                UpdateDeformationGradientForTet(ele_id);

                tMatrix3d F_new = mF[ele_id];
                tMatrix3d dFdx_num = (F_new - F_old) / eps;
                tMatrix3d dFdx_ana = dfdx_ana_lst[3 * j + comp];
                tMatrix3d dFdx_diff = dFdx_num - dFdx_ana;
                float dFdx_diff_norm = dFdx_diff.norm();
                if (dFdx_diff_norm > 10 * eps)
                {
                    std::cout << "[error] failed, dFdx diff norm = " << dFdx_diff_norm << std::endl;
                    std::cout << "dFdx num = \n"
                              << dFdx_num << std::endl;
                    std::cout << "dFdx ana = \n"
                              << dFdx_ana << std::endl;
                    std::cout << "dFdx diff = \n"
                              << dFdx_diff << std::endl;
                }

                mXcur[3 * cur_v + comp] -= eps;
                SyncPosToVectorArray();
                UpdateDeformationGradientForTet(ele_id);
            }
        }
    }

    SIM_INFO("check df/dx succ");
}

#include "sim/softbody/FourOrderTensor.h"

// forward calculate function

void cSoftBodyImplicit::CheckDFDF()
{
    tMatrix3d F = tMatrix3d::Random();
    tMatrixXd DFDF_ana = CalcDFDF(F).ExpandToMatrix();
    auto CalcItself = [](const tMatrix3d &val) -> tMatrix3d
    {
        return val;
    };
    tMatrixXd DFDF_num = CalcNumericalDerivaitve_MatrixWRTMatrix(CalcItself, F, tVector2i(3, 3)).ExpandToMatrix();
    tMatrixXd DFDF_diff = DFDF_ana - DFDF_num;
    double DFDF_diff_norm = DFDF_diff.norm();
    std::cout << "dFdF diff norm = " << DFDF_diff_norm << std::endl;
    if (DFDF_diff_norm > 1e-4)
    {
        std::cout << "[error] dFdF ana = \n"
                  << DFDF_ana << std::endl;
        std::cout << "[error] DFDF num = \n"
                  << DFDF_num << std::endl;
        exit(1);
    }
    // else{
    //     std::cout << "dFdF = \n"
    //               << DFDF_ana << std::endl;
    // }
}
extern tMatrix3d CalcGreenStrain(const tMatrix3d &F);
void cSoftBodyImplicit::CheckDEDF()
{
    // CalcDPDF();
    // 1. check dEdF (E is green strain)
    tMatrix3d F = tMatrix3d::Random();
    auto dEdF_num = CalcNumericalDerivaitve_MatrixWRTMatrix(CalcGreenStrain, F, tVector2i(3, 3));
    std::cout << "dEdF shape = " << dEdF_num.GetShape().transpose() << std::endl;
    auto dEdF_ana = CalcDEDF(F);
    {
        auto dEdF_num_mat = dEdF_num.ExpandToMatrix();
        auto dEdF_ana_mat = dEdF_ana.ExpandToMatrix();
        auto dEdF_diff = dEdF_ana_mat - dEdF_num_mat;
        double dEdF_diff_norm = dEdF_diff.norm();
        std::cout << "dEdF diff norm = " << dEdF_diff_norm << std::endl;
        if (dEdF_diff_norm > 1e-4)
        {
            std::cout << "[error] dEdF ana = \n"
                      << dEdF_ana_mat << std::endl;
            std::cout << "[error] dEdF num = \n"
                      << dEdF_num_mat << std::endl;
            exit(1);
        }
    }
}

// void cSoftBodyImplicit::CheckDPDF_part2()
// {
// tMatrix3d F = tMatrix3d::Random();
// // 1. check dPdF (P is PK1)
// auto dPdF_part2_num = CalcNumericalDerivaitve_MatrixWRTMatrix(CalcPK1_part2, F, tVector2i(3, 3)).ExpandToMatrix();

// auto dPdF_part2_ana = CalcDPDF_part2(F).ExpandToMatrix();
// auto dPdF_part2_diff = (dPdF_part2_ana - dPdF_part2_num) / (2 * gMu);
// double dPdF_part2_diff_norm = dPdF_part2_diff.norm();
// std::cout << "dPdF_part2_diff_norm = " << dPdF_part2_diff_norm << std::endl;
// if (dPdF_part2_diff_norm > 1e-4)
// {
//     std::cout << "dPdF_part2_num = \n"
//               << dPdF_part2_num << std::endl;
//     std::cout << "dPdF_part2_ana = \n"
//               << dPdF_part2_ana << std::endl;
//     std::cout << "dPdF_part2_diff = \n"
//               << dPdF_part2_diff << std::endl;
//     exit(1);
// }
// }
// void cSoftBodyImplicit::CheckDPDF() {}
// void cSoftBodyImplicit::CheckDPDF_part1()
// {
//     tMatrix3d F = tMatrix3d::Random();
//     // 1. check dPdF (P is PK1)
//     auto dPdF_part1_num = CalcNumericalDerivaitve_MatrixWRTMatrix(CalcPK1_part1, F, tVector2i(3, 3)).ExpandToMatrix();

//     auto dPdF_part1_ana = CalcDPDF_part1(F).ExpandToMatrix();
//     auto dPdF_part1_diff = (dPdF_part1_ana - dPdF_part1_num) / (2 * gMu);
//     double dPdF_part1_diff_norm = dPdF_part1_diff.norm();
//     std::cout << "dPdF_part1_diff_norm = " << dPdF_part1_diff_norm << std::endl;
//     if (dPdF_part1_diff_norm > 1e-4)
//     {
//         std::cout << "dPdF_part1_num = \n"
//                   << dPdF_part1_num << std::endl;
//         std::cout << "dPdF_part1_ana = \n"
//                   << dPdF_part1_ana << std::endl;
//         std::cout << "dPdF_part1_diff = \n"
//                   << dPdF_part1_diff << std::endl;
//         exit(1);
//     }
//     // std::cout << "dPdF_part1 shape = " << dPdF_part1_num.GetShape().transpose() << std::endl;
// }
void cSoftBodyImplicit::InitDDsDxTensor()
{
    mDDsDxTensor_const.resize(12, tMatrix3d::Zero());
    for (int i = 0; i < 12; i++)
    {
        // for the first 3 points
        if (i < 9)
        {
            int row_id = i % 3;
            int col_id = std::floor(i / 3);
            mDDsDxTensor_const[i](row_id, col_id) = 1;
        }
        else
        {
            // for the last point
            mDDsDxTensor_const[i].row(i - 9).noalias() = tVector3d::Ones() * -1;
        }
        // printf("d(Ds)/dx[%d] = ", i);
        // std::cout << mDDsDxTensor_const[i] << std::endl;
    }
}

// cFourOrderTensor CalcDEDF(const tMatrix3d &F)
// {
//     tMatrix3d FT = F.transpose();
//     cFourOrderTensor deriv(3, 3, 3, 3);
//     for (int i = 0; i < 3; i++)
//         for (int j = 0; j < 3; j++)
//             for (int k = 0; k < 3; k++)
//                 for (int l = 0; l < 3; l++)
//                 {
//                     if (l == j)
//                     {
//                         deriv(k, l)(i, j) += 0.5 * FT(k, i);
//                     }
//                     if (k == j)
//                     {
//                         deriv(k, l)(i, j) += 0.5 * F(i, l);
//                     }
//                 }
//     return deriv;
// }

cFourOrderTensor CalcDFDF(const tMatrix3d &F)
{
    cFourOrderTensor deriv(3, 3, 3, 3);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
                for (int l = 0; l < 3; l++)
                {
                    if (i == k && j == l)
                    {
                        deriv(i, j)(k, l) = 1;
                    }
                }
    return deriv;
}

// cFourOrderTensor CalcDPDF_part1(const tMatrix3d &F)
// {
//     cFourOrderTensor DFDF = CalcDFDF(F);
//     tMatrix3d E = CalcGreenStrain(F);
//     cFourOrderTensor DEDF = CalcDEDF(F);
//     DFDF.TensorMatrixProductWithoutCopy(0, 1, E);
//     DEDF.MatrixTensorProductWithoutCopy(0, 1, F);
//     cFourOrderTensor deriv = (DFDF + DEDF) * 2 * gMu;
//     return deriv;
// }

// cFourOrderTensor CalcDPDF_part2(const tMatrix3d &F)
// {
//     // tMatrix3d E = CalcGreenStrain(F);
//     // cFourOrderTensor DTrEIDF = CalcDTrEIDF(F),
//     //                  trE_DFDF = CalcDFDF(F) * E.trace();
//     // DTrEIDF.TensorMatrixProductWithoutCopy(0, 1, F);
//     // return (DTrEIDF + trE_DFDF) * gLambda;
// }

cFourOrderTensor CalcDTrEIDF(const tMatrix3d &F)
{
    cFourOrderTensor dEdF = CalcDEDF(F);
    tMatrix3d sum = dEdF(0, 0) + dEdF(1, 1) + dEdF(2, 2);
    cFourOrderTensor DTrEIDF(3, 3, 3, 3);
    DTrEIDF(0, 0) = sum;
    DTrEIDF(1, 1) = sum;
    DTrEIDF(2, 2) = sum;
    return DTrEIDF;
}

void cSoftBodyImplicit::CheckElementStiffnessMat(int tet_id)
{
    UpdateDeformationGradient();
    tMatrixXd K_ana = CalcElementStiffnessMatrix(tet_id);
    // std::cout << "K_ana = \n"
    //           << K_ana << std::endl;
    // 1. get old force
    cSoftBody::UpdateIntForce();
    tVectorXd old_force_total = mIntForce;

    // 2. change x, get new force, get deriv
    tVectorXd old_force = GetTetForce(tet_id, old_force_total);
    tVector4i tet_vertices_id = this->mTetArrayShared[tet_id]->mVertexId;
    std::vector<int> v_lst = {};
    for (size_t i = 0; i < 4; i++)
    {
        v_lst.push_back(3 * tet_vertices_id[i] + 0);
        v_lst.push_back(3 * tet_vertices_id[i] + 1);
        v_lst.push_back(3 * tet_vertices_id[i] + 2);
    }
    float eps = 1e-9;
    tMatrixXd K_num = tMatrixXd::Zero(12, 12);
    // std::cout << "old force = " << old_force.transpose() << std::endl;
    // std::cout << "old F = \n"
    //           << mF[tet_id] << std::endl;
    // std::cout << "old P = \n"
    //           << CalcPK1(mF[tet_id]) << std::endl;
    // std::cout << "old force recalc = " << CalcTetIntForce(tet_id).transpose() << std::endl;
    // exit(1);
    for (size_t _idx = 0; _idx < 12; _idx++)
    {
        size_t pos = v_lst[_idx];
        mXcur[pos] += eps;
        SyncPosToVectorArray();
        UpdateDeformationGradientForTet(tet_id);
        cSoftBody::UpdateIntForce();
        tVectorXd new_force = GetTetForce(tet_id, mIntForce);
        K_num.col(_idx) = (new_force - old_force) / eps;
        std::cout << "new force " << _idx << " = " << new_force.transpose() << std::endl;
        mXcur[pos] -= eps;
        SyncPosToVectorArray();
        UpdateDeformationGradientForTet(tet_id);
        cSoftBody::UpdateIntForce();
    }
    // std::cout << "K num = \n"
    //           << K_num << std::endl;
    // std::cout << "K ana = \n"
    //           << K_ana << std::endl;
    tMatrixXd K_diff = K_num - K_ana;
    double K_diff_norm = K_diff.norm();
    std::cout << "[info] K diff norm = " << K_diff_norm << std::endl;
    if (K_diff_norm > 1)
    {
        std::cout << "[error] K_ana = " << K_ana << std::endl;
        std::cout << "[error] K_num = " << K_num << std::endl;
        std::cout << "[error] K_diff = " << K_diff << std::endl;
    }
}

void cSoftBodyImplicit::CheckGlobalStiffnessMat()
{
    // 1. calculate global stiffness
    UpdateDeformationGradient();
    int global_dof = mXcur.size();
    tMatrixXd globalK_ana = CalcGlobalStiffnessMatrix();
    tMatrixXd globalK_num = tMatrixXd::Zero(global_dof, global_dof);

    // 2. calculate global force
    UpdateIntForce();
    tVectorXd globalF_old = mIntForce;

    double eps = 1e-9;
    for (size_t i = 0; i < mXcur.size(); i++)
    {
        mXcur[i] += eps;
        SyncPosToVectorArray();
        UpdateDeformationGradient();
        UpdateIntForce();
        tVectorXd globalK_num_comp = (mIntForce - globalF_old) / eps;
        globalK_num.row(i) = globalK_num_comp;

        mXcur[i] -= eps;
        SyncPosToVectorArray();
        UpdateDeformationGradient();
        UpdateIntForce();
    }
    tMatrixXd globalK_diff = globalK_ana - globalK_num;
    double globalK_diff_norm = globalK_diff.norm();
    std::cout << "globalK_diff_norm = " << globalK_diff_norm << std::endl;
    if (globalK_diff_norm > 1)
    {
        std::cout << "[error] check failed: globalK_diff = " << globalK_diff << std::endl;
        std::cout << "[error] check failed: globalK_ana = " << globalK_ana << std::endl;
        std::cout << "[error] check failed: globalK_num = " << globalK_num << std::endl;
        exit(1);
    }

    return;
}

tEigenArr<tMatrix3d> GetTTensor()
{
    tEigenArr<tMatrix3d> T(12, tMatrix3d::Zero());
    for (int i = 0; i < 9; i++)
    {
        T[i](i % 3, i / 3) = 1;
    }
    T[9].row(0) = tVector3d::Ones() * -1;
    T[10].row(1) = tVector3d::Ones() * -1;
    T[11].row(2) = tVector3d::Ones() * -1;
    return T;
}
/**
 * \brief       calculate element stiffness matrix
*/
tMatrixXd cSoftBodyImplicit::CalcElementStiffnessMatrix(int tet_id)
{
    cThreeOrderTensor dPdx = CalcDPDx(mMaterial->CalcDPDF(mF[tet_id]), CalcDFDx(tet_id));

    dPdx.TensorMatrixProductWithoutCopy(0, 1, mInvDm[tet_id].transpose());

    // begin to do selection
    auto T_lst = GetTTensor();
    tMatrixXd K = tMatrixXd::Zero(12, 12);
    for (size_t i = 0; i < 12; i++)
        for (size_t j = 0; j < 12; j++)
        {
            K(i, j) = -mInitTetVolume[tet_id] * T_lst[i].cwiseProduct(dPdx(j)).sum();
        }
    // std::cout << "ele volume " << tet_id << " = " << mInitTetVolume[tet_id] << std::endl;
    // std::cout << "K_ele for tet " << tet_id << " = \n" << K << std::endl;
    return K;
}

void cSoftBodyImplicit::CheckDPDx()
{
    UpdateDeformationGradient();
    int tet_id = 0;

    cThreeOrderTensor dPdx_ana = CalcDPDx(mMaterial->CalcDPDF(mF[tet_id]), CalcDFDx(tet_id));

    // 1. calculate P_old
    tVectorXd tet_pos = GetTetVerticesPos(tet_id, mXcur);
    UpdateDeformationGradientForTet(tet_id);

    tMatrix3d P_old = mMaterial->CalcP(mF[tet_id]);

    double eps = 1e-9;
    for (size_t i = 0; i < 12; i++)
    {
        // 1. add eps over the displacement
        tet_pos[i] += eps;
        SetTetVerticesPos(tet_id, tet_pos);

        // 2. update F
        UpdateDeformationGradientForTet(tet_id);

        // 3. calculate P new
        tMatrix3d P_new = mMaterial->CalcP(mF[tet_id]);

        // 2. calculate dPdx_num_component
        tMatrix3d dPdx_num_comp = (P_new - P_old) / eps;
        tMatrix3d dPdx_ana_comp = dPdx_ana.GetInnerArray()[i];
        tMatrix3d dPdx_diff = (dPdx_num_comp - dPdx_ana_comp);
        double dPdx_diff_norm = dPdx_diff.norm();
        std::cout << "dPdx_diff_norm = " << dPdx_diff_norm << std::endl;
        if (dPdx_diff_norm > 1)
        {
            std::cout << "dPdx_ana = \n"
                      << dPdx_ana_comp << std::endl;
            std::cout << "dPdx_num = \n"
                      << dPdx_num_comp << std::endl;
            std::cout << "dPdx_diff = \n"
                      << dPdx_diff << std::endl;
            exit(1);
        }

        // 3. add it into the final result

        // 4. restore F, restore displacement
        tet_pos[i] -= eps;
        SetTetVerticesPos(tet_id, tet_pos);
        UpdateDeformationGradientForTet(tet_id);
    }
}

void cSoftBodyImplicit::SetTetVerticesPos(size_t tet_id, const tVectorXd &tet_vertices_pos)
{
    for (size_t i = 0; i < 4; i++)
    {
        size_t v_id = mTetArrayShared[tet_id]->mVertexId[i];
        mXcur.segment(3 * v_id, 3) = tet_vertices_pos.segment(3 * i, 3);
        mVertexArray[v_id]->mPos.segment(0, 3) = tet_vertices_pos.segment(3 * i, 3);
    }
}

// void cSoftBodyImplicit::CheckSelectionMat()
// {
//     UpdateDeformationGradient();
//     size_t tet_id = 0;
//     // 1. calculate old tet internal force
//     tVectorXd tet_int_force_old = CalcTetIntForce(tet_id);
//     // 2. calculate tet internal force by selection mat
//     tVectorXd tet_int_force_sel = CalcTetIntForceBySelectionMatrix(tet_id);
//     // 3. compare
//     tVectorXd diff = tet_int_force_sel - tet_int_force_old;
//     std::cout << "diff = " << diff.transpose() << std::endl;
//     std::cout << "tet_int_force_sel = " << tet_int_force_sel.transpose() << std::endl;
//     std::cout << "tet_int_force_old = " << tet_int_force_old.transpose() << std::endl;
//     exit(1);
// }

/**
 * \brief               global stiffness matrix assembly
*/
tMatrixXd cSoftBodyImplicit::CalcGlobalStiffnessMatrix()
{
    // matrix assmebly
    // dense storage now
    int global_dof = this->mXcur.size();
    tMatrixXd global_K = tMatrixXd::Zero(global_dof, global_dof);
    for (size_t tet_id = 0; tet_id < GetNumOfTets(); tet_id++)
    {
        auto tet = mTetArrayShared[tet_id];
        tMatrixXd ele_K = CalcElementStiffnessMatrix(tet_id);
        for (size_t i = 0; i < 4; i++)
        {
            size_t global_vi_idx = tet->mVertexId[i];
            for (size_t j = 0; j < 4; j++)
            {
                size_t global_vj_idx = tet->mVertexId[j];

                global_K.block(3 * global_vi_idx, 3 * global_vj_idx, 3, 3) += ele_K.block(3 * i, 3 * j, 3, 3);
            }
        }
    }
    return global_K;
}

tSparseMatd cSoftBodyImplicit::CalcGlobalStiffnessSparseMatrix()
{
    size_t num_of_v = this->mVertexArray.size();
    std::vector<tTriplet> tripletList = {};
    tSparseMatd mat(3 * num_of_v, 3 * num_of_v);
    tripletList.reserve(GetNumOfTets() * 4 * 4 * 9);
    OMP_PARALLEL_FOR
    for (int tet_id = 0; tet_id < GetNumOfTets(); tet_id++)
    {
        std::vector<tTriplet> sub_triples = {};
        // std::cout << "omp num thread = " << omp_get_num_threads() << std::endl;
        auto tet = mTetArrayShared[tet_id];
        tMatrixXd ele_K = CalcElementStiffnessMatrix(tet_id);
        for (size_t i = 0; i < 4; i++)
        {
            size_t global_vi_idx = tet->mVertexId[i];
            for (size_t j = 0; j < 4; j++)
            {
                size_t global_vj_idx = tet->mVertexId[j];

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

#pragma omp critical
        tripletList.insert(tripletList.end(), sub_triples.begin(), sub_triples.end());
    }

    mat.setFromTriplets(tripletList.begin(), tripletList.end());
    return mat;
}
#include <imgui.h>
void cSoftBodyImplicit::UpdateImGui()
{
    cSoftBody::UpdateImGui();
    for (auto &x : gProfRec)
    {
        ImGui::Text("%s %d ms", x.first.c_str(), int(x.second));
    }
}

tSparseMatd cSoftBodyImplicit::GetGlobalStiffnessMatrix() const
{
    return mGlobalStiffnessSparseMatrix;
}
tVectorXd cSoftBodyImplicit::GetMassMatrixDiag()
{
    return mInvLumpedMassMatrixDiag.cwiseInverse();
}
tVector2f cSoftBodyImplicit::GetRayleightDamping()
{
    return tVector2f(mRayleighDamplingA, mRayleighDamplingB);
}