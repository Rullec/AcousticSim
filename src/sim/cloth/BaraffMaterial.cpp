#include "BaraffMaterial.h"
#include "geometries/Primitives.h"
#include "sim/BaseObject.h"
#include "utils/RotUtil.h"
#include <iostream>
tMatrix9d KroneckerProduct(const tMatrix3d &a, const tMatrix3d &b)
{
    tMatrix9d res = tMatrix9d::Zero();
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
        {
            res.block(3 * i, j * 3, 3, 3) = a(i, j) * b;
        }
    return res;
}

cBaraffMaterial::cBaraffMaterial()
{
    mS.setZero();
    mS(0, 0) = -1;
    mS(0, 1) = -1;
    mS(1, 0) = 1;
    mS(2, 1) = 1;

    mKwarpweftshear[0] = 1e4;
    mKwarpweftshear[1] = 1e4;
    mKwarpweftshear[2] = 1e4;
}

void cBaraffMaterial::Init(cBaseObjectPtr object,
                           const tVector3d &Kwarpweftbias)
{
    mObject = object;
    mNumOfVertices = mObject->GetNumOfVertices();
    mNumOfTriangles = mObject->GetNumOfTriangles();
    mKwarpweftshear = Kwarpweftbias;
    Allocate();
    InitN();
}
void cBaraffMaterial::Allocate()
{
    mNLst.resize(mNumOfTriangles);
    mNprimeLst.resize(mNumOfTriangles);
    mFLst.resize(mNumOfTriangles);
    mFprimeLst.resize(mNumOfTriangles);
    mCLst.resize(mNumOfTriangles);
    mCprimeLst.resize(mNumOfTriangles);

    mnLst.resize(mNumOfTriangles);
    mnprimeLst.resize(mNumOfTriangles);
    mgLst.resize(mNumOfTriangles);
    mgprimeLst.resize(mNumOfTriangles);

    mIntForceLst.resize(mNumOfTriangles);
    mKLst.resize(mNumOfTriangles);
    mELst.resize(mNumOfTriangles);
    global_K_buf.resize(3 * mNumOfVertices, 3 * mNumOfVertices);
    global_K_buf.reserve(mNumOfTriangles * 3 * 3 * 9);
    total_triplets_buf.reserve(mNumOfTriangles * 3 * 3 * 9);
}
void cBaraffMaterial::InitN()
{
    auto &tri_lst = mObject->GetTriangleArray();
    auto &v_lst = mObject->GetVertexArray();
    tMatrix2d R_neg45 = cRotUtil::RotMat2D(-M_PI / 4);
    for (int i = 0; i < mNumOfTriangles; i++)
    {
        auto t = tri_lst[i];

        tMatrix2d DmInv = tMatrix2d::Zero();
        DmInv.col(0) =
            (v_lst[t->mId1]->muv - v_lst[t->mId0]->muv).cast<double>();
        DmInv.col(1) =
            (v_lst[t->mId2]->muv - v_lst[t->mId0]->muv).cast<double>();

        DmInv = DmInv.inverse().eval();

        // 2. calculate N' = S * Dminv * R^{-45}
        mNLst[i].noalias() = mS * DmInv;
        mNprimeLst[i].noalias() = mNLst[i] * R_neg45;
    }
}
/**
 * \brief           f = -dE/dx = - k_i * C_i * g_i
 * where g_i = N_i \otimes n_i
 * For more details, please check the note "Baraff
 * 内力和刚度矩阵计算(利用Kronecker简化).md"
 */
tVector9d
cBaraffMaterial::CalcStretchForce(const tMatrix3d &pos,
                                  const tMatrix32d &rest_texture_coords) const
{
    // 1. calculate F and N
    // std::cout << "--------[calc force]------\n";
    tMatrix32d F, N;
    // std::cout << "pos = \n"
    //   << pos << std::endl;
    // std::cout << "uv = \n"
    //   << rest_texture_coords << std::endl;
    CalcFAndN(pos, rest_texture_coords, F, N);
    // std::cout << "F = \n"
    //   << F << std::endl;
    // 2. calculate C
    tVector2f C = CalcC(F);
    // std::cout << "C = "
    //   << C.transpose() << std::endl;
    tMatrix32d n = Calcn(pos, N);
    tMatrixXd g = Calcg(N, n);
    tVector9d int_force = (-mKwarpweftshear[0] * C[0]) * g.col(0) +
                          (-mKwarpweftshear[1] * C[1]) * g.col(1);
    // std::cout << "force = " << int_force.transpose() << std::endl;
    return int_force;
}

void cBaraffMaterial::Update(bool calc_energy /*= true*/,
                             bool calc_fint /*= true*/,
                             bool calc_stiffnessmatrix /*= true*/)
{
    // 1. update value
    auto &tri_lst = mObject->GetTriangleArray();
    auto &v_lst = mObject->GetVertexArray();

    // #pragma omp parallel for num_threads(15)
    for (int i = 0; i < mNumOfTriangles; i++)
    {
        tMatrix3d X;
        auto t = tri_lst[i];
        X.col(0) = v_lst[t->mId0]->mPos.segment(0, 3);
        X.col(1) = v_lst[t->mId1]->mPos.segment(0, 3);
        X.col(2) = v_lst[t->mId2]->mPos.segment(0, 3);

        const tMatrix32d &F = X * mNLst[i];
        mnLst[i].noalias() = F.colwise().normalized();
        mCLst[i].noalias() =
            tVector2d(F.col(0).norm() - 1, F.col(1).norm() - 1);

        mFLst[i].noalias() = F;

        const tMatrix32d &Fprime = X * mNprimeLst[i];
        mnprimeLst[i].noalias() = Fprime.colwise().normalized();
        mCprimeLst[i].noalias() =
            tVector2d(Fprime.col(0).norm() - 1, Fprime.col(1).norm() - 1);
        // std::cout << "tri " << i << " stretch C = " << mCLst[i].transpose()
        //           << " C' = " << mCprimeLst[i].transpose() << std::endl;
        mFprimeLst[i].noalias() = Fprime;

        mgLst[i].noalias() = cBaraffMaterial::Calcg(mNLst[i], mnLst[i]);
        mgprimeLst[i].noalias() =
            cBaraffMaterial::Calcg(mNprimeLst[i], mnprimeLst[i]);

        // 1.0 calculate energy
        mELst[i] = 0;
        if (calc_energy)
        {
            mELst[i] =
                0.5 * (mCLst[i][0] * mCLst[i][0] * mKwarpweftshear[0] +
                       mCLst[i][1] * mCLst[i][1] * mKwarpweftshear[1]) +
                0.5 *
                    (mCprimeLst[i][0] * mCprimeLst[i][0] * mKwarpweftshear[2] +
                     mCprimeLst[i][1] * mCprimeLst[i][1] * mKwarpweftshear[2]);
        }
        // 1.1. calculate fint: stretch and shearing
        mIntForceLst[i].setZero();
        if (calc_fint)
        {
            tVector9d fstretch =
                (-mKwarpweftshear[0] * mCLst[i][0]) * mgLst[i].col(0) +
                (-mKwarpweftshear[1] * mCLst[i][1]) * mgLst[i].col(1);

            tVector9d fshear =
                -mKwarpweftshear[2] * (mCprimeLst[i][0] * mgprimeLst[i].col(0) +
                                       mCprimeLst[i][1] * mgprimeLst[i].col(1));

            mIntForceLst[i].noalias() = fstretch + fshear;
            // std::cout << "tri " << i << " f cwiseabs max= "
            //           << mIntForceLst[i].cwiseAbs().maxCoeff() << std::endl;
            // std::cout << "[] triangle " << i << " force = " << (fstretch +
            // fshear).transpose() << std::endl; std::cout << "C = " <<
            // mCLst[i].transpose() << std::endl; std::cout << "g = \n" <<
            // mgLst[i] << std::endl; std::cout << "Cprime = " <<
            // mCprimeLst[i].transpose() << std::endl; std::cout << "gprime = "
            // << mgprimeLst[i].transpose() << std::endl;
            // mIntForceLst[i].noalias() = fstretch;
        }

        // 1.2 calculate K: stretch and shearing
        mKLst[i].setZero();
        if (calc_stiffnessmatrix)
        {

            // stretch part
            tMatrix3d P[2];
            {
                tVector2d X_Ni_norm = F.colwise().norm();

                // 3. calculate P_i

                CalcPi(mnLst[i], P[0], P[1]);

                // 4. add
                tMatrix9d K = tMatrix9d::Zero();
                for (size_t j = 0; j < 2; j++)
                {
                    K.noalias() +=
                        mKwarpweftshear[j] *
                        (-mgLst[i].col(j) * mgLst[i].col(j).transpose() -
                         mCLst[i][j] / X_Ni_norm[j] *
                             KroneckerProduct(mNLst[i].col(j) *
                                                  mNLst[i].col(j).transpose(),
                                              P[j])

                        );
                }
                mKLst[i].noalias() = K;
            }

            {
                tVector2d X_Ni_norm = Fprime.colwise().norm();
                CalcPi(mnprimeLst[i], P[0], P[1]);

                // 4. add
                tMatrix9d K = tMatrix9d::Zero();
                for (size_t j = 0; j < 2; j++)
                {
                    K.noalias() +=
                        mKwarpweftshear[2] *
                        (-mgprimeLst[i].col(j) *
                             mgprimeLst[i].col(j).transpose() -
                         mCprimeLst[i][j] / X_Ni_norm[j] *

                             KroneckerProduct(
                                 mNprimeLst[i].col(j) *
                                     mNprimeLst[i].col(j).transpose(),
                                 P[j])

                        );
                }

                mKLst[i] += K;
            }
            // std::cout << "[] triangle " << i << " K = \n" << mKLst[i] <<
            // std::endl;
        }
    }

    // 3. update internal force
}

void cBaraffMaterial::UpdateFnC() {}

// void cBaraffMaterial::UpdateC()
// {
//     for (int i = 0; i < mNumOfTriangles; i++)
//     {
//         mCLst[i].noalias() = CalcC(mFLst[i]);
//         mCprimeLst[i].noalias() = CalcC(mFprimeLst[i]);
//     }
// }
// tMatrix9d cBaraffMaterial::CalcStiffMatrix(const tMatrix3d &pos, const
// tMatrix32d &rest_texture_coords) const
// {
//     return tMatrix9d::Zero();
// }

/**
 * \brief           E = 1/2 * C^T * K * C
 *      let F = [Wu, Wv] \in R^{3 \times 2}
 *      C is the "baraff stretch condition" = [|Wu| - Bu, |Wv| - Bv]
 *      K is the stiffness mat = diag(Kwarp, Kweft)
 */
double
cBaraffMaterial::CalcStretchEnergy(const tMatrix3d &pos,
                                   const tMatrix32d &rest_texture_coords) const
{
    // 1. calculate F and N
    tMatrix32d F, N;
    CalcFAndN(pos, rest_texture_coords, F, N);

    // 2. calculate C
    tVector2f C = CalcC(F);
    return 0.5 * (C[0] * C[0] * mKwarpweftshear[0] +
                  C[1] * C[1] * mKwarpweftshear[1]);
}

/**
 * \brief       F = X * S * D_m^{-1}
 */
void cBaraffMaterial::CalcFAndN(const tMatrix3d &pos,
                                const tMatrix32d &uv_coords, tMatrix32d &F,
                                tMatrix32d &N) const
{
    // 1. calculate Dminv
    tMatrix2d DmInv = tMatrix2d::Zero();
    DmInv.col(0) = uv_coords.row(1) - uv_coords.row(0);
    DmInv.col(1) = uv_coords.row(2) - uv_coords.row(0);

    DmInv = DmInv.inverse().eval();

    // 2. calculate N = S * Dminv
    N.noalias() = mS * DmInv;

    // 3. calculate F = X * N
    F.noalias() = pos * N;
}

// void cBaraffMaterial::InitN()
// {
//     int num_of_tri = mObject->GetNumOfTriangles();
//     mNLst.reserve(num_of_tri);
// }
/**
 * \brief       C = [|F0| - 1, |F1| - 1]
 */
tVector2f cBaraffMaterial::CalcC(tMatrix32d &F)
{
    return tVector2f(F.col(0).norm() - 1, F.col(1).norm() - 1);
}

/**
 * \brief           gi = Ni \otimes ni
 */
tMatrix92d cBaraffMaterial::Calcg(const tMatrix32d &N, const tMatrix32d &n)
{
    tMatrix92d g = tMatrix92d::Zero(9, 2);
    for (size_t i = 0; i < 2; i++)
    {
        auto Ni = N.col(i);
        auto ni = n.col(i);
        for (int j = 0; j < 3; j++)
        {
            g.col(i).segment(3 * j, 3) = Ni[j] * ni;
        }
    }
    return g;
}

tMatrix3d VecToMat(const tVectorXd &vec)
{
    tMatrix3d mat;
    mat.col(0) = vec.segment(0, 3);
    mat.col(1) = vec.segment(3, 3);
    mat.col(2) = vec.segment(6, 3);
    return mat;
}
tVectorXd MatToVec(const tMatrix3d &mat)
{
    tVectorXd vec = tVectorXd::Zero(9);
    vec.segment(0, 3) = mat.col(0);
    vec.segment(3, 3) = mat.col(1);
    vec.segment(6, 3) = mat.col(2);
    return vec;
}
/**
 * \breif           check the derivative of energy
 */
void cBaraffMaterial::CheckStretchForce()
{
    tMatrix32d uv_coords = tMatrix32d::Random();
    tVectorXd pos_vec = tVectorXd::Random(9);

    double e_old = CalcStretchEnergy(VecToMat(pos_vec), uv_coords);
    double eps = 1e-3;
    tVectorXd f_num = tVectorXd::Zero(pos_vec.size());
    tVectorXd f_ana = CalcStretchForce(VecToMat(pos_vec), uv_coords);
    for (int i = 0; i < pos_vec.size(); i++)

    {
        pos_vec[i] += eps;
        double e_ = CalcStretchEnergy(VecToMat(pos_vec), uv_coords);
        f_num[i] = -(e_ - e_old) / eps;
        pos_vec[i] -= eps;
    }
    tVectorXd ana_diff = (f_num - f_ana) / mKwarpweftshear.segment(0, 2).mean();
    double ana_diff_norm = ana_diff.norm();
    // std::cout << "[stretch] f num = " << f_num.transpose() << std::endl;
    // std::cout << "[stretch] f ana = " << f_ana.transpose() << std::endl;
    std::cout << "[stretch] ana_diff_norm = " << ana_diff_norm << std::endl;
}

void cBaraffMaterial::CheckStretchStiffnessMatrix()
{
    tMatrix32d uv_coords = tMatrix32d::Random();
    tVectorXd pos_vec = tVectorXd::Random(9);

    tVectorXd force_old = CalcStretchForce(VecToMat(pos_vec), uv_coords);
    double eps = 1e-3;
    tMatrixXd K_num = tMatrixXd::Zero(pos_vec.size(), pos_vec.size());
    tMatrixXd K_ana = CalcStretchStiffMatrix(VecToMat(pos_vec), uv_coords);
    for (int i = 0; i < pos_vec.size(); i++)

    {
        pos_vec[i] += eps;
        tVectorXd force_ = CalcStretchForce(VecToMat(pos_vec), uv_coords);
        K_num.col(i) = (force_ - force_old) / eps;
        pos_vec[i] -= eps;
    }
    tMatrixXd K_diff = (K_num - K_ana);

    double K_diff_norm = K_diff.cwiseAbs().maxCoeff();
    // std::cout << "K num = " << K_num.transpose() << std::endl;
    // std::cout << "K ana = " << K_ana.transpose() << std::endl;
    // std::cout << "K diff = \n"
    //           << K_diff << std::endl;
    std::cout << "[stretch stiffness] K diff max coef = " << K_diff_norm
              << std::endl;
    // std::cout << "ana_diff_norm = " << ana_diff_norm << std::endl;
}

/**
 * \brief       K =
 *      k_i * (
 *          g_i * g_i^T
 *          - C_i / |X * N_i|
 *              * [ (N_i * N_i^T ) \otimes P_i  ]
 *          )
 */
tMatrix9d
cBaraffMaterial::CalcStretchStiffMatrix(const tMatrix3d &pos,
                                        const tMatrix32d &uv_coords) const
{
    tMatrix32d F, N;
    CalcFAndN(pos, uv_coords, F, N);

    tMatrix32d n = Calcn(pos, N);
    // 1. calculate gi
    tMatrixXd g = Calcg(N, n);

    // 2. calculate Ci, |X * N_i|
    tVector2f C = CalcC(F);
    tVector2d X_Ni_norm = (pos * N).colwise().norm();

    // 3. calculate P_i
    tMatrix3d P[2];
    CalcPi(n, P[0], P[1]);

    // 4. add
    tMatrix9d K = tMatrix9d::Zero();
    for (size_t i = 0; i < 2; i++)
    {
        K.noalias() +=
            mKwarpweftshear[i] *
            (-g.col(i) * g.col(i).transpose() -
             C[i] / X_Ni_norm[i] *

                 KroneckerProduct(N.col(i) * N.col(i).transpose(), P[i])

            );
    }
    return K;
}

tMatrix32d cBaraffMaterial::Calcn(const tMatrix3d &pos, tMatrix32d &N) const
{
    return (pos * N).colwise().normalized();
}

void cBaraffMaterial::CalcPi(const tMatrix32d &n, tMatrix3d &P0, tMatrix3d &P1)
{
    P0.noalias() = tMatrix3d::Identity() - n.col(0) * n.col(0).transpose();
    P1.noalias() = tMatrix3d::Identity() - n.col(1) * n.col(1).transpose();
}

void cBaraffMaterial::CheckShearingForce()
{
    tMatrix32d uv_coords = tMatrix32d::Random();
    tVectorXd pos_vec = tVectorXd::Random(9);

    double e_old = CalcShearingEnergy(VecToMat(pos_vec), uv_coords);
    double eps = 1e-5;
    tVectorXd f_num = tVectorXd::Zero(pos_vec.size());
    tVectorXd f_ana = CalcShearingForce(VecToMat(pos_vec), uv_coords);
    for (int i = 0; i < pos_vec.size(); i++)

    {
        pos_vec[i] += eps;
        double e_ = CalcShearingEnergy(VecToMat(pos_vec), uv_coords);
        f_num[i] = -(e_ - e_old) / eps;
        pos_vec[i] -= eps;
    }
    tVectorXd ana_diff = (f_num - f_ana) / mKwarpweftshear.segment(0, 2).mean();
    double ana_diff_norm = ana_diff.norm();
    std::cout << "[shear] f num = " << f_num.transpose() << std::endl;
    std::cout << "[shear] f ana = " << f_ana.transpose() << std::endl;
    std::cout << "[shear] ana_diff_norm = " << ana_diff_norm << std::endl;
}

void cBaraffMaterial::CheckShearingStiffnessMatrix()
{

    tMatrix32d uv_coords = tMatrix32d::Random();
    tVectorXd pos_vec = tVectorXd::Random(9);

    tVectorXd force_old = CalcShearingForce(VecToMat(pos_vec), uv_coords);
    double eps = 1e-3;
    tMatrixXd K_num = tMatrixXd::Zero(pos_vec.size(), pos_vec.size());
    tMatrixXd K_ana = CalcShearingStiffMatrix(VecToMat(pos_vec), uv_coords);
    for (int i = 0; i < pos_vec.size(); i++)

    {
        pos_vec[i] += eps;
        tVectorXd force_ = CalcShearingForce(VecToMat(pos_vec), uv_coords);
        K_num.col(i) = (force_ - force_old) / eps;
        pos_vec[i] -= eps;
    }
    tMatrixXd K_diff = (K_num - K_ana);

    double K_diff_norm = K_diff.cwiseAbs().maxCoeff();
    // std::cout << "K num = " << K_num.transpose() << std::endl;
    // std::cout << "K ana = " << K_ana.transpose() << std::endl;
    // std::cout << "K diff = \n"
    //           << K_diff << std::endl;
    std::cout << "[stretch stiffness] K diff norm = " << K_diff_norm
              << std::endl;
}

// void cBaraffMaterial::SetSheaingK(double Kshear)
// {
//     mKwarpweftshear[2] = Kshear;
// }

/**
 * \biref               f = -k_i C_i' g_i'
 */
tVector9d cBaraffMaterial::CalcShearingForce(const tMatrix3d &pos,
                                             const tMatrix32d &uv_coords) const
{
    // 1. calculate F‘ and N’
    tMatrix32d F_prime, N_prime;
    CalcFAndN_shearing(pos, uv_coords, F_prime, N_prime);
    // 2. calculate C_prime
    tVector2f C_prime = CalcC(F_prime);
    tMatrix32d n_prime = Calcn(pos, N_prime);
    tMatrixXd g_prime = Calcg(N_prime, n_prime);
    // tVector9d int_force = (-mKwarpweftshear[0] * C[0]) * g.col(0) +
    // (-mKwarpweftshear[1] * C[1]) * g.col(1); std::cout << "F' = \n"
    //   << F_prime << std::endl;
    // std::cout << "C' = \n"
    //   << C_prime << std::endl;
    // std::cout << "n' = \n"
    //   << n_prime << std::endl;
    // std::cout << "g' = \n"
    //   << g_prime << std::endl;
    tVector9d int_force = -mKwarpweftshear[2] * (C_prime[0] * g_prime.col(0) +
                                                 C_prime[1] * g_prime.col(1));
    return int_force;
}

tMatrix9d
cBaraffMaterial::CalcShearingStiffMatrix(const tMatrix3d &pos,
                                         const tMatrix32d &uv_coords) const
{
    tMatrix32d F, N;
    CalcFAndN_shearing(pos, uv_coords, F, N);

    tMatrix32d n = Calcn(pos, N);
    // 1. calculate gi
    tMatrixXd g = Calcg(N, n);

    // 2. calculate Ci, |X * N_i|
    tVector2f C = CalcC(F);
    tVector2d X_Ni_norm = (pos * N).colwise().norm();

    // 3. calculate P_i
    tMatrix3d P[2];
    CalcPi(n, P[0], P[1]);

    // 4. add
    tMatrix9d K = tMatrix9d::Zero();
    for (size_t i = 0; i < 2; i++)
    {
        K.noalias() +=
            mKwarpweftshear[i] *
            (-g.col(i) * g.col(i).transpose() -
             C[i] / X_Ni_norm[i] *

                 KroneckerProduct(N.col(i) * N.col(i).transpose(), P[i])

            );
    }
    return K;
}

void cBaraffMaterial::CalcFAndN_shearing(const tMatrix3d &pos,
                                         const tMatrix32d &uv_coords,
                                         tMatrix32d &F_prime,
                                         tMatrix32d &N_prime) const
{
    tMatrix2d DmInv = tMatrix2d::Zero();
    DmInv.col(0) = uv_coords.row(1) - uv_coords.row(0);
    DmInv.col(1) = uv_coords.row(2) - uv_coords.row(0);

    DmInv = DmInv.inverse().eval();

    tMatrix2d R_neg45 = cRotUtil::RotMat2D(-M_PI / 4);
    // 2. calculate N' = S * Dminv * R^{-45}
    N_prime.noalias() = mS * DmInv * R_neg45;

    // 3. calculate F' = X * N'
    F_prime.noalias() = pos * N_prime;
}

double cBaraffMaterial::CalcShearingEnergy(const tMatrix3d &pos,
                                           const tMatrix32d &uv_coords) const
{
    tMatrix32d F, N;
    CalcFAndN_shearing(pos, uv_coords, F, N);

    // 2. calculate C
    tVector2f C = CalcC(F);
    return 0.5 * (C[0] * C[0] * mKwarpweftshear[2] +
                  C[1] * C[1] * mKwarpweftshear[2]);
}

double cBaraffMaterial::GetEnergy(int tri_id) const { return mELst[tri_id]; }
tVector9d cBaraffMaterial::GetForce(int tri_id) const
{
    return mIntForceLst[tri_id];
};
tMatrix9d cBaraffMaterial::GetStiffMatrix(int tri_id) const
{
    return mKLst[tri_id];
}

void cBaraffMaterial::CheckForce()
{
    // 1. calculate old energy and force
    // for (int i = 0; i < 3 * mNumOfVertices; i++)
    // {
    //     mObject->GetVertexArrayRef()[int(i / 3)]->mPos[i % 3] +=
    //     cMathUtil::RandDouble(-0.01, 0.01);
    // }
    this->Update(true, true, false);
    double old_e = CalcTotalEnergy();
    tVectorXd force_ana = CalcTotalForce();
    tVectorXd force_num = tVectorXd::Zero(force_ana.size());

    // 2. update pos, update energy, get energy, get num ana
    double eps = 1e-5;
    for (int i = 0; i < 3 * mNumOfVertices; i++)
    {
        mObject->GetVertexArrayRef()[int(i / 3)]->mPos[i % 3] += eps;
        Update(true, false, false);
        force_num[i] = -(CalcTotalEnergy() - old_e) / eps;
        mObject->GetVertexArrayRef()[int(i / 3)]->mPos[i % 3] -= eps;
    }
    Update(true, true, false);

    tVectorXd force_diff = force_ana - force_num;
    std::cout << "force ana = " << force_ana.segment(0, 12).transpose()
              << std::endl;
    std::cout << "force num = " << force_num.segment(0, 12).transpose()
              << std::endl;
    std::cout << "force diff = " << force_diff.cwiseAbs().maxCoeff()
              << std::endl;
    // exit(1);
}
void cBaraffMaterial::CheckStiffnessMatrix()
{
    //  for (int i = 0; i < 3 * mNumOfVertices; i++)
    // {
    //     mObject->GetVertexArrayRef()[int(i / 3)]->mPos[i % 3] +=
    //     cMathUtil::RandDouble(-0.01, 0.01);
    // }
    this->Update(true, true, true);
    tVectorXd old_f = CalcTotalForce();
    tMatrixXd K_ana = CalcTotalStiffnessMatrix().toDense();
    tMatrixXd K_num = tMatrixXd::Zero(old_f.size(), old_f.size());

    // 2. update pos, update energy, get energy, get num ana
    double eps = 1e-5;
    for (int i = 0; i < 3 * mNumOfVertices; i++)
    {
        mObject->GetVertexArrayRef()[int(i / 3)]->mPos[i % 3] += eps;
        Update(true, true, false);
        K_num.col(i) = (CalcTotalForce() - old_f) / eps;
        mObject->GetVertexArrayRef()[int(i / 3)]->mPos[i % 3] -= eps;
    }
    Update(true, true, true);

    tMatrixXd K_diff = K_ana - K_num;
    // std::cout << "K ana = " << K_ana.transpose() << std::endl;
    // std::cout << "K num = " << K_num.transpose() << std::endl;
    std::cout << "K diff = " << K_diff.cwiseAbs().maxCoeff() << std::endl;
    // exit(1);
}

void cBaraffMaterial::SetK(const tVector3d &Kwarpweftbias) {}
double cBaraffMaterial::CalcTotalEnergy() const
{
    double sum = 0;
    for (auto &e : mELst)
        sum += e;
    return sum;
}
tVectorXd cBaraffMaterial::CalcTotalForce() const
{
    tVectorXd total_force = tVectorXd::Zero(3 * mNumOfVertices);
    auto &tri_lst = mObject->GetTriangleArray();
    for (int i = 0; i < mNumOfTriangles; i++)
    {
        auto t = tri_lst[i];
        auto &ele_f = mIntForceLst[i];
        total_force.segment(3 * t->mId0, 3) += ele_f.segment(0, 3);
        total_force.segment(3 * t->mId1, 3) += ele_f.segment(3, 3);
        total_force.segment(3 * t->mId2, 3) += ele_f.segment(6, 3);
    }
    // std::cout << "[mat] stretch total force max = " << total_force.cwiseAbs().maxCoeff() << std::endl;
    return total_force;
}
#include "utils/DefUtil.h"
tSparseMatd cBaraffMaterial::CalcTotalStiffnessMatrix()
{
#define num_divide 15
    auto &tri_lst = mObject->GetTriangleArray();
    global_K_buf.setZero();
    total_triplets_buf.clear();
    std::vector<tTriplet> sub_triples_set[num_divide];
    for (int i = 0; i < num_divide; i++)
    {
        sub_triples_set[i].reserve(3 * 3 * 9);
    }
    // #pragma omp parallel for num_threads(num_divide)
    for (int t = 0; t < mNumOfTriangles; t++)
    {
        int thread_num = omp_get_thread_num();
        auto &sub_triples = sub_triples_set[thread_num];
        sub_triples.clear();
        auto cur_tri = tri_lst[t];
        int v_id[3] = {cur_tri->mId0, cur_tri->mId1, cur_tri->mId2};
        for (size_t i = 0; i < 3; i++)
        {
            size_t global_vi_idx = v_id[i];
            for (size_t j = 0; j < 3; j++)
            {
                size_t global_vj_idx = v_id[j];

                tMatrix3d ele_K_part = mKLst[t].block(3 * i, 3 * j, 3, 3);
                size_t i3 = 3 * global_vi_idx, j3 = 3 * global_vj_idx;
                sub_triples.emplace_back(i3, j3, ele_K_part(0, 0));
                sub_triples.emplace_back(i3, j3 + 1, ele_K_part(0, 1));
                sub_triples.emplace_back(i3, j3 + 2, ele_K_part(0, 2));

                sub_triples.emplace_back(i3 + 1, j3, ele_K_part(1, 0));
                sub_triples.emplace_back(i3 + 1, j3 + 1, ele_K_part(1, 1));
                sub_triples.emplace_back(i3 + 1, j3 + 2, ele_K_part(1, 2));

                sub_triples.emplace_back(i3 + 2, j3, ele_K_part(2, 0));
                sub_triples.emplace_back(i3 + 2, j3 + 1, ele_K_part(2, 1));
                sub_triples.emplace_back(i3 + 2, j3 + 2, ele_K_part(2, 2));
            }
        }
        // #pragma omp critical
        total_triplets_buf.insert(total_triplets_buf.end(), sub_triples.begin(),
                                  sub_triples.end());
    }
    global_K_buf.setFromTriplets(total_triplets_buf.begin(),
                                 total_triplets_buf.end());
    return global_K_buf;
}

template <typename S, typename T>
void Fetch(std::vector<S> &dest, tEigenArr<T> &ori)
{
    dest.resize(ori.size());
    for (size_t i = 0; i < ori.size(); i++)
    {
        dest[i].noalias() = ori[i].cast<float>();
    }
}

void cBaraffMaterial::GetFLst(std::vector<tMatrix32f> &vec)
{
    Fetch(vec, this->mFLst);
}
void cBaraffMaterial::GetFprimeLst(std::vector<tMatrix32f> &vec)
{
    Fetch(vec, this->mFprimeLst);
}

void cBaraffMaterial::GetnLst(std::vector<tMatrix32f> &vec)
{
    Fetch(vec, this->mnLst);
}
void cBaraffMaterial::GetgLst(std::vector<tMatrix92f> &vec)
{
    Fetch(vec, this->mgLst);
}
void cBaraffMaterial::GetCLst(std::vector<tVector2f> &vec)
{
    Fetch(vec, this->mCLst);
}

void cBaraffMaterial::GetnprimeLst(std::vector<tMatrix32f> &vec)
{
    Fetch(vec, this->mnprimeLst);
}
void cBaraffMaterial::GetgprimeLst(std::vector<tMatrix92f> &vec)
{
    Fetch(vec, this->mgprimeLst);
}
void cBaraffMaterial::GetCprimeLst(std::vector<tVector2f> &vec)
{
    Fetch(vec, this->mCprimeLst);
}

void cBaraffMaterial::GetEleKLst(std::vector<tMatrix9f> &mat)
{
    Fetch(mat, this->mKLst);
}

void cBaraffMaterial::GetEleFintLst(std::vector<tVector9f> &mat)
{
    Fetch(mat, this->mIntForceLst);
}