#include "BaraffMaterial.h"
#include <iostream>
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

/**
 * \brief           f = -dE/dx = - k_i * C_i * g_i
 * where g_i = N_i \otimes n_i
 * For more details, please check the note "Baraff 内力和刚度矩阵计算(利用Kronecker简化).md"
*/
tVector9d cBaraffMaterial::CalcStretchForce(const tMatrix3d &pos, const tMatrix32d &rest_texture_coords) const
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
    // std::cout << "C = " << C.transpose() << std::endl;
    // std::cout << "g = \n" << g.transpose() << std::endl;
    tVector9d int_force = (-mKwarpweftshear[0] * C[0]) * g.col(0) + (-mKwarpweftshear[1] * C[1]) * g.col(1);
    // std::cout << "force = " << int_force.transpose() << std::endl;
    return int_force;
}

void cBaraffMaterial::SetStretchK(double Kwarp, double Kweft)
{
    this->mKwarpweftshear[0] = Kwarp;
    this->mKwarpweftshear[1] = Kweft;
}
// tMatrix9d cBaraffMaterial::CalcStiffMatrix(const tMatrix3d &pos, const tMatrix32d &rest_texture_coords) const
// {
//     return tMatrix9d::Zero();
// }

/**
 * \brief           E = 1/2 * C^T * K * C
 *      let F = [Wu, Wv] \in R^{3 \times 2}
 *      C is the "baraff stretch condition" = [|Wu| - Bu, |Wv| - Bv]
 *      K is the stiffness mat = diag(Kwarp, Kweft)
*/
double cBaraffMaterial::CalcStretchEnergy(const tMatrix3d &pos, const tMatrix32d &rest_texture_coords) const
{
    // 1. calculate F and N
    tMatrix32d F, N;
    CalcFAndN(pos, rest_texture_coords, F, N);

    // 2. calculate C
    tVector2f C = CalcC(F);
    return 0.5 * (C[0] * C[0] * mKwarpweftshear[0] + C[1] * C[1] * mKwarpweftshear[1]);
}

/**
 * \brief       F = X * S * D_m^{-1}
*/
void cBaraffMaterial::CalcFAndN(const tMatrix3d &pos, const tMatrix32d &uv_coords, tMatrix32d &F, tMatrix32d &N) const
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

/**
 * \brief       C = [|F0| - 1, |F1| - 1]
*/
tVector2f cBaraffMaterial::CalcC(tMatrix32d &F) const
{
    return tVector2f(F.col(0).norm() - 1, F.col(1).norm() - 1);
}

/**
 * \brief           gi = Ni \otimes ni
*/
tMatrixXd cBaraffMaterial::Calcg(const tMatrix32d &N, const tMatrix32d &n) const
{
    tMatrixXd g = tMatrixXd::Zero(9, 2);
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

extern tMatrix3d VecToMat(const tVectorXd &vec);
// {
//     tMatrix3d mat;
//     mat.col(0) = vec.segment(0, 3);
//     mat.col(1) = vec.segment(3, 3);
//     mat.col(2) = vec.segment(6, 3);
//     return mat;
// }
extern tVectorXd MatToVec(const tMatrix3d &mat);
// {
//     tVectorXd vec = tVectorXd::Zero(9);
//     vec.segment(0, 3) = mat.col(0);
//     vec.segment(3, 3) = mat.col(1);
//     vec.segment(6, 3) = mat.col(2);
//     return vec;
// }
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
        double e_new = CalcStretchEnergy(VecToMat(pos_vec), uv_coords);
        f_num[i] = -(e_new - e_old) / eps;
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
        tVectorXd force_new = CalcStretchForce(VecToMat(pos_vec), uv_coords);
        K_num.col(i) = (force_new - force_old) / eps;
        pos_vec[i] -= eps;
    }
    tMatrixXd K_diff = (K_num - K_ana);

    double K_diff_norm = K_diff.cwiseAbs().maxCoeff();
    // std::cout << "K num = " << K_num.transpose() << std::endl;
    // std::cout << "K ana = " << K_ana.transpose() << std::endl;
    // std::cout << "K diff = \n"
    //           << K_diff << std::endl;
    std::cout << "[stretch stiffness] K diff max coef = " << K_diff_norm << std::endl;
    // std::cout << "ana_diff_norm = " << ana_diff_norm << std::endl;
}

extern tMatrix9d KroneckerProduct(const tMatrix3d &a, const tMatrix3d &b);
// {
//     tMatrix9d res = tMatrix9d::Zero();
//     for (int i = 0; i < 3; i++)
//         for (int j = 0; j < 3; j++)
//         {
//             res.block(3 * i, j * 3, 3, 3) = a(i, j) * b;
//         }
//     return res;
// }

/**
 * \brief       K = 
 *      k_i * (
 *          g_i * g_i^T 
 *          - C_i / |X * N_i| 
 *              * [ (N_i * N_i^T ) \otimes P_i  ]
 *          )
*/
tMatrix9d cBaraffMaterial::CalcStretchStiffMatrix(const tMatrix3d &pos, const tMatrix32d &uv_coords) const
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
        K.noalias() += mKwarpweftshear[i] * (-g.col(i) * g.col(i).transpose() - C[i] / X_Ni_norm[i] *

                                                                                    KroneckerProduct(N.col(i) * N.col(i).transpose(), P[i])

                                            );
    }
    return K;
}

tMatrix32d cBaraffMaterial::Calcn(const tMatrix3d &pos, tMatrix32d &N) const
{

    return (pos * N).colwise().normalized();
}

void cBaraffMaterial::CalcPi(const tMatrix32d &n, tMatrix3d &P0, tMatrix3d &P1) const
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
        double e_new = CalcShearingEnergy(VecToMat(pos_vec), uv_coords);
        f_num[i] = -(e_new - e_old) / eps;
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
        tVectorXd force_new = CalcShearingForce(VecToMat(pos_vec), uv_coords);
        K_num.col(i) = (force_new - force_old) / eps;
        pos_vec[i] -= eps;
    }
    tMatrixXd K_diff = (K_num - K_ana);

    double K_diff_norm = K_diff.cwiseAbs().maxCoeff();
    // std::cout << "K num = " << K_num.transpose() << std::endl;
    // std::cout << "K ana = " << K_ana.transpose() << std::endl;
    // std::cout << "K diff = \n"
    //           << K_diff << std::endl;
    std::cout << "[stretch stiffness] K diff norm = " << K_diff_norm << std::endl;
}

void cBaraffMaterial::SetSheaingK(double Kshear)
{
    mKwarpweftshear[2] = Kshear;
}

/**
 * \biref               f = -k_i C_i' g_i'
*/
tVector9d cBaraffMaterial::CalcShearingForce(const tMatrix3d &pos, const tMatrix32d &uv_coords) const
{
    // 1. calculate F‘ and N’
    tMatrix32d F_prime, N_prime;
    CalcFAndN_shearing(pos, uv_coords, F_prime, N_prime);
    // 2. calculate C_prime
    tVector2f C_prime = CalcC(F_prime);
    tMatrix32d n_prime = Calcn(pos, N_prime);
    tMatrixXd g_prime = Calcg(N_prime, n_prime);
    // tVector9d int_force = (-mKwarpweftshear[0] * C[0]) * g.col(0) + (-mKwarpweftshear[1] * C[1]) * g.col(1);
    // std::cout << "F' = \n"
    //   << F_prime << std::endl;
    // std::cout << "C' = \n"
    //   << C_prime << std::endl;
    // std::cout << "n' = \n"
    //   << n_prime << std::endl;
    // std::cout << "g' = \n"
    //   << g_prime << std::endl;
    // std::cout << "Cprime = " << C_prime.transpose() << std::endl;
    // std::cout << "gprime = \n" << g_prime.transpose() << std::endl;
    tVector9d int_force = -mKwarpweftshear[2] * (C_prime[0] * g_prime.col(0) + C_prime[1] * g_prime.col(1));
    return int_force;
}

tMatrix9d cBaraffMaterial::CalcShearingStiffMatrix(const tMatrix3d &pos, const tMatrix32d &uv_coords) const
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
        K.noalias() += mKwarpweftshear[i] * (-g.col(i) * g.col(i).transpose() - C[i] / X_Ni_norm[i] *

                                                                                    KroneckerProduct(N.col(i) * N.col(i).transpose(), P[i])

                                            );
    }
    return K;
}

void cBaraffMaterial::CalcFAndN_shearing(const tMatrix3d &pos, const tMatrix32d &uv_coords, tMatrix32d &F_prime, tMatrix32d &N_prime) const
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

double cBaraffMaterial::CalcShearingEnergy(const tMatrix3d &pos, const tMatrix32d &uv_coords) const
{
    tMatrix32d F, N;
    CalcFAndN_shearing(pos, uv_coords, F, N);

    // 2. calculate C
    tVector2f C = CalcC(F);
    return 0.5 * (C[0] * C[0] * mKwarpweftshear[2] + C[1] * C[1] * mKwarpweftshear[2]);
}

double cBaraffMaterial::CalcEnergy(const tMatrix3d &pos, const tMatrix32d &uv_coords) const
{
    return this->CalcStretchEnergy(pos, uv_coords) +
           this->CalcShearingEnergy(pos, uv_coords);
}
tVector9d cBaraffMaterial::CalcForce(const tMatrix3d &pos, const tMatrix32d &uv_coords) const
{
    return this->CalcStretchForce(pos, uv_coords) +
           this->CalcShearingForce(pos, uv_coords);
};
tMatrix9d cBaraffMaterial::CalcStiffMatrix(const tMatrix3d &pos, const tMatrix32d &uv_coords) const
{
    return this->CalcStretchStiffMatrix(pos, uv_coords) +
           this->CalcShearingStiffMatrix(pos, uv_coords);
}

void cBaraffMaterial::CheckForce()
{
    CheckStretchForce();
    CheckShearingForce();
}
void cBaraffMaterial::CheckStiffnessMatrix()
{
    CheckStretchStiffnessMatrix();
    CheckShearingStiffnessMatrix();
}