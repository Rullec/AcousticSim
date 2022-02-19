#include "NeoHookeanMaterial.h"
#include "utils/JsonUtil.h"
#include <iostream>
extern cFourOrderTensor CalcDFDF(const tMatrix3d &F);
float CalcJ(const tMatrix3d &F);
cFourOrderTensor CalcDFinvTDF(const tMatrix3d &F);
void CheckDFinvTDF(const tMatrix3d &F);

cFourOrderTensor CalcDJDF(const tMatrix3d &F);
void CheckDJDF(const tMatrix3d &F);

typedef std::function<tMatrixXd(const tMatrixXd &)> tForwardCalcFunc;

cFourOrderTensor CalcNumericalDerivaitve_MatrixWRTMatrix(tForwardCalcFunc func, const tMatrixXd &X_raw, const tVector2i output_dims, float eps = 1e-5)
{
    tMatrixXd X = X_raw;
    tMatrixXd old_value = func(X);
    cFourOrderTensor res(output_dims[0], output_dims[1], X.rows(), X.cols());
    tVector4i derive_dims = res.GetShape();
    SIM_ASSERT(old_value.hasNaN() == false);
    for (size_t k = 0; k < X.rows(); k++)
        for (size_t l = 0; l < X.cols(); l++)
        {
            X(k, l) += eps;
            tMatrixXd new_value = func(X);
            SIM_ASSERT(new_value.hasNaN() == false);
            tMatrixXd slice = (new_value - old_value) / eps;
            // std::cout << "slice = \n"
            //           << slice << std::endl;
            // std::cout << "new value = \n"
            //           << new_value << std::endl;
            res.SetLastTwoComp(k, l, slice);
            X(k, l) -= eps;
        }
    return res;
}

void cNeoHookeanMaterial::Init(const Json::Value &conf)
{
    mMu = cJsonUtil::ParseAsDouble("youngs", conf);
    mLambda = cJsonUtil::ParseAsDouble("poisson_ratio", conf);
    SIM_INFO("parse mu {} lambda {}", mMu, mLambda);
}

tMatrix3d cNeoHookeanMaterial::CalcP(const tMatrix3d &F) const
{
    // tMatrix3d FinvT = F.transpose().inverse();
    // double J = CalcJ(F);
    return CalcP_part1(F) + CalcP_part2(F);
} // PK1

tMatrix3d cNeoHookeanMaterial::CalcP_part1(const tMatrix3d &F) const
{
    tMatrix3d FinvT = F.transpose().inverse();
    return mMu * (F - FinvT);
}
tMatrix3d cNeoHookeanMaterial::CalcP_part2(const tMatrix3d &F) const
{
    tMatrix3d FinvT = F.transpose().inverse();
    double J = CalcJ(F);
    return mLambda * std::log(J) * FinvT;
}
cFourOrderTensor cNeoHookeanMaterial::CalcDPDF(const tMatrix3d &F) const
{
    return CalcDPDF_part1(F) + CalcDPDF_part2(F);
}

cFourOrderTensor cNeoHookeanMaterial::CalcDPDF_part1(const tMatrix3d &F) const
{
    // mu * (Q - R)
    cFourOrderTensor DFDF_Q = CalcDFDF(F);
    cFourOrderTensor DFinvTDF_R = CalcDFinvTDF(F);
    return (DFDF_Q - DFinvTDF_R) * mMu;
}
cFourOrderTensor cNeoHookeanMaterial::CalcDPDF_part2(const tMatrix3d &F) const
{
    float J = CalcJ(F);
    cFourOrderTensor DJIDF = CalcDJDF(F);
    DJIDF.TensorMatrixProductWithoutCopy(0, 1, F.inverse().transpose());
    DJIDF = (DJIDF * (1.0 / J) + CalcDFinvTDF(F) * std::log(J)) * mLambda;
    return DJIDF;
}

void cNeoHookeanMaterial::CheckDPDF(const tMatrix3d &F) const
{
    // CheckDPDF_part1(F);
    // CheckDPDF_part2(F);
    // CheckDJDF(F);
    // exit(1);
    tMatrixXd deriv_ana_mat = CalcDPDF(F).ExpandToMatrix();
    tMatrixXd X = F;
    tMatrixXd old_value = CalcP(X);
    double eps = 1e-5;
    cFourOrderTensor deriv_num(3, 3, X.rows(), X.cols());
    tVector4i derive_dims = deriv_num.GetShape();
    std::cout << "old_value = \n"
              << old_value << std::endl;
    SIM_ASSERT(old_value.hasNaN() == false);
    for (size_t k = 0; k < X.rows(); k++)
        for (size_t l = 0; l < X.cols(); l++)
        {
            X(k, l) += eps;
            tMatrixXd new_value = CalcP(X);
            SIM_ASSERT(new_value.hasNaN() == false);
            tMatrixXd slice = (new_value - old_value) / eps;
            // std::cout << "slice = \n"
            //           << slice << std::endl;
            // std::cout << "new value = \n"
            //           << new_value << std::endl;
            deriv_num.SetLastTwoComp(k, l, slice);
            X(k, l) -= eps;
        }
    tMatrixXd deriv_num_mat = deriv_num.ExpandToMatrix();
    tMatrixXd deriv_diff = deriv_num_mat - deriv_ana_mat;
    double deriv_diff_norm = deriv_diff.norm();
    std::cout << "deriv_diff_norm = " << deriv_diff_norm << std::endl;
    std::cout << "deriv_num = " << deriv_num_mat << std::endl;
    std::cout << "deriv_ana = " << deriv_ana_mat << std::endl;
    std::cout << "deriv_diff = " << deriv_diff << std::endl;
    // return deriv_num;
}

cFourOrderTensor CalcDFinvTDF(const tMatrix3d &F)
{
    tMatrix3d Finv = F.inverse();
    cFourOrderTensor DFinvTDF(3, 3, 3, 3);
    for (size_t i = 0; i < 3; i++)
        for (size_t j = 0; j < 3; j++)
        {
            tMatrixXd &comp = DFinvTDF(i, j);
            for (size_t k = 0; k < 3; k++)
                for (size_t l = 0; l < 3; l++)
                {
                    comp(k, l) = -Finv(j, k) * Finv(l, i);
                }
        }
    return DFinvTDF;
}
void CheckDFinvTDF(const tMatrix3d &F)
{
    const cFourOrderTensor &deriv_ana = CalcDFinvTDF(F);
    auto CalcFinvT = [](const tMatrixXd &F)
    {
        return F.inverse().transpose();
    };
    tMatrixXd deriv_ana_mat = deriv_ana.ExpandToMatrix();
    tMatrixXd deriv_num_mat = CalcNumericalDerivaitve_MatrixWRTMatrix(CalcFinvT, F, tVector2i(3, 3)).ExpandToMatrix();
    tMatrixXd diff = deriv_num_mat - deriv_ana_mat;
    double diff_norm = diff.norm();
    std::cout << "CheckDFinvTDF, diff norm = " << diff_norm << std::endl;

    // std::cout << "deriv ana = \n"
    //           << deriv_ana_mat << std::endl;
    // std::cout << "deriv num = \n"
    //           << deriv_num_mat << std::endl;
}

float CalcJ(const tMatrix3d &F)
{
    return std::fabs(F.determinant());
}
cFourOrderTensor CalcDJDF(const tMatrix3d &F)
{
    tMatrix3d comp = (F.determinant() > 0 ? 1 : -1) * F.determinant() * F.inverse().transpose();
    cFourOrderTensor res(3, 3, 3, 3);
    res(0, 0) = comp;
    res(1, 1) = comp;
    res(2, 2) = comp;
    return res;
}
void CheckDJDF(const tMatrix3d &F)
{
    auto deriv_ana_mat = CalcDJDF(F).ExpandToMatrix();
    auto calc_JI = [](const tMatrixXd &F)
    {
        return CalcJ(F) * tMatrix3d::Identity();
    };

    auto deriv_num_mat = CalcNumericalDerivaitve_MatrixWRTMatrix(calc_JI, F, tVector2i(3, 3), 1e-5).ExpandToMatrix();
    tMatrixXd deriv_diff = deriv_num_mat - deriv_ana_mat;
    double deriv_diff_norm = deriv_diff.norm();
    std::cout << "DJDF ana = \n"
              << deriv_ana_mat << std::endl;
    std::cout << "DJDF num = \n"
              << deriv_num_mat << std::endl;
    std::cout << "DJDF diff = \n"
              << deriv_diff << std::endl;
    std::cout << "DJDF diff norm = " << deriv_diff_norm << std::endl;
}

void cNeoHookeanMaterial::CheckDPDF_part1(const tMatrix3d &F) const
{
    tMatrixXd deriv_ana_mat = CalcDPDF_part1(F).ExpandToMatrix();
    tMatrixXd X = F;
    tMatrixXd old_value = CalcP_part1(X);
    double eps = 1e-5;
    cFourOrderTensor deriv_num(3, 3, X.rows(), X.cols());
    tVector4i derive_dims = deriv_num.GetShape();
    std::cout << "old_value = \n"
              << old_value << std::endl;
    SIM_ASSERT(old_value.hasNaN() == false);
    for (size_t k = 0; k < X.rows(); k++)
        for (size_t l = 0; l < X.cols(); l++)
        {
            X(k, l) += eps;
            tMatrixXd new_value = CalcP_part1(X);
            SIM_ASSERT(new_value.hasNaN() == false);
            tMatrixXd slice = (new_value - old_value) / eps;
            // std::cout << "slice = \n"
            //           << slice << std::endl;
            // std::cout << "new value = \n"
            //           << new_value << std::endl;
            deriv_num.SetLastTwoComp(k, l, slice);
            X(k, l) -= eps;
        }
    tMatrixXd deriv_num_mat = deriv_num.ExpandToMatrix();
    tMatrixXd deriv_diff = deriv_num_mat - deriv_ana_mat;
    double deriv_diff_norm = deriv_diff.norm();
    std::cout << "[part1] deriv_diff_norm = " << deriv_diff_norm << std::endl;
    // std::cout << "[part1] deriv_num = " << deriv_num_mat << std::endl;
    // std::cout << "[part1] deriv_ana = " << deriv_ana_mat << std::endl;
    // std::cout << "[part1] deriv_diff = " << deriv_diff << std::endl;
}
void cNeoHookeanMaterial::CheckDPDF_part2(const tMatrix3d &F) const
{
    tMatrixXd deriv_ana_mat = CalcDPDF_part2(F).ExpandToMatrix();
    tMatrixXd X = F;
    tMatrixXd old_value = CalcP_part2(X);
    double eps = 1e-4;
    cFourOrderTensor deriv_num(3, 3, X.rows(), X.cols());
    tVector4i derive_dims = deriv_num.GetShape();
    std::cout << "old_value = \n"
              << old_value << std::endl;
    SIM_ASSERT(old_value.hasNaN() == false);
    for (size_t k = 0; k < X.rows(); k++)
        for (size_t l = 0; l < X.cols(); l++)
        {
            X(k, l) += eps;
            tMatrixXd new_value = CalcP_part2(X);
            SIM_ASSERT(new_value.hasNaN() == false);
            tMatrixXd slice = (new_value - old_value) / eps;
            // std::cout << "slice = \n"
            //           << slice << std::endl;
            // std::cout << "new value = \n"
            //           << new_value << std::endl;
            deriv_num.SetLastTwoComp(k, l, slice);
            X(k, l) -= eps;
        }
    tMatrixXd deriv_num_mat = deriv_num.ExpandToMatrix();
    tMatrixXd deriv_diff = deriv_num_mat - deriv_ana_mat;
    double deriv_diff_norm = deriv_diff.norm();
    std::cout << "[part2] deriv_ana_mat = " << deriv_ana_mat << std::endl;
    std::cout << "[part2] deriv_num_mat = " << deriv_num_mat << std::endl;
    std::cout << "[part2] deriv_diff = " << deriv_diff << std::endl;
    std::cout << "[part2] deriv_diff_norm = " << deriv_diff_norm << std::endl;
}