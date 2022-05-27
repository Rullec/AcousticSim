#include "StvkMaterial.h"
#include "utils/JsonUtil.h"

tMatrix3d CalcGreenStrain(const tMatrix3d &F)
{
    tMatrix3d I = tMatrix3d::Identity();
    return 0.5 * (F.transpose() * F - I);
}

/**
 * \brief           calc dE/dF
 *      For more details, please check note "矩阵对矩阵求导" 中二次型的情况
 */
cFourOrderTensor CalcDEDF(const tMatrix3d &F)
{
    tMatrix3d FT = F.transpose();
    cFourOrderTensor deriv(3, 3, 3, 3);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
                for (int l = 0; l < 3; l++)
                {
                    if (l == j)
                    {
                        deriv(k, l)(i, j) += 0.5 * FT(k, i);
                    }
                    if (k == j)
                    {
                        deriv(k, l)(i, j) += 0.5 * F(i, l);
                    }
                }
    return deriv;
}

cStvkMaterial::cStvkMaterial() : cBaseMaterial(eMaterialType::STVK) {}
void cStvkMaterial::Init(const Json::Value &conf) { cBaseMaterial::Init(conf); }
tMatrix3d cStvkMaterial::CalcP(const tMatrix3d &F) const
{
    return CalcP_part1(F) + CalcP_part2(F);
} // PK1

/**
 * \brief               calculate dPdF part 1
 *      For more details, please check note "FEM Course 第三部分 离散化
 * 刚度矩阵计算"
 */
cFourOrderTensor cStvkMaterial::CalcDPDF(const tMatrix3d &F) const
{
    return CalcDPDF_part1(F) + CalcDPDF_part2(F);
}

void cStvkMaterial::CheckDPDF(const tMatrix3d &F) const {}
extern cFourOrderTensor CalcDFDF(const tMatrix3d &F);
extern cFourOrderTensor CalcDTrEIDF(const tMatrix3d &F);

/**
 * \brief               calculate dPdF part 1
 *      For more details, please check note "FEM Course 第三部分 离散化
 * 刚度矩阵计算"
 */
cFourOrderTensor cStvkMaterial::CalcDPDF_part1(const tMatrix3d &F) const
{
    cFourOrderTensor DFDF = CalcDFDF(F);
    tMatrix3d E = CalcGreenStrain(F);
    cFourOrderTensor DEDF = CalcDEDF(F);
    DFDF.TensorMatrixProductWithoutCopy(0, 1, E);
    DEDF.MatrixTensorProductWithoutCopy(0, 1, F);
    cFourOrderTensor deriv = (DFDF + DEDF) * 2 * mMu;
    return deriv;
}

cFourOrderTensor cStvkMaterial::CalcDPDF_part2(const tMatrix3d &F) const
{
    tMatrix3d E = CalcGreenStrain(F);
    cFourOrderTensor DTrEIDF = CalcDTrEIDF(F),
                     trE_DFDF = CalcDFDF(F) * E.trace();
    DTrEIDF.TensorMatrixProductWithoutCopy(0, 1, F);
    return (DTrEIDF + trE_DFDF) * mLambda;
}
void cStvkMaterial::CheckDPDF_part1(const tMatrix3d &F) const {}

void cStvkMaterial::CheckDPDF_part2(const tMatrix3d &F) const {}
tMatrix3d cStvkMaterial::CalcP_part1(const tMatrix3d &F) const
{

    tMatrix3d E = CalcGreenStrain(F);
    return 2 * mMu * F * E;
}
tMatrix3d cStvkMaterial::CalcP_part2(const tMatrix3d &F) const
{

    tMatrix3d E = CalcGreenStrain(F);
    return mLambda * E.trace() * F;
}