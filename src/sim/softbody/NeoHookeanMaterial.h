#include "BaseMaterial.h"

class cNeoHookeanMaterial : public cBaseMaterial
{
public:
    virtual void Init(const Json::Value &conf) override;
    virtual tMatrix3d CalcP(const tMatrix3d &F) const override; // PK1
    virtual cFourOrderTensor CalcDPDF(const tMatrix3d &F) const override;
    virtual void CheckDPDF(const tMatrix3d &F) const override;

protected:
    double mMu;     // young's modulus
    double mLambda; // Poisson's ratio
    virtual cFourOrderTensor CalcDPDF_part1(const tMatrix3d &F) const;
    virtual cFourOrderTensor CalcDPDF_part2(const tMatrix3d &F) const;
    virtual void CheckDPDF_part1(const tMatrix3d &F) const;
    virtual void CheckDPDF_part2(const tMatrix3d &F) const;
    tMatrix3d CalcP_part1(const tMatrix3d & F) const;
    tMatrix3d CalcP_part2(const tMatrix3d & F) const;
};

void CheckDFinvTDF(const tMatrix3d &F);
// void CheckDdetFIDF(const tMatrix3d &F);
