#pragma once
#include "sim/softbody/FourOrderTensor.h"
#include "sim/softbody/ThreeOrderTensor.h"
#include "utils/DefUtil.h"
#include "utils/MathUtil.h"
namespace Json
{
class Value;
};

enum eMaterialType
{
    LINEAR_ELASTICITY = 0,
    COROTATED,
    FIX_COROTATED,
    STVK,
    NEO_HOOKEAN,
    NUM_OF_MATERIAL_MODEL
};

class cBaseMaterial : std::enable_shared_from_this<cBaseMaterial>
{
public:
    explicit cBaseMaterial(eMaterialType type);
    eMaterialType GetType() const;
    virtual void Init(const Json::Value &conf);
    virtual tMatrix3d CalcP(const tMatrix3d &F) const = 0; // PK1
    virtual cFourOrderTensor CalcDPDF(const tMatrix3d &F) const = 0;
    virtual void CheckDPDF(const tMatrix3d &F) const;
    virtual double GetPoissonRatio() const;
    virtual double GetYoungsModulus() const;
    virtual double GetRho() const;

protected:
    double mMu;     // young's modulus
    double mLambda; // Poisson's ratio
    double mRho;    // density kg.m-3
    eMaterialType mType;
};

std::string BuildMaterialTypeStr(eMaterialType type);
eMaterialType BuildMaterialTypeFromStr(std::string);
extern std::string gMaterialModelTypeStr[];