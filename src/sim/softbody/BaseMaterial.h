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

struct cBaseMaterial : std::enable_shared_from_this<cBaseMaterial>
{
public:
    explicit cBaseMaterial(eMaterialType type);
    eMaterialType GetType() const;
    virtual void Init(std::string math_path);
    // virtual void Init(const Json::Value &conf);
    // virtual void Init(std::string math_path, double youngs_modulus, double poisson_ratio, double rho, double damping_a, double damping_b);
    virtual tMatrix3d CalcP(const tMatrix3d &F) const = 0; // PK1
    virtual cFourOrderTensor CalcDPDF(const tMatrix3d &F) const = 0;
    virtual void CheckDPDF(const tMatrix3d &F) const;
    virtual double GetLameFirstCoefLambda() const;
    virtual double GetLameSecondCoefMu() const;
    virtual double GetRho() const;
    virtual std::string GetMaterialPath() const;
    double mYoungsModulusNew;     // young's modulus
    double mPoissonRatioNew; // Poisson's ratio
    double mRho;    // density kg.m-3
    float mRayleighDamplingA,
        mRayleighDamplingB; // rayleigh damping for mass mat and stiffness mat
    std::string mMatPath;   // material path
    std::string mName;
    eMaterialType mType;
};

std::string BuildMaterialTypeStr(eMaterialType type);
eMaterialType BuildMaterialTypeFromStr(std::string);
extern std::string gMaterialModelTypeStr[];