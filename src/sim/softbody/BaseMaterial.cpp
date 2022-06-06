#include "BaseMaterial.h"
#include "utils/JsonUtil.h"
#include "utils/LogUtil.h"
std::string gMaterialModelTypeStr[eMaterialType::NUM_OF_MATERIAL_MODEL] = {
    "LINEAR_ELASTICITY", "COROTATED", "FIX_COROTATED", "STVK", "NEO_HOOKEAN"};
std::string BuildMaterialTypeStr(eMaterialType type)
{
    return gMaterialModelTypeStr[type];
};

cBaseMaterial::cBaseMaterial(eMaterialType type) : mType(type) {}

void cBaseMaterial::CheckDPDF(const tMatrix3d &F) const
{
    SIM_ERROR("hasn't been impled");
    exit(1);
}

eMaterialType BuildMaterialTypeFromStr(std::string name)
{
    for (size_t i = 0; i < eMaterialType::NUM_OF_MATERIAL_MODEL; i++)
    {
        if (gMaterialModelTypeStr[i] == name)
        {
            return static_cast<eMaterialType>(i);
        }
    }
    return eMaterialType::NUM_OF_MATERIAL_MODEL;
}

eMaterialType cBaseMaterial::GetType() const { return mType; }

// double cBaseMaterial::GetPoissonRatio() const { return this->mLambda; }
// double cBaseMaterial::GetYoungsModulus() const { return this->mMu; }

/*
    \lambda = E * poisson / ( (1 + poisson) * (1 - 2 * poisson) )
*/
double cBaseMaterial::GetLameFirstCoefLambda() const
{
    return mYoungsModulusNew * mPoissonRatioNew / ((1 + mPoissonRatioNew) * (1 - 2 * mPoissonRatioNew));
}
/*
    \mu = E / (2 * (1 + poisson))
*/
double cBaseMaterial::GetLameSecondCoefMu() const
{
    return mYoungsModulusNew / ( 2 * (1 + mPoissonRatioNew));
}
double cBaseMaterial::GetRho() const { return this->mRho; }

void cBaseMaterial::Init(const Json::Value &conf)
{

    if (conf.isMember("material"))
    {
        mMatPath = cJsonUtil::ParseAsString("material", conf);
        Json::Value material_json;
        SIM_ASSERT(cJsonUtil::LoadJson(mMatPath, material_json))
        mYoungsModulusNew = cJsonUtil::ParseAsDouble("youngs", material_json);
        mPoissonRatioNew = cJsonUtil::ParseAsDouble("poisson_ratio", material_json);
        mRho = cJsonUtil::ParseAsDouble("rho", material_json);
        mRayleighDamplingA =
            cJsonUtil::ParseAsFloat("rayleigh_damping_a", material_json);
        mRayleighDamplingB =
            cJsonUtil::ParseAsFloat("rayleigh_damping_b", material_json);

        SIM_INFO("rayleight damping a {} b {}", mRayleighDamplingA,
                 mRayleighDamplingB);
    }
    else
    {
        SIM_ERROR("no material key")
    }
}
std::string cBaseMaterial::GetMaterialPath() const { return this->mMatPath; }