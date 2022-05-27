#include "BaseMaterial.h"
#include "utils/LogUtil.h"
#include "utils/JsonUtil.h"
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

double cBaseMaterial::GetPoissonRatio() const { return this->mLambda; }
double cBaseMaterial::GetYoungsModulus() const { return this->mMu; }

double cBaseMaterial::GetRho() const { return this->mRho; }

void cBaseMaterial::Init(const Json::Value &conf)
{

    if (conf.isMember("material"))
    {
        std::string material_path = cJsonUtil::ParseAsString("material", conf);
        Json::Value material_json;
        SIM_ASSERT(cJsonUtil::LoadJson(material_path, material_json))
        mMu = cJsonUtil::ParseAsDouble("youngs", material_json);
        mLambda = cJsonUtil::ParseAsDouble("poisson_ratio", material_json);
        mRho = cJsonUtil::ParseAsDouble("rho", material_json);
    }
    else
    {
        SIM_ERROR("no material key")
    }
}