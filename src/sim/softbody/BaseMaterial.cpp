#include "BaseMaterial.h"
#include "utils/LogUtil.h"

std::string gMaterialModelTypeStr[eMaterialType::NUM_OF_MATERIAL_MODEL] =
    {
        "LINEAR_ELASTICITY",
        "COROTATED",
        "FIX_COROTATED",
        "STVK",
        "NEO_HOOKEAN"};
std::string BuildMaterialTypeStr(eMaterialType type)
{
    return gMaterialModelTypeStr[type];
};

cBaseMaterial::cBaseMaterial(eMaterialType type) : mType(type)
{
}

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

eMaterialType cBaseMaterial::GetType() const
{
    return mType;
}