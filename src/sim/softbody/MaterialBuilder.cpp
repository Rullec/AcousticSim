#include "MaterialBuilder.h"
#include "utils/JsonUtil.h"
#include "sim/softbody/BaseMaterial.h"
#include "sim/softbody/NeoHookeanMaterial.h"
#include "sim/softbody/StvkMaterial.h"

cBaseMaterialPtr BuildMaterial(const Json::Value &conf)
{
    std::string material_type_str = cJsonUtil::ParseAsString("material_type", conf);
    cBaseMaterialPtr mat = nullptr;

    switch (BuildMaterialTypeFromStr(material_type_str))
    {
    case eMaterialType::STVK:
    {
        mat = std::make_shared<cStvkMaterial>();
        printf("[log] build stvk material\n");
        break;
    }
    case eMaterialType::NEO_HOOKEAN:
    {
        mat = std::make_shared<cNeoHookeanMaterial>();
        printf("[log] build neo-hookean material\n");
        break;
    }
    default:
    {
        SIM_ERROR("unsupported material type {}", material_type_str);
        break;
    }
    }
    mat->Init(conf);
    return mat;
}

cBaseMaterialPtr BuildDefaultMaterial(eMaterialType type)
{
    float mMu = 2e4;
    float mLambda = 0.7;
    Json::Value conf;
    conf["youngs"] = mMu;
    conf["poisson_ratio"] = mLambda;
    cBaseMaterialPtr mat = nullptr;
    switch (type)
    {
    case eMaterialType::STVK:
    {
        mat = std::make_shared<cStvkMaterial>();
        break;
    }
    case eMaterialType::NEO_HOOKEAN:
    {
        mat = std::make_shared<cNeoHookeanMaterial>();
        break;
    }
    default:
    {
        SIM_ERROR("unsupported material type {}", type);
        break;
    }
    }
    mat->Init(conf);
    return mat;
}