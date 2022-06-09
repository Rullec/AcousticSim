#include "MaterialBuilder.h"
#include "sim/softbody/BaseMaterial.h"
#include "sim/softbody/NeoHookeanMaterial.h"
#include "sim/softbody/StvkMaterial.h"
#include "utils/JsonUtil.h"

cBaseMaterialPtr BuildMaterial(std::string mat_path, eMaterialType mat_type)
{
    cBaseMaterialPtr mat = nullptr;

    switch (mat_type)
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
        SIM_ERROR("unsupported material type {}", mat_type);
        break;
    }
    }
    mat->Init(mat_path);

    return mat;
}
