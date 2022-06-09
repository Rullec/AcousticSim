#include "sim/softbody/BaseMaterial.h"
namespace Json
{
    class Value;
};
SIM_DECLARE_PTR(cBaseMaterial);
cBaseMaterialPtr BuildMaterial(std::string mat_path, eMaterialType mat_type);