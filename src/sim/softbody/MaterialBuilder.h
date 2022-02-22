#include "sim/softbody/BaseMaterial.h"
namespace Json
{
    class Value;
};
SIM_DECLARE_PTR(cBaseMaterial);
cBaseMaterialPtr BuildMaterial(const Json::Value &);
cBaseMaterialPtr BuildDefaultMaterial(eMaterialType);
