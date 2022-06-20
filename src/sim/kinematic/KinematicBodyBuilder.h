#include "utils/DefUtil.h"
#include <string>
SIM_DECLARE_CLASS_AND_PTR(cBaseObject);
namespace Json
{
class Value;
};
cBaseObjectPtr BuildKinematicBody(const Json::Value &conf, int id_);
cBaseObjectPtr BuildKinematicBodyFromObjPath(std::string name, std::string obj_path, int id_);