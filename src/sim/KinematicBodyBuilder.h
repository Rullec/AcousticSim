#include "utils/DefUtil.h"
#include <string>
SIM_DECLARE_CLASS_AND_PTR(cKinematicBody);
namespace Json
{
class Value;
};
cKinematicBodyPtr BuildKinematicBody(const Json::Value &conf, int id_);
cKinematicBodyPtr BuildKinematicBodyFromObjPath(std::string name, std::string obj_path, int id_);