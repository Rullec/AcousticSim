#include "utils/DefUtil.h"
#include <string>
SIM_DECLARE_CLASS_AND_PTR(cKinematicBody);
namespace Json
{
    class Value;
};
cBaseObjectPtr BuildSimObj(const Json::Value &conf, int id_);