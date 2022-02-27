#pragma once
#include "utils/DefUtil.h"
SIM_DECLARE_CLASS_AND_PTR(cBaseObject)
namespace Json
{
class Value;
};
cBaseObjectPtr BuildCloth(Json::Value conf, int obj_id);
