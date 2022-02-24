#pragma once
#include "utils/DefUtil.h"
SIM_DECLARE_CLASS_AND_PTR(cBaseCloth)
SIM_DECLARE_CLASS_AND_PTR(cSimScene)
namespace Json
{
class Value;
};
cBaseClothPtr BuildCloth(Json::Value conf, int obj_id);
