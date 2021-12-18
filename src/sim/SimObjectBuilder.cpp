#include "SimObjectBuilder.h"
#include "utils/JsonUtil.h"
#include "sim/BaseObject.h"
cBaseObjectPtr BuildSimObj(const Json::Value &conf, int id_)
{
    eObjectType type = cBaseObject::BuildObjectType(cJsonUtil::ParseAsString("object_type", conf));
    switch (type)
    {
    case eObjectType::ACOUSTIC_TYPE:
    {
    }
    break;
    }
}