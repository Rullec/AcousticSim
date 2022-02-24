#include "SimObjectBuilder.h"
#include "utils/JsonUtil.h"
#include "sim/BaseObject.h"
#include "sim/softbody/SoftBody.h"
#include "sim/softbody/SoftBodyImplicit.h"
#include "sim/AcousticSoftBody.h"
#include "sim/cloth/ClothBuilder.h"
cBaseObjectPtr BuildSimObj(const Json::Value &conf, int id_)
{
    eObjectType type = cBaseObject::BuildObjectType(cJsonUtil::ParseAsString("object_type", conf));
    cBaseObjectPtr object = nullptr;
    switch (type)
    {
    case eObjectType::ACOUSTIC_TYPE:
    {
        SIM_ERROR("do not support acoustic type");
        break;
    }
    case eObjectType::SOFTBODY_TYPE:
    {
        // object = std::make_shared<cSoftBody>(id_);
        object = std::make_shared<cSoftBodyImplicit>(id_);
        // object = std::make_shared<cAcousticSoftBody>(id_);
        break;
    }
    case eObjectType::CLOTH_TYPE:
    {
        auto result = BuildCloth(conf, id_);
        // object = std::dynamic_pointer_cast<cBaseObject>(BuildCloth(conf, id_));
        break;
    }
    default:
        SIM_ERROR("unrecognized object type {}", type);
        break;
    };
    if (object)
        object->Init(conf);
    return object;
}