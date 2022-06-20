#include "SimObjectBuilder.h"
#include "sim/ViscoMassSpring.h"
#include "sim/acoustic/ModalSoftBody.h"
#include "sim/acoustic/TransferSoftBody.h"
#include "sim/BaseObject.h"
#include "sim/kinematic/KinematicBodyBuilder.h"
#include "sim/cloth/ClothBuilder.h"
#include "sim/softbody/SoftBody.h"
#include "sim/softbody/SoftBodyImplicit.h"
#include "utils/JsonUtil.h"
cBaseObjectPtr BuildSimObj(const Json::Value &conf, int id_)
{
    eObjectType type = cBaseObject::BuildObjectType(
        cJsonUtil::ParseAsString("object_type", conf));
    cBaseObjectPtr object = nullptr;
    switch (type)
    {
    case eObjectType::MODAL_ANALYSIS_TYPE:
    {
        object = std::make_shared<cModalSoftBody>(id_);
        break;
    }
    case eObjectType::ACOUSTIC_TRANSFER_TYPE:
    {
        object = std::make_shared<cTransferSoftBody>(id_);
        break;
    }
    case eObjectType::SOFTBODY_TYPE:
    {
        // object = std::make_shared<cSoftBody>(id_);
        object = std::make_shared<cSoftBodyImplicit>(id_);

        break;
    }
    case eObjectType::KINEMATICBODY_TYPE:
    {
        object = BuildKinematicBody(conf, id_);
        break;
    }
    case eObjectType::CLOTH_TYPE:
    {
        object = BuildCloth(conf, id_);
        break;
    }
    case eObjectType::VISCOSITY_MASS_SPRING_TYPE:{
        object = std::make_shared<cViscoMassSpring>(id_);
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