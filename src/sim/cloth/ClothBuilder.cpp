#include "sim/cloth/ClothBuilder.h"
#include "sim/cloth/BaseCloth.h"
#include "sim/cloth/BaraffCloth.h"
// #include "sim/cloth/ImplicitCloth.h"
// #include "sim/cloth/PBDCloth.h"
// #include "sim/cloth/PDCloth.h"
// #include "sim/cloth/SemiCloth.h"
#include "utils/JsonUtil.h"
#include <iostream>
#include <string>

const std::string CLOTH_TYPE_KEY = "cloth_type";
cBaseObjectPtr BuildCloth(Json::Value conf, int obj_id)
{
    std::string cloth_type_str = cJsonUtil::ParseAsString(CLOTH_TYPE_KEY, conf);
    eClothType type = cBaseCloth::BuildClothType(cloth_type_str);
    cBaseObjectPtr ptr = nullptr;
    switch (type)
    {
    // case eClothType::IMPLICIT_CLOTH :
    //     S
    //     break;
    // case eClothType::SEMI_IMPLICIT_CLOTH:
    //     ptr = std::make_shared<cSemiCloth>(obj_id);
    //     break;
    // case eClothType::PBD_CLOTH:
    //     ptr = std::make_shared<cPBDCloth>(obj_id);
    //     break;
    // case eClothType::IMPLICIT_CLOTH:
    //     ptr = std::make_shared<cImplicitCloth>(obj_id);
    //     break;
    // case eClothType::PD_CLOTH:
    //     ptr = std::make_shared<cPDCloth>(obj_id);
    //     break;
    case eClothType::FEM_CLOTH:
        ptr = std::make_shared<cBaraffCloth>(obj_id);
        break;
    default:
        // std::cout << "pbd type = " << eClothType::PBD_CLOTH
        //           << " cur type = " << type << std::endl;
        SIM_ERROR("unsupported cloth type enum {}", type);
        break;
    }
    return ptr;
}