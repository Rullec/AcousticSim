#include "KinematicBodyBuilder.h"
#include "KinematicBody.h"
#include "utils/JsonUtil.h"
cKinematicBodyPtr BuildKinematicBody(const Json::Value &conf, int id)
{
    cKinematicBodyPtr ptr = std::make_shared<cKinematicBody>(id);
    ptr->Init(conf);
    return ptr;
}

cKinematicBodyPtr BuildKinematicBodyFromObjPath(std::string name,
                                                std::string obj_path, int id_)
{
    Json::Value root;
    root["object_name"] = name;
    root["type"] = "custom";
    root["mesh_path"] = obj_path;
    root["target_aabb"] = Json::arrayValue;
    for (int i = 0; i < 3; i++)
    {
        root["target_aabb"].append(0);
    }
    root["scale"] = Json::arrayValue;
    for (int i = 0; i < 3; i++)
    {
        root["scale"].append(1);
    }

    root["translation"] = Json::arrayValue;
    for (int i = 0; i < 3; i++)
    {
        root["translation"].append(0);
    }
    root["orientation"] = Json::arrayValue;
    for (int i = 0; i < 3; i++)
    {
        root["orientation"].append(0);
    }
    return BuildKinematicBody(root, id_);
}