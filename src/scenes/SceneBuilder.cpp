#include "SceneBuilder.h"
#include "DrawScene.h"
#include "DrawSceneImGUI.h"
#include "utils/JsonUtil.h"
#include "utils/LogUtil.h"
std::shared_ptr<cScene>
cSceneBuilder::BuildScene(const std::string type, bool enable_draw /*= true*/,
                          bool enable_draw_imgui /*= false */)
{
    if (enable_draw == false && enable_draw_imgui == true)
    {
        SIM_ERROR("enable draw is {} but enable draw imgui is {}: illegal",
                  enable_draw, enable_draw_imgui);
        exit(1);
    }
    std::shared_ptr<cScene> scene = nullptr;
    if (enable_draw == true)
    {
        if (enable_draw_imgui == false)
        {
            scene = std::dynamic_pointer_cast<cScene>(
                std::make_shared<cDrawScene>());
        }
        else
        {
            scene = std::dynamic_pointer_cast<cScene>(
                std::make_shared<cDrawSceneImGui>());
        }
    }
    else
    {
        SIM_ASSERT(false);
    }
    return scene;
}
std::shared_ptr<cSimScene>
cSceneBuilder::BuildSimScene(const std::string config_file)
{
    Json::Value root;
    cJsonUtil::LoadJson(config_file, root);
    std::string type = cJsonUtil::ParseAsString("scene_type", root);
    eSceneType scheme = cSimScene::BuildSceneType(type);
    std::shared_ptr<cSimScene> scene = nullptr;
    switch (scheme)
    {
    case eSceneType::SCENE_SIM:
        scene = std::make_shared<cSimScene>();
        break;
    default:
        SIM_ERROR("unsupported sim scene {}", type);
        break;
    }
    return scene;
}