#ifdef _WIN32
#define VK_USE_PLATFORM_WIN32_KHR
#define GLFW_INCLUDE_VULKAN
#include "GLFW/glfw3.h"
#define GLFW_EXPOSE_NATIVE_WIN32
#include <GLFW/glfw3native.h>
#endif

#ifdef __linux__
#define VK_USE_PLATFORM_XCB_KHR
#include <GLFW/glfw3.h>
// #define GLFW_EXPOSE_NATIVE_
#include <GLFW/glfw3native.h>
#endif

#ifdef __APPLE__
#include <GLFW/glfw3.h>
#endif

#include "scenes/SceneBuilder.h"
#include <iostream>
#include <memory>

#include "utils/FileUtil.h"
#include "utils/JsonUtil.h"
#include "utils/LogUtil.h"
#include <cmath>

GLFWwindow *window = nullptr;
std::shared_ptr<cDrawScene> draw_scene = nullptr;
std::shared_ptr<cScene> scene = nullptr;
bool esc_pushed = false;
bool gPause = true;
int gWindowWidth, gWindowHeight;

static void ResizeCallback(GLFWwindow *window, int w, int h)
{
    draw_scene->Resize(w, h);
}

static void CursorPositionCallback(GLFWwindow *window, double xpos, double ypos)
{
    draw_scene->CursorMove(xpos, ypos);
}

void MouseButtonCallback(GLFWwindow *window, int button, int action, int mods)
{
    draw_scene->MouseButton(button, action, mods);
}

void KeyCallback(GLFWwindow *window, int key, int scancode, int action,
                 int mods)
{
    draw_scene->Key(key, scancode, action, mods);
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
    {
        esc_pushed = true;
    }
    else if (key == GLFW_KEY_I && action == GLFW_PRESS)
    {
        // gPause = !gPause;
        // std::cout << "[log] simulation paused\n";
    }
}

void ScrollCallback(GLFWwindow *window, double xoffset, double yoffset)
{
    draw_scene->Scroll(xoffset, yoffset);
}
void InitGlfw()
{
    glfwInit();
    glfwWindowHint(GLFW_CLIENT_API, GLFW_NO_API);
    window = glfwCreateWindow(gWindowWidth, gWindowHeight, "Acoustic Simulator",
                              nullptr, nullptr);
    glfwSetFramebufferSizeCallback(window, ResizeCallback);
    glfwSetCursorPosCallback(window, CursorPositionCallback);
    glfwSetMouseButtonCallback(window, MouseButtonCallback);
    glfwSetKeyCallback(window, KeyCallback);
    glfwSetScrollCallback(window, ScrollCallback);
}

#include "utils/LogUtil.h"
#include "utils/TimeUtil.hpp"
bool gEnableDraw = true;
// bool gEnableDraw = false;
void SimDraw(const std::string &conf_path, bool disable_imgui);
void SimNoDraw(const std::string &conf_path);
void ParseConfig(std::string conf);

#include "cxxopts.hpp"
void ParseArg(int argc, char *argv[], std::string &config_path,
              bool &disable_imgui)
{
    try
    {
        cxxopts::Options options(argv[0], " - simple acoustic simulator");
        options.positional_help("[optional args]").show_positional_help();

        options.add_options()("conf", "config path",
                              cxxopts::value<std::string>())(
            "d,disable_imgui", "enalbe imgui rendering",
            cxxopts::value<bool>()->default_value("false"));

        options.parse_positional({"conf"});
        auto result = options.parse(argc, argv);

        if (result.count("conf"))
        {
            config_path = result["conf"].as<std::string>();
            std::cout << "saw param config = " << config_path << std::endl;
        }
        if (result.count("d"))
        {
            disable_imgui = result["disable_imgui"].as<bool>();
            std::cout << "saw param disable_imgui = " << disable_imgui
                      << std::endl;
        }
    }
    catch (const cxxopts::OptionException &e)
    {
        std::cout << "[error] when parsing, " << e.what() << std::endl;
        exit(1);
    }
    SIM_INFO("conf path {}, enable imgui {}", config_path, disable_imgui);
}

#include "sim/AudioOutput.h"
cAudioOutputPtr gAudioOutput;
int main(int argc, char **argv)
{
    std::string conf = "";
    bool disable_imgui = false;
    ParseArg(argc, argv, conf, disable_imgui);

    ParseConfig(conf);

    if (gEnableDraw == true)
    {
        SimDraw(conf, disable_imgui);
    }
    else
    {
        SimNoDraw(conf);
    }
    return 0;
}
#include "scenes/DrawSceneImGUI.h"
void SimDraw(const std::string &conf, bool disable_imgui)
{
    InitGlfw();
    scene = cSceneBuilder::BuildScene("sim_draw", true, !disable_imgui);
    draw_scene = std::dynamic_pointer_cast<cDrawSceneImGui>(scene);
    draw_scene->Init(conf);
    auto last = cTimeUtil::GetCurrentTime_chrono();
    while (glfwWindowShouldClose(window) == false && esc_pushed == false)
    {
        glfwPollEvents();

        // 1. calc delta time for real time simulation
        auto cur = cTimeUtil::GetCurrentTime_chrono();
        float delta_time = cTimeUtil::CalcTimeElaspedms(last, cur) * 1e-3;

        // 2. update
        // delta_time = 1e-3;
        // delta_time /= 4;
        float limit = 1.0 / 30;
        // float limit = 1e-4;
#ifdef _WIN32
        delta_time = min(delta_time, limit);
#else
        delta_time = std::min(delta_time, limit);
#endif

        // delta_time = 1e-4;
        if (gPause == false)
        {
            // cTimeUtil::Begin("scene_update");
            draw_scene->Update(delta_time);
            // cTimeUtil::End("scene_update");
        }

        last = cur;
    }

    glfwDestroyWindow(window);

    glfwTerminate();
}

void SimNoDraw(const std::string &conf_path)
{
    // InitGlfw();
    scene = cSceneBuilder::BuildSimScene(conf_path);
    scene->Init(conf_path);

    int max_iters = 1e3;
    int cur_iter = 0;
    float dt = 1e-2;
    while (++cur_iter < max_iters)
    {
        scene->Update(dt);
        printf("[debug] iters %d/%d\n", cur_iter, max_iters);
    }
}

void ParseConfig(std::string conf)
{
    SIM_ASSERT(cFileUtil::ExistsFile(conf) == true);
    Json::Value root;
    cJsonUtil::LoadJson(conf, root);
    gPause = cJsonUtil::ParseAsBool("pause_at_first", root);
    gEnableDraw = cJsonUtil::ParseAsBool("enable_draw", root);
    if (gEnableDraw == true)
    {
        gWindowWidth = cJsonUtil::ParseAsInt("window_width", root);
        gWindowHeight = cJsonUtil::ParseAsInt("window_height", root);
    }
    SIM_INFO("pause at first = {}", gPause);
}
