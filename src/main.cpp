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
#include "sim/AudioOutput.h"
extern cAudioOutputPtr gAudioOutput;
#include "sim/softbody/FourOrderTensorTest.hpp"

GLFWwindow *window = nullptr;
std::shared_ptr<cDrawScene> draw_scene = nullptr;
std::shared_ptr<cScene> scene = nullptr;
bool esc_pushed = false;
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
// bool gEnableDraw = true;
// bool gEnableDraw = false;
void SimDraw(const std::string &conf_path);
void ParseConfig(std::string conf);

#include "cxxopts.hpp"
void ParseArg(int argc, char *argv[], std::string &config_path)
{
    try
    {
        cxxopts::Options options(argv[0], " - simple acoustic simulator");
        options.positional_help("[optional args]").show_positional_help();

        options.add_options()("conf", "config path",
                              cxxopts::value<std::string>());
        //                   (
        // "d,disable_imgui", "enalbe imgui rendering",
        // cxxopts::value<bool>()->default_value("false"));

        options.parse_positional({"conf"});
        auto result = options.parse(argc, argv);

        if (result.count("conf"))
        {
            config_path = result["conf"].as<std::string>();
            std::cout << "saw param config = " << config_path << std::endl;
        }
        // if (result.count("d"))
        // {
        //     disable_imgui = result["disable_imgui"].as<bool>();
        //     std::cout << "saw param disable_imgui = " << disable_imgui
        //               << std::endl;
        // }
    }
    catch (const cxxopts::OptionException &e)
    {
        std::cout << "[error] when parsing, " << e.what() << std::endl;
        exit(1);
    }
    // SIM_INFO("conf path {}, enable imgui {}", config_path, disable_imgui);
    SIM_INFO("conf path {}", config_path);
}

#include "CGSolver.h"
void test_cg()
{
    int dims = 1000;
    tMatrixXd A = tMatrixXd::Random(dims, dims);
    A = A.transpose() * A;
    tVectorXd b = tVectorXd::Random(dims);
    cTimeUtil::Begin("cg");
    tVectorXd x = cCGSolver::Solve(A, b);
    cTimeUtil::End("cg");
    tVectorXd residual = A * x - b;
    std::cout << "CG residual = " << residual.norm() << std::endl;
    // std::cout << "CG x = " << x.transpose() << std::endl;
    cTimeUtil::Begin("inv");
    tVectorXd inv_x = A.inverse() * b;
    cTimeUtil::End("inv");
    tVectorXd inv_residual = A * inv_x - b;
    std::cout << "inv residual = " << inv_residual.norm() << std::endl;
    // std::cout << "inv x = " << inv_x.transpose() << std::endl;
    exit(1);
}
#include "sim/cloth/BaraffMaterial.h"

#include "sim/cloth/QBendingMaterial.h"
int main(int argc, char **argv)
{
    auto mat = std::make_shared<cBaraffMaterial>();
    mat->CheckShearingForce();
    // mat->CheckStretchForce();
    exit(1);
    // cQBendingMaterial::CheckStretchForce();
    // exit(1);
    // gAudioOutput =  std::make_shared<cAudioOutput>();
    // gAudioOutput->Init();
    std::string conf = "";
    ParseArg(argc, argv, conf);

    ParseConfig(conf);

    SimDraw(conf);

    return 0;
}
#include "scenes/DrawSceneImGUI.h"
void SimDraw(const std::string &conf)
{
    InitGlfw();
    scene = cSceneBuilder::BuildScene("sim_draw", true, true);
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
        float limit = 1.0 / 30;
#ifdef _WIN32
        delta_time = min(delta_time, limit);
#else
        delta_time = std::min(delta_time, limit);
#endif

        // delta_time = 1e-4;
        // if (gPause == false)
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
void ParseConfig(std::string conf)
{
    if (cFileUtil::ExistsFile(conf) == false)
    {
        if (conf.size() == 0)
        {
            SIM_ERROR("Please specify the config file!");
        }
        else
        {
            SIM_ERROR("current config file {} doesn't exist", conf);
        }
    }
    Json::Value root;
    cJsonUtil::LoadJson(conf, root);
    // gPause = cJsonUtil::ParseAsBool("pause_at_first", root);
    gWindowWidth = cJsonUtil::ParseAsInt("window_width", root);
    gWindowHeight = cJsonUtil::ParseAsInt("window_height", root);
    // SIM_INFO("pause at first = {}", gPause);
}
