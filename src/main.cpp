
// #ifdef _WIN32
// #define VK_USE_PLATFORM_WIN32_KHR
// #define GLFW_INCLUDE_VULKAN
// #include "GLFW/glfw3.h"
// #define GLFW_EXPOSE_NATIVE_WIN32
// #include <GLFW/glfw3native.h>
// #endif

// #ifdef __linux__
// #define VK_USE_PLATFORM_XCB_KHR
// #include <GLFW/glfw3.h>
// // #define GLFW_EXPOSE_NATIVE_
// #include <GLFW/glfw3native.h>
// #endif

// #ifdef __APPLE__
// #include <GLFW/glfw3.h>
// #endif

// #include "scenes/SceneBuilder.h"
// #include <iostream>
// #include <memory>

// #include "utils/FileUtil.h"
// #include "utils/JsonUtil.h"
// #include "utils/LogUtil.h"
// #include <cmath>

// GLFWwindow *window = nullptr;
// std::shared_ptr<cDrawScene> draw_scene = nullptr;
// std::shared_ptr<cScene> scene = nullptr;
// bool esc_pushed = false;
// bool gPause = true;
// int gWindowWidth, gWindowHeight;

// static void ResizeCallback(GLFWwindow *window, int w, int h)
// {
//     draw_scene->Resize(w, h);
// }

// static void CursorPositionCallback(GLFWwindow *window, double xpos, double ypos)
// {
//     draw_scene->CursorMove(xpos, ypos);
// }

// void MouseButtonCallback(GLFWwindow *window, int button, int action, int mods)
// {
//     draw_scene->MouseButton(button, action, mods);
// }

// void KeyCallback(GLFWwindow *window, int key, int scancode, int action,
//                  int mods)
// {
//     draw_scene->Key(key, scancode, action, mods);
//     if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
//     {
//         esc_pushed = true;
//     }
//     else if (key == GLFW_KEY_I && action == GLFW_PRESS)
//     {
//         // gPause = !gPause;
//         // std::cout << "[log] simulation paused\n";
//     }
// }

// void ScrollCallback(GLFWwindow *window, double xoffset, double yoffset)
// {
//     draw_scene->Scroll(xoffset, yoffset);
// }
// void InitGlfw()
// {
//     glfwInit();
//     glfwWindowHint(GLFW_CLIENT_API, GLFW_NO_API);
//     window = glfwCreateWindow(gWindowWidth, gWindowHeight, "Acoustic Simulator",
//                               nullptr, nullptr);
//     glfwSetFramebufferSizeCallback(window, ResizeCallback);
//     glfwSetCursorPosCallback(window, CursorPositionCallback);
//     glfwSetMouseButtonCallback(window, MouseButtonCallback);
//     glfwSetKeyCallback(window, KeyCallback);
//     glfwSetScrollCallback(window, ScrollCallback);
// }

// #include "utils/LogUtil.h"
// #include "utils/TimeUtil.hpp"
// bool gEnableDraw = true;
// // bool gEnableDraw = false;
// void SimDraw(const std::string &conf_path, bool disable_imgui);
// void SimNoDraw(const std::string &conf_path);
// void ParseConfig(std::string conf);

// #include "cxxopts.hpp"
// void ParseArg(int argc, char *argv[], std::string &config_path,
//               bool &disable_imgui)
// {
//     try
//     {
//         cxxopts::Options options(argv[0], " - simple acoustic simulator");
//         options.positional_help("[optional args]").show_positional_help();

//         options.add_options()("conf", "config path",
//                               cxxopts::value<std::string>())(
//             "d,disable_imgui", "enalbe imgui rendering",
//             cxxopts::value<bool>()->default_value("false"));

//         options.parse_positional({"conf"});
//         auto result = options.parse(argc, argv);

//         if (result.count("conf"))
//         {
//             config_path = result["conf"].as<std::string>();
//             std::cout << "saw param config = " << config_path << std::endl;
//         }
//         if (result.count("d"))
//         {
//             disable_imgui = result["disable_imgui"].as<bool>();
//             std::cout << "saw param disable_imgui = " << disable_imgui
//                       << std::endl;
//         }
//     }
//     catch (const cxxopts::OptionException &e)
//     {
//         std::cout << "[error] when parsing, " << e.what() << std::endl;
//         exit(1);
//     }
//     SIM_INFO("conf path {}, enable imgui {}", config_path, disable_imgui);
// }

// int main(int argc, char **argv)
// {
//     // for (int i = 0; i < 100; i++)
//     // {
//     //     std::cout << cMathUtil::RandInt(0, 10) << std::endl;
//     // }
//     // exit(1);
//     std::string conf = "";
//     bool disable_imgui = false;
//     ParseArg(argc, argv, conf, disable_imgui);

//     ParseConfig(conf);

//     if (gEnableDraw == true)
//     {
//         SimDraw(conf, disable_imgui);
//     }
//     else
//     {
//         SimNoDraw(conf);
//     }
//     return 0;
// }

// void SimDraw(const std::string &conf, bool disable_imgui)
// {
//     InitGlfw();
//     scene = cSceneBuilder::BuildScene("sim_draw", true, !disable_imgui);
//     draw_scene = std::dynamic_pointer_cast<cDrawScene>(scene);
//     draw_scene->Init(conf);
//     auto last = cTimeUtil::GetCurrentTime_chrono();
//     while (glfwWindowShouldClose(window) == false && esc_pushed == false)
//     {
//         glfwPollEvents();

//         // 1. calc delta time for real time simulation
//         auto cur = cTimeUtil::GetCurrentTime_chrono();
//         double delta_time = cTimeUtil::CalcTimeElaspedms(last, cur) * 1e-3;

//         // 2. update
//         // delta_time = 1e-3;
//         // delta_time /= 4;
//         double limit = 1.0 / 30;
//         // double limit = 1e-4;
// #ifdef _WIN32
//         delta_time = std::min(delta_time, limit);
// #else
//         delta_time = std::min(delta_time, limit);
// #endif

//         // delta_time = 1e-4;
//         if (gPause == false)
//         {
//             // cTimeUtil::Begin("scene_update");
//             draw_scene->Update(delta_time);
//             // cTimeUtil::End("scene_update");
//         }

//         last = cur;
//     }

//     glfwDestroyWindow(window);

//     glfwTerminate();
// }

// void SimNoDraw(const std::string &conf_path)
// {
//     // InitGlfw();
//     scene = cSceneBuilder::BuildSimScene(conf_path);
//     scene->Init(conf_path);

//     int max_iters = 1e3;
//     int cur_iter = 0;
//     double dt = 1e-2;
//     while (++cur_iter < max_iters)
//     {
//         scene->Update(dt);
//         printf("[debug] iters %d/%d\n", cur_iter, max_iters);
//     }
// }

// void ParseConfig(std::string conf)
// {
//     SIM_ASSERT(cFileUtil::ExistsFile(conf) == true);
//     Json::Value root;
//     cJsonUtil::LoadJson(conf, root);
//     gPause = cJsonUtil::ParseAsBool("pause_at_first", root);
//     gEnableDraw = cJsonUtil::ParseAsBool("enable_draw", root);
//     if (gEnableDraw == true)
//     {
//         gWindowWidth = cJsonUtil::ParseAsInt("window_width", root);
//         gWindowHeight = cJsonUtil::ParseAsInt("window_height", root);
//     }
//     SIM_INFO("pause at first = {}", gPause);
// }

/*
Demonstrates playback of a sine wave.
Since all this example is doing is playing back a sine wave, we can disable decoding (and encoding) which will slightly
reduce the size of the executable. This is done with the `MA_NO_DECODING` and `MA_NO_ENCODING` options.
The generation of sine wave is achieved via the `ma_waveform` API. A waveform is a data source which means it can be
seamlessly plugged into the `ma_data_source_*()` family of APIs as well.
A waveform is initialized using the standard config/init pattern used throughout all of miniaudio. Frames are read via
the `ma_waveform_read_pcm_frames()` API.
This example works with Emscripten.
*/
#include <iostream>
#define MA_NO_DECODING
#define MA_NO_ENCODING
#define MINIAUDIO_IMPLEMENTATION
#include "./miniaudio.h"

#include <stdio.h>

#ifdef __EMSCRIPTEN__
#include <emscripten.h>

void main_loop__em()
{
}
#endif
#include <cmath>
#define DEVICE_FORMAT ma_format_f32
#define DEVICE_CHANNELS 1
#define DEVICE_SAMPLE_RATE 48000
using uint = unsigned int;
int gCurFrame = 0;
void data_callback(ma_device *pDevice, void *pOutput, const void *pInput, ma_uint32 frameCount)
{
    /*
            float* pFramesOutF32 = (float*)pFramesOut;
        for (iFrame = 0; iFrame < frameCount; iFrame += 1) {
            float s = ma_waveform_sine_f32(pWaveform->time, pWaveform->config.amplitude);
            pWaveform->time += pWaveform->advance;

            for (iChannel = 0; iChannel < pWaveform->config.channels; iChannel += 1) {
                pFramesOutF32[iFrame*pWaveform->config.channels + iChannel] = s;
            }
        }
    */
    float *real_output = static_cast<float *>(pOutput);
    std::cout << "cur frame = " << gCurFrame << std::endl;
    for (uint cur_frame = 0; cur_frame < frameCount; cur_frame++, gCurFrame++)
    {

        float cur_value = std::sin(gCurFrame / 1000);
        real_output[cur_frame * 1 + 0] = cur_value;
    }
    // ma_waveform *pSineWave;

    // MA_ASSERT(pDevice->playback.channels == DEVICE_CHANNELS);

    // pSineWave = (ma_waveform *)pDevice->pUserData;
    // MA_ASSERT(pSineWave != NULL);

    // ma_waveform_read_pcm_frames(pSineWave, pOutput, frameCount);
    // std::cout << "frame count = " << frameCount << std::endl;
    (void)pInput; /* Unused. */
}
/*

The use of the decoder API might be unnecessary. So you're decoding audio via FFMpeg before anything even touches miniaudio? In that case, I would consider just bypassing ma_decoder completely and just write straight to the output buffer in the device's data callback. Decoders and devices are 100% decoupled - if you already have the raw decoded data, you can just pump it straight through to the device.

If that's the approach you want to take, take a look at the simple_playback example for how to implement the data callback. No need for the decoder stuff so just remove all that. Instead just replace that with your FFMpeg code. Make sure you set the device's data format to that of your raw audio data.
*/
int main(int argc, char **argv)
{
    ma_waveform sineWave;
    ma_device_config deviceConfig;
    ma_device device;
    ma_waveform_config sineWaveConfig;

    sineWaveConfig = ma_waveform_config_init(DEVICE_FORMAT, DEVICE_CHANNELS, DEVICE_SAMPLE_RATE, ma_waveform_type_sine, 0.2, 220);
    ma_waveform_init(&sineWaveConfig, &sineWave);

    deviceConfig = ma_device_config_init(ma_device_type_playback);
    deviceConfig.playback.format = DEVICE_FORMAT;
    deviceConfig.playback.channels = DEVICE_CHANNELS;
    deviceConfig.sampleRate = DEVICE_SAMPLE_RATE;
    deviceConfig.dataCallback = data_callback;
    deviceConfig.pUserData = &sineWave;

    if (ma_device_init(NULL, &deviceConfig, &device) != MA_SUCCESS)
    {
        printf("Failed to open playback device.\n");
        return -4;
    }

    printf("Device Name: %s\n", device.playback.name);

    if (ma_device_start(&device) != MA_SUCCESS)
    {
        printf("Failed to start playback device.\n");
        ma_device_uninit(&device);
        return -5;
    }

#ifdef __EMSCRIPTEN__
    emscripten_set_main_loop(main_loop__em, 0, 1);
#else
    printf("Press Enter to quit...\n");
    getchar();
#endif

    ma_device_uninit(&device);

    (void)argc;
    (void)argv;
    return 0;
}