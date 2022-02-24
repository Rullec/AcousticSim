#include "AudioOutput.h"
#include "utils/MathUtil.h"
#include <cmath>
#include <iostream>
#define MA_NO_DECODING
#define MA_NO_ENCODING
#define MINIAUDIO_IMPLEMENTATION
#include "miniaudio.h"
using uint = unsigned int;
#define DEVICE_FORMAT ma_format_f32
#define DEVICE_CHANNELS 1
#define DEVICE_SAMPLE_RATE 48000 // frames per second

ma_device_config deviceConfig;
ma_device device;
cAudioOutputPtr gAudioOutput = nullptr;
void data_callback(ma_device *pDevice, void *pOutput, const void *pInput, ma_uint32 frameCount)
{
    if (gAudioOutput != nullptr)
    {
        float *buf = static_cast<float *>(pOutput);
        gAudioOutput->SetContent(frameCount, buf);
        (void)pInput;
    }
}

static int gCurFrame = 0;
cAudioOutput::cAudioOutput()
{
}
cTimePoint prev = cTimeUtil::GetCurrentTime_chrono();
void cAudioOutput::Init()
{
    deviceConfig = ma_device_config_init(ma_device_type_playback);
    deviceConfig.playback.format = DEVICE_FORMAT;
    deviceConfig.playback.channels = DEVICE_CHANNELS;
    deviceConfig.sampleRate = DEVICE_SAMPLE_RATE;
    deviceConfig.dataCallback = data_callback;
    deviceConfig.pUserData = nullptr;

    if (ma_device_init(NULL, &deviceConfig, &device) != MA_SUCCESS)
    {
        printf("Failed to open playback device.\n");
        exit(4);
    }

    printf("Device Name: %s\n", device.playback.name);

    if (ma_device_start(&device) != MA_SUCCESS)
    {
        printf("Failed to start playback device.\n");
        ma_device_uninit(&device);
        exit(5);
    }
}
extern tVectorXd sum_wave;
tVectorXd Interpolate(const tVectorXd &old, int fps)
{
    tVectorXd new_vec = tVectorXd::Zero(fps);
    size_t old_size = old.size();
    for (size_t i = 0; i < fps; i++)
    {
        // 1. find the gap
        float cur_perc = i * 1.0 / fps;
        size_t old_st = min(old.size() - 1, int(cur_perc * old.size()));
        size_t old_ed = min(old.size() - 1, old_st + 1);
        float st_perc = (cur_perc - old_st * 1.0 / old.size()) / (1.0 / old.size());
        float ed_perc = 1.0 - st_perc;
        new_vec[i] = (st_perc * old[old_st] +
                      ed_perc * old[old_ed]) *
                     1e3;
        // std::cout << "new vec " << i << " = " << new_vec[i] << std::endl;
    }
    return new_vec;
}
void cAudioOutput::SetContent(unsigned int frame_count, float *buf)
{
    {
        // cTimePoint cur_time = cTimeUtil::GetCurrentTime_chrono();
        // float *real_output = static_cast<float *>(buf);
        // float elasped_time = cTimeUtil::CalcTimeElaspedms(prev, cur_time);
        // std::cout << "cur frame = " << gCurFrame << " frame_count = " << frame_count << ", interval = " << elasped_time << " ms, frenquencies = " << 1.0 / (1e-3 * elasped_time) << "HZ\n";

        // for (uint cur_frame = 0; cur_frame < frame_count; cur_frame++, gCurFrame++)
        // {
        //     float cur_value = std::sin(gCurFrame / 10) * 0.1;

        //     real_output[cur_frame * 1 + 0] = cur_value;
        // }

        // prev = cur_time;
    }
    if (sum_wave.size() != 0)
    {
        // std::cout << "sum_wave = " << sum_wave.size() << std::endl;
        float *real_output = static_cast<float *>(buf);
        int one_sec_cur_frame = gCurFrame % DEVICE_SAMPLE_RATE;
        tVectorXd new_wave = Interpolate(sum_wave, DEVICE_SAMPLE_RATE);
        // std::cout << "new_wave = " << new_wave.transpose() << std::endl;
        float max_abs = 0;
        // if (int(gCurFrame / DEVICE_SAMPLE_RATE) % 2 != 0)
        //     return;
        for (uint cur_frame = 0; cur_frame < frame_count; cur_frame++, gCurFrame++)
        {
            real_output[cur_frame] = new_wave[(gCurFrame % DEVICE_SAMPLE_RATE)];
            // max_abs = max(std::fabs(real_output[cur_frame]), max_abs);
            float cur_value = std::sin(gCurFrame / 10) * 0.1;

            // real_output[cur_frame * 1 + 0] = cur_value;
        }
        // std::cout << "max_abs = " << max_abs << std::endl;
    }
}

cAudioOutput::~cAudioOutput()
{
    ma_device_uninit(&device);
}