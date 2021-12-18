#include "AudioOutput.h"
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
    float *buf = static_cast<float *>(pOutput);
    gAudioOutput->SetContent(frameCount, buf);
    (void)pInput;
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
void cAudioOutput::SetContent(unsigned int frame_count, float *buf)
{
    cTimePoint cur_time = cTimeUtil::GetCurrentTime_chrono();
    float *real_output = static_cast<float *>(buf);
    float elasped_time = cTimeUtil::CalcTimeElaspedms(prev, cur_time);
    std::cout << "cur frame = " << gCurFrame << " frame_count = " << frame_count << ", interval = " << elasped_time << " ms, frenquencies = " << 1.0 / (1e-3 * elasped_time) << "HZ\n";

    for (uint cur_frame = 0; cur_frame < frame_count; cur_frame++, gCurFrame++)
    {
        float cur_value = std::sin(gCurFrame / 10) * 0.1;
        real_output[cur_frame * 1 + 0] = cur_value;
    }

    prev = cur_time;
}

cAudioOutput::~cAudioOutput()
{
    ma_device_uninit(&device);
}