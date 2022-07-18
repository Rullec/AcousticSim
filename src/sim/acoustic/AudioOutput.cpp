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

void data_callback(ma_device *pDevice, void *pOutput, const void *pInput,
                   ma_uint32 frameCount)
{
    cAudioOutputPtr mAudioOutput = cAudioOutput::getInstance();
    float *target_buf = static_cast<float *>(pOutput);
    if (mAudioOutput != nullptr && true == mAudioOutput->GetAudioFlag())
        mAudioOutput->SetContent(frameCount, target_buf);
    else
    {
        printf("didn't set any data into audio channel\n");
    }
    (void)pInput;
}

static int gCurFrame = 0;
//下面这个静态成员变量在类加载的时候就已经初始化好了
cAudioOutputPtr cAudioOutput::instance(new cAudioOutput(),
                                       cAudioOutput::DestroyInstance);

cAudioOutput::cAudioOutput()
{
    printf("AudioOutput constructed\n");
    mEnableAudio = true;
    this->Init();
}
cAudioOutput::~cAudioOutput()
{
    printf("AudioOutput deconstructed\n");
    ma_device_uninit(&device);
}
bool cAudioOutput::GetAudioFlag() const { return mEnableAudio; }
void cAudioOutput::SetAudioFlag(bool val) { mEnableAudio = val; }
cAudioOutputPtr cAudioOutput::getInstance() { return instance; }
void cAudioOutput::DestroyInstance(cAudioOutput *obj) { delete obj; }
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
void cAudioOutput::SetContent(unsigned int frame_count, float *tar_buf)
{
    if (mCurWave != nullptr)
    {
        // std::cout << "sum_wave = " << sum_wave.size() << std::endl;
        float *real_output = static_cast<float *>(tar_buf);

        float max_abs = 0;
        int num_of_data = mCurWave->GetNumOfData();
        for (uint cur_frame = 0; cur_frame < frame_count;
             cur_frame++, gCurFrame++)
        {
            real_output[cur_frame] = mCurWave->data[gCurFrame % num_of_data];
            // sine wave
            if (gCurFrame > num_of_data)
            {
                gCurFrame %= num_of_data;
                // printf("after a cycle\n");
            }
            // float cur_value = std::sin(gCurFrame / 10) * 0.1;
            // real_output[cur_frame * 1 + 0] = cur_value;
        }
        // std::cout << "max_abs = " << max_abs << std::endl;
    }
}

#include <fstream>

void cAudioOutput::SetWave(const tDiscretedWavePtr &wave)
{
    mCurWave = wave;
    // mCurWave->ChangeFrequency(DEVICE_SAMPLE_RATE);

    printf("[set wave] cur new freq = %d\n", mCurWave->GetFrequency());
    // std::ofstream fout("log/wave_play.txt");
    // fout << mCurWave->GetFrequency() << " HZ\n";
    // fout << mCurWave->data.transpose() << std::endl;
    // std::cout << "output to wave_play.txt\n";
}