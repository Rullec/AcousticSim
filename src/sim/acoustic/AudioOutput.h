#ifndef AUDIO_OUTPUT_H_
#pragma once
#include "sim/acoustic/AudioWave.h"
#include "utils/DefUtil.h"
#include "utils/TimeUtil.hpp"

SIM_DECLARE_CLASS_AND_PTR(cAudioOutput);
class cAudioOutput : std::enable_shared_from_this<cAudioOutput>
{

public:
    static std::shared_ptr<cAudioOutput> getInstance();
    void SetContent(unsigned int frame_count, float *buf);
    void SetWave(const tDiscretedWavePtr &wave);
    bool GetAudioFlag() const;
    void SetAudioFlag(bool);

private:
    bool mEnableAudio;
    virtual void Init();
    cAudioOutput();
    cAudioOutput(const cAudioOutput &);
    cAudioOutput &operator=(const cAudioOutput &);
    ~cAudioOutput();
    static void DestroyInstance(cAudioOutput *); // define a deleter

    static std::shared_ptr<cAudioOutput> instance; // singleton, a shared ptr
    tDiscretedWavePtr mCurWave;
};

#endif