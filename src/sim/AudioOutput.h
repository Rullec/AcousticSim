#ifndef AUDIO_OUTPUT_H_
#pragma once
#include "sim/AudioWave.h"
#include "utils/DefUtil.h"
#include "utils/TimeUtil.hpp"

class cAudioOutput : std::enable_shared_from_this<cAudioOutput>
{
public:
    explicit cAudioOutput();
    virtual ~cAudioOutput();
    void Init();
    void SetContent(unsigned int frame_count, float *buf);
    void SetWave(const tDiscretedWavePtr &wave);

protected:
    tDiscretedWavePtr mCurWave;
};
SIM_DECLARE_PTR(cAudioOutput);
#endif