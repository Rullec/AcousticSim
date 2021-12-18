#ifndef AUDIO_OUTPUT_H_
#pragma once
#include "utils/TimeUtil.hpp"
#include "utils/DefUtil.h"

class cAudioOutput : std::enable_shared_from_this<cAudioOutput>
{
public:
    explicit cAudioOutput();
    virtual ~cAudioOutput();
    void Init();
    void SetContent(unsigned int frame_count, float *buf);
protected:
};
SIM_DECLARE_PTR(cAudioOutput);
#endif