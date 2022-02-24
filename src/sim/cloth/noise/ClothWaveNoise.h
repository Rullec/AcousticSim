#pragma once
#include "utils/DefUtil.h"
#include "utils/JsonUtil.h"
#include "ClothBaseNoise.h"

struct tClothWaveNoise : tClothNoiseBase
{
    tClothWaveNoise(const Json::Value &value);
    bool mEnableWaveNoise;
    double mMaxWaveAmp;
    int mMinWaveNum, mMaxWaveNum;
    int mNumOfSamplesPerProp;
    bool mEnableNoise;
};

SIM_DECLARE_PTR(tClothWaveNoise);