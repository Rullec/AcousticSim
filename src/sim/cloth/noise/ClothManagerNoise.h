#pragma once
#include "utils/JsonUtil.h"

// apply noise on it!
class cClothNoiseManager
{
public:
    cClothNoiseManager();
    virtual void Init(const Json::Value &value);

    int mNumOfSamplesPerProp;
};
SIM_DECLARE_PTR(cClothNoiseManager);