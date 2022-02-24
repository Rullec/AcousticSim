#pragma once
#include <string>
enum eClothNoiseType
{
    WAVE_NOISE_TYPE = 0,
    HAND_NOISE_TYPE,
    NUM_NOISE_TYPE
};

extern std::string gClothNoiseTypeStr[];
struct tClothNoiseBase
{
public:
    tClothNoiseBase();
    virtual double GetProbabilty();

    double mProbability;
};