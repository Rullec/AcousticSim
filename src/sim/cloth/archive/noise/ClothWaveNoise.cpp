#include "ClothWaveNoise.h"
#include <iostream>
/**
 * \brief       init the data augmentation strucutre
 */
#define _USE_MATH_DEFINES
#include <math.h>
tClothWaveNoise::tClothWaveNoise(const Json::Value &conf)
{
    mEnableNoise = cJsonUtil::ParseAsBool("enable_noise", conf);
    mNumOfSamplesPerProp = cJsonUtil::ParseAsInt("samples_per_prop", conf);
    // mEnableInitYRotation = cJsonUtil::ParseAsBool("enable_init_rotation",
    // conf); mEnableFoldNoise = cJsonUtil::ParseAsBool("enable_fold_noise",
    // conf); mEnableInitYPosNoise =
    //     cJsonUtil::ParseAsBool("enable_gaussian_pos_noise", conf);
    // mInitYPosNoiseStd = cJsonUtil::ParseAsDouble("gaussian_std", conf);
    // mFoldCoef = cJsonUtil::ParseAsDouble("fold_coef", conf);
    mEnableLowFreqNoise =
        cJsonUtil::ParseAsDouble("enable_low_freq_noise", conf);
    mMaxFoldAmp = cJsonUtil::ParseAsDouble("max_fold_amp", conf);
    mMinFoldNum = cJsonUtil::ParseAsInt("min_fold_num", conf);
    mMaxFoldNum = cJsonUtil::ParseAsInt("max_fold_num", conf);
    // SIM_ASSERT(mEnableInitYRotation == false);
    // SIM_ASSERT(mEnableFoldNoise == true);
    // std::cout << mNumOfSamplesPerProp << " " << mEnableInitYRotation << " "
    // << mEnableInitYPosNoise << " " << this->mInitYPosNoiseStd << std::endl;
    // exit(0);
}

