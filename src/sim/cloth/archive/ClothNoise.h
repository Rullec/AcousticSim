#pragma once
#include "utils/DefUtil.h"
namespace Json{
    class Value;
}
enum eClothNoiseType
{
    WAVE_NOISE_TYPE = 0,     // wave noise
    MANUAL_NOISE_TYPE = 1,   // manually handing noise
    SPANNING_NOISE_TYPE = 2, // spanning result from given meshes
    EMPTY_NOISE_TYPE = 3,    // empty noise (no noise)
    NUM_NOISE_TYPE
};

extern std::string BuildClothNoiseNameFromType(eClothNoiseType);
extern eClothNoiseType BuildClothNoiseTypeFromName(std::string);
SIM_DECLARE_CLASS_AND_PTR(cLinctexCloth);

struct tClothNoise
{
    explicit tClothNoise(const Json::Value &value);
    eClothNoiseType RandomTheNoise() const;
    virtual bool GetEnableNoise() const;
    virtual int GetNumOfSamplesPerProp() const;
    virtual void GetWaveNoiseParam(double &_mMaxFoldAmp, int &_mMinFoldNum,
                                   int &_mMaxFoldNum) const;
    virtual void ApplySpanningNoise(cLinctexClothPtr cloth_ptr);
    virtual void Reset();

protected:
    double mClothProb[NUM_NOISE_TYPE];
    bool mEnableNoise;
    int mNumOfSamplesPerProp;
    double mMaxFoldAmp;
    int mMinFoldNum, mMaxFoldNum;

    // mesh data for
    std::string mSpanningMeshDataDir;           //
    std::vector<std::string> mSpanningMeshData; // data path
    std::vector<bool> mSpannedVisited;
    virtual void TestRandom() const;
    virtual void InitSpannedData();
    virtual int GetNumCurrentVisited();
    virtual int SelectUnvisitedMeshRandom() const;
};
SIM_DECLARE_PTR(tClothNoise);