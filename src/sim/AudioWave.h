#pragma once
#include "utils/DefUtil.h"
#include "utils/MathUtil.h"

class tDiscretedWave
{
public:
    tDiscretedWave(float dt, float duration);

    tVectorXf data;
    void Allocate();
    float GetDuration() const;
    int GetNumOfData() const;
    void SetData(const tVectorXf &data);
    void ChangeFrequency(int tar_freq);
    int GetFrequency() const;
    bool LoadFromFile(std::string path);
    void DumpToFile(std::string path);

protected:
    float dt;
    float duration;
};
SIM_DECLARE_PTR(tDiscretedWave);