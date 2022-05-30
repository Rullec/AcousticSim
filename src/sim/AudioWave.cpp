#include "AudioWave.h"
#include "utils/JsonUtil.h"
#include "utils/LogUtil.h"
tDiscretedWave::tDiscretedWave(float dt_, float duration_)
{
    dt = dt_;
    duration = duration_;
    Allocate();
}
void tDiscretedWave::Allocate()
{
    SIM_ASSERT(dt > 1e-6);
    SIM_ASSERT(duration > 0);

    data.resize(duration / dt);
}
float tDiscretedWave::GetDuration() const { return duration; }

int tDiscretedWave::GetNumOfData() const { return data.size(); }

void tDiscretedWave::SetData(const tVectorXf &data_)
{
    // SIM_ASSERT(data_.size() == GetNumOfData());
    data = data_;
    duration = dt * data.size();
}

tVectorXf Interpolate(const tVectorXf &old, int fps, double duration)
{
    int tar_num = int(fps * duration);
    tVectorXf new_vec = tVectorXf::Zero(tar_num);
    size_t old_size = old.size();
    double new_dt = 1.0 / (fps * 1.0);
    double old_dt = duration / (old.size() - 1);

    for (size_t i = 0; i < tar_num; i++)
    {
        // 1. find the gap
        double cur_time = i * new_dt;
        double cur_perc = i * 1.0 / tar_num;
        int old_st_idx = SIM_MIN(old.size() - 1, int(cur_perc * old.size()));
        size_t old_ed_idx = SIM_MIN(old.size() - 1, old_st_idx + 1);

        double st_perc = (cur_perc - old_st_idx * old_dt) / old_dt;
        double ed_perc = 1.0 - st_perc;
        new_vec[i] = (st_perc * old[old_st_idx] + ed_perc * old[old_ed_idx]);
        // std::cout << "new vec " << i << " = " << new_vec[i] << std::endl;
    }
    return new_vec;
}

void tDiscretedWave::ChangeFrequency(int tar_freq)
{
    data = Interpolate(data, tar_freq, GetDuration());
}
int tDiscretedWave::GetFrequency() const { return int(1.0 / this->dt); }

bool tDiscretedWave::LoadFromFile(std::string path)
{
    Json::Value root;
    if (false == cJsonUtil::LoadJson(path, root))
    {
        return false;
    }
    else
    {

        int freq = cJsonUtil::ParseAsInt("freq", root);
        dt = 1.0 / (freq * 1.0);
        duration = cJsonUtil::ParseAsDouble("duration", root);

        data = cJsonUtil::ReadVectorJson(cJsonUtil::ParseAsValue("data", root))
                   .cast<float>();
    }
    return true;
}

#include "utils/JsonUtil.h"

void tDiscretedWave::DumpToFile(std::string path)
{
    
    Json::Value root;
    root["freq"] = int(1.0 / this->dt);
    root["duration"] = this->duration;
    root["data"] = cJsonUtil::BuildVectorJson(data.cast<double>());

    cJsonUtil::WriteJson(path, root, true);
    printf("[log] dump to file %s\n", path);
}