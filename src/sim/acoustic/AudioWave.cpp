#include "AudioWave.h"
#include "utils/JsonUtil.h"
#include "utils/LogUtil.h"
tDiscretedWave::tDiscretedWave(float dt_)
{
    dt = dt_;
    data.resize(1);
    duration = dt;
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

tVectorXf Interpolate(const tVectorXf &old_data, double old_dt, double new_dt)
{
    // 1. get total duration
    double duration = old_data.size() * old_dt;
    // 2. get new data size
    int num_of_data_new = duration / new_dt;
    tVectorXf data = tVectorXf::Zero(num_of_data_new);
    for (int i = 0; i < num_of_data_new; i++)
    {
        // 1. get interval in old data
        int old_data_st = -1, old_data_ed = -1;
        double new_perc = i * 1.0 / num_of_data_new;
        {

            old_data_st = int(new_perc * old_data.size());
            old_data_ed = SIM_MIN(old_data_st + 1, old_data.size() - 1);
        }

        // 2. interpolate inside the interval
        {
            double old_st_perc = 1.0 * old_data_st / old_data.size();
            double old_ed_perc = 1.0 * old_data_ed / old_data.size();
            double gap = old_ed_perc - old_st_perc;
            double ed_weight = (new_perc - old_st_perc) / gap;
            double st_weight = (old_ed_perc - new_perc) / gap;
            data[i] = st_weight * old_data[old_data_st] +
                      ed_weight * old_data[old_data_ed];
        }
    }
    return data;
}

void tDiscretedWave::ChangeFrequency(int tar_freq)
{
    double new_dt = 1.0 / (tar_freq * 1.0);
    data = Interpolate(data, dt, new_dt);

    // new dt
    dt = new_dt;
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
    printf("[log] dump to file %s\n", path.c_str());
}
#include "AudioFile.h"

bool tDiscretedWave::LoadFromWAV(std::string path)
{
    AudioFile<double> audioFile;
    audioFile.load(path);
    audioFile.printSummary();
    uint32_t sampling_rate = audioFile.getSampleRate();
    int num_of_channels = audioFile.getNumChannels();
    int num_of_samples = audioFile.getNumSamplesPerChannel();
    duration = audioFile.getLengthInSeconds();
    dt = 1.0 / (sampling_rate * 1.0);
    data.resize(num_of_samples);
    for (int i = 0; i < num_of_samples; i++)
        data[i] = audioFile.samples[0][i];
    return true;
}
#include "utils/FileUtil.h"
void tDiscretedWave::DumpToWAV(std::string path)
{
    std::string tar_dir = cFileUtil::GetDir(path);
    if (cFileUtil::ExistsDir(tar_dir) == false)
    {
        cFileUtil::CreateDir(tar_dir.c_str());
    }
    // cFileUtil::CreateDir( .c_str());
    AudioFile<double> audioFile;
    uint32_t rate = int(1.0 / dt);
    audioFile.setNumChannels(1);
    audioFile.setSampleRate(rate);

    // set info
    audioFile.setAudioBufferSize(1, data.size());

    // set data
    for (int i = 0; i < data.size(); i++)
    {
        //
        audioFile.samples[0][i] = data[i];
    }

    // set result
    audioFile.save(path, AudioFileFormat::Wave);
    std::cout << "[log] save wav to " << path << std::endl;
}

tAnalyticWave::tAnalyticWave() { Clear(); }
void tAnalyticWave::Clear()
{
    mStrength.clear();
    mFreqHZ.clear();
}
void tAnalyticWave::AddWave(double strength, double freq)
{
    mStrength.push_back(strength);
    mFreqHZ.push_back(freq);
}
double tAnalyticWave::Evaluate(double t) const
{
    double total_amp = 0;
    for (int i = 0; i < mStrength.size(); i++)
    {
        double amp = mStrength[i];
        double freq_hz = mFreqHZ[i];
        total_amp += amp * std::cos(2 * M_PI * freq_hz * t);
    }
    return total_amp;
}

tDiscretedWavePtr DiscretizeWave(tAnalyticWavePtr ana_wave, double duration,
                                 double dt)
{
    tDiscretedWavePtr d_wave = std::make_shared<tDiscretedWave>(dt);
    // create data
    int num_of_samples = int(duration / dt);
    tVectorXf data = tVectorXf::Zero(num_of_samples);

    for (int i = 0; i < num_of_samples; i++)
    {
        data[i] = ana_wave->Evaluate(i * dt);
    }
    d_wave->SetData(data);
    return d_wave;
}
