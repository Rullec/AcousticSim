#include "ClothNoise.h"
#include <iostream>
/**
 * \brief       init the data augmentation strucutre
 */
#define _USE_MATH_DEFINES
#include <cmath>
tClothNoise::tClothNoise(const Json::Value &conf)
{
    mEnableNoise = cJsonUtil::ParseAsBool("enable_noise", conf);
    mNumOfSamplesPerProp = cJsonUtil::ParseAsInt("samples_per_prop", conf);
    // mEnableLowFreqNoise =
    //     cJsonUtil::ParseAsDouble("enable_low_freq_noise", conf);
    mMaxFoldAmp = cJsonUtil::ParseAsDouble("max_fold_amp", conf);
    mMinFoldNum = cJsonUtil::ParseAsInt("min_fold_num", conf);
    mMaxFoldNum = cJsonUtil::ParseAsInt("max_fold_num", conf);

    // begin to init the probs
    mSpanningMeshDataDir =
        cJsonUtil::ParseAsString("feature_mesh_datadir", conf);
    // std::cout << "mSpanningMeshDataDir = " << mSpanningMeshDataDir <<
    // std::endl; exit(1);
    Json::Value prob_conf = cJsonUtil::ParseAsValue("noise_prob", conf);
    float sum_prob = 0;

    for (int i = 0; i < NUM_NOISE_TYPE; i++)
    {
        mClothProb[i] = cJsonUtil::ParseAsDouble(
            BuildClothNoiseNameFromType(static_cast<eClothNoiseType>(i)),
            prob_conf);
        sum_prob += mClothProb[i];
        std::cout << "[log] init noise for type "
                  << BuildClothNoiseNameFromType(
                         static_cast<eClothNoiseType>(i))
                  << " prob = " << mClothProb[i] << std::endl;
    }
    for (int i = 0; i < NUM_NOISE_TYPE; i++)
        mClothProb[i] /= sum_prob;

    InitSpannedData();
    Reset();
}

static std::string gClothNoiseStr[eClothNoiseType::NUM_NOISE_TYPE] = {
    "wave_noise", "manual_noise", "spanning_noise", "empty_noise"};
/**
 * \brief       given the cloth noise type, return the namae (string)
 */
std::string BuildClothNoiseNameFromType(eClothNoiseType type)
{
    return gClothNoiseStr[static_cast<int>(type)];
}

/**
 * \brief       given the cloth noise name, return the type
 */
eClothNoiseType BuildClothNoiseTypeFromName(std::string name)
{
    for (int i = 0; i < eClothNoiseType::NUM_NOISE_TYPE; i++)
    {
        if (gClothNoiseStr[i] == name)
        {
            return static_cast<eClothNoiseType>(i);
        }
    }
    SIM_ERROR("fail to determine the cloth noise type for string {}", name);
    exit(1);
    return eClothNoiseType::NUM_NOISE_TYPE;
}

/**
 * \brief           Cloth Noise
 */
eClothNoiseType tClothNoise::RandomTheNoise() const
{
    std::vector<double> probs = {};
    for (int i = 0; i < eClothNoiseType::NUM_NOISE_TYPE; i++)
    {

        probs.push_back(mClothProb[i]);
    }
    int cur_selection = cMathUtil::RandIntCategorical(probs);
    std::cout << "[debug] select cloth noise type: "
              << BuildClothNoiseNameFromType(
                     static_cast<eClothNoiseType>(cur_selection))
              << std::endl;
    return static_cast<eClothNoiseType>(cur_selection);
}

/**
 * \brief           Test random function
 */
void tClothNoise::TestRandom() const
{
    int times[eClothNoiseType::NUM_NOISE_TYPE] = {};
    int iters = 10000;
    for (int i = 0; i < iters; i++)
    {
        times[static_cast<int>(RandomTheNoise())] += 1;
    }

    std::cout << "times[0] = " << times[0]
              << " ideal = " << iters * mClothProb[0] << std::endl;
    std::cout << "times[1] = " << times[1]
              << " ideal = " << iters * mClothProb[1] << std::endl;
    std::cout << "times[2] = " << times[2]
              << " ideal = " << iters * mClothProb[2] << std::endl;
    exit(1);
}

bool tClothNoise::GetEnableNoise() const { return this->mEnableNoise; }

int tClothNoise::GetNumOfSamplesPerProp() const
{
    return this->mNumOfSamplesPerProp;
}
void tClothNoise::GetWaveNoiseParam(double &_mMaxFoldAmp, int &_mMinFoldNum,
                                    int &_mMaxFoldNum) const
{
    _mMaxFoldAmp = mMaxFoldAmp;
    _mMinFoldNum = mMinFoldNum;
    _mMaxFoldNum = mMaxFoldNum;
}
#include "geometries/ObjExport.h"
#include "sim/cloth/LinctexCloth.h"
/**
 * \brief           apply spanning noise (load mesh data, and set pose)
 */
void tClothNoise::ApplySpanningNoise(cLinctexClothPtr cloth_ptr)
{
    // 1. get an index which hasn't been visited
    // std::cout << "---begin to apply spanning noise---\n";
    int num_cur_visited = GetNumCurrentVisited();
    int cur_selected = -1;
    if (num_cur_visited >= mSpanningMeshData.size())
    {
        cur_selected = cMathUtil::RandInt(0, mSpanningMeshData.size());
        SIM_WARN("cur visited {} >= total spanning mesh data {}",
                 num_cur_visited, mSpanningMeshData.size());
        // exit(1);
    }
    else
    {
        cur_selected = SelectUnvisitedMeshRandom();
    }
    std::cout << "[debug] now unvisited "
              << mSpanningMeshData.size() - GetNumCurrentVisited()
              << " total num " << mSpanningMeshData.size() << " selected  "
              << cur_selected << std::endl;
    mSpannedVisited[cur_selected] = true;

    // 2. begin to load mesh data, and set pos
    std::string mesh_data_path = mSpanningMeshData[cur_selected];
    // std::cout << "[debug] begin to load " << mesh_data_path << std::endl;
    Json::Value mesh_json_root;
    if (false == cJsonUtil::LoadJson(mesh_data_path, mesh_json_root))
    {
        SIM_ERROR("fail to load the mesh data from {}", mesh_data_path);
    }
    tVectorXd new_mesh_pos = cJsonUtil::ReadVectorJson(mesh_json_root["input"]);
    if (cloth_ptr->GetPos().size() != new_mesh_pos.size())
    {
        SIM_ERROR("loaded mesh length {} , but current length is {}",
                  new_mesh_pos.size(), cloth_ptr->GetPos().size());
    }

    cloth_ptr->SetPos(new_mesh_pos);
    // cObjExporter::ExportObj("tmp.obj", cloth_ptr->GetVertexArray(),
    //                         cloth_ptr->GetTriangleArray(), false);
    // std::cout << "[debug] set mesh pos succ!\n";
}

/**
 * \brief           Reset cloth noise
 */
void tClothNoise::Reset()
{
    // reset the final result
    mSpannedVisited.resize(mSpanningMeshData.size(), 0);
    for (auto &x : mSpannedVisited)
        x = false;
    SIM_WARN("cloth noise reset, visited set to zero! current visited num "
             "should be {} = 0",
             this->GetNumCurrentVisited());
}

/**
 * \brief           init spanned data given cloth
 */
#include "utils/FileUtil.h"
void tClothNoise::InitSpannedData()
{
    SIM_ASSERT(cFileUtil::ExistsDir(mSpanningMeshDataDir) == true);
    std::cout << "init spanned data from data dir " << mSpanningMeshDataDir
              << std::endl;
    mSpanningMeshData.clear();
    for (auto &file : cFileUtil::ListDir(this->mSpanningMeshDataDir))
    {
        mSpanningMeshData.push_back(file);
    }
}

/**
 * \brief           Get num of current visited elements
 */
int tClothNoise::GetNumCurrentVisited()
{
    int times = 0;
    for (auto &x : mSpannedVisited)
    {
        times += (x == true);
    }
    return times;
}

/**
 * \brief           random select a unvisited mesh
 */
int tClothNoise::SelectUnvisitedMeshRandom() const
{
    std::vector<int> unvisited = {};
    for (int idx = 0; idx < mSpanningMeshData.size(); idx++)
    {
        if (mSpannedVisited[idx] == false)
        {
            unvisited.push_back(idx);
        }
    }
    if (unvisited.size() == 0)
    {
        SIM_ERROR("We have no unvisited mesh! total mesh num {}",
                  mSpanningMeshData.size());
    }

    int selected_idx = unvisited[cMathUtil::RandInt(0, unvisited.size())];
    SIM_ASSERT(mSpannedVisited[selected_idx] == false);
    return selected_idx;
}