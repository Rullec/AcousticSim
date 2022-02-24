#ifdef _WIN32
#include "ClothPropertyManager.h"
#include "ClothProperty.h"
#include "utils/FileUtil.h"
#include "utils/JsonUtil.h"
#include <iostream>

tPhyPropertyManager::tPhyPropertyManager(const Json::Value &conf)
{
    mPropertyPath = cJsonUtil::ParseAsString(PROPERTY_SAMPLES_PATH_KEY, conf);

    Json::Value props = cJsonUtil::ParseAsValue("properties", conf);
    SIM_ASSERT(props.size() == tPhyProperty::mNumOfProperties);
    mVisibilities.resize(tPhyProperty::mNumOfProperties);
    mVisibleIndex = tVectorXi::Ones(tPhyProperty::mNumOfProperties) * -1;
    for (int i = 0; i < tPhyProperty::mNumOfProperties; i++)
    {
        Json::Value sub_value =
            cJsonUtil::ParseAsValue(tPhyProperty::mPropertiesName[i], props);
        mVisibilities[i] = cJsonUtil::ParseAsBool("visible", sub_value);
        if (mVisibilities[i] == true)
        {
            mVisibleIndex[i] = cJsonUtil::ParseAsInt("index", sub_value);
        }
    }
    InitFeatures();
}

tPhyPropertyPtr tPhyPropertyManager::GetProperty(int i)
{
    auto ptr = std::make_shared<tBatchProperty>();
    ptr->ReadFeatureVector(mAllPropertyFeatures.row(i));
    ptr->SetVisilibities(this->mVisibilities, this->mVisibleIndex);

    return ptr;
}

tVectorXd vec2eigen(std::vector<double> res)
{
    tVectorXd vec = tVectorXd::Zero(res.size());
    for (int i = 0; i < res.size(); i++)
    {
        vec[i] = res[i];
    }
    return vec;
}

int tPhyPropertyManager::GetNumOfProperties() const
{
    return mAllPropertyFeatures.rows();
}

/**
 * \brief               init properties from reading a given config
 */
void tPhyPropertyManager::InitFeatures()
{
    if (cFileUtil::ExistsFile(mPropertyPath) == false)
    {
        SIM_ERROR("mPropertyPath {} doesn't exist", mPropertyPath);
        exit(1);
    }
    Json::Value root;
    SIM_ASSERT(cJsonUtil::LoadJson(mPropertyPath, root) == true);

    // 1. validate the property names
    auto prop_name_lst = cJsonUtil::ParseAsValue("prop_name_list", root);
    SIM_ASSERT(prop_name_lst.size() == tPhyProperty::mNumOfProperties);
    for (int id = 0; id < prop_name_lst.size(); id++)
    {
        SIM_ASSERT(tPhyProperty::mPropertiesName[id] ==
                   prop_name_lst[id].asString());
    }

    // begin to load property names
    auto prop_samples = cJsonUtil::ParseAsValue("prop_samples", root);
    int num_of_samples = prop_samples.size();
    SIM_ASSERT(
        cJsonUtil::ReadMatrixJson(prop_samples, this->mAllPropertyFeatures));
    mPropertyStartId =
        cJsonUtil::ParseAsInt(PROPERTY_SAMPLES_START_ID_KEY, root);
    printf("[log] load %d properties from %s\n", num_of_samples,
           mPropertyPath.c_str());
    std::cout << "[log] begin feature = " << mAllPropertyFeatures.row(0)
              << std::endl;
    std::cout << "[log] end feature = "
              << mAllPropertyFeatures.row(mAllPropertyFeatures.rows() - 1)
              << std::endl;
    for (int i = 0; i < mAllPropertyFeatures.rows(); i++)
    {
        float linear_warp = mAllPropertyFeatures(i, 0);
        float linear_weft = mAllPropertyFeatures(i, 1);
        if (linear_warp < linear_weft)
        {
            SIM_ERROR("linear warp {} < linear weft, illegal, please check the "
                      "data {}",
                      linear_warp, linear_weft, mPropertyPath);
            exit(1);
        }
    }
}

int tPhyPropertyManager::GetCurrentPropertyStartId() const
{
    return mPropertyStartId;
}

#endif