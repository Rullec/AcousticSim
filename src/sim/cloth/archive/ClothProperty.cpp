#ifdef _WIN32
#include "ClothProperty.h"
#include "utils/JsonUtil.h"
#include "utils/LogUtil.h"
#include <iostream>
void tPhyProperty::Init(const Json::Value &root)
{
    Json::Value conf = cJsonUtil::ParseAsValue("cloth_property", root);
    // std::cout << "[debug] init physical prop = \n " << conf << std::endl;
    mStretchWarp = cJsonUtil::ParseAsDouble(mPropertiesName[0], conf);
    mStretchWeft = cJsonUtil::ParseAsDouble(mPropertiesName[1], conf);
    mStretchBias = cJsonUtil::ParseAsDouble(mPropertiesName[2], conf);
    mLinearBendingWarp = cJsonUtil::ParseAsDouble(mPropertiesName[3], conf);
    mLinearBendingWeft = cJsonUtil::ParseAsDouble(mPropertiesName[4], conf);
    mLinearBendingBias = cJsonUtil::ParseAsDouble(mPropertiesName[5], conf);

#ifdef SE_HAVE_NONLINEAR_PROPERTIES
    mNonlinearBendingWarp = cJsonUtil::ParseAsDouble(mPropertiesName[6], conf);
    mNonlinearBendingWeft = cJsonUtil::ParseAsDouble(mPropertiesName[7], conf);
    mNonlinearBendingBias = cJsonUtil::ParseAsDouble(mPropertiesName[8], conf);
#endif
    SIM_ASSERT(conf.size() == mNumOfProperties);
}
tVectorXd tPhyProperty::BuildFullFeatureVector() const
{
    tVectorXd feature = tVectorXd::Zero(mNumOfProperties);
    feature[0] = mStretchWarp;
    feature[1] = mStretchWeft;
    feature[2] = mStretchBias;
    feature[3] = mLinearBendingWarp;
    feature[4] = mLinearBendingWeft;
    feature[5] = mLinearBendingBias;
#ifdef SE_HAVE_NONLINEAR_PROPERTIES
    feature[6] = mNonlinearBendingWarp;
    feature[7] = mNonlinearBendingWeft;
    feature[8] = mNonlinearBendingBias;
#endif
    return feature;
}

/**
 * \brief           Given a name, return the feature value
 */
double tPhyProperty::GetFeature(std::string name) const
{
    if (name == "stretch_warp")
    {
        return mStretchWarp;
    }
    else if (name == "stretch_weft")
    {
        return mStretchWeft;
    }
    else if (name == "stretch_bias")
    {
        return mStretchBias;
    }
    else if (name == "bending_warp")
    {
        return mLinearBendingWarp;
    }
    else if (name == "bending_weft")
    {
        return mLinearBendingWeft;
    }
    else if (name == "bending_bias")
    {
        return mLinearBendingBias;
    }
#ifdef SE_HAVE_NONLINEAR_PROPERTIES
    else if (name == "bending_warp_nonlinear")
    {
        return mNonlinearBendingWarp;
    }
    else if (name == "bending_weft_nonlinear")
    {
        return mNonlinearBendingWeft;
    }
    else if (name == "bending_bias_nonlinear")
    {
        return mNonlinearBendingBias;
    }
#endif
    else
    {
        SIM_ERROR("unrecognized feature name {}", name);
        exit(0);
    }
    return 0;
}

/**
 * \brief           Given a name, get feature index
 */
int tPhyProperty::GetFeatureIdx(std::string name)
{
    if (name == "stretch_warp")
    {
        return 0;
    }
    else if (name == "stretch_weft")
    {
        return 1;
    }
    else if (name == "stretch_bias")
    {
        return 2;
    }
    else if (name == "bending_warp")
    {
        return 3;
    }
    else if (name == "bending_weft")
    {
        return 4;
    }
    else if (name == "bending_bias")
    {
        return 5;
    }
#ifdef SE_HAVE_NONLINEAR_PROPERTIES
    else if (name == "nonlinear_bending_warp")
    {
        return 6;
    }
    else if (name == "nonlinear_bending_weft")
    {
        return 7;
    }
    else if (name == "nonlinear_bending_bias")
    {
        return 8;
    }
#endif
    else
    {
        SIM_ERROR("unrecognized feature name {}", name);
        exit(0);
    }
    return -1;
}
/**
 * \brief           Given a full feature vector, load its value from vector to
 * discrete values
 */
void tPhyProperty::ReadFeatureVector(const tVectorXd &vec)
{
    SIM_ASSERT(mNumOfProperties == vec.size());
    mStretchWarp = vec[0];
    mStretchWeft = vec[1];
    mStretchBias = vec[2];
    mLinearBendingWarp = vec[3];
    mLinearBendingWeft = vec[4];
    mLinearBendingBias = vec[5];
#ifdef SE_HAVE_NONLINEAR_PROPERTIES
    mNonlinearBendingWarp = vec[6];
    mNonlinearBendingWeft = vec[7];
    mNonlinearBendingBias = vec[8];
#endif
}

/**
 * \brief           Given a json value and a bool vector, init the property
 */
void tBatchProperty::SetVisilibities(std::vector<bool> visibilities,
                                     const tVectorXi &visible_faeture_index)
{
    SIM_ASSERT(visibilities.size() == mNumOfProperties);
    for (int i = 0; i < mNumOfProperties; i++)
    {
        mVisibility[i] = visibilities[i];
    }

    mVisibleFeatureIndex = visible_faeture_index;
    SIM_ASSERT(mVisibleFeatureIndex.size() == mNumOfProperties);
    mNumOfVisibleFeature = 0;
    for (int i = 0; i < mNumOfProperties; i++)
    {
        if (mVisibleFeatureIndex[i] != -1)
        {
            mNumOfVisibleFeature++;
        }
    }

    // verify times
    tVectorXi times = tVectorXi::Zero(mNumOfVisibleFeature);
    for (int i = 0; i < mNumOfProperties; i++)
    {
        if (mVisibleFeatureIndex[i] != -1)
        {
            times[mVisibleFeatureIndex[i]] += 1;
        }
    }

    for (int i = 0; i < mNumOfVisibleFeature; i++)
    {
        SIM_ASSERT(times[i] == 1);
    }
}

/**
 * \brief           Build visible feature vector
 */
tVectorXd tPhyProperty::BuildVisibleFeatureVector() const
{
    SIM_ASSERT(false && "unsupported vec");
    exit(0);
    return tVectorXd::Zero(0);
}
tVectorXd tBatchProperty::BuildVisibleFeatureVector() const
{
    tVectorXd visible_feature = tVectorXd::Zero(mNumOfVisibleFeature);
    for (int i = 0; i < mNumOfProperties; i++)
    {
        int vis_id = mVisibleFeatureIndex[i];
        if (vis_id != -1)
        {
            visible_feature[vis_id] = GetFeature(mPropertiesName[i]);
        }
    }
    return visible_feature;
}
// extern std::string GetCurrentLinctexHash();

std::vector<std::pair<float, float>>
    stretch_map_guitosim = {{0.0f, 1e3f * STRETCH_SCALE},
                            {9.0f, 1e4f * STRETCH_SCALE},
                            {27.0f, 1e5f * STRETCH_SCALE},
                            {57.0f, 4e5f * STRETCH_SCALE},
                            {93.0f, 4e6f * STRETCH_SCALE},
                            {100.0f, 1e7f * STRETCH_SCALE}},
    bending_map_guitosim = {
        {0.0f, 0e0f * BENDING_SCALE},  {10.0f, 1e2f * BENDING_SCALE},
        {28.0f, 1e3f * BENDING_SCALE}, {48.0f, 3e3f * BENDING_SCALE},
        {65.0f, 2e4f * BENDING_SCALE}, {83.0f, 2e5f * BENDING_SCALE},
        {100.0f, 2e6f * BENDING_SCALE}};

double gNonlinearBendingScale = 0.1;
std::vector<std::pair<float, float>> nonlinear_bending_map_guitosim = {
    {0.0f, 0e0f * BENDING_SCALE *gNonlinearBendingScale},
    {10.0f, 1e2f * BENDING_SCALE *gNonlinearBendingScale},
    {28.0f, 1e3f * BENDING_SCALE *gNonlinearBendingScale},
    {48.0f, 3e3f * BENDING_SCALE *gNonlinearBendingScale},
    {65.0f, 2e4f * BENDING_SCALE *gNonlinearBendingScale},
    {83.0f, 2e5f * BENDING_SCALE *gNonlinearBendingScale},
    {100.0f, 2e6f * BENDING_SCALE *gNonlinearBendingScale}};

float CalculateSimValueFromGUI(
    const std::vector<std::pair<float, float>> &arrays, float gui_value)
{
    // std::cout << "bending_scale = " << bending_scale << std::endl;
    // std::cout << "stretch_scale = " << stretch_scale << std::endl;
    // exit(1);
    // 1. exceed the boundary
    int last_id = arrays.size() - 1;
    float threshold = 1e-6;
    if (gui_value >= arrays[last_id].first - threshold)
    {
        return arrays[last_id].second;
    }
    if (gui_value <= arrays[0].first + threshold)
    {
        return arrays[0].second;
    }

    // 2. find the interval
    for (int st = 0; st < arrays.size() - 1; st++)
    {
        float st_key = arrays[st].first;
        float ed_key = arrays[st + 1].first;
        float st_value = arrays[st].second;
        float ed_value = arrays[st + 1].second;
        float gap = ed_key - st_key;
        if ((st_key <= gui_value) && (ed_key >= gui_value))
        {
            // do interplotion and return
            return (1.0 - (gui_value - st_key) / gap) * st_value +
                   (1.0 - (ed_key - gui_value) / gap) * ed_value;
        }
    }
    return std::nan("");
}

float CalculateGUIValueFromSim(
    const std::vector<std::pair<float, float>> &arrays, float sim_value)
{
    // 1. exceed the boundary
    int last_id = arrays.size() - 1;
    float threshold = 1e-6;
    if (sim_value < arrays[0].second)
    {
        return arrays[0].first;
    }
    if (arrays[last_id].second < sim_value)
    {
        return arrays[last_id].first;
    }

    // 2. find the interval
    for (int st = 0; st < arrays.size() - 1; st++)
    {
        float st_key = arrays[st].second, ed_key = arrays[st + 1].second,
              st_value = arrays[st].first, ed_value = arrays[st + 1].first;
        float gap = ed_key - st_key;
        if ((st_key <= sim_value) && (ed_key >= sim_value))
        {
            return (1.0 - (sim_value - st_key) / gap) * st_value +
                   (1.0 - (ed_key - sim_value) / gap) * ed_value;
        }
    }
    return std::nan("");
}
/**
 * \brief           convert gui bending coef to simulation coef
 */
float tPhyProperty::ConvertLinearBendingCoefFromGUIToSim(float gui_value)
{
    // 1. judge illegal
    float val = CalculateSimValueFromGUI(bending_map_guitosim, gui_value);
    // {
    // printf("[linear_bending] gui value %f, sim value %f\n", gui_value, val);
    // }
    // std::cout << "convert bending gui " << gui_value << " to sim " << val
    //           << std::endl;
    // std::cout << " BENDING_SCALE = " << BENDING_SCALE << std::endl;
    return val;
}
/**
 * \brief           convert simulation bending coef to gui coef
 */
float tPhyProperty::ConvertLinearBendingCoefFromSimToGUI(float sim_value)
{
    float val = CalculateGUIValueFromSim(bending_map_guitosim, sim_value);
    // std::cout << "convert bending sim " << sim_value << " to gui " << val
    //           << std::endl;
    // std::cout << " BENDING_SCALE = " << BENDING_SCALE << std::endl;
    return val;
}
/**
 * \brief           convert stretch gui coef to simulation coef
 */
float tPhyProperty::ConvertStretchCoefFromGUIToSim(float gui_value)
{
    float val = CalculateSimValueFromGUI(stretch_map_guitosim, gui_value);
    // std::cout << "convert stretch gui " << gui_value << " to sim " << val
    //           << std::endl;
    // std::cout << " STRETCH_SCALE = " << STRETCH_SCALE << std::endl;
    return val;
}
/**
 * \brief           convert stretch simulation coef to gui coef
 */
float tPhyProperty::ConvertStretchCoefFromSimToGUI(float sim_value)
{
    float val = CalculateGUIValueFromSim(stretch_map_guitosim, sim_value);
    // std::cout << "convert stretch sim " << sim_value << " to gui " << val
    //           << std::endl;
    // std::cout << " STRETCH_SCALE = " << STRETCH_SCALE << std::endl;
    return val;
}

float &tPhyProperty::GetStretchWarpRef() { return mStretchWarp; }
float &tPhyProperty::GetStretchWeftRef() { return mStretchWeft; }
float &tPhyProperty::GetStretchBiasRef() { return mStretchBias; }

float tPhyProperty::GetLinearBendingWarp() { return mLinearBendingWarp; }
float tPhyProperty::GetLinearBendingWeft() { return mLinearBendingWeft; }
float tPhyProperty::GetLinearBendingBias() { return mLinearBendingBias; }

float tPhyProperty::GetStretchWarp() { return mStretchWarp; }
float tPhyProperty::GetStretchWeft() { return mStretchWeft; }
float tPhyProperty::GetStretchBias() { return mStretchBias; }

void tPhyProperty::SetLinearBendingWarp(float val) { mLinearBendingWarp = val; }
void tPhyProperty::SetLinearBendingWeft(float val) { mLinearBendingWeft = val; }
void tPhyProperty::SetLinearBendingBias(float val) { mLinearBendingBias = val; }
void tPhyProperty::SetStretchWarp(float val) { mStretchWarp = val; }
void tPhyProperty::SetStretchWeft(float val) { mStretchWeft = val; }
void tPhyProperty::SetStretchBias(float val) { mStretchBias = val; }

float &tPhyProperty::GetLinearBendingWarpRef() { return mLinearBendingWarp; }
float &tPhyProperty::GetLinearBendingWeftRef() { return mLinearBendingWeft; }
float &tPhyProperty::GetLinearBendingBiasRef() { return mLinearBendingBias; }

#ifdef SE_HAVE_NONLINEAR_PROPERTIES
float &tPhyProperty::GetNonLinearBendingWarpRef()
{
    return mNonlinearBendingWarp;
}
float &tPhyProperty::GetNonLinearBendingWeftRef()
{
    return mNonlinearBendingWeft;
}
float &tPhyProperty::GetNonLinearBendingBiasRef()
{
    return mNonlinearBendingBias;
}
float tPhyProperty::GetNonLinearBendingWarp() { return mNonlinearBendingWarp; }
float tPhyProperty::GetNonLinearBendingWeft() { return mNonlinearBendingWeft; }
float tPhyProperty::GetNonLinearBendingBias() { return mNonlinearBendingBias; }

// set value
void tPhyProperty::SetNonLinearBendingWarp(float val)
{
    mNonlinearBendingWarp = val;
}
void tPhyProperty::SetNonLinearBendingWeft(float val)
{
    mNonlinearBendingWeft = val;
}
void tPhyProperty::SetNonLinearBendingBias(float val)
{
    mNonlinearBendingBias = val;
}

float tPhyProperty::ConvertNonLinearBendingCoefFromGUIToSim(float gui_value)
{
    // 1. get the abs value
    float abs_gui_value = std::fabs(gui_value);
    float sign = gui_value > 0 ? 1 : -1;
    // 2. do the same convertion with the linear item
    float sim_value = CalculateSimValueFromGUI(nonlinear_bending_map_guitosim,
                                               abs_gui_value) *
                      sign;
    // 3. apply the sign
    // if (gui_value < 0)
    // {
    // printf("[nonlinear_bending] gui value %f, abs gui value %f, sign value "
    //        "%f, sim value %f\n",
    //        gui_value, abs_gui_value, sign, sim_value);
    // }
    return sim_value;
}
float tPhyProperty::ConvertNonLinearBendingCoefFromSimToGUI(float sim_value)
{
    // 1. get the abs value
    float abs_sim_value = std::fabs(sim_value);
    float sign = sim_value > 0 ? 1 : -1;

    // 2. do the same convertion with the linear item
    float gui_value = CalculateGUIValueFromSim(nonlinear_bending_map_guitosim,
                                               abs_sim_value) *
                      sign;
    // 3. apply the sign
    return gui_value;
}
#endif

/**
 * \brief       Given scale factor (for bending stiffnes [SI], instead of [GUI])
 *          Given base gui warp, weft, bias
 * 1. calculate the base sim [SI] value
 * 2. calculate the new sim [SI] value, convert to sim [GUI]
 * 3. set the new sim [GUI] value
 */
void tPhyProperty::ScaleAndSetTheLinearBendingStiffness(
    float scale_factor, float bending_linear_gui_warp,
    float bending_linear_gui_weft, float bending_linear_gui_bias)
{
    SIM_ASSERT(scale_factor > 0);
    // 1. step 1 calc the sim value from
    // printf("[scale] cur bending gui %.1f, %.1f, %.1f\n", mLinearBendingWarp,
    //        mLinearBendingWeft, mLinearBendingBias);
    float base_warp_si = tPhyProperty::ConvertLinearBendingCoefFromGUIToSim(
              bending_linear_gui_warp),
          base_weft_si = tPhyProperty::ConvertLinearBendingCoefFromGUIToSim(
              bending_linear_gui_weft),
          base_bias_si = tPhyProperty::ConvertLinearBendingCoefFromGUIToSim(
              bending_linear_gui_bias);

    // 2. calculate new sim value, convert to GUI
    base_warp_si *= scale_factor;
    base_weft_si *= scale_factor;
    base_bias_si *= scale_factor;

    float new_warp_gui =
        tPhyProperty::ConvertLinearBendingCoefFromSimToGUI(base_warp_si);
    float new_weft_gui =
        tPhyProperty::ConvertLinearBendingCoefFromSimToGUI(base_weft_si);
    float new_bias_gui =
        tPhyProperty::ConvertLinearBendingCoefFromSimToGUI(base_bias_si);

    // 3. set new value
    SetLinearBendingWarp(new_warp_gui);
    SetLinearBendingWeft(new_weft_gui);
    SetLinearBendingBias(new_bias_gui);

    // printf("[scale] new bending gui %.1f, %.1f, %.1f\n", mLinearBendingWarp,
    //        mLinearBendingWeft, mLinearBendingBias);
}
#endif
