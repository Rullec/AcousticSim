#ifdef _WIN32
#pragma once
// #include "SePhysicalProperties.h"
#include "utils/MathUtil.h"

namespace Json
{
class Value;
};

class tPhyProperty
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    virtual void Init(const Json::Value &conf); // normal init
    virtual tVectorXd BuildFullFeatureVector() const;
    virtual tVectorXd BuildVisibleFeatureVector() const;
    virtual void ReadFeatureVector(const tVectorXd &vec);
    static int GetFeatureIdx(std::string name);
    double GetFeature(std::string name) const;
    static float ConvertLinearBendingCoefFromGUIToSim(float gui_value);
    static float ConvertLinearBendingCoefFromSimToGUI(float sim_value);
    static float ConvertStretchCoefFromGUIToSim(float gui_value);
    static float ConvertStretchCoefFromSimToGUI(float sim_value);

    // get vlaue
    float GetLinearBendingWarp();
    float GetLinearBendingWeft();
    float GetLinearBendingBias();
    float GetStretchWarp();
    float GetStretchWeft();
    float GetStretchBias();

    // get value ref (for imgui purpose)
    float &GetLinearBendingWarpRef();
    float &GetLinearBendingWeftRef();
    float &GetLinearBendingBiasRef();
    float &GetStretchWarpRef();
    float &GetStretchWeftRef();
    float &GetStretchBiasRef();

    // set value
    void SetLinearBendingWarp(float val);
    void SetLinearBendingWeft(float val);
    void SetLinearBendingBias(float val);
    void SetStretchWarp(float val);
    void SetStretchWeft(float val);
    void SetStretchBias(float val);

    //
    void ScaleAndSetTheLinearBendingStiffness(float scale_factor,
                                              float bending_linear_gui_warp,
                                              float bending_linear_gui_weft,
                                              float bending_linear_gui_bias);
protected:
    float mStretchWarp;       // GUI value
    float mStretchWeft;       // GUI value
    float mStretchBias;       // GUI value
    float mLinearBendingWarp; // GUI value
    float mLinearBendingWeft; // GUI value
    float mLinearBendingBias; // GUI value
    
#ifdef SE_HAVE_NONLINEAR_PROPERTIES
public:
    inline static const int mNumOfProperties = 9;
    inline static const std::string
        mPropertiesName[tPhyProperty::mNumOfProperties] = {
            "stretch_warp",           "stretch_weft",
            "stretch_bias",           "bending_warp",
            "bending_weft",           "bending_bias",
            "bending_warp_nonlinear", "bending_weft_nonlinear",
            "bending_bias_nonlinear"};

    float &GetNonLinearBendingWarpRef();
    float &GetNonLinearBendingWeftRef();
    float &GetNonLinearBendingBiasRef();
    float GetNonLinearBendingWarp();
    float GetNonLinearBendingWeft();
    float GetNonLinearBendingBias();

    // set value
    void SetNonLinearBendingWarp(float val);
    void SetNonLinearBendingWeft(float val);
    void SetNonLinearBendingBias(float val);
    // convert method
    static float ConvertNonLinearBendingCoefFromGUIToSim(float gui_value);
    static float ConvertNonLinearBendingCoefFromSimToGUI(float sim_value);

protected:
    float mNonlinearBendingWarp;
    float mNonlinearBendingWeft;
    float mNonlinearBendingBias;

#else
public:
    inline static const int mNumOfProperties = 6;
    inline static const std::string
        mPropertiesName[tPhyProperty::mNumOfProperties] = {
            "stretch_warp", "stretch_weft", "stretch_bias",
            "bending_warp", "bending_weft", "bending_bias"};
#endif
};

struct tBatchProperty : public tPhyProperty
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    bool mVisibility[mNumOfProperties] = {};
    tVectorXi mVisibleFeatureIndex;
    int mNumOfVisibleFeature;
    void SetVisilibities(std::vector<bool> visibilities,
                         const tVectorXi &visible_faeture_index);
    virtual tVectorXd BuildVisibleFeatureVector() const override;
};
#endif