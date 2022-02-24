#ifdef _WIN32
#pragma once
#include "utils/DefUtil.h"
#include "utils/MathUtil.h"
#include <utility>
namespace Json
{
class Value;
};
SIM_DECLARE_CLASS_AND_PTR(tPhyProperty);
class tPhyPropertyManager
{
public:
    inline static const std::string PROPERTY_SAMPLES_PATH_KEY =
                                        "property_samples_path",
                                    PROPERTY_SAMPLES_START_ID_KEY =
                                        "property_samples_start_id";
    explicit tPhyPropertyManager(const Json::Value &conf);
    tPhyPropertyPtr GetProperty(int idx);
    int GetNumOfProperties() const;
    int GetCurrentPropertyStartId() const;

protected:
    tMatrixXd mAllPropertyFeatures;
    std::string mPropertyPath;
    int mPropertyStartId;
    std::vector<bool> mVisibilities;
    tVectorXi mVisibleIndex;
   
    void InitFeatures();
};
#endif