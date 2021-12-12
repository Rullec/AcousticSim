#pragma once

#include "LogUtil.h"
#include "MathUtil.h"
#include "json/json.h"
#include <string>

namespace Json
{
class Value;
};
class cJsonUtil
{
public:
    // static Json::Value BuildVectorJson(const tVector &vec);
    static tVectorXd ReadVectorJson(const Json::Value &root);
    static Json::Value BuildVectorJson(const Eigen::VectorXd &vec);
    static std::string BuildVectorString(const Eigen::VectorXd &vec);
    static bool ReadVectorJson(const Json::Value &root,
                               Eigen::VectorXd &out_vec);
    static bool ReadVectorJson(const Json::Value &root, tVector &out_vec);
    static bool ReadMatrixJson(const Json::Value &root, tMatrixXd &out_mat);
    static bool LoadJson(const std::string &path, Json::Value &value);
    static bool WriteJson(const std::string &path, Json::Value &value,
                          bool indent = true);

    static int ParseAsInt(const std::string &data_field_name,
                          const Json::Value &root);
    static std::string ParseAsString(const std::string &data_field_name,
                                     const Json::Value &root);
    static double ParseAsDouble(const std::string &data_field_name,
                                const Json::Value &root);
    static float ParseAsFloat(const std::string &data_field_name,
                              const Json::Value &root);
    static bool ParseAsBool(const std::string &data_field_name,
                            const Json::Value &root);
    static Json::Value ParseAsValue(const std::string &data_field_name,
                                    const Json::Value &root);
    static bool HasValue(const std::string & name, const Json::Value & root);
};
