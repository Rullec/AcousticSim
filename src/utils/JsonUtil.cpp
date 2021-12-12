#include "JsonUtil.h"
#include "FileUtil.h"
#include <fstream>
#include <iostream>
#include <memory>

// tLogger cJsonUtil::mLogger = cLogUtil::CreateLogger("cJsonUtil");
Json::Value cJsonUtil::BuildVectorJson(const tVectorXd &vec)
{
    Json::Value json = Json::arrayValue;

    for (int i = 0; i < vec.size(); ++i)
    {
        json.append(vec[i]);
    }
    return json;
}

bool cJsonUtil::ReadVectorJson(const Json::Value &root, tVector &out_vec)
{
    bool succ = false;
    int num_vals = root.size();
    assert(num_vals <= 4);
    num_vals = std::min(num_vals, static_cast<int>(out_vec.size()));

    if (root.isArray())
    {
        out_vec.setZero();
        for (int i = 0; i < num_vals; ++i)
        {
            Json::Value json_elem = root.get(i, 0);
            out_vec[i] = json_elem.asDouble();
        }
        succ = true;
    }

    return succ;
}

// std::string cJsonUtil::BuildVectorJson(const Eigen::VectorXd &vec)
// {
//     std::string json = BuildVectorString(vec);
//     json = "[" + json + "]";
//     return json;
// }

std::string cJsonUtil::BuildVectorString(const Eigen::VectorXd &vec)
{
    std::string str = "";
    char str_buffer[32];
    for (int i = 0; i < vec.size(); ++i)
    {
        if (i != 0)
        {
            str += ",";
        }
        sprintf(str_buffer, "%20.10f", vec[i]);
        str += std::string(str_buffer);
    }
    return str;
}

tVectorXd cJsonUtil::ReadVectorJson(const Json::Value &root)
{
    tVectorXd xd;
    cJsonUtil::ReadVectorJson(root, xd);
    return xd;
}
bool cJsonUtil::ReadVectorJson(const Json::Value &root,
                               Eigen::VectorXd &out_vec)
{
    bool succ = false;
    int num_vals = root.size();

    if (root.isArray())
    {
        out_vec.resize(num_vals);
        for (int i = 0; i < num_vals; ++i)
        {
            Json::Value json_elem = root.get(i, 0);
            out_vec[i] = json_elem.asDouble();
        }
        succ = true;
    }

    return succ;
}

bool cJsonUtil::LoadJson(const std::string &path, Json::Value &value)
{
    // cFileUtil::AddLock(path);
    // std::cout <<"parsing " << path << " begin \n";
    std::ifstream fin(path);
    if (fin.fail() == true)
    {
        std::cout << "[error] cJsonUtil::LoadJson file " << path
                  << " doesn't exist\n";
        return false;
    }
    Json::CharReaderBuilder rbuilder;
    std::string errs;
    bool parsingSuccessful =
        Json::parseFromStream(rbuilder, fin, &value, &errs);
    if (!parsingSuccessful)
    {
        // report to the user the failure and their locations in the
        // document.
        std::cout << "[error] cJsonUtil::LoadJson: Failed to parse json\n"
                  << errs << ", file path " << path << std::endl;
        return false;
    }
    // std::cout <<"parsing " << path << " end \n";
    // cFileUtil::DeleteLock(path);
    return true;
}

bool cJsonUtil::WriteJson(const std::string &path, Json::Value &value,
                          bool indent /* = true*/)
{
    // cFileUtil::AddLock(path);
    Json::StreamWriterBuilder builder;
    if (indent == false)
        builder.settings_["indentation"] = "";
    std::unique_ptr<Json::StreamWriter> writer(builder.newStreamWriter());
    std::ofstream fout(path);
    if (fout.fail() == true)
    {
        SIM_ERROR("WriteJson open {} failed", path);
        exit(1);
    }
    writer->write(value, &fout);
    fout.close();
    // cFileUtil::DeleteLock(path);
    return fout.fail() == false;
}

#define JSONUTIL_ASSERT_NULL(root, data) (root.isMember(data))

int cJsonUtil::ParseAsInt(const std::string &data_field_name,
                          const Json::Value &root)
{
    if (false == JSONUTIL_ASSERT_NULL(root, data_field_name))
    {
        SIM_ERROR("ParseAsInt {} failed", data_field_name.c_str());
        exit(0);
    }
    return root[data_field_name].asInt();
}

std::string cJsonUtil::ParseAsString(const std::string &data_field_name,
                                     const Json::Value &root)
{
    if (false == JSONUTIL_ASSERT_NULL(root, data_field_name))
    {
        SIM_ERROR("ParseAsString {} failed", data_field_name.c_str());
        exit(0);
    }
    return root[data_field_name].asString();
}

double cJsonUtil::ParseAsDouble(const std::string &data_field_name,
                                const Json::Value &root)
{
    if (false == JSONUTIL_ASSERT_NULL(root, data_field_name))
    {
        SIM_ERROR("ParseAsDouble {} failed", data_field_name.c_str());
        exit(0);
    }
    return root[data_field_name].asDouble();
}

float cJsonUtil::ParseAsFloat(const std::string &data_field_name,
                              const Json::Value &root)
{
    if (false == JSONUTIL_ASSERT_NULL(root, data_field_name))
    {
        SIM_ERROR("ParseAsFloat {} failed", data_field_name.c_str());
        exit(0);
    }
    return root[data_field_name].asFloat();
}

bool cJsonUtil::ParseAsBool(const std::string &data_field_name,
                            const Json::Value &root)
{
    if (false == JSONUTIL_ASSERT_NULL(root, data_field_name))
    {
        SIM_ERROR("ParseAsBool {} failed", data_field_name.c_str());
        exit(0);
    }
    return root[data_field_name].asBool();
}

Json::Value cJsonUtil::ParseAsValue(const std::string &data_field_name,
                                    const Json::Value &root)
{
    if (false == JSONUTIL_ASSERT_NULL(root, data_field_name))
    {
        SIM_ERROR("ParseAsValue {} failed", data_field_name.c_str());
        exit(0);
    }
    return root[data_field_name];
}

/**
 * \brief           read matrix json
 */
bool cJsonUtil::ReadMatrixJson(const Json::Value &root, tMatrixXd &out_mat)
{
    out_mat.resize(0, 0);
    tEigenArr<tVectorXd> mat(0);
    int num_of_cols = -1;
    for (int i = 0; i < root.size(); i++)
    {
        tVectorXd res;
        bool succ = cJsonUtil::ReadVectorJson(root[i], res);
        if (succ == false)
        {
            // fail to load the row vector
            return false;
        }
        if (num_of_cols == -1)
        {
            num_of_cols = res.size();
            mat.push_back(res);
        }
        else
        {
            // the dimension doesn't meet
            if (num_of_cols != res.size())
            {
                return false;
            }
            else
            {
                mat.push_back(res);
            }
        }
    }
    out_mat.noalias() = tMatrixXd::Zero(mat.size(), num_of_cols);
    for (int i = 0; i < mat.size(); i++)
    {
        out_mat.row(i) = mat[i].transpose();
    }
    return true;
}

bool cJsonUtil::HasValue(const std::string &name, const Json::Value &root)
{
    return JSONUTIL_ASSERT_NULL(root, name);
}