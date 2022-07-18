#include "TransferSoftBody.h"
#include "AudioOutput.h"
#include "AudioWave.h"
#include "ModeVibration.h"
#include "Monopole.h"
#include "geometries/Primitives.h"
#include "imgui.h"
#include "utils/FileUtil.h"
#include "utils/JsonUtil.h"
#include <iostream>

tVector RandomPointInAABB(const tVector &aabb_min, const tVector &aabb_max)
{
    tVector pos = tVector::Zero();
    for (int i = 0; i < 3; i++)
    {
        pos[i] = cMathUtil::RandDouble(aabb_min[i], aabb_max[i]);
    }
    return pos;
}

cTransferSoftBody::cTransferSoftBody(int id_) : cKinematicBody(id_)
{
    mType = eObjectType::ACOUSTIC_TRANSFER_TYPE;
}

/*
a_n(t) = - A * w * w
*/
void PresetVertexAccel(tVectorXd &vertex_accel, double omega,
                       const std::vector<tVertexPtr> &v_array)
{
    double A = 1e-6;
    vertex_accel.resize(v_array.size());
    for (int i = 0; i < vertex_accel.size(); i++)
    {
        vertex_accel[i] = -A * omega * omega;
    }
}
void cTransferSoftBody::Init(const Json::Value &conf)
{
    // load modal analysis and init kinematic body
    mNumOfPolesPerMode = cJsonUtil::ParseAsInt("num_of_poles_per_mode", conf);
    LoadModalAnalysisResult(conf);

    // calc AABB for pole result
    CalcAABB(mAABBMin, mAABBMax);

    int num_of_mode = mModalVibrationWaveArray.size();
    // begin to solve acoustic transfer
    mPolesArrayArray.resize(num_of_mode, {});
    mPolesWeightArray.clear();
    for (int i = 0; i < num_of_mode; i++)
    {
        // init pole pos
        InitPole(mPolesArrayArray[i], i);
        // calculate pole amplitue
        tVectorXd weight = SolveForTargetMode(i);
        mPolesWeightArray.push_back(weight);
    }

    // generate sound from sources
    GenerateSound(mPolesArrayArray, mPolesWeightArray);
}

void cTransferSoftBody::InitPole(cPoleArray &pole_array, int mode_idx)
{
    double wd = mModalVibrationWaveArray[mode_idx]->mWd;
    pole_array.clear();
    for (int i = 0; i < mNumOfPolesPerMode; i++)
    {
        tVector pole_pos;
        if (i == 0)
        {
            pole_pos = (mAABBMin + mAABBMax) / 2;
        }
        else
        {
            pole_pos = RandomPointInAABB(mAABBMin, mAABBMax);
        }
        auto pole = std::make_shared<cMonopole>(0, wd, pole_pos.segment(0, 3));
        // std::cout << "pole pos = " << pole_pos.transpose() << std::endl;
        pole_array.push_back(pole);
    }
}

void cTransferSoftBody::LoadModalAnalysisResult(const Json::Value &root)
{
    Json::Value conf_ = root;
    mModalAnalysisPath =
        cJsonUtil::ParseAsString("modal_analysis_result", conf_);

    Json::Value modal_analysis_conf;
    SIM_ASSERT(cJsonUtil::LoadJson(mModalAnalysisPath, modal_analysis_conf));

    // set kinematic mesh path by modal anlaysis result
    std::string mesh_path =
        cJsonUtil::ParseAsString("surface_obj_path", modal_analysis_conf);
    conf_["mesh_path"] = mesh_path;
    cKinematicBody::Init(conf_);
    // get the result of modal analysis
    Json::Value modes_array_json, coef_array_json;
    modes_array_json = cJsonUtil::ParseAsValue("modes", modal_analysis_conf);
    coef_array_json = cJsonUtil::ParseAsValue("coef", modal_analysis_conf);
    // modal vibration wave
    mModalVibrationWaveArray.clear();
    // restore modal analysis result from modes and coef
    int num_of_modes = modes_array_json.size();
    for (auto &coef_json : modes_array_json)
    {
        tVector coef = cJsonUtil::ReadVectorJson(coef_json).segment(0, 4);
        mModalVibrationWaveArray.push_back(
            std::make_shared<tModeVibration>(coef));
    }

    // vertex combination coef for each mode
    int num_of_v = mVertexArray.size();
    mVertexModesCoef.resize(3 * num_of_v, num_of_modes);

    for (int i = 0; i < 3 * num_of_v; i++)
    {
        mVertexModesCoef.row(i) = cJsonUtil::ReadVectorJson(coef_array_json[i]);
    }
}

void cTransferSoftBody::Update(float dt)
{
    // return;
    if ((mCurSolvedCamPos - mCamPos).norm() > 0.1)
    {
        GenerateSound(mPolesArrayArray, mPolesWeightArray);
        mCurSolvedCamPos = mCamPos;
    }
}
void cTransferSoftBody::UpdateImGui()
{
    ImGui::Text("cam pos %.2f %.2f %.2f\n", mCamPos[0], mCamPos[1], mCamPos[2]);
}

void cTransferSoftBody::SetCameraPos(const tVector3d &pos)
{
    this->mCamPos = pos;
}

tVectorXd cTransferSoftBody::CalcVertexAccel(int mode_idx)
{
    int num_of_v = GetNumOfVertices();
    tVectorXd an = tVectorXd::Zero(num_of_v);
    tVectorXd vertex_mode_coef = mVertexModesCoef.col(mode_idx);
    double c = mModalVibrationWaveArray[mode_idx]->mCoef;
    double wd = mModalVibrationWaveArray[mode_idx]->mWd;
    for (int i = 0; i < num_of_v; i++)
    {
        // 1. get normal
        tVector3d n = mVertexArray[i]->mNormal.segment(0, 3);
        double part1 = vertex_mode_coef.segment(3 * i, 3).dot(n);

        // 2. an = a.dot(n)
        // a = coef * c * wd * e^{iwt}
        // -> an = [coef].dot(n) * c * d
        an[i] = std::fabs(part1 * c * wd);
    }
    return an;
}

tMatrixXd
cTransferSoftBody::CalcSoundPressureGradient(const cPoleArray &pole_array,
                                             int mode_idx)
{
    int num_of_pole = pole_array.size();
    int num_of_v = mVertexArray.size();

    tMatrixXd A = tMatrixXd::Zero(num_of_v, num_of_pole);

    for (int i = 0; i < num_of_v; i++)
    {
        // 1. get monopole normal gradient
        tVector3d n = mVertexArray[i]->mNormal.segment(0, 3);
        tVector3d pos = mVertexArray[i]->mPos.segment(0, 3);
        for (int j = 0; j < num_of_pole; j++)
        {
            A(i, j) = pole_array[j]->CalcdPdn_Re(pos, n);
        }
    }
    return A;
}

tVectorXd cTransferSoftBody::SolveForTargetMode(int mode_idx)
{

    // 2. calculate the normal gradient of monopole on each point
    tVectorXd an = CalcVertexAccel(mode_idx);
    tMatrixXd A =
        CalcSoundPressureGradient(mPolesArrayArray[mode_idx], mode_idx);

    int num_of_pole = 1;
    // init pole pos

    // 3. combine the linear system and solve
    double air_density = 1.13;
    tVectorXd b = -air_density * an;
    tVectorXd C = (A.transpose() * A).inverse() * A.transpose() * b;
    double residual = (A * C - b).norm();
    std::cout << "mode " << mode_idx << " residual = " << residual << std::endl;
    return C;
}

void cTransferSoftBody::GenerateSound(
    const std::vector<std::vector<cMonopolePtr>> &pole_array,
    const tEigenArr<tVectorXd> &mPolesWeightArray)
{
    int num_of_mode = pole_array.size();
    double duration = 1.0;
    double freq = 48000;
    tVectorXd data = tVectorXd::Zero(int(duration * freq));
    for (int i = 0; i < num_of_mode; i++)
    {
        tVectorXd cur_mode_weight = mPolesWeightArray[i];
        auto pole_per_mode_array = pole_array[i];
        double total_pressure = 0;
        int offset = 0;
        int num_of_pole = pole_per_mode_array.size();
        for (int j = 0; j < num_of_pole; j++)
        {
            auto pole = pole_per_mode_array[j];
            int cur_dof = pole->GetNumOfDof();
            tVectorXd weight = cur_mode_weight.segment(offset, cur_dof);

            double pressure =
                pole->CalcPressureForSoundSynthesis(this->mCamPos, weight);
            total_pressure += pressure;
            offset += cur_dof;
        }
        // 1.
        data += total_pressure *
                mModalVibrationWaveArray[i]->GetWaveExpAndSin(duration, freq);
    }
    auto wave = std::make_shared<tDiscretedWave>(1.0 / freq);
    wave->SetData(data.cast<float>());
    auto output = cAudioOutput::getInstance();
    output->SetWave(wave);
}