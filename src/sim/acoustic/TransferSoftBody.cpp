#include "TransferSoftBody.h"
#include "AudioOutput.h"
#include "AudioWave.h"
#include "Monopole.h"
#include "geometries/Primitives.h"
#include "imgui.h"
#include "utils/FileUtil.h"
#include "utils/JsonUtil.h"
#include <iostream>
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
void cTransferSoftBody::LoadModalAnalysisResult()
{
    // SIM_ASSERT(cFileUtil::ExistsFile(mModalAnalysisResult));
    // Json::Value root;
    // SIM_ASSERT(cJsonUtil::LoadJson(mModalAnalysisResult, root));

    // Json::Value modes, coef;
    // std::string mSurfaceObjPath = cJsonUtil::ParseAsValue("surface_obj_path",
    // root); modes = cJsonUtil::ParseAsValue("modes", root); coef =
    // cJsonUtil::ParseAsValue("coef", root);

    // // 1. parse tet path, judge equality
    // {
    //     std::string tet_path = cJsonUtil::ParseAsString("tet_path",
    //     geometry); SIM_ASSERT(tet_path == mTetMeshPath);
    // }
    // int num_of_surface_v = mSurfaceVertexIdArray.size();
    // std::cout << "num_of_surface_v = " << num_of_surface_v << std::endl;
    // mModalVibrationArray.clear();
    // // 2. parse modes
    // {
    //     for (auto &info : modes)
    //     {
    //         tVector mode_vib_info = cJsonUtil::ReadVectorJson(info);

    //         auto ptr = std::make_shared<tModeVibration>(mode_vib_info);
    //         mModalVibrationArray.push_back(ptr);
    //     }
    //     std::cout << "num of valid modes = " << mModalVibrationArray.size()
    //               << std::endl;
    // }
    // // 3. parse weight
    // {
    //     mSurfaceVertexDOFCoef.resize(3 * num_of_surface_v,
    //                                  mModalVibrationArray.size());

    //     for (int i = 0; i < coef.size(); i++)
    //     {
    //         mSurfaceVertexDOFCoef.row(i) =
    //         cJsonUtil::ReadVectorJson(coef[i]);
    //     }
    // }
}
#include "ModeVibration.h"
#include "utils/JsonUtil.h"
void cTransferSoftBody::Init(const Json::Value &conf)
{
    Json::Value conf_ = conf;
    mModalAnalysisResult =
        cJsonUtil::ParseAsString("modal_analysis_result", conf_);
    {
        Json::Value modal_analysis_conf;
        SIM_ASSERT(
            cJsonUtil::LoadJson(mModalAnalysisResult, modal_analysis_conf));

        std::string mesh_path =
            cJsonUtil::ParseAsString("surface_obj_path", modal_analysis_conf);
        conf_["mesh_path"] = mesh_path;
        cKinematicBody::Init(conf_);

        // 1. get the result of modal analysis
        Json::Value modes_array_json, coef_array_json;
        modes_array_json =
            cJsonUtil::ParseAsValue("modes", modal_analysis_conf);
        coef_array_json = cJsonUtil::ParseAsValue("coef", modal_analysis_conf);
        mModalVibrationWaveArray.clear();
        // restore modal analysis result from modes and coef
        int num_of_modes = modes_array_json.size();
        for (auto &coef_json : modes_array_json)
        {
            tVector coef = cJsonUtil::ReadVectorJson(coef_json).segment(0, 4);
            mModalVibrationWaveArray.push_back(
                std::make_shared<tModeVibration>(coef));
        }

        // restore vertex coef
        int num_of_v = mVertexArray.size();
        mVertexModesCoef.resize(3 * num_of_v, num_of_modes);

        for (int i = 0; i < 3 * num_of_v; i++)
        {
            mVertexModesCoef.row(i) =
                cJsonUtil::ReadVectorJson(coef_array_json[i]);
        }
        printf("init done\n");
    }
    SetModalAnalysisSound();
    int num_of_mode = mModalVibrationWaveArray.size();

    for (int i = 0; i < num_of_mode; i++)
    {
        tVectorXd weight;
        auto poles = SolveForTargetMode(i, weight);
        pole_array_array.push_back(poles);
        pole_weight_array.push_back(weight);
    }
    SolveSoundForPoles(pole_array_array, pole_weight_array);
    return;
    int num_of_v = mVertexArray.size();
    double omega = 2000;
    // 1. preset a sound pressure on each surface
    tVectorXd vertex_an;
    PresetVertexAccel(vertex_an, omega, mVertexArray);

    // 2. set the monopole place (origin)
    tVector3d monopole_pos = tVector3d::Zero();
    double rho = 1.13;
    // 3. solve the sound pressure

    // 3.1 calculate matrix A
    mPole = std::make_shared<cMonopole>(0, omega, tVector3d::Zero(3));
    tVectorXd A = tVectorXd::Zero(num_of_v);
    for (int i = 0; i < num_of_v; i++)
    {
        auto v = mVertexArray[i];
        // only a single pole
        A(i, 0) =
            mPole->CalcdPdn(v->mPos.segment(0, 3), v->mNormal.segment(0, 3));
    }

    tVectorXd b = -rho * vertex_an;
    double C = 1.0 / A.dot(A) * A.dot(b);
    std::cout << "C = " << C << std::endl;
    std::cout << "pole pos = " << mPole->mPos.transpose() << std::endl;
    {
        for (int i = 0; i < num_of_v; i++)
        {
            auto v = mVertexArray[i];
            // 1. get its an
            double an = vertex_an[i];
            double RHS = -rho * an;
            tVector3d pos = v->mPos.segment(0, 3);
            tVector3d n = v->mNormal.segment(0, 3);

            // 2. get its sound pressure grad
            double LHS = mPole->CalcdPdn(pos, n) * C;
            printf("v %d LHS %.3f RHS %.3f\n", i, LHS, RHS);
            // std::cout << "pos = " << pos.transpose() << " n = " <<
            // n.transpose() << std::endl;
        }
    }
    // 3.2 calculate sound pressure
    mOmega = omega;
    mC = C;
    mCamPos.setOnes();
    mAnaWave = std::make_shared<tAnalyticWave>();
    ResolveSoundByCamPos();
    mCurSolvedCamPos = mCamPos;
}

void cTransferSoftBody::Update(float dt)
{
    // return;
    if ((mCurSolvedCamPos - mCamPos).norm() > 0.1)
    {
        // ResolveSoundByCamPos();
        SolveSoundForPoles(pole_array_array, pole_weight_array);
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

/*
p(x, t) = C * p(x) * cos (w (t - r / c))
*/
void cTransferSoftBody::ResolveSoundByCamPos()
{
    // 1. calculate analytic wave
    mAnaWave->Clear();

    // double strength = C * p(x)
    double strength = mC * mPole->EvaluatePressure(mCamPos);

    // omega =
    double f_hz = mOmega / (2 * M_PI);
    printf("strength %.1f frequench %.1f hz\n", strength, f_hz);
    std::cout << "cam pos = " << mCamPos.transpose() << std::endl;
    mAnaWave->AddWave(strength, f_hz);
    // 2. convert to discrete wave
    double duration = 2.0;
    double dt = 1.0 / (44000);
    auto discrete_wave = DiscretizeWave(mAnaWave, duration, dt);
    // 3. set the audio manager
    auto output = cAudioOutput::getInstance();
    output->SetWave(discrete_wave);
}

void cTransferSoftBody::SetModalAnalysisSound()
{
    int tar_vid = 0;
    double sec = 1.0;
    tVectorXf data = tVectorXf::Zero(0);
    double freq = 48000;
    tVectorXd cur_coef = mVertexModesCoef.row(tar_vid);
    std::cout << "combination coef = " << cur_coef.transpose() << std::endl;
    for (int i = 0; i < mModalVibrationWaveArray.size(); i++)
    {
        auto mode = mModalVibrationWaveArray[i];

        double coef = cur_coef[i];

        printf("mode %d coef %.3f info = ", i, coef);
        std::cout << mode->GetCoefVec().transpose() << std::endl;
        if (data.size() == 0)
            data = coef * mode->GetWave(sec, freq).cast<float>();
        else
            data += coef * mode->GetWave(sec, freq).cast<float>();
    }
    // scale data
    data = 1.0 / data.cwiseAbs().maxCoeff() * data;
    auto wave = std::make_shared<tDiscretedWave>(1.0 / freq);
    wave->SetData(data);
    auto output = cAudioOutput::getInstance();
    output->SetWave(wave);
}
tVector RandomPointInAABB(const tVector &aabb_min, const tVector &aabb_max)
{
    tVector pos = tVector::Zero();
    for (int i = 0; i < 3; i++)
    {
        pos[i] = cMathUtil::RandDouble(aabb_min[i], aabb_max[i]);
    }
    return pos;
}
std::vector<cMonopolePtr>
cTransferSoftBody::SolveForTargetMode(int tar_mode, tVectorXd &weight)
{
    std::vector<cMonopolePtr> pole_array = {};
    // 1. for this mode, get the vertex accel on each vertex
    int num_of_v = GetNumOfVertices();
    tVectorXd an = tVectorXd::Zero(num_of_v);
    tVectorXd vertex_mode_coef = mVertexModesCoef.col(tar_mode);
    double air_density = 1.113;
    double c = mModalVibrationWaveArray[tar_mode]->mCoef;
    double wd = mModalVibrationWaveArray[tar_mode]->mWd;
    for (int i = 0; i < num_of_v; i++)
    {
        // 1. get normal
        tVector3d n = mVertexArray[i]->mNormal.segment(0, 3);
        double part1 = vertex_mode_coef.segment(3 * i, 3).dot(n);

        an[i] = std::fabs(part1 * c * wd);
        // printf("v %d an = %.3f\n", i, an[i]);
    }

    // 2. calculate the normal gradient of monopole on each point

    // 2.1 confirm monopole pos
    tVector aabb_min, aabb_max;
    CalcAABB(aabb_min, aabb_max);
    int num_of_pole = 1;
    // init pole pos
    pole_array.clear();
    for (int i = 0; i < num_of_pole; i++)
    {
        tVector pole_pos;
        if (i == 0)
        {
            pole_pos = (aabb_max + aabb_min) / 2;
        }
        else
        {
            pole_pos = RandomPointInAABB(aabb_min, aabb_max);
        }
        auto pole = std::make_shared<cMonopole>(0, wd, pole_pos.segment(0, 3));
        pole_array.push_back(pole);
    }

    tMatrixXd A = tMatrixXd::Zero(num_of_v, num_of_pole);

    for (int i = 0; i < num_of_v; i++)
    {
        // 1. get monopole normal gradient
        tVector3d n = mVertexArray[i]->mNormal.segment(0, 3);
        tVector3d pos = mVertexArray[i]->mPos.segment(0, 3);
        for (int j = 0; j < num_of_pole; j++)
        {
            A(i, j) = pole_array[j]->CalcdPdn(pos, n);
        }
    }

    // 3. combine the linear system and solve
    tVectorXd b = -air_density * an;
    tVectorXd C = (A.transpose() * A).inverse() * A.transpose() * b;
    double residual = (A * C - b).norm();
    // std::cout << "C = " << C.transpose() << " residual = " << residual
    //           << std::endl;
    // 4. check the BCs diff
    // tVectorXd pred_vec = A * C;
    // for (int i = 0; i < num_of_v; i++)
    // {
    //     printf("v %d LHS %.3f RHS %.3f\n", i, pred_vec[i], b[i]);
    // }
    // exit(1);
    weight = C;
    return pole_array;
}

void cTransferSoftBody::SolveSoundForPoles(
    const std::vector<std::vector<cMonopolePtr>> &pole_array,
    const tEigenArr<tVectorXd> &pole_weight_array)
{
    int num_of_mode = pole_array.size();
    double duration = 1.0;
    double freq = 48000;
    tVectorXd data = tVectorXd::Zero(int(duration * freq));
    for (int i = 0; i < num_of_mode; i++)
    {
        auto pole_per_mode_array = pole_array[i];
        double total_pressure = 0;
        for (auto pole : pole_per_mode_array)
        {
            double pressure = pole->EvaluatePressure(this->mCamPos);
            total_pressure += pressure;
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