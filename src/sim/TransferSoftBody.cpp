#include "sim/TransferSoftBody.h"
#include "geometries/Primitives.h"
#include "sim/AudioOutput.h"
#include "sim/AudioWave.h"
#include "sim/ModeVibration.h"
#include "sim/Monopole.h"
#include "utils/FileUtil.h"
#include "utils/JsonUtil.h"

#include <iostream>
const double gAirDensity = 1.225; // kg.m-3
extern cAudioOutputPtr gAudioOutput;
extern tDiscretedWavePtr
CalculateVertexVibration(int v_id, int num_of_selected_modes,
                         const std::vector<tDiscretedWavePtr> &mModalWaves,
                         const tMatrixXd &mVertexModesCoef, double sampling_hz);
cTransferSoftBody::cTransferSoftBody(int id_) : cSoftBodyImplicit(id_) {}

void cTransferSoftBody::Init(const Json::Value &conf)
{
    cSoftBodyImplicit::Init(conf);
    mModalAnlysisResultPath =
        cJsonUtil::ParseAsString("modal_analysis_result", conf);
    mEnableDumpedModalSolution = false;
    mNumOfMonopolesPerFreq = 1;
    InitSurfaceNormal();

    LoadModalAnalysisResult();

    int num_of_modes = mSurfaceVertexDOFCoef.cols();
    InitPole();

    for (int i = 0; i < num_of_modes; i++)
    {
        printf("------mode %d ------\n", i);
        CalcSoundPressure(i);
        SolveMonopole(i);
    }
    SetModalVibrationSound();
}
void cTransferSoftBody::InitSurfaceNormal()
{
    mSurfaceVertexNormalArray.clear();
    for (auto &vid : this->mSurfaceVertexIdArray)
    {
        mSurfaceVertexNormalArray.push_back(mVertexArray[vid]->mNormal);
    }
}

void cTransferSoftBody::Update(float dt) { cSoftBodyImplicit::Update(dt); }

void cTransferSoftBody::UpdateImGui() { cSoftBodyImplicit::UpdateImGui(); }

int cTransferSoftBody::GetNumOfDrawTriangles() const
{
    return cSoftBodyImplicit::GetNumOfDrawTriangles();
}

int cTransferSoftBody::GetNumOfDrawEdges() const
{
    return cSoftBodyImplicit::GetNumOfDrawEdges();
}

int cTransferSoftBody::GetNumOfDrawVertices() const
{
    return cSoftBodyImplicit::GetNumOfDrawVertices();
}

void cTransferSoftBody::CalcTriangleDrawBuffer(Eigen::Map<tVectorXf> &res,
                                               int &st) const
{
    cSoftBodyImplicit::CalcTriangleDrawBuffer(res, st);
}

void cTransferSoftBody::CalcEdgeDrawBuffer(Eigen::Map<tVectorXf> &res,
                                           int &st) const
{
    cSoftBodyImplicit::CalcEdgeDrawBuffer(res, st);
}

/**
 * \brief           solve monopole positions and strengths for each mode
 */
void cTransferSoftBody::SolveMonopole(int mode_idx)
{
    // int num_of_modes = mModalVibrationsInfo.size();
    // int num_of_surface_id = mSurfaceVertexDOFCoef.size() / 3;
    // tVector mode_info = mModalVibrationsInfo[mode_idx];
    // tVectorXd dof_coef_for_modes = mSurfaceVertexDOFCoef.col(mode_idx);
}

void cTransferSoftBody::LoadModalAnalysisResult()
{
    SIM_ASSERT(cFileUtil::ExistsFile(mModalAnlysisResultPath));
    Json::Value root;
    SIM_ASSERT(cJsonUtil::LoadJson(mModalAnlysisResultPath, root));

    Json::Value geometry, modes, coef;
    geometry = cJsonUtil::ParseAsValue("geometry", root);
    modes = cJsonUtil::ParseAsValue("modes", root);
    coef = cJsonUtil::ParseAsValue("coef", root);

    // 1. parse tet path, judge equality
    {
        std::string tet_path = cJsonUtil::ParseAsString("tet_path", geometry);
        SIM_ASSERT(tet_path == mTetMeshPath);
    }
    int num_of_surface_v = mSurfaceVertexIdArray.size();
    std::cout << "num_of_surface_v = " << num_of_surface_v << std::endl;
    mModalVibrationArray.clear();
    // 2. parse modes
    {
        for (auto &info : modes)
        {
            tVector mode_vib_info = cJsonUtil::ReadVectorJson(info);

            auto ptr = std::make_shared<tModeVibration>(mode_vib_info);
            mModalVibrationArray.push_back(ptr);
        }
        std::cout << "num of valid modes = " << mModalVibrationArray.size()
                  << std::endl;
    }
    // 3. parse weight
    {
        mSurfaceVertexDOFCoef.resize(3 * num_of_surface_v,
                                     mModalVibrationArray.size());

        for (int i = 0; i < coef.size(); i++)
        {
            mSurfaceVertexDOFCoef.row(i) = cJsonUtil::ReadVectorJson(coef[i]);
        }
    }
}

/**
 * \brief       calculate sound pressure (vector) for each vertex
 *
 *      In mode j, the displacement wave of vertex i is:
 *      a_i^j(t) = - dt * U_{ij} * f_j^U  * w_j * sin(wj * t)
 *               = - U_{ij} * coef_j * wj * sin(wj * t)
 *      coef_j = dt * f_j^U
 *
 *      p_i^j = rho_air * wj * coef_j * U_{ij}
 */
tVectorXd cTransferSoftBody::CalcSoundPressure(int mode_idx)
{
    auto mode_vib = mModalVibrationArray[mode_idx];
    int num_of_surface_v = this->mSurfaceVertexIdArray.size();
    SIM_ASSERT(mSurfaceVertexDOFCoef.rows() == 3 * num_of_surface_v);
    int num_of_modes = mModalVibrationArray.size();
    printf("w = %.3f\n", mode_vib->mW);
    tVectorXd pressure_amp_vec = tVectorXd::Zero(num_of_surface_v);
    for (int i = 0; i < num_of_surface_v; i++)
    {
        tVector3d normal = mSurfaceVertexNormalArray[i].segment(0, 3);
        // 1. get the value of this sound pressure
        tVector3d vertex_modes_coef =
            mSurfaceVertexDOFCoef.block(3 * i, mode_idx, 3, 1);

        // 2. calculate pressure (directional)
        tVector3d pressure =
            gAirDensity * mode_vib->mW * mode_vib->mCoef * vertex_modes_coef;

        // 3. normalize to vertex normal direction

        double pressure_amp = std::fabs(pressure.dot(normal));
        // std::cout << "v " << i << " pressure = " << pressure.transpose()
        //           << " normal = " << normal.transpose()
        //           << " amp = " << pressure_amp << std::endl;
        pressure_amp_vec[i] = pressure_amp;
    }
    return pressure_amp_vec;
}

void cTransferSoftBody::SetModalVibrationSound()
{
    // for surface vertex 0, combine the modal vibration, then set the audio
    double sampling_hz = 48000;
    int num_of_modes = mModalVibrationArray.size();
    double dt = 1.0 / sampling_hz;
    double duration = 1.0;

    std::vector<tDiscretedWavePtr> modal_waves = {};
    for (auto &x : this->mModalVibrationArray)
    {
        auto mode_wave = std::make_shared<tDiscretedWave>(dt);
        mode_wave->SetData(x->GetWave(duration, sampling_hz).cast<float>());
        modal_waves.push_back(mode_wave);
    }

    auto play_wave = CalculateVertexVibration(
        0, num_of_modes, modal_waves, mSurfaceVertexDOFCoef, sampling_hz);

    gAudioOutput->SetWave(play_wave);
}

#include "sim/Monopole.h"

// init create pole
void cTransferSoftBody::InitPole()
{
    mPolesArray.clear();
    int num_of_modes = mModalVibrationArray.size();
    int cur_id = 0;
    mPolesArray = {};
    for (int i = 0; i < num_of_modes; i++)
    {
        auto modal = mModalVibrationArray[i];
        tPolesPerFreq data;
        // std::vector<cMonopolePtr> modes = {};
        data.mPoles = {};
        for (int j = 0; j < mNumOfMonopolesPerFreq; j++)
        {
            auto pole = std::make_shared<cMonopole>(cur_id++, modal->mW);
            data.mPoles.push_back(pole);
        }
        data.mOmega = mModalVibrationArray[i]->mW;
        // mPolesArray.push_back(modes);
    }
}

/**
 * \brief           get energy
 */
double cTransferSoftBody::GetEnergy(int mode_idx) {
    // calcualte residual 
}

/**
 * \brief           get grad
 */
tVectorXd cTransferSoftBody::GetGrad(int mode_idx) {}

/**
 * \brief           get solution X
 * X = [c_0, o_0, c_1, o_1]
 */
tVectorXd cTransferSoftBody::GetX(int mode_idx)
{
    tVectorXd x = tVectorXd::Zero(4 * mNumOfMonopolesPerFreq);

    for (int i = 0; i < mNumOfMonopolesPerFreq; i++)
    {
        auto cur_pole = mPolesArray[mode_idx].mPoles[i];
        // coef
        x[4 * i] = cur_pole->mStrength;
        x.segment(4 * i + 1, 3) = cur_pole->mCenterPos;
    }
    return x;
}