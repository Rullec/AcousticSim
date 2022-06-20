#include "sim/acoustic/TransferSoftBody.h"
#include "geometries/Primitives.h"
#include "sim/acoustic/AudioOutput.h"
#include "sim/acoustic/AudioWave.h"
#include "sim/acoustic/ModeVibration.h"
#include "sim/acoustic/Monopole.h"
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
static int gGDIters, gGDPrintOuts;
static double gGDStepsize;
void cTransferSoftBody::Init(const Json::Value &conf)
{
    cSoftBodyImplicit::Init(conf);
    mModalAnlysisResultPath =
        cJsonUtil::ParseAsString("modal_analysis_result", conf);
    mEnableDumpedModalSolution = false;
    mNumOfMonopolesPerFreq = cJsonUtil::ParseAsInt("num_of_poles", conf);
    gGDIters = cJsonUtil::ParseAsInt("gd_iters", conf);
    gGDStepsize = cJsonUtil::ParseAsDouble("gd_stepsize", conf);
    gGDPrintOuts = cJsonUtil::ParseAsInt("gd_printgap", conf);
    CalcAABB(this->mAABBMin, this->mAABBMax);
    mFixPolePosition = cJsonUtil::ParseAsBool("fix_pole_position", conf);

    InitSurfaceNormal();

    LoadModalAnalysisResult();

    int num_of_modes = mSurfaceVertexDOFCoef.cols();
    InitPole();

    for (int i = 0; i < num_of_modes; i++)
    {
        // printf("------mode %d ------\n", i);
        if (mFixPolePosition == true)
        {
            SolveSoundFixPole(i);
        }
        else
        {
            SIM_ERROR("for moving pole, hasn't been impled (Gauss-Newton)");
        }
        // SolveMonopole(i, sound_pressure);
    }
    SetModalVibrationSound();
}

void cTransferSoftBody::SolveSoundFixPole(int mode_idx)
{
    // 1. get the coefficient matrix C
    /*
        C * A = b
        where c_{ij} = gradient of p_j on surface vertex i
        b_i = - rho * ai.dot(ni)
    */
    // PrintPoleInfo(mode_idx);
    int num_of_surface_v = this->mSurfaceVertexIdArray.size();
    int num_of_poles = this->mNumOfMonopolesPerFreq;
    tVectorXd surface_an = CalcSurfaceVertexAccelNormal(mode_idx);
    tVectorXd b = -gAirDensity * surface_an;
    tMatrixXd C = tMatrixXd::Zero(num_of_surface_v, num_of_poles);
    for (int vid = 0; vid < num_of_surface_v; vid++)
    {
        int vid_global = mSurfaceVertexIdArray[vid];
        auto v = mVertexArray[vid_global];
        tVector3d v_pos = v->mPos.segment(0, 3);
        tVector3d v_normal = v->mNormal.segment(0, 3);

        // 1. calculate cij
        for (int pid = 0; pid < num_of_poles; pid++)
        {
            auto cur_pole = mPolesArray[mode_idx].mPoles[pid];
            double cij = cur_pole->CalcCoef(v_pos, v_normal);
            C(vid, pid) = cij;
        }
    }
    // std::cout << "C = \n" << C << std::endl;
    // std::cout << "b = " << b.transpose() << std::endl;
    /*
        C_{nm} * A_{m} = b_{n}
        C^T * C A = CT * b
        A = (CTC).inv() * CT * b
    */
    tVectorXd pole_weight = (C.transpose() * C).inverse() * C.transpose() * b;

    // std::cout << "mode " << mode_idx
    //           << "pole weight = " << pole_weight.transpose() << std::endl;
    for (int i = 0; i < num_of_poles; i++)
    {
        mPolesArray[mode_idx].mPoles[i]->mA = pole_weight[i];
    }

    // PrintPoleInfo(mode_idx);
    // tVectorXd dPdn_pred = tVectorXd::Zero(num_of_surface_v);
    // for (int vid = 0; vid < num_of_surface_v; vid++)
    // {
    //     int vid_global = mSurfaceVertexIdArray[vid];
    //     auto v = mVertexArray[vid_global];
    //     tVector3d v_pos = v->mPos.segment(0, 3);
    //     tVector3d v_normal = v->mNormal.segment(0, 3);

    //     // 1. calculate cij
    //     for (int pid = 0; pid < num_of_poles; pid++)
    //     {
    //         auto cur_pole = mPolesArray[mode_idx].mPoles[pid];
    //         dPdn_pred[vid] += cur_pole->CalcdPdn(v_pos, v_normal);
    //     }
    // }
    // tVectorXd dpdn_diff = dPdn_pred - b;
    // std::cout << "dpdn pred = " << dPdn_pred.transpose() << std::endl;
    // std::cout << "dpdn real = " << b.transpose() << std::endl;
    // std::cout << "dpdn diff = " << dpdn_diff.transpose() << std::endl;
    // std::cout << "dpdn diff norm = " << dpdn_diff.norm() << std::endl;
    // exit(1);
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
void cTransferSoftBody::SolveMonopole(int mode_idx,
                                      const tVectorXd &sound_pressure)
{
    // int max_iters = gGDIters;
    // double stepsize = gGDStepsize;
    // tVectorXd cur_x = GetX(mode_idx);
    // for (int iters = 0; iters < max_iters; iters++)
    // {
    //     // this->GetX(mode_idx)

    //     double e = GetEnergy(mode_idx, sound_pressure);
    //     if (iters % gGDPrintOuts == 0)
    //     {
    //         // std::cout << "iter " << iters << " x = " << cur_x.transpose()
    //         //           << std::endl;

    //         printf("iter %d ener %.3f\n", iters, e);
    //         std::cout << "AABB info: " << mAABBMin.transpose()
    //                   << mAABBMax.transpose() << std::endl;
    //         PrintPoleInfo(mode_idx);
    //     }
    //     auto sound_diff = GetSoundPressureDiff(mode_idx, sound_pressure);
    //     tVectorXd grad = GetGrad(mode_idx, sound_diff);
    //     cur_x += -stepsize * grad;
    //     ConstrainX(cur_x);
    //     SetX(mode_idx, cur_x);
    // }
    // PrintSoundPressureDiffInfo(mode_idx, sound_pressure);
    // auto diff = GetSoundPressureDiff(mode_idx, sound_pressure);
    // for (int i = 0; i < this->mSurfaceVertexIdArray.size(); i++)
    // {
    //     std::cout << "vertex " << i << " diff = " << diff[i].transpose()
    //               << std::endl;
    // }
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
 */
tVectorXd cTransferSoftBody::CalcSurfaceVertexAccelNormal(int mode_idx)
{
    auto mode_vib = mModalVibrationArray[mode_idx];
    int num_of_surface_v = this->mSurfaceVertexIdArray.size();
    SIM_ASSERT(mSurfaceVertexDOFCoef.rows() == 3 * num_of_surface_v);
    int num_of_modes = mModalVibrationArray.size();
    // printf("w = %.3f\n", mode_vib->mW);
    tVectorXd vertex_acc_normal = tVectorXd::Zero(num_of_surface_v);
    for (int i = 0; i < num_of_surface_v; i++)
    {
        tVector3d normal = mSurfaceVertexNormalArray[i].segment(0, 3);
        // 1. get the value of this sound pressure
        tVector3d vertex_modes_coef =
            mSurfaceVertexDOFCoef.block(3 * i, mode_idx, 3, 1);

        // 2. calculate pressure (directional)
        tVector3d pressure =
            gAirDensity * mode_vib->mWd * mode_vib->mCoef * vertex_modes_coef;

        // 3. normalize to vertex normal direction

        double pressure_amp = std::fabs(pressure.dot(normal));
        // std::cout << "v " << i << " pressure = " << pressure.transpose()
        //           << " normal = " << normal.transpose()
        //           << " amp = " << pressure_amp << std::endl;
        vertex_acc_normal[i] = pressure_amp;
    }
    return vertex_acc_normal;
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
            pole->mPos =
                tVector3d(cMathUtil::RandDouble(mAABBMin[0], mAABBMax[0]),
                          cMathUtil::RandDouble(mAABBMin[1], mAABBMax[1]),
                          cMathUtil::RandDouble(mAABBMin[2], mAABBMax[2]));
            pole->mA = cMathUtil::RandDouble(0, 1);
            data.mPoles.push_back(pole);
        }
        data.mOmega = mModalVibrationArray[i]->mW;

        mPolesArray.push_back(data);
    }
}

// /**
//  * \brief           get energy
//  */
// double cTransferSoftBody::GetEnergy(int mode_idx,
//                                     const tVectorXd &sound_pressure)
// {
//     // calcualte residual
//     auto diff = this->GetSoundPressureDiff(mode_idx, sound_pressure);
//     double e = 0;
//     for (auto &x : diff)
//     {
//         e += x.squaredNorm();
//     }
//     return e;
// }

// /**
//  * \brief           get grad
//  */
// tVectorXd
// cTransferSoftBody::GetGrad(int mode_idx,
//                            const tEigenArr<tVector3d> &sound_pressure_diff)
// {
//     tVectorXd grad = tVectorXd::Zero(4 * this->mNumOfMonopolesPerFreq);

//     for (int pole_id = 0; pole_id < this->mNumOfMonopolesPerFreq; pole_id++)
//     {
//         auto pole = mPolesArray[mode_idx].mPoles[pole_id];

//         // for each boundary very
//         for (int v_id = 0; v_id < this->mSurfaceVertexIdArray.size(); v_id++)
//         {
//             int v_gid = mSurfaceVertexIdArray[v_id];
//             tVector3d v_global_pos = mVertexArray[v_gid]->mPos.segment(0, 3);
//             // 1. calcualte d
//             tVector3d d = sound_pressure_diff[v_id];
//             // 2. calculate dPdocef and dPdcom
//             tVector3d dPdcoef = pole->EvaluatedPdcoef(v_global_pos);
//             tMatrix3d dPdcenter = pole->EvaulatedPdcenter(v_global_pos);
//             // 3. multiply
//             grad[4 * pole_id] += -2 * d.dot(dPdcoef);
//             grad.segment(4 * pole_id + 1, 3) += -2 * d.transpose() *
//             dPdcenter; if (grad.hasNaN())
//             {
//                 std::cout << "grad hasNan, dpdcoef = " << dPdcoef.transpose()
//                           << " dpdc = " << dPdcenter << std::endl;
//                 std::cout << "d = " << d.transpose() << std::endl;
//                 exit(1);
//             }
//             // 4. write down
//         }

//         // gradient to coef

//         // grad to center
//     }

//     return grad;
// }

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
        x[4 * i] = cur_pole->mA;
        x.segment(4 * i + 1, 3) = cur_pole->mPos;
    }
    return x;
}

// tEigenArr<tVector3d>
// cTransferSoftBody::GetSoundPressureDiff(int mode_idx,
//                                         const tVectorXd &sound_pressure)
// {
//     tEigenArr<tVector3d> pressure_diff = {};
//     auto pred_pressure = GetPredSoundPressure(mode_idx, sound_pressure);
//     for (int v_id = 0; v_id < this->mSurfaceVertexIdArray.size(); v_id++)
//     {
//         pressure_diff.push_back(mSurfaceVertexNormalArray[v_id].segment(0, 3)
//         *
//                                     sound_pressure[v_id] -
//                                 pred_pressure[v_id]);
//     }
//     return pressure_diff;
// }

// tEigenArr<tVector3d>
// cTransferSoftBody::GetPredSoundPressure(int mode_idx,
//                                         const tVectorXd &sound_pressure)
// {
//     tEigenArr<tVector3d> pressure = {};
//     for (int v_id = 0; v_id < this->mSurfaceVertexIdArray.size(); v_id++)
//     {
//         tVector3d v_pos =
//             mVertexArray[mSurfaceVertexIdArray[v_id]]->mPos.segment(0, 3);
//         // prediction
//         tVector3d pred = tVector3d::Zero();
//         for (int m_id = 0; m_id < mNumOfMonopolesPerFreq; m_id++)
//         {
//             pred +=
//             mPolesArray[mode_idx].mPoles[m_id]->EvaluatePressure(v_pos);
//         }
//         pressure.push_back(pred);
//     }
//     return pressure;
// }

// void cTransferSoftBody::CheckGradient()
// {
//     // 1. get old energy
//     int mode_idx = 0;
//     tVectorXd sound_pressure_gt = CalcSoundPressure(mode_idx);
//     double old_e = GetEnergy(mode_idx, sound_pressure_gt);

//     // 2. get grad ana
//     auto sound_pressure_diff =
//         GetSoundPressureDiff(mode_idx, sound_pressure_gt);
//     tVectorXd grad_ana = GetGrad(mode_idx, sound_pressure_diff);
//     tVectorXd grad_num = tVectorXd::Zero(grad_ana.size());
//     tVectorXd curX = GetX(mode_idx);
//     // 3. for each dof, add and set; calculate grad num
//     double eps = 1e-6;
//     for (int i = 0; i < grad_num.size(); i++)
//     {
//         curX[i] += eps;
//         SetX(mode_idx, curX);
//         double new_e = this->GetEnergy(mode_idx, sound_pressure_gt);
//         grad_num[i] = (new_e - old_e) / eps;
//         curX[i] -= eps;
//         SetX(mode_idx, curX);
//     }
//     std::cout << "ana = " << grad_ana.transpose() << std::endl;
//     std::cout << "num = " << grad_num.transpose() << std::endl;
//     auto diff = grad_ana - grad_num;
//     std::cout << "diff = " << diff.transpose() << std::endl;
//     // exit(1);
// }

void cTransferSoftBody::SetX(int mode_idx, const tVectorXd &sol)
{
    //  tVectorXd x = tVectorXd::Zero(4 * mNumOfMonopolesPerFreq);

    for (int i = 0; i < mNumOfMonopolesPerFreq; i++)
    {
        auto cur_pole = mPolesArray[mode_idx].mPoles[i];
        // coef
        cur_pole->mA = sol[4 * i];
        cur_pole->mPos = sol.segment(4 * i + 1, 3);
    }
}

void cTransferSoftBody::ConstrainX(tVectorXd &sol)
{
    for (int i = 0; i < this->mNumOfMonopolesPerFreq; i++)
    {
        double &strength = sol[4 * i];
        tVector3d pos = sol.segment(4 * i + 1, 3);
        // 1. constrain strength
        strength = SIM_MAX(0.0, strength);
        for (int j = 0; j < 3; j++)
        {
            pos[j] = cMathUtil::Clamp(pos[j], mAABBMin[j] + 1e-3,
                                      mAABBMax[j] - 1e-3);
        }
        sol.segment(4 * i + 1, 3) = pos;
        sol[4 * i] = strength;
    }
}

void cTransferSoftBody::PrintPoleInfo(int mode_idx)
{
    for (int i = 0; i < this->mNumOfMonopolesPerFreq; i++)
    {
        auto pole = mPolesArray[mode_idx].mPoles[i];
        printf("pole %d pos %.3f %.3f %.3f, strength %.3f\n", i, pole->mPos[0],
               pole->mPos[1], pole->mPos[2], pole->mA);
    }
}

// void cTransferSoftBody::PrintSoundPressureDiffInfo(
//     int mode_idx, const tVectorXd &sound_pressure)
// {
//     tEigenArr<tVector3d> pressure_diff = {};
//     auto pred_pressure = GetPredSoundPressure(mode_idx, sound_pressure);
//     for (int v_id = 0; v_id < this->mSurfaceVertexIdArray.size(); v_id++)
//     {
//         tVector3d real = mSurfaceVertexNormalArray[v_id].segment(0, 3) *
//                          sound_pressure[v_id];
//         tVector3d pred = pred_pressure[v_id];
//         tVector3d diff = real - pred;
//         std::cout << "v " << v_id << " real = " << real.transpose()
//                   << " pred = " << pred.transpose()
//                   << " diff = " << diff.transpose() << std::endl;
//     }
// }