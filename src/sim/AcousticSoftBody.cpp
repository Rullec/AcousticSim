#include "sim/AcousticSoftBody.h"
#include "geometries/Primitives.h"
#include "imgui.h"
#include "sim/AudioOutput.h"
#include "sim/AudioWave.h"
#include "sim/softbody/AcousticLinearElasticity.h"
#include "sim/softbody/BaseMaterial.h"
#include "utils/FileUtil.h"
#include "utils/JsonUtil.h"
#include <fstream>
#include <iostream>
extern cAudioOutputPtr gAudioOutput;
eMaterialType GetHyperMaterialFromAcousticMaterial(std::string name);
std::string GetAcousticMaterialNameFromHyperMaterial(eMaterialType type);
int GetAcousticMaterialFromName(std::string name)
{
    if (name == "stvk")
    {
        return 0;
    }
    else if (name == "neohookean")
    {
        return 1;
    }
    else if (name == "linear_elasticity")
    {
        return 2;
    }
    else
    {
        SIM_ERROR("no material for name {}", name);
    }
}

std::string GetAcousticMaterialNameFromIdx(int idx)
{
    std::string name = "";
    switch (idx)
    {
    case 0:
        name = "stvk";
        break;
    case 1:
        name = "neohookean";
        break;
    case 2:
        name = "linear_elasticity";
        break;
    default:
        SIM_ERROR("material name from idx {}", name);
    }
    return name;
}

cAcousticSoftBody::cAcousticSoftBody(int id_) : cSoftBodyImplicit(id_) {}

void cAcousticSoftBody::Init(const Json::Value &conf)
{
    cSoftBodyImplicit::Init(conf);
    mAcousticMaterialModelIdx = GetAcousticMaterialFromName(
        GetAcousticMaterialNameFromHyperMaterial(this->mMaterial->GetType()));
    mAcousticSamplingFreqHZ =
        cJsonUtil::ParseAsInt("acoustic_sampling_freq_HZ", conf);
    bool always_resolve =
        cJsonUtil::ParseAsBool("always_resolve_modal_analysis", conf);
    mAcousticDuration = cJsonUtil::ParseAsFloat("acoustic_duration", conf);
    printf("[log] acoustic sampling freq %d HZ, duration %.1f s\n",
           mAcousticSamplingFreqHZ, mAcousticDuration);
    mHammerForce = cJsonUtil::ParseAsDouble("hammer_force", conf);
    mNumOfMaxModes = cJsonUtil::ParseAsInt("num_of_max_modes", conf);

    // load many material types ... for selection
    LoadMaterialParams();

    // init linear elasticity stiffness matrix
    mLinearElasticityStiffnessMatrix =
        cAcousticLinearElasticity::CalcGlobalStiffness(
            mVertexArray, mTetArrayShared, this->mMaterial->mYoungsModulusNew,
            this->mMaterial->mPoissonRatioNew,
            this->mMaterial->mRayleighDamplingA,
            this->mMaterial->mRayleighDamplingB);
    // std::cout << "linear K = \n"
    //           << mLinearElasticityStiffnessMatrix.toDense() << std::endl;
    ResolveVibration();
}
void cAcousticSoftBody::ResolveVibration()
{
    printf("[info] ------------resolve vibration-------------\n");
    // Init noise
    UpdateDeformationGradient();
    mGlobalStiffnessSparseMatrix = CalcGlobalStiffnessSparseMatrix();

    EigenDecompose(mEigenValues, mEigenVecs);

    CalculateModalVibration();
    auto play_wave = CalculateVertexVibration(0);
    gAudioOutput->SetWave(play_wave);

    // save
    std::string mat_name = cFileUtil::RemoveExtension(
        cFileUtil::GetFilename(this->mMaterial->mMatPath));

    play_wave->DumpToWAV(cFileUtil::ConcatFilename("log/", mat_name + ".wav"));
    std::cout << "total mass = " << CalcTotalMass() << std::endl;
}
tVectorXd CalcWave(double coef_j, double w, double xi, double wd, double dt,
                   int num_of_steps)
{
    // q(t) = coef_j / w_d e^{- \xi * w * t} sin(w_d * t)
    tVectorXd data = tVectorXd::Zero(num_of_steps);
    for (int i = 0; i < num_of_steps; i++)
    {
        double t = i * dt;
        data[i] = coef_j / wd * std::exp(-xi * w * t) * std::sin(wd * t);
    }
    return data;
}
int cAcousticSoftBody::GetNumOfModes() const
{
    if (mNumOfMaxModes <= 0)
    {
        return mEigenValues.size();
    }
    else
        return SIM_MIN(mNumOfMaxModes, mEigenValues.size());
}
/**
 * \brief   calcualte modal vibration
 */
void cAcousticSoftBody::CalculateModalVibration()
{
    int num_of_modes = GetNumOfModes();
    // solve mode vibration
    tVector3d hammer_force = tVector3d::Ones() * this->mHammerForce;
    tVectorXd UTFinit =
        mEigenVecs.transpose().block(0, 0, num_of_modes, 3) * hammer_force;
    /*
        modal vibration equation:
        qddot + (a  + b * eigen_val) qdot + eigen_val * q = UTf_init

        for mode j, solution:
        q_j(t)  = dt * Utf_init_j / (w_d) e^{- \xi * w * t} sin(w_d * t)
                = coef_j / w_d e^{- \xi * w * t} sin(w_d * t)
        we need to calculate
            1. coef_j
            3. w = \sqrt(eigen_val)
            4. xi = (a + b * eigen_val) / (2 * w)
            2. w_d = w * (1 - xi^2)
    */
    double dt = 1.0 / mAcousticSamplingFreqHZ;
    mModalVibrationsInfo.clear();
    mModalWaves.clear();
    mIsModeValid.clear();
    for (int mode_id = 0; mode_id < num_of_modes; mode_id++)
    {
        double coef_j = dt * UTFinit[mode_id];
        double eigen_val = mEigenValues[mode_id];
        double w = std::sqrt(eigen_val);
        bool failed_or_invisibled = false;
        if (w < 2 * M_PI * 20 || w > 2 * M_PI * 20000)
        {
            printf("[debug] mode %d w %.1f is invisible, ignore\n", mode_id, w);
            failed_or_invisibled = true;
        }
        if (std::isnan(w))
            failed_or_invisibled = true;
        double xi =
            (this->mMaterial->mRayleighDamplingA +
             this->mMaterial->mRayleighDamplingB * mEigenValues[mode_id]) /
            (2 * w);
        if (!failed_or_invisibled && xi >= 1)
        {
            SIM_WARN("mode {} xi {} >=1!, system is overdamped set to 1",
                     mode_id, xi);
            xi = 1;
            failed_or_invisibled = true;
        }

        double wd = w * std::sqrt(1 - xi * xi);

        printf("mode %d w %.1f xi %.2f wd %.1f\n", mode_id, w, xi, wd);
        if (failed_or_invisibled == true)
        {
            // mModalVibrationsInfo.push_back(tVector(0, 0, 0, 0));
            printf("invalid mode, set to zero\n");
            coef_j = 0, w = 0, xi = 0, wd = 0;
            mIsModeValid.push_back(false);
        }
        else
        {
            mIsModeValid.push_back(true);
        }

        mModalVibrationsInfo.push_back(tVector(coef_j, w, xi, wd));

        tVectorXd data = CalcWave(coef_j, w, xi, wd, dt,
                                  mAcousticDuration * mAcousticSamplingFreqHZ);
        auto mode_wave = std::make_shared<tDiscretedWave>(dt);
        mode_wave->SetData(data.cast<float>());
        mModalWaves.push_back(mode_wave);
    }
}

void cAcousticSoftBody::Update(float dt)
{
    // cSoftBodyImplicit::Update(dt);

    // update second
}

void cAcousticSoftBody::SolveForMonopole()
{
    // 1. get outer surface vertex, and accel
    // 2. set a monopole place, calculate the residual
    // 3. output
}
#include "utils/FileUtil.h"
#include "utils/StringUtil.h"
std::string cAcousticSoftBody::GetWaveName() const
{
    // 1. mesh name
    auto mesh_name = cFileUtil::GetFilename(mTetMeshPath);
    // damping_a_b

    auto mat = GetMaterial();
    // material name
    std::string mat_path = mat->GetMaterialPath();
    std::cout << "mat path = " << mat_path << std::endl;
    std::string mat_name =
        cFileUtil::RemoveExtension(cFileUtil::GetFilename(mat_path));
    std::cout << "mat name = " << mat_name << std::endl;
    // sampling freq
    std::string name = "data/wave_" + mesh_name + "_" + mat_name;
    return name;
}

tVectorXd CalcAccel(const tVectorXd &x, float dt)
{
    tVectorXd accel(x.size() - 2);
    for (int i = 1; i < x.size() - 1; i++)
    {
        accel[i - 1] = (x[i + 1] + x[i - 1] - 2 * x[i]) / (dt * dt);
    }
    return accel;
}
tVectorXd SolveModeITimeSeries(double m, double c, double k, double fprev,
                               double dt, double steps)
{

    if (k < 1e-6)
        return tVectorXd::Zero(steps);

    /*
            w = sqrt(k / m)
            xi = c / (2 * m * w)
            wd = w * sqrt(1 - xi^2)
            eps = e^{- xi * w * h}
            theta = w_d * h
            gamma = arcsin(xi)
        */
    double w = std::sqrt(k / m);
    double xi = c / (2 * m * w);
    double wd = w * std::sqrt(1 - xi * xi);
    // std::cout << "mode freq = " << w << std::endl;
    // std::cout << "wd = " << wd << std::endl;

    double qprev = 0;

    tVectorXd q_vec = tVectorXd::Zero(steps);
    if (std::isnan(wd) == false)
    {

        double eps = std::exp(-xi * w * dt);
        double theta = wd * dt;
        double gamma = std::asin(xi);
        double qpprev = 0;

        for (size_t t = 0; t < steps - 1; t++)
        {
            /*
            xcur = 2 * eps * cos(theta) * qprev
                    - eps^2 * qpprev
                    + 2 * fprev *
                        (
                            eps * cos(theta + gamma)
                            - eps^2 * cos(2 * theta + gamma)
                        )
                        /
                        (3 * m * w * wd)
            */
            double qcur = 2 * eps * std::cos(theta) * qprev -
                          eps * eps * qpprev +
                          2 * fprev *
                              (eps * std::cos(theta + gamma) -
                               eps * eps * std::cos(2 * theta + gamma)) /
                              (3 * m * w * wd);

            // update qprev and qcur, fprev
            qpprev = qprev;
            qprev = qcur;

            // the force is applied at first
            fprev = 0;

            q_vec[t + 1] = qcur;
            if (std::isnan(qcur))
            {
                std::cout << "m " << m << std::endl;
                std::cout << "w " << w << std::endl;
                std::cout << "wd " << wd << std::endl;
                std::cout << "xi " << xi << std::endl;
                std::cout << "c " << c << std::endl;
                exit(1);
            }
        }
    }
    else
    {
        printf("[warn] w %.1f wd %.1f nan, return\n", w, wd);
    }

    return q_vec;
}
void cAcousticSoftBody::EigenDecompose(tVectorXd &eigenValues,
                                       tMatrixXd &eigenVecs) const
{
    {
        Eigen::GeneralizedEigenSolver<tMatrixXd> ges;

        int num_of_dof = mInvLumpedMassMatrixDiag.size();
        // tMatrixXd dense_M = tMatrixXd::Zero(num_of_dof, num_of_dof);
        // dense_M.diagonal() = mInvLumpedMassMatrixDiag.cwiseInverse();
        tMatrixXd dense_M = mRawMassMatrix.toDense();

        tMatrixXd dense_K;
        if (this->mAcousticMaterialModelIdx == 2)
        {
            printf("[info] take linear elasticity's K matrix\n");
            dense_K = mLinearElasticityStiffnessMatrix;
        }
        else
        {
            printf("[info] take hyper elasticity's K matrix\n");
            dense_K = -mGlobalStiffnessSparseMatrix.toDense();
        }

        ges.compute(dense_K, dense_M);
        eigenValues = ges.eigenvalues().real();
        eigenVecs = ges.eigenvectors().real().transpose();
        // std::cout << "eigenvalues = " << eigenvalues.transpose() <<
        // std::endl;
    }
}

tDiscretedWavePtr cAcousticSoftBody::CalculateVertexVibration(int v_id)
{
    int num_of_modes = this->GetNumOfModes();
    tVectorXd weight = mEigenVecs.row(v_id).segment(0, num_of_modes);
    tVectorXf new_data = tVectorXf::Zero(mModalWaves[0]->data.size());
    for (int i = 0; i < num_of_modes; i++)
    {
        if (mIsModeValid[i] == true)
            new_data += weight[i] * mModalWaves[i]->data;
        else
        {
            printf("mode %d is invalid\n");
        }
    }

    // scale sound
    new_data = 1.0 / new_data.cwiseAbs().maxCoeff() * new_data;
    auto mode_wave = std::make_shared<tDiscretedWave>(GetDt());
    mode_wave->SetData(new_data);
    return mode_wave;
}

double cAcousticSoftBody::GetDt() const
{
    return 1.0 / mAcousticSamplingFreqHZ;
}

double cAcousticSoftBody::CalcTotalMass() const
{
    std::cout << mInvLumpedMassMatrixDiag.cwiseInverse().transpose()
              << std::endl;
    return mInvLumpedMassMatrixDiag.cwiseInverse().sum();
}

#include "sim/softbody/MaterialBuilder.h"

void cAcousticSoftBody::UpdateImGui()
{
    // material model selection
    if (ImGui::BeginCombo(
            "material model idx",
            GetAcousticMaterialNameFromIdx(mAcousticMaterialModelIdx).c_str()))
    {
        for (int i = 0; i < 3; i++)
        {
            std::string name = GetAcousticMaterialNameFromIdx(i);
            bool is_selected = (i == mAcousticMaterialModelIdx);

            if (ImGui::Selectable(name.c_str(), is_selected))
            {
                if (mAcousticMaterialModelIdx != i)
                {
                    // change material
                    ChangeMaterial(mAcousticMaterialModelIdx, i);
                }
            }
            if (is_selected)
            {
                ImGui::SetItemDefaultFocus();
            }
        }
        ImGui::EndCombo();
    }

    // material model selection
    if (ImGui::BeginCombo(
            "material param idx",
            this->mMaterialParamNameLst[this->mMaterialParamIdx].c_str()))
    {
        for (int i = 0; i < mMaterialParamNameLst.size(); i++)
        {
            std::string name = mMaterialParamNameLst[i];
            bool is_selected = (i == mMaterialParamIdx);

            if (ImGui::Selectable(name.c_str(), is_selected))
            {
                if (mMaterialParamIdx != i)
                {
                    // change material
                    mMaterialParamIdx = i;
                    mMaterial->Init(mMaterialParamPathLst[mMaterialParamIdx]);
                    ResolveVibration();
                }
            }
            if (is_selected)
            {
                ImGui::SetItemDefaultFocus();
            }
        }
        ImGui::EndCombo();
    }
    if (ImGui::Button("resolve vibration\n"))
    {
        ResolveVibration();
    }
}

std::string GetAcousticMaterialNameFromHyperMaterial(eMaterialType type)
{
    if (type == eMaterialType::STVK)
        return "stvk";
    else if (type == eMaterialType::NEO_HOOKEAN)
        return "neohookean";
    return "";
}
eMaterialType GetHyperMaterialFromAcousticMaterial(std::string name)
{
    if (name == "stvk")
        return eMaterialType::STVK;
    else if (name == "neohookean")
        return eMaterialType::NEO_HOOKEAN;
}
#include "sim/softbody/NeoHookeanMaterial.h"
#include "sim/softbody/StvkMaterial.h"
void cAcousticSoftBody::ChangeMaterial(int old_idx, int new_idx)
{
    std::string old_name = GetAcousticMaterialNameFromIdx(old_idx);
    std::string new_name = GetAcousticMaterialNameFromIdx(new_idx);
    bool old_is_hyper = (old_name == "stvk") || (old_name == "neohookean");
    bool new_is_hyper = (new_name == "stvk") || (new_name == "neohookean");

    bool enable_hyper_convert = false;
    if (old_is_hyper == true)
    {
        if (new_is_hyper)
        {
            // hyper convert
            enable_hyper_convert = true;
        }
        else
        {
            // new is linear
            // do nothing
        }
    }
    else
    {
        // old is linear
        // new must be hyper
        if (new_name != GetAcousticMaterialNameFromHyperMaterial(
                            this->mMaterial->GetType()))
        {
            // hyper convert
            enable_hyper_convert = true;
        }
        else
        {
            // do nothing
        }
    }

    if (true == enable_hyper_convert)
    {
        eMaterialType new_hyper_type =
            GetHyperMaterialFromAcousticMaterial(new_name);
        ;
        cBaseMaterialPtr new_mat = nullptr;
        if (new_hyper_type == eMaterialType::STVK)
        {
            new_mat = std::make_shared<cStvkMaterial>();
        }
        else if (new_hyper_type == eMaterialType::NEO_HOOKEAN)
        {
            new_mat = std::make_shared<cNeoHookeanMaterial>();
        }
        new_mat->Init(mMaterial->mMatPath);
        mMaterial = new_mat;
    }
    mAcousticMaterialModelIdx = new_idx;

    // recalculate vibration
    ResolveVibration();
}

/**
 * \brief       load material params
 */
void cAcousticSoftBody::LoadMaterialParams()
{
    SIM_ASSERT(mMaterial != nullptr);
    // 1. get root dir
    std::string material_params_rootdir =
        cFileUtil::GetDir(mMaterial->GetMaterialPath());
    mMaterialParamNameLst.clear();
    mMaterialParamPathLst.clear();
    mMaterialParamIdx = -1;
    auto mat_paths = cFileUtil::ListDir(material_params_rootdir);
    for (int i = 0; i < mat_paths.size(); i++)
    {
        // 2. list dir, load json; load their name
        Json::Value root;
        auto mat_path = mat_paths[i];
        cJsonUtil::LoadJson(mat_path, root);
        std::string mat_name = root["name"].asString();
        mMaterialParamNameLst.push_back(mat_name);
        mMaterialParamPathLst.push_back(mat_path);
        if (mMaterial->mName == mat_name)
        {
            mMaterialParamIdx = i;
        }
    }

    // 3. judge current
}