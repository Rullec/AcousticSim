#include "sim/acoustic/ModalSoftBody.h"
#include "ModeVibration.h"
#include "geometries/Primitives.h"
#include "imgui.h"
#include "sim/acoustic/AudioOutput.h"
#include "sim/acoustic/AudioWave.h"
#include "sim/softbody/AcousticLinearElasticity.h"
#include "sim/softbody/BaseMaterial.h"
#include "utils/FileUtil.h"
#include "utils/JsonUtil.h"
#include <fstream>
#include <iostream>
#include <set>
extern cAudioOutputPtr gAudioOutput;
eMaterialType GetHyperMaterialFromAcousticMaterial(std::string name);
std::string GetAcousticMaterialNameFromHyperMaterial(eMaterialType type);

tDiscretedWavePtr
CalculateVertexVibration(int v_id, int num_of_selected_modes,
                         const std::vector<tDiscretedWavePtr> &mModalWaves,
                         const tMatrixXd &mVertexModesCoef, double sampling_hz);
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

cModalSoftBody::cModalSoftBody(int id_) : cSoftBodyImplicit(id_) {}

void cModalSoftBody::Init(const Json::Value &conf)
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
    mSelectedVertexForHearing = 0;
    mEnableDrawVertexNormal = false;
    mEnableDrawTriNormal = true;

    // confirm normal
    InitArrowFromNormal();

    // load many material types ... for selection
    LoadMaterialParams();
    // std::cout << "linear K = \n"
    //           << mLinearElasticityStiffnessMatrix.toDense() << std::endl;
    ResolveModalVibration();
    mNumOfSelectedModes = this->mModalWaves.size() / 2;
    GenerateSound();
}
void cModalSoftBody::ResolveModalVibration()
{
    printf("[info] ------------resolve vibration-------------\n");
    // Init noise
    UpdateDeformationGradient();
    mGlobalStiffnessSparseMatrix = CalcGlobalStiffnessSparseMatrix();

    // init linear elasticity stiffness matrix
    mLinearElasticityStiffnessMatrix =
        cAcousticLinearElasticity::CalcGlobalStiffness(
            mVertexArray, mTetArrayShared, this->mMaterial->mYoungsModulusNew,
            this->mMaterial->mPoissonRatioNew,
            this->mMaterial->mRayleighDamplingA,
            this->mMaterial->mRayleighDamplingB);

    tVectorXd eigen_vals;
    tMatrixXd eigen_vecs;
    EigenDecompose(eigen_vals, eigen_vecs);

    CalculateModalVibration(eigen_vals, eigen_vecs);
}
void cModalSoftBody::GenerateSound()
{
    auto play_wave = CalculateVertexVibration(
        mSelectedVertexForHearing, mNumOfSelectedModes, mModalWaves,
        mVertexModesCoef, mAcousticSamplingFreqHZ);
    gAudioOutput->SetWave(play_wave);

    // save
    std::string mat_name = cFileUtil::RemoveExtension(
        cFileUtil::GetFilename(this->mMaterial->mMatPath));

    play_wave->DumpToWAV(cFileUtil::ConcatFilename("log/", mat_name + ".wav"));
}

/**
 * \brief   calcualte modal vibration
 */
void cModalSoftBody::CalculateModalVibration(const tVectorXd &eigen_vals,
                                             const tMatrixXd &eigen_vecs)
{
    int num_of_dof = eigen_vals.size();
    // solve mode vibration
    tVector3d hammer_force = tVector3d::Ones() * this->mHammerForce;
    tVectorXd UTFinit =
        eigen_vecs.transpose().block(0, 0, num_of_dof, 3) * hammer_force;

    // 1. remove surface rows in eigen vectors (prepare for coef extraction)
    int num_of_surface_v = mSurfaceVertexIdArray.size();
    tMatrixXd surface_vertex_vecs =
        tMatrixXd::Zero(3 * num_of_surface_v, eigen_vecs.cols());
    for (int i = 0; i < num_of_surface_v; i++)
    {
        for (int j = 0; j < 3; j++)
            surface_vertex_vecs.row(3 * i + j) =
                eigen_vecs.row(3 * mSurfaceVertexIdArray[i] + j);
    }

    double dt = 1.0 / mAcousticSamplingFreqHZ;
    mModalVibrationsInfo.clear();
    mModalWaves.clear();

    // 2. solve the solution for each mode
    std::vector<int> valid_dof_array = {};
    for (int mode_id = 0; mode_id < num_of_dof; mode_id++)
    {
        double coef_j = dt * UTFinit[mode_id];
        double eigen_val = eigen_vals[mode_id];
        double w = std::sqrt(eigen_val);
        bool failed_or_invisibled = false;
        if (w < 2 * M_PI * 20 || w > 2 * M_PI * 20000)
        {
            printf("[debug] mode %d w %.1f is invisible, ignore\n", mode_id, w);
            failed_or_invisibled = true;
        }
        if (std::isnan(w))
            failed_or_invisibled = true;
        double xi = (this->mMaterial->mRayleighDamplingA +
                     this->mMaterial->mRayleighDamplingB * eigen_val) /
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
            continue;
        }
        else
        {
            // 3. get the valid & visible modes
            auto cur_vib = std::make_shared<tModeVibration>(coef_j, w, xi, wd);
            mModalVibrationsInfo.push_back(cur_vib);

            valid_dof_array.push_back(mode_id);
            tVectorXd data =
                cur_vib->GetWave(mAcousticDuration, mAcousticSamplingFreqHZ);
            auto mode_wave = std::make_shared<tDiscretedWave>(dt);
            mode_wave->SetData(data.cast<float>());
            mModalWaves.push_back(mode_wave);
        }
    }

    // 4. get the coef for valid id
    int num_of_valid_dof = valid_dof_array.size();
    mVertexModesCoef.resize(3 * num_of_surface_v, num_of_valid_dof);
    for (int i = 0; i < num_of_valid_dof; i++)
    {
        mVertexModesCoef.col(i) = surface_vertex_vecs.col(valid_dof_array[i]);
    }
    printf("[info] get surface vertex %d/%d, valid modes %d/%d\n",
           num_of_surface_v, mVertexArray.size(), num_of_valid_dof, num_of_dof);
}

void cModalSoftBody::Update(float dt)
{
    // cSoftBodyImplicit::Update(dt);

    // update second
}

void cModalSoftBody::SolveForMonopole()
{
    // 1. get outer surface vertex, and accel
    // 2. set a monopole place, calculate the residual
    // 3. output
}
#include "utils/FileUtil.h"
#include "utils/StringUtil.h"
std::string cModalSoftBody::GetWaveName() const
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
void cModalSoftBody::EigenDecompose(tVectorXd &eigenValues,
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

tDiscretedWavePtr
CalculateVertexVibration(int v_id, int num_of_selected_modes,
                         const std::vector<tDiscretedWavePtr> &mModalWaves,
                         const tMatrixXd &mVertexModesCoef, double sampling_hz)
{
    // total num of modes
    int num_of_modes = mModalWaves.size();

    tVectorXd weight = mVertexModesCoef.row(v_id);
    tVectorXf new_data = tVectorXf::Zero(mModalWaves[0]->data.size());

    // 1. begin to select modes
    std::set<int> selected_modes = {};
    if (num_of_selected_modes > num_of_modes)
    {
        printf("[log] use full %d modes to build sound\n", num_of_modes);

        for (int i = 0; i < num_of_modes; i++)
            selected_modes.insert(i);
    }
    else
    {
        SIM_ASSERT(num_of_selected_modes >= 2);
        double stepsize = num_of_modes * 1.0 / (num_of_selected_modes - 1);
        int st = 0;
        for (int i = 0; i < num_of_selected_modes; i++)
        {
            int cur_selection = int(st + i * stepsize);
            cur_selection =
                cMathUtil::Clamp(cur_selection, 0, num_of_modes - 1);
            selected_modes.insert(cur_selection);
        }
        printf("[info] get %d modes where target %d modes\n",
               selected_modes.size(), num_of_modes);
        while (selected_modes.size() < num_of_selected_modes)
        {
            selected_modes.insert(
                cMathUtil::RandInt(0, num_of_selected_modes - 1));
        }

        // generate sound
    }
    printf("[info] generate sound for vertex %d\n", v_id);
    for (auto &i : selected_modes)
    {
        // std::cout << "mode " << i << " weight = " << weight[i] << " data = "
        //           << mModalWaves[i]->data.segment(0, 2).transpose()
        //           << std::endl;
        new_data += weight[i] * mModalWaves[i]->data;
    }
    printf("[info] generate sound from %d modes\n", selected_modes.size());
    // 2. get result

    // scale sound
    new_data = 1.0 / new_data.cwiseAbs().maxCoeff() * new_data;
    auto mode_wave = std::make_shared<tDiscretedWave>(1.0 / sampling_hz);
    mode_wave->SetData(new_data);
    return mode_wave;
}

double cModalSoftBody::GetDt() const { return 1.0 / mAcousticSamplingFreqHZ; }

double cModalSoftBody::CalcTotalMass() const
{
    std::cout << mInvLumpedMassMatrixDiag.cwiseInverse().transpose()
              << std::endl;
    return mInvLumpedMassMatrixDiag.cwiseInverse().sum();
}

#include "sim/softbody/MaterialBuilder.h"

void cModalSoftBody::UpdateImGui()
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
                    ResolveModalVibration();
                    GenerateSound();
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
        ResolveModalVibration();
        GenerateSound();
    }
    if (ImGui::Button("dump vibration"))
    {
        DumpResultForTransfer();
    }
    // select vid
    {
        int old_v = mSelectedVertexForHearing;
        ImGui::SliderInt("hearing vertex", &mSelectedVertexForHearing, 0,
                         mVertexArray.size() - 1);
        if (mSelectedVertexForHearing != old_v)
        {
            GenerateSound();
        }
    }

    // option
    ImGui::Checkbox("draw tri normal", &mEnableDrawTriNormal);
    ImGui::Checkbox("draw vertex normal", &mEnableDrawVertexNormal);

    int num_of_modes = this->mModalWaves.size();
    ImGui::Text("num of modes (solved): %d", num_of_modes);

    {
        int old_val = mNumOfSelectedModes;
        ImGui::DragInt("num of selected modes", &old_val, 1.0, 2, num_of_modes);

        old_val = SIM_MAX(old_val, 2);
        if (old_val != mNumOfSelectedModes)
        {
            mNumOfSelectedModes = old_val;
            GenerateSound();
        }
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
void cModalSoftBody::ChangeMaterial(int old_idx, int new_idx)
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
    ResolveModalVibration();
    GenerateSound();
}

/**
 * \brief       load material params
 */
void cModalSoftBody::LoadMaterialParams()
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

void cModalSoftBody::DumpResultForTransfer()
{
    Json::Value geometry, modes, coef;
    // 1. geo info: surface vertex pos, id, normal vector
    {
        int num_of_surface_v = mSurfaceVertexIdArray.size();
        Json::Value v_lst = Json::arrayValue, pos_lst = Json::arrayValue,
                    normal_lst = Json::arrayValue;
        for (int i = 0; i < num_of_surface_v; i++)
        {
            int vid = mSurfaceVertexIdArray[i];
            auto v = mVertexArray[vid];

            tVector3d pos = v->mPos.segment(0, 3);
            tVector3d normal = v->mNormal.segment(0, 3);
            v_lst.append(vid);
            pos_lst.append(cJsonUtil::BuildVectorJson(pos));
            normal_lst.append(cJsonUtil::BuildVectorJson(normal));
        }
        geometry["tet_path"] = mTetMeshPath;
        geometry["vertex_id_lst"] = v_lst;
        geometry["pos_lst"] = pos_lst;
        geometry["normal_lst"] = normal_lst;
    }

    // 2. modes

    {
        modes = Json::arrayValue;
        int num_of_modes = this->mModalVibrationsInfo.size();
        for (auto &info : mModalVibrationsInfo)
        {
            modes.append(cJsonUtil::BuildVectorJson(info->GetCoefVec()));
        }
    }

    // 3. coef of the modes for each info
    {
        coef = Json::arrayValue;
        for (auto &surface_v_id : mSurfaceVertexIdArray)
        {
            // ! need to get all eigen vecs, it seems current dump is lacked of
            // DOF!
            coef.append(cJsonUtil::BuildVectorJson(
                mVertexModesCoef.row(3 * surface_v_id + 0)));
            coef.append(cJsonUtil::BuildVectorJson(
                mVertexModesCoef.row(3 * surface_v_id + 1)));
            coef.append(cJsonUtil::BuildVectorJson(
                mVertexModesCoef.row(3 * surface_v_id + 2)));
        }
    }
    // export name: "vib_" + material_param_name + material_name
    std::string output =
        "data/vib_" + this->mMaterialParamNameLst[this->mMaterialParamIdx] +
        "_" + GetAcousticMaterialNameFromIdx(this->mAcousticMaterialModelIdx) +
        ".json";
    Json::Value root;
    root["geometry"] = geometry;
    root["modes"] = modes;
    root["coef"] = coef;
    cJsonUtil::WriteJson(output, root);
    std::cout << "dump to " << output << std::endl;
}
#include "geometries/Arrow.h"
void cModalSoftBody::InitArrowFromNormal()
{
    // only surface vertex
    tVector min, max;
    CalcAABB(min, max);
    tVector size = max - min;
    double obj_scale = size.segment(0, 3).mean() / 20;

    mDrawVertexNormalArray.clear();
    for (auto &vid : this->mSurfaceVertexIdArray)
    {
        auto arrow = std::make_shared<cArrow>();
        arrow->SetStEd(mVertexArray[vid]->mPos,
                       mVertexArray[vid]->mPos +
                           mVertexArray[vid]->mNormal * obj_scale);
        mDrawVertexNormalArray.push_back(arrow);
    }

    // only surface triangles
    mDrawTriNormalArray.clear();
    for (auto &tid : this->mSurfaceTriangleIdArray)
    {
        tVector normal = mTriangleArray[tid]->mNormal;
        tVector pos = (mVertexArray[mTriangleArray[tid]->mId0]->mPos +
                       mVertexArray[mTriangleArray[tid]->mId1]->mPos +
                       mVertexArray[mTriangleArray[tid]->mId2]->mPos) /
                      3;

        auto arrow = std::make_shared<cArrow>();
        arrow->SetStEd(pos, pos + normal * obj_scale);
        mDrawTriNormalArray.push_back(arrow);
    }
}

int cModalSoftBody::GetNumOfDrawTriangles() const
{
    int raw = cSoftBodyImplicit::GetNumOfDrawTriangles();
    int num_of_arrow_tri = mDrawVertexNormalArray[0]->GetNumOfDrawTriangles();
    if (mEnableDrawVertexNormal)
    {
        raw += mDrawVertexNormalArray.size() * num_of_arrow_tri;
    }
    if (mEnableDrawTriNormal)
    {
        raw += mDrawTriNormalArray.size() * num_of_arrow_tri;
    }
    return raw;
}
int cModalSoftBody::GetNumOfDrawEdges() const
{
    int raw = cSoftBodyImplicit::GetNumOfDrawEdges();
    int num_of_arrow_edge = mDrawVertexNormalArray[0]->GetNumOfDrawEdges();
    if (mEnableDrawVertexNormal)
    {
        raw += mDrawVertexNormalArray.size() * num_of_arrow_edge;
    }
    if (mEnableDrawTriNormal)
    {
        raw += mDrawTriNormalArray.size() * num_of_arrow_edge;
    }
    return raw;
}
int cModalSoftBody::GetNumOfDrawVertices() const
{
    int raw = cSoftBodyImplicit::GetNumOfDrawVertices();
    return raw;
}

void cModalSoftBody::CalcTriangleDrawBuffer(Eigen::Map<tVectorXf> &res,
                                            int &st) const
{
    cSoftBodyImplicit::CalcTriangleDrawBuffer(res, st);
    if (mEnableDrawVertexNormal)
    {
        for (auto &v : this->mDrawVertexNormalArray)
        {
            v->CalcTriangleDrawBuffer(res, st);
        }
    }
    if (mEnableDrawTriNormal)
    {
        for (auto &v : this->mDrawTriNormalArray)
        {
            v->CalcTriangleDrawBuffer(res, st);
        }
    }
}

void cModalSoftBody::CalcEdgeDrawBuffer(Eigen::Map<tVectorXf> &res,
                                        int &st) const
{
    cSoftBodyImplicit::CalcEdgeDrawBuffer(res, st);
    if (mEnableDrawVertexNormal)
    {
        for (auto &v : this->mDrawVertexNormalArray)
        {
            v->CalcEdgeDrawBuffer(res, st);
        }
    }
    if (mEnableDrawTriNormal)
    {
        for (auto &v : this->mDrawTriNormalArray)
        {
            v->CalcEdgeDrawBuffer(res, st);
        }
    }
}