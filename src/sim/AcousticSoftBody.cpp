#include "sim/AcousticSoftBody.h"
#include "geometries/Primitives.h"
#include "sim/AudioOutput.h"
#include "utils/FileUtil.h"
#include "utils/JsonUtil.h"
#include <iostream>
extern cAudioOutputPtr gAudioOutput;
cAcousticSoftBody::cAcousticSoftBody(int id_) : cSoftBodyImplicit(id_) {}

void cAcousticSoftBody::Init(const Json::Value &conf)
{
    cSoftBodyImplicit::Init(conf);
    mAcousticSamplingFreq =
        cJsonUtil::ParseAsInt("acoustic_sampling_freq", conf);
    mAcousticDuration = cJsonUtil::ParseAsFloat("acoustic_duration", conf);
    printf("[log] acoustic sampling freq %d HZ, duration %.1f s\n",
           mAcousticSamplingFreq, mAcousticDuration);

    // Init noise
    UpdateDeformationGradient();
    mGlobalStiffnessSparseMatrix = CalcGlobalStiffnessSparseMatrix();

    // check whether this wave exists or not. if exist, load; else calculate and
    // dump
    std::string wave_name = GetWaveName();
    tDiscretedWavePtr wave = nullptr;
    // if (cFileUtil::ExistsFile(wave_name) == false)
    // {
    // for (int i = 0; i < this->GetNumOfTets(); i++)
    // {
    //     CheckElementStiffnessMat(i);
    // }
    // exit(1);
    tEigenArr<tVector> normal_array = {};
    for (auto &v : this->mVertexArray)
    {
        normal_array.push_back(v->mNormal);
    }

    wave =
        SolveVibration(
            normal_array,
            mInvLumpedMassMatrixDiag.cwiseInverse(),
                       mGlobalStiffnessSparseMatrix, GetRayleightDamping(),
                       mXcur, mXprev, mAcousticSamplingFreq, mAcousticDuration);

    //     wave->DumpToFile(wave_name);
    // }
    // else
    // {
    //     wave = std::make_shared<tDiscretedWave>(1e-3, 1);
    //     wave->LoadFromFile(wave_name);
    // }

    gAudioOutput->SetWave(wave);
    // SolveForMonopole();
    // exit(1);
    // sum each vertex's displacement
}

void cAcousticSoftBody::Update(float dt)
{
    cSoftBodyImplicit::Update(dt);

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
    // sampling freq
    std::string name = "wave_" + mesh_name + "_" +
                       std::to_string(mRayleighDamplingA) + "_" +
                       std::to_string(mRayleighDamplingB) + "_" +
                       std::to_string(mat->GetPoissonRatio()) + "_ " +
                       std::to_string(mat->GetYoungsModulus()) + "_" +
                       std::to_string(this->mRho) + "_" +
                       std::to_string(mAcousticSamplingFreq) + "_" +
                       std::to_string(this->mAcousticDuration);
    return name;
}