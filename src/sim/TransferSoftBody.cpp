#include "sim/TransferSoftBody.h"
#include "geometries/Primitives.h"
#include "utils/FileUtil.h"
#include "utils/JsonUtil.h"
#include <iostream>

cTransferSoftBody::cTransferSoftBody(int id_) : cSoftBodyImplicit(id_) {}

void cTransferSoftBody::Init(const Json::Value &conf)
{
    cSoftBodyImplicit::Init(conf);
    mModalAnlysisResultPath =
        cJsonUtil::ParseAsString("modal_analysis_result", conf);
    mEnableDumpedModalSolution = false;
    mNumOfMonopoles = 1;
    InitSurfaceNormal();

    LoadModalAnalysisResult();

    int num_of_modes = mSurfaceVertexDOFCoef.cols();

    for (int i = 0; i < num_of_modes; i++)
    {
        SolveMonopole(i);
    }
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
    int num_of_modes = mModalVibrationsInfo.size();
    int num_of_surface_id = mSurfaceVertexDOFCoef.size() / 3;
    tVector mode_info = mModalVibrationsInfo[mode_idx];
    tVectorXd dof_coef_for_modes = mSurfaceVertexDOFCoef.col(mode_idx);
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
    mModalVibrationsInfo.clear();
    // 2. parse modes
    {
        for (auto &info : modes)
        {
            tVector mode_vib_info = cJsonUtil::ReadVectorJson(info);
            mModalVibrationsInfo.push_back(mode_vib_info);
        }
        std::cout << "num of valid modes = " << mModalVibrationsInfo.size()
                  << std::endl;
    }
    // 3. parse weight
    {
        mSurfaceVertexDOFCoef.resize(3 * num_of_surface_v,
                                     mModalVibrationsInfo.size());

        for (int i = 0; i < coef.size(); i++)
        {
            mSurfaceVertexDOFCoef.row(i) = cJsonUtil::ReadVectorJson(coef[i]);
        }
    }
}

/**
 * \brief       calculate sound pressure
 */
void cTransferSoftBody::CalcSoundPressure() 
{

}