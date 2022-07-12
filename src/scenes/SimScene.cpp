#include "SimScene.h"
// #include "sim/collision/CollisionDetecter.h"
#include "geometries/Primitives.h"
#include "geometries/Triangulator.h"
#include "imgui.h"
#include "scenes/SimStateMachine.h"
#include "sim/kinematic/kinematicBody.h"
#include "sim/Perturb.h"
#include "sim/SimObjectBuilder.h"
#include "sim/collision/BVHCollisionDetecter.h"
#include "utils/ColorUtil.h"
#include "utils/JsonUtil.h"
#include <iostream>
// std::string gSceneTypeStr[eSceneType::NUM_OF_SCENE_TYPES] = {"sim",
// "acoustic"};
std::string gSceneTypeStr[eSceneType::NUM_OF_SCENE_TYPES] = {"sim"};

eSceneType cSimScene::BuildSceneType(const std::string &str)
{
    int i = 0;
    for (i = 0; i < eSceneType::NUM_OF_SCENE_TYPES; i++)
    {
        // std::cout << gSceneTypeStr[i] << std::endl;
        if (str == gSceneTypeStr[i])
        {
            break;
        }
    }

    SIM_ASSERT(i != eSceneType::NUM_OF_SCENE_TYPES);
    return static_cast<eSceneType>(i);
}

cSimScene::cSimScene()
{
    // mTriangleArray.clear();
    // mEdgeArray.clear();
    // mVertexArray.clear();

    mPerturb = nullptr;
    mSimStateMachine = std::make_shared<tSimStateMachine>();
    mSceneType = eSceneType::SCENE_SIM;
    mCamPos.setZero();
    // mColDetecter = nullptr;
}

eSceneType cSimScene::GetSceneType() const { return this->mSceneType; }
#include "sim/acoustic/AudioOutput.h"

void cSimScene::Init(const std::string &conf_path)
{
    // 1. load config
    Json::Value root;
    cJsonUtil::LoadJson(conf_path, root);

    mEnableProfiling =
        cJsonUtil::ParseAsBool(cSimScene::ENABLE_PROFLINE_KEY, root);

    mEnableCollisionDetection =
        cJsonUtil::ParseAsBool(cSimScene::ENABLE_COLLISION_DETECTION_KEY, root);
    BuildObjects(cJsonUtil::ParseAsValue(cSimScene::OBJECT_LIST_KEY, root));

    CreateCollisionDetecter();

    UpdateRenderingResource();
    mIsRenderingResourceSizeUpdated = false;
    InitRaycaster(root);

    mSimStateMachine->SimulatorInitDone(
        cJsonUtil::ParseAsBool("pause_at_first", root) == true
            ? eSimState::SIMSTATE_PAUSE
            : eSimState::SIMSTATE_RUN);
}
void cSimScene::BuildObjects(const Json::Value &obj_conf_)
{
    // 1. parse the number of obstacles
    Json::Value obj_conf = obj_conf_;
    int num_of_obstacles = obj_conf.size();
    SIM_ASSERT(num_of_obstacles == obj_conf.size());
    for (int i = 0; i < num_of_obstacles; i++)
    {
        auto obs = BuildSimObj(obj_conf[i], mObjectList.size());
        // auto obs = BuildKinematicBody(obj_conf[i], GetNumOfObjects());
        mObjectList.push_back(obs);
    }
}

/**
 * \breif       save current scene (obstacles to objes)
 */
#include "geometries/ObjExport.h"
void cSimScene::SaveCurrentScene()
{
    for (auto &x : this->mObjectList)
    {
        cObjExporter::ExportObj(x->GetObjName() + ".obj", x->GetVertexArray(),
                                x->GetTriangleArray());
    }
}

/**
 * \brief           Init the raycasting strucutre
 */
#include "geometries/Raycaster.h"
void cSimScene::InitRaycaster(const Json::Value &conf)
{
    mRaycaster = std::make_shared<cRaycaster>();
    mRaycaster->Init(conf);
    for (auto &x : mObjectList)
    {
        mRaycaster->AddResources(x);
    }
    std::cout << "[debug] add resources to raycaster done, num of objects = "
              << mObjectList.size() << std::endl;
}
/**
 * \brief           Update the simulation procedure
 */
#include "utils/TimeUtil.hpp"
static int frame = 0;
void cSimScene::Update(double delta_time)
{

    if (mSimStateMachine->IsRunning() == true)
    {
        // std::cout << "--------------frame " << frame++ << "-----------\n";
        // double default_dt = mIdealDefaultTimestep;
        // if (delta_time < default_dt)
        //     default_dt = delta_time;

        // printf("[debug] sim scene update cur time = %.4f\n", mCurTime);
        cScene::Update(delta_time);

        UpdateObjects();
        PerformCollisionDetection();
        // clear force
        // apply ext force
        // update position
    }
    mSimStateMachine->SimulationForwardOneFrameDone();
}

/**
 * \brief           update obstacles
 */
void cSimScene::UpdateObjects()
{

    // 1. update perturb on objects if possible
    if (mPerturb != nullptr)
    {
        // ApplyUserPerturbForceOnce
        mPerturb->mObject->ApplyUserPerturbForceOnce(mPerturb);
    }

    // 2. update objects

    for (auto &obs : this->mObjectList)
    {

        obs->Update(mCurdt);
    }
}
/**
 * \brief       do (discrete) collision detection
 */
void cSimScene::PerformCollisionDetection()
{
    if (mEnableCollisionDetection == true)
    {
        mColDetecter->Update();
        mColDetecter->PerformCD();
        // auto pts = mColDetecter->GetContactPoints();
        // std::cout << "[debug] num of contacts = " << pts.size() << std::endl;
        for (auto &obs : this->mObjectList)
        {
            obs->SetPointTriangleCollisionInfo(
                mColDetecter->GetObjPointTriangleCollisionInfo(
                    obs->GetObjId()));
            obs->SetEdgeEdgeCollisionInfo(
                mColDetecter->GetObjEdgeEdgeCollisionInfo(obs->GetObjId()));
        }
    }
}
/**
 * \brief           Reset the whole scene
 */
void cSimScene::Reset()
{
    cScene::Reset();
    ReleasePerturb();
    ClearForce();
    for (auto &x : mObjectList)
        x->Reset();
}

/**
 * \brief           Get number of vertices
 */
int cSimScene::GetNumOfVertices() const
{
    int num_of_vertices = 0;
    for (auto &x : mObjectList)
    {
        num_of_vertices += x->GetNumOfVertices();
    }
    return num_of_vertices;
}

/**
 * \brief       clear all forces
 */
void cSimScene::ClearForce() {}

int cSimScene::GetNumOfFreedom() const { return GetNumOfVertices() * 3; }

int cSimScene::GetNumOfDrawEdges() const
{
    int num_of_edges = 0;
    for (auto &x : mObjectList)
    {
        num_of_edges += x->GetNumOfDrawEdges();
    }
    return num_of_edges;
}

int cSimScene::GetNumOfDrawTriangles() const
{
    int num_of_triangles = 0;
    for (auto &x : mObjectList)
    {
        num_of_triangles += x->GetNumOfDrawTriangles();
    }
    return num_of_triangles;
}
/**
 * \brief       external force
 */
const tVectorXf &cSimScene::GetTriangleDrawBuffer()
{
    return mTriangleDrawBuffer;
}
/**
 * \brief           Calculate vertex rendering data
 */
void cSimScene::CalcTriangleDrawBuffer()
{
    mTriangleDrawBuffer.fill(std::nan(""));
    // 1. calculate for sim triangle
    int st = 0;
    Eigen::Map<tVectorXf> ref(mTriangleDrawBuffer.data(),
                              mTriangleDrawBuffer.size());

    // 2. calculate for obstacle triangle

    for (auto &x : mObjectList)
    {
        x->CalcTriangleDrawBuffer(ref, st);
    }
}

const tVectorXf &cSimScene::GetEdgesDrawBuffer() { return mEdgesDrawBuffer; }
bool cSimScene::IsRenderingResourceSizeUpdated() const
{
    return mIsRenderingResourceSizeUpdated;
}
void cSimScene::UpdateRenderingResource()
{
    mIsRenderingResourceSizeUpdated = false;
    {
        int num_of_triangles = GetNumOfDrawTriangles();
        int num_of_vertices = num_of_triangles * 3;
        int size_per_vertices = RENDERING_SIZE_PER_VERTICE;
        int tar_size = num_of_vertices * size_per_vertices;
        if (mTriangleDrawBuffer.size() != tar_size)
        {
            mIsRenderingResourceSizeUpdated = true;
            mTriangleDrawBuffer.resize(tar_size);
        }
        // std::cout << "triangle draw buffer size = " <<
        // mTriangleDrawBuffer.size() << std::endl; exit(0);
    }
    {

        int size_per_edge = 2 * RENDERING_SIZE_PER_VERTICE;
        int tar_size = GetNumOfDrawEdges() * size_per_edge;
        if (mEdgesDrawBuffer.size() != tar_size)
        {
            mIsRenderingResourceSizeUpdated = true;
            mEdgesDrawBuffer.resize(tar_size);
        }
    }
    {
        int num_of_v = GetNumOfVertices();
        int tar_size = num_of_v * RENDERING_SIZE_PER_VERTICE;
        if (tar_size != mPointDrawBuffer.size())
        {
            mIsRenderingResourceSizeUpdated = true;
            mPointDrawBuffer.resize(tar_size);
        }
    }
    CalcEdgesDrawBuffer();
    CalcTriangleDrawBuffer();
    CalcPointDrawBuffer();
}

int cSimScene::CalcEdgesDrawBuffer(int st /* = 0 */)
{
    SIM_ASSERT(st == 0);
    mEdgesDrawBuffer.fill(std::nan(""));

    Eigen::Map<tVectorXf> render_ref(mEdgesDrawBuffer.data() + st,
                                     mEdgesDrawBuffer.size() - st);

    // 2. for draw buffer

    for (auto &x : mObjectList)
    {
        x->CalcEdgeDrawBuffer(render_ref, st);
    }

    return st;
}

cSimScene::~cSimScene()
{
    // for (auto x : mVertexArray)
    //     delete x;
    // mVertexArray.clear();
    // for (auto &x : mTriangleArray)
    //     delete x;
    // mTriangleArray.clear();
    // for (auto &x : mEdgeArray)
    //     delete x;
    // mEdgeArray.clear();
}

/**
 * \brief               Event response (add perturb)
 */
void cSimScene::CursorMove(int xpos, int ypos) {}

void cSimScene::UpdatePerturbPos(const tVector &camera_pos, const tVector &dir)
{
    if (mPerturb != nullptr)
    {
        mPerturb->UpdatePerturbPos(camera_pos, dir);
    }
}
/**
 * \brief               Event response (add perturb)
 */

void cSimScene::MouseButton(int button, int action, int mods) {}

#include "GLFW/glfw3.h"
void cSimScene::Key(int key, int scancode, int action, int mods)
{
    // std::cout << "[sim scene] key = " << key << std::endl;
    mSimStateMachine->Key(key, scancode, action, mods);
    switch (key)
    {
    case GLFW_KEY_S:
        std::cout << "[draw scene] key S, save now\n";
        SaveCurrentScene();
        break;
    }
}

bool cSimScene::CreatePerturb(tRay *ray)
{
    SIM_ASSERT(mRaycaster != nullptr);

    cRaycaster::tRaycastResult res = mRaycaster->RayCast(ray);
    if (res.mObject == nullptr)
    {
        return false;
    }
    else
    {
        std::cout << "[debug] add perturb on triangle " << res.mLocalTriangleId
                  << std::endl;
    }

    // 2. we have a triangle to track
    SIM_ASSERT(mPerturb == nullptr);

    mPerturb = new tPerturb();

    mPerturb->mObject = res.mObject;
    mPerturb->mAffectedTriId = res.mLocalTriangleId;
    // std::cout << "[debug] affect id = " << res.mLocalTriangleId << std::endl;
    const auto &ver_array = mPerturb->mObject->GetVertexArray();
    const auto &tri_array = mPerturb->mObject->GetTriangleArray();

    mPerturb->mBarycentricCoords =
        cMathUtil::CalcBarycentric(
            res.mIntersectionPoint,
            ver_array[tri_array[res.mLocalTriangleId]->mId0]->mPos,
            ver_array[tri_array[res.mLocalTriangleId]->mId1]->mPos,
            ver_array[tri_array[res.mLocalTriangleId]->mId2]->mPos)
            .segment(0, 3);
    // std::cout << "[perturb] intersection pt (from raycast) = "
    //           << res.mIntersectionPoint.transpose() << std::endl;
    // std::cout << "[perturb] bary = " <<
    // mPerturb->mBarycentricCoords.transpose()
    //           << std::endl;

    tVector restore_intersection_pt =
        ver_array[tri_array[res.mLocalTriangleId]->mId0]->mPos *
            mPerturb->mBarycentricCoords[0] +
        ver_array[tri_array[res.mLocalTriangleId]->mId1]->mPos *
            mPerturb->mBarycentricCoords[1] +
        ver_array[tri_array[res.mLocalTriangleId]->mId2]->mPos *
            mPerturb->mBarycentricCoords[2];
    // std::cout << "[perturb] restore_intersection_pt = "
    //           << restore_intersection_pt.transpose() << std::endl;

    // std::cout
    //     << "uv = "
    //     << ver_array[tri_array[res.mLocalTriangleId]->mId0]->muv.transpose()
    //     << " vid = " << tri_array[res.mLocalTriangleId]->mId0 << std::endl;
    SIM_ASSERT(mPerturb->mBarycentricCoords.hasNaN() == false);
    mPerturb->InitTangentRect(-1 * ray->mDir);
    mPerturb->UpdatePerturbPos(ray->mOrigin, ray->mDir);

    // // change the color
    mPerturb->mObject->ChangeTriangleColor(
        res.mLocalTriangleId, ColorShoJoHi.segment(0, 3).cast<float>());
    return true;
}
void cSimScene::ReleasePerturb()
{
    if (mPerturb != nullptr)
    {
        // restore the color

        mPerturb->mObject->ChangeTriangleColor(
            mPerturb->mAffectedTriId, ColorBlue.segment(0, 3).cast<float>());
        // 1, 0);
        delete mPerturb;
        mPerturb = nullptr;
    }
}

#include "sim/kinematic/kinematicBodyBuilder.h"

void cSimScene::CreateObstacle(const Json::Value &conf)
{

    // printf("[debug] create %d obstacle(s) done\n", mObstacleList.size());
    // exit(0);
}

// /**
//  * \brief                   Raycast the whole scene
//  * @param ray:              the given ray
//  * @param selected_tri:     a reference to selected triangle pointer
//  * @param selected_tri_id:  a reference to selected triangle id
//  * @param raycast_point:    a reference to intersection point
//  */
// void cSimScene::RayCastScene(const tRay *ray, cBaseObjectPtr casted_obj,
//                              int obj_triangle_id, tVector inter_point) const
// {
//     SIM_ASSERT(mRaycaster != nullptr);
//     cRaycaster::tRaycastResult res = mRaycaster->RayCast(ray);
//     casted_obj = res.mObject;
//     obj_triangle_id
// }

/**
 * \brief                   Collision Detection
 */
void cSimScene::CreateCollisionDetecter()
{
    if (mEnableCollisionDetection)
    {
        mColDetecter = std::make_shared<cBVHCollisionDetecter>();
        // add resources into the collision detecter now
        for (auto &x : this->mObjectList)
        {
            bool enable_self_collision =
                (x->GetObjectType() == eObjectType::CLOTH_TYPE);
            mColDetecter->AddObject(x, enable_self_collision);
        }
        mColDetecter->Init();
    }
}

/**
 * \brief                   Get number of objects
 */
int cSimScene::GetNumOfObjects() const
{
    int num_of_objects = 0;
    // num_of_objects += mObstacleList.size();
    return num_of_objects;
}

/**
 * \brief                   Get obstacle list
 */
std::vector<cKinematicBodyPtr> cSimScene::GetObstacleList()
{
    // return this->mObstacleList;
    return {};
}

/**
 * \brief                   Is sim paused
 */
bool cSimScene::IsSimPaused() const
{
    return mSimStateMachine->IsRunning() == false;
}

/**
 * \brief                   Update imgui for simulation scene
 */
void cSimScene::UpdateImGui()
{
    auto cur_state = mSimStateMachine->GetCurrentState();
    auto name = tSimStateMachine::BuildStateStr(cur_state);
    // std::cout << "cSimScene::UpdateImGui, cur state = " << cur_state << "
    // name = " << name << std::endl;
    ImGui::Text("simulation state: %s", name.c_str());
    for (auto &obj : this->mObjectList)
    {
        obj->UpdateImGui();
    }
    // update vertices and triangle number
    int v_total = 0, t_total = 0;
    for (auto &obj : this->mObjectList)
    {
        int num_of_v = obj->GetNumOfVertices();
        int num_of_t = obj->GetNumOfTriangles();
        v_total += num_of_v;
        t_total += num_of_t;
        ImGui::Text("%s v %d t %d", obj->GetObjName().c_str(), num_of_v,
                    num_of_t);
    }
    ImGui::Text("total v %d t %d", v_total, t_total);
}

const tVectorXf &cSimScene::GetPointDrawBuffer() { return mPointDrawBuffer; }

void cSimScene::CalcPointDrawBuffer()
{
    mPointDrawBuffer.fill(std::nan(""));

    Eigen::Map<tVectorXf> render_ref(mPointDrawBuffer.data(),
                                     mPointDrawBuffer.size());

    // 2. for draw buffer
    int st = 0;
    for (auto &x : mObjectList)
    {
        x->CalcPointDrawBuffer(render_ref, st);
    }
}

void cSimScene::SetCameraPos(const tVector3d &cam_pos) { mCamPos = cam_pos; }