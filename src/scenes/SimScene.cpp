#include <iostream>
#include "SimScene.h"
#include "utils/ColorUtil.h"
#include "utils/JsonUtil.h"
#include "geometries/CollisionDetecter.h"
#include "geometries/Primitives.h"
#include "geometries/Triangulator.h"
#include "sim/KinematicBody.h"
#include "sim/Perturb.h"
#include "sim/SimObjectBuilder.h"
// std::string gSceneTypeStr[eSceneType::NUM_OF_SCENE_TYPES] = {"sim", "acoustic"};
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
    mPauseSim = false;
    mSceneType = eSceneType::SCENE_SIM;
    // mColDetecter = nullptr;
}

eSceneType cSimScene::GetSceneType() const { return this->mSceneType; }
#include "sim/AudioOutput.h"
extern cAudioOutputPtr gAudioOutput;

void cSimScene::Init(const std::string &conf_path)
{
    // 1. load config
    Json::Value root;
    cJsonUtil::LoadJson(conf_path, root);

    mEnableProfiling =
        cJsonUtil::ParseAsBool(cSimScene::ENABLE_PROFLINE_KEY, root);

    mEnableCollisionDetection =
        cJsonUtil::ParseAsBool(cSimScene::ENABLE_COLLISION_DETECTION_KEY, root);
    // gAudioOutput->Init();
    BuildObjects(
        cJsonUtil::ParseAsValue(cSimScene::OBJECT_LIST_KEY, root));

    CreateCollisionDetecter();

    InitDrawBuffer();
    InitRaycaster(root);
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
void cSimScene::PauseSim() { mPauseSim = !mPauseSim; }

/**
 * \breif       save current scene (obstacles to objes)
 */
#include "geometries/ObjExport.h"
void cSimScene::SaveCurrentScene()
{
    // for (auto &x : mObstacleList)
    // {
    //     cObjExporter::ExportObj(x->GetObjName() + ".obj", x->GetVertexArray(),
    //                             x->GetTriangleArray());
    // }
}

void cSimScene::InitDrawBuffer()
{
    // 2. build arrays
    // init the buffer
    {
        int num_of_triangles = GetNumOfTriangles();
        int num_of_vertices = num_of_triangles * 3;
        int size_per_vertices = RENDERING_SIZE_PER_VERTICE;
        int trinalge_size = num_of_vertices * size_per_vertices;

        mTriangleDrawBuffer.resize(trinalge_size);
        // std::cout << "triangle draw buffer size = " <<
        // mTriangleDrawBuffer.size() << std::endl; exit(0);
    }
    {

        int size_per_edge = 2 * RENDERING_SIZE_PER_VERTICE;
        mEdgesDrawBuffer.resize(GetNumOfDrawEdges() * size_per_edge);
    }

    UpdateRenderingResource();
}

/**
 * \brief           Init the raycasting strucutre
 */
#include "geometries/Raycaster.h"
void cSimScene::InitRaycaster(const Json::Value &conf)
{
    // auto total_triangle_array = mTriangleArray;
    // auto total_vertex_array = mVertexArray;
    // std::cout << "begin to add obstacle data array\n";
    // for (auto &x : mObstacleList)
    // {
    //     auto obstacle_v_array =  x->GetVertexArray();
    //     auto obstacle_triangle_array =  x->GetVertexArray();

    // }
    // for (int i = 0; i < this->)
    mRaycaster = std::make_shared<cRaycaster>();
    mRaycaster->Init(conf);
    // for (auto &x : mObstacleList)
    // {
    //     // auto obstacle_v_array = x->GetVertexArray();
    //     // auto obstacle_triangle_array = x->GetTriangleArray();
    //     mRaycaster->AddResources(x);
    // }
    // std::cout << "[debug] add resources to raycaster done, num of obstacles =
    // "
    //           << mObstacleList.size() << std::endl;
}
/**
 * \brief           Update the simulation procedure
 */
#include "utils/TimeUtil.hpp"
void cSimScene::Update(double delta_time)
{
    // double default_dt = mIdealDefaultTimestep;
    // if (delta_time < default_dt)
    //     default_dt = delta_time;
    // printf("[debug] sim scene update cur time = %.4f\n", mCurTime);
    cScene::Update(delta_time);

    UpdateObstacles();
    // clear force
    // apply ext force
    // update position
}

/**
 * \brief           update obstacles
 */
void cSimScene::UpdateObstacles()
{
    for (auto &obs : this->mObjectList)
    {
        
        // if (false == obs->IsStatic())
        {
            // std::cout << " obstacle " << obs->GetObjName()
            //           << " is not static, need to update\n";
            obs->Update(mCurdt);
        }
    }
}
/**
 * \brief       do (discrete) collision detection
 */
void cSimScene::PerformCollisionDetection()
{
    if (mEnableCollisionDetection == true)
    {
        mColDetecter->PerformCD();
        // auto pts = mColDetecter->GetContactPoints();
        // std::cout << "[debug] num of contacts = " << pts.size() << std::endl;
    }
}
/**
 * \brief           Reset the whole scene
 */
void cSimScene::Reset()
{
    cScene::Reset();
    ClearForce();
}

/**
 * \brief           Get number of vertices
 */
int cSimScene::GetNumOfVertices() const
{
    int num_of_vertices = 0;
    // for (auto &x : mObstacleList)
    // {
    //     num_of_vertices += x->GetNumOfVertices();
    // }
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
        num_of_edges += x->GetNumOfEdges();
    }
    return num_of_edges;
}

int cSimScene::GetNumOfTriangles() const
{
    int num_of_triangles = 0;
    for (auto &x : mObjectList)
    {
        num_of_triangles += x->GetNumOfTriangles();
    }
    return num_of_triangles;
}
/**
 * \brief       external force
 */
extern const tVector gGravity;

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

void cSimScene::UpdateRenderingResource()
{
    CalcEdgesDrawBuffer();
    CalcTriangleDrawBuffer();
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
#include "scenes/DrawScene.h"
void cSimScene::MouseButton(int button, int action, int mods)
{
}
#include "GLFW/glfw3.h"
void cSimScene::Key(int key, int scancode, int action, int mods)
{
    if (key == GLFW_KEY_I && action == GLFW_PRESS)
    {
        PauseSim();
    }
    else if (key == GLFW_KEY_S && action == GLFW_PRESS)
    {
        SaveCurrentScene();
    }
}

bool cSimScene::CreatePerturb(tRay *ray)
{
    return true;
}
void cSimScene::ReleasePerturb()
{

    if (mPerturb != nullptr)
    {
        // restore the color

        mPerturb->mObject->ChangeTriangleColor(
            mPerturb->mAffectedTriId,
            ColorBlue.segment(0, 3).cast<float>());
        // 1, 0);
        delete mPerturb;
        mPerturb = nullptr;
    }
}

#include "sim/KinematicBodyBuilder.h"

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

        mColDetecter = std::make_shared<cCollisionDetecter>();
        // add resources into the collision detecter now
        // for (auto &x : this->mObstacleList)
        // {
        //     mColDetecter->AddObject(x, false);
        // }
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
bool cSimScene::IsSimPaused() const { return this->mPauseSim; }