#pragma once
#include "Scene.h"
#include "utils/DefUtil.h"
#include "utils/MathUtil.h"

namespace Json
{
    class Value;
};

enum eSceneType
{
    SCENE_SIM = 0,  // default simulation scene
    SCENE_ACOUSTIC, // acoustic simulation scene
    NUM_OF_SCENE_TYPES
};

struct tVertex;
struct tEdge;
struct tTriangle;
struct tRay;
struct tPerturb;
// class cDrawScene;
SIM_DECLARE_CLASS_AND_PTR(cKinematicBody)
SIM_DECLARE_CLASS_AND_PTR(cBaseObject)
SIM_DECLARE_CLASS_AND_PTR(cRaycaster)
SIM_DECLARE_CLASS_AND_PTR(cCollisionDetecter)
class cSimScene : public cScene
{
public:
    inline static const std::string ENABLE_PROFLINE_KEY = "enable_profiling",
                                    ENABLE_OBSTACLE_KEY = "enable_obstacle",
                                    OBSTACLE_CONF_KEY = "obstacle_conf",
                                    ENABLE_COLLISION_DETECTION_KEY =
                                        "enable_collision_detection";

    cSimScene();
    ~cSimScene();
    virtual void Init(const std::string &conf_path) override;
    virtual void Update(double dt) override;
    virtual void UpdateRenderingResource();
    virtual void Reset() override;
    virtual const tVectorXf &GetTriangleDrawBuffer();
    virtual const tVectorXf &GetEdgesDrawBuffer();
    static eSceneType BuildSceneType(const std::string &str);
    eSceneType GetSceneType() const;
    virtual bool CreatePerturb(tRay *ray);
    virtual void ReleasePerturb();
    virtual void UpdatePerturbPos(const tVector &camera_pos,
                                  const tVector &dir);
    virtual void CursorMove(int xpos, int ypos);
    virtual void MouseButton(int button, int action, int mods);
    virtual void Key(int key, int scancode, int action, int mods);
    void RayCastScene(const tRay *ray, tTriangle **selected_triangle,
                      int &selected_triangle_id,
                      tVector &ray_cast_position) const;
    virtual int GetNumOfObjects() const;
    virtual std::vector<cKinematicBodyPtr> GetObstacleList();
    virtual bool IsSimPaused() const;

protected:
    tPerturb *mPerturb;
    eSceneType mSceneType;
    bool mEnableProfiling;
    bool mEnableObstacle; // using obstacle?
    bool mEnableCollisionDetection;
    std::vector<cBaseObjectPtr> mObjectList;
    cRaycasterPtr mRaycaster; // raycaster

    tVectorXf mTriangleDrawBuffer,
        mEdgesDrawBuffer; // buffer to triangle buffer drawing (should use index
                          // buffer to improve the velocity)
    std::vector<tRay *> mRayArray;
    // std::vector<tVertex *> mVertexArray;     // total vertices
    // std::vector<tEdge *> mEdgeArray;         // total edges
    // std::vector<tTriangle *> mTriangleArray; // total triangles

    cCollisionDetecterPtr mColDetecter; // collision detecter

    // base methods
    void CalcDampingForce(const tVectorXd &vel, tVectorXd &damping) const;
    virtual void InitDrawBuffer();
    virtual void InitRaycaster(const Json::Value &conf);

    void ClearForce(); // clear all forces
    void SaveCurrentScene();

    virtual void BuildObjects(const Json::Value &obj_conf_path) const;
    virtual void CalcTriangleDrawBuffer();       //
    virtual int CalcEdgesDrawBuffer(int st = 0); //
    virtual int GetNumOfVertices() const;
    virtual int GetNumOfFreedom() const;
    virtual int GetNumOfDrawEdges() const;
    virtual int GetNumOfTriangles() const;
    void CalcNodePositionVector(tVectorXd &pos) const;
    virtual void UpdateObstacles();
    virtual void CreateObstacle(const Json::Value &conf);
    virtual void CreateCollisionDetecter();
    bool mPauseSim;
    virtual void PauseSim();
    virtual void PerformCollisionDetection();
};