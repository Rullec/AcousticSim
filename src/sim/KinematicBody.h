#pragma once
#include "sim/BaseObject.h"
#include "utils/MathUtil.h"

struct tTriangle;
struct tEdge;
struct tVertex;
enum eKinematicBodyShape
{
    KINEMATIC_PLANE = 0,
    KINEMATIC_CUBE,
    KINEMATIC_SPHERE,
    KINEMATIC_CAPSULE,
    KINEMATIC_CUSTOM,
    NUM_OF_KINEMATIC_SHAPE,
    KINEMATIC_INVALID
};

class cKinematicBody : public cBaseObject
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    inline const static std::string
        TYPE_KEY = "type",
        MESH_PATH_KEY = "mesh_path", TARGET_AABB_KEY = "target_aabb",
        SCALE_KEY = "scale", TRANSLATION_KEY = "translation",
        ORIENTATION_KEY = "orientation", PLANE_EQUATION_KEY = "equation",
        PLANE_SCALE_KEY = "plane_scale",
        // the below setting only work for moving object
        IS_STATIC_KEY = "is_static",
        TARGET_TRANSLATION_KEY = "target_translation",
        TARGET_ORIENTATION_KEY = "target_orientation",
        ELASPED_TIME_SEC_KEY = "elasped_time_sec";
    cKinematicBody(int id_);
    virtual ~cKinematicBody();
    virtual void Init(const Json::Value &conf) override;
    static eKinematicBodyShape BuildKinematicBodyShape(std::string type_str);
    bool IsStatic() const;
    eKinematicBodyShape GetBodyShape() const;
    
    virtual void Update(float dt) override;
    virtual void ApplyUserPerturbForceOnce(tPerturb *) override;
    // virtual void UpdatePos(double dt) override final;
    // virtual void UpdateRenderingResource() override final;
    virtual tMatrix GetCurWorldTransform() const; //
    void Reset() override;
    virtual tVector CalcCOM() const;
    virtual void MoveTranslation(const tVector &shift);
    virtual void ApplyScale(float scale);

protected:
    float mCurTime;
    eKinematicBodyShape mBodyShape;
    std::string mCustomMeshPath;
    tVectorXd
        mScaledMeshVertices; // the loaded obj mesh pos after AABB rescaled
    tVector mTargetAABBDontUseDirectly;     // scale setting
    tVector mScaleDontUseDirectly;          // scale setting explictly
    tVector mInitPos;                       // init position for kinematic body
    tVector mInitOrientation;               // init orietnation for kinect body,
    tVector mTargetPos, mTargetOrientation; // target pos and orientation
    bool mIsStatic;
    tVector mPlaneEquation;
    double mPlaneScale;
    double mMovingElaspedTimeSec; // how many seconds does it cost for moving
                                  // object?
    tMatrix mCurWorldTransform;   //
    // methods
    void BuildCustomKinematicBody();
    void SetMeshPos();
    void BuildPlane();
    void UpdateCurWorldTransformByTime();
    tVector GetScaleVec() const;
    // virtual void InitDrawBuffer() override final;
};