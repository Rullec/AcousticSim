#pragma once
#include "utils/DefUtil.h"
#include "utils/MathUtil.h"
#include <memory>
#include <string>

/**
 * \brief           The ULTIMATE object type collections
 */
enum eObjectType
{
    KINEMATICBODY_TYPE,
    RIGIDBODY_TYPE,
    CLOTH_TYPE,
    FLUID_TYPE,
    SOFTBODY_TYPE,
    ACOUSTIC_TYPE,
    NUM_OBJ_TYPES,
    INVALID_OBJ_TYPE
};

/**
 * \brief           base object class
 *
 */
namespace Json
{
class Value;
};

struct tVertex;
struct tEdge;
struct tTriangle;
struct tPerturb;
SIM_DECLARE_PTR(tVertex);
SIM_DECLARE_PTR(tEdge);
SIM_DECLARE_PTR(tTriangle);
SIM_DECLARE_STRUCT_AND_PTR(tPointTriangleCollisionInfo);
class cBaseObject : public std::enable_shared_from_this<cBaseObject>
{
public:
    inline static const std::string OBJECT_NAME_KEY = "object_name";
    explicit cBaseObject(eObjectType type, int obj_id);
    virtual ~cBaseObject();
    virtual int GetObjId() const;
    virtual void SetObjName(std::string);
    virtual std::string GetObjName() const;
    virtual void Init(const Json::Value &conf);
    static eObjectType BuildObjectType(std::string type);
    eObjectType GetObjectType() const;
    virtual void CalcTriangleDrawBuffer(Eigen::Map<tVectorXf> &res,
                                        int &st) const;
    virtual void CalcEdgeDrawBuffer(Eigen::Map<tVectorXf> &res,
                                    int &st) const = 0;
    virtual void CalcPointDrawBuffer(Eigen::Map<tVectorXf> &res, int &st) const;
    virtual void Update(float dt) = 0;
    virtual void ApplyUserPerturbForceOnce(tPerturb *) = 0;
    virtual void SetGravity(const tVector3d &g);
    // triangularize methods to visit the mesh data
    virtual int GetNumOfTriangles() const;
    virtual int GetNumOfEdges() const;
    virtual int GetNumOfVertices() const;
    void SetVertexColorAlpha(float val);
    float GetVertexColorAlpha() const;

    const std::vector<tVertexPtr> &GetVertexArray() const;
    const std::vector<tEdgePtr> &GetEdgeArray() const;
    const std::vector<tTrianglePtr> &GetTriangleArray() const;

    std::vector<tVertexPtr> &GetVertexArrayRef();
    std::vector<tEdgePtr> &GetEdgeArrayRef();
    std::vector<tTrianglePtr> &GetTriangleArrayRef();
    void SetPointTriangleCollisionInfo(
        const std::vector<tPointTriangleCollisionInfoPtr> &info);
    void ChangeTriangleColor(int tri_id, const tVector3f &color);
    virtual void CalcAABB(tVector &min, tVector &max) const;
    double CalcTotalArea() const;
    virtual void UpdateImGui();
    virtual void Reset() = 0;

protected:
    float mColorAlpha = 1.0;
    int mObjId;
    std::string mObjName;
    eObjectType mType;
    tVector3d mGravity;
    bool mEnableDrawBuffer; // enable to open draw buffer
    // std::vector<tVertexPtr > mVertexArray;
    // std::vector<tEdge *> mEdgeArray;
    // std::vector<tTriangle *> mTriangleArray;
    std::vector<std::vector<int>> mVertexConnectedTriangles;
    std::vector<tPointTriangleCollisionInfoPtr> mPointTriangleCollisionInfo;
    std::vector<tVertexPtr> mVertexArray;
    std::vector<tEdgePtr> mEdgeArray;
    std::vector<tTrianglePtr> mTriangleArray;
    virtual void UpdateTriangleNormal();
    virtual void UpdateVertexNormalFromTriangleNormal();
};

SIM_DECLARE_PTR(cBaseObject);