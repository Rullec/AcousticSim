#pragma once
#include "utils/MathUtil.h"
#include "utils/DefUtil.h"
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
class cBaseObject : std::enable_shared_from_this<cBaseObject>
{
public:
    inline static const std::string OBJECT_NAME_KEY = "object_name";
    explicit cBaseObject(eObjectType type, int obj_id);
    virtual ~cBaseObject();
    virtual void SetObjName(std::string);
    virtual std::string GetObjName() const;
    virtual void Init(const Json::Value &conf);
    static eObjectType BuildObjectType(std::string type);
    eObjectType GetObjectType() const;
    virtual void CalcTriangleDrawBuffer(Eigen::Map<tVectorXf> &res,
                                        int &st) const = 0;
    virtual void CalcEdgeDrawBuffer(Eigen::Map<tVectorXf> &res,
                                    int &st) const = 0;
    // virtual void Update(double dt) = 0;

    // triangularize methods to visit the mesh data
    virtual int GetNumOfTriangles() const;
    virtual int GetNumOfEdges() const;
    virtual int GetNumOfVertices() const;
    void SetVertexColorAlpha(float val);
    float GetVertexColorAlpha() const;

    const std::vector<tVertex *> &GetVertexArray() const;
    const std::vector<tEdge *> &GetEdgeArray() const;
    const std::vector<tTriangle *> &GetTriangleArray() const;

    std::vector<tVertex *> &GetVertexArrayRef();
    std::vector<tEdge *> &GetEdgeArrayRef();
    std::vector<tTriangle *> &GetTriangleArrayRef();

    void ChangeTriangleColor(int tri_id, const tVector3f &color);
    virtual void CalcAABB(tVector &min, tVector &max) const;
    double CalcTotalArea() const;
protected:
    float mColorAlpha = 1.0;
    int mObjId;
    std::string mObjName;
    eObjectType mType;
    bool mEnableDrawBuffer; // enable to open draw buffer
    std::vector<tVertex *> mVertexArray;
    std::vector<tEdge *> mEdgeArray;
    std::vector<tTriangle *> mTriangleArray;
    virtual void UpdateTriangleNormal();
    virtual void UpdateVertexNormalFromTriangleNormal();
};

SIM_DECLARE_PTR(cBaseObject);