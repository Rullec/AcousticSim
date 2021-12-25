#pragma once
#include "sim/BaseObject.h"

/**
 * \brief           softbody object (mostly based on 3D FEM)
*/

struct tTriangle;
struct tEdge;
struct tVertex;
struct tTet;
SIM_DECLARE_PTR(tTriangle);
SIM_DECLARE_PTR(tEdge);
SIM_DECLARE_PTR(tVertex);
SIM_DECLARE_PTR(tTet);
class cSoftBody : public cBaseObject
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    explicit cSoftBody(int id);
    virtual ~cSoftBody();

    virtual void Init(const Json::Value &conf) override;
    virtual void CalcTriangleDrawBuffer(Eigen::Map<tVectorXf> &res,
                                        int &st) const override;
    virtual void CalcEdgeDrawBuffer(Eigen::Map<tVectorXf> &res,
                                    int &st) const override;
    virtual int GetNumOfTriangles() const override;
    virtual int GetNumOfEdges() const override;
    virtual int GetNumOfVertices() const override;

protected:
    std::string mTetMeshPath = "";
    std::vector<tVertexPtr> mVertexArrayShared = {};
    std::vector<tEdgePtr> mEdgeArrayShared = {};
    std::vector<tTrianglePtr> mTriangleArrayShared = {};
    std::vector<tTetPtr> mTetArrayShared = {};
    virtual void UpdateTriangleNormal() override;

    virtual void UpdateVertexNormalFromTriangleNormal() override;
};

SIM_DECLARE_PTR(cSoftBody);