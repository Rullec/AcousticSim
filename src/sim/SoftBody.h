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

    virtual void Update(float dt) override;
    virtual void Init(const Json::Value &conf) override;
    virtual void CalcTriangleDrawBuffer(Eigen::Map<tVectorXf> &res,
                                        int &st) const override;
    virtual void CalcEdgeDrawBuffer(Eigen::Map<tVectorXf> &res,
                                    int &st) const override;
    virtual int GetNumOfTriangles() const override;
    virtual int GetNumOfEdges() const override;
    virtual int GetNumOfVertices() const override;
    virtual void SetVerticesPos(const tVectorXf & pos);
protected:
    std::string mTetMeshPath = "";
    std::vector<tVertexPtr> mVertexArrayShared = {};
    std::vector<tEdgePtr> mEdgeArrayShared = {};
    std::vector<tTrianglePtr> mTriangleArrayShared = {};
    std::vector<tTetPtr> mTetArrayShared = {};
    tEigenArr<tMatrix3f> mF;        // deformation gradient, F
    tEigenArr<tMatrix3f> mInvDm;    // DmInv, used to calculate F, please read the SIGGRAPH 2012 course for more details
    tVectorXf mGravityForce, mExtForce,  mIntForce, mUserForce; // internal & external force vector on each node, \in R^{3n}
    tVectorXf mXcur, mXprev;        // timestep previous and current
    tVectorXf mInvMassMatrixDiag;
    
    virtual void InitInvDm();
    virtual void InitPos();

    void UpdateIntForce();
    void UpdateExtForce();
    float CalcEnergy();
    virtual void UpdateTriangleNormal() override;
    virtual void UpdateVertexNormalFromTriangleNormal() override;
    virtual void UpdateDeformationGradient();
    void SolveForNextPos(float dt);
    void SyncPosToVectorArray();
};

SIM_DECLARE_PTR(cSoftBody);