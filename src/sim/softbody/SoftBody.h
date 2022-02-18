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
enum eMaterialModelType
{
    LINEAR_ELASTICITY = 0,
    COROTATED,
    FIX_COROTATED,
    STVK,
    NEO_HOOKEAN,
    NUM_OF_MATERIAL_MODEL
};

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
    virtual void SetVerticesPos(const tVectorXd &pos);
    virtual void ApplyUserPerturbForceOnce(tPerturb *) override;
    virtual eMaterialModelType GetMaterial() const;
    virtual void UpdateImGUi() override;
    virtual void Reset() override;

protected:
    std::string mTetMeshPath = "";
    // std::vector<tVertexPtr> mVertexArrayShared = {};
    // std::vector<tEdgePtr> mEdgeArrayShared = {};
    // std::vector<tTrianglePtr> mTriangleArrayShared = {};
    std::vector<tTetPtr> mTetArrayShared = {};
    tEigenArr<tMatrix3d> mF;                                   // deformation gradient, F
    tEigenArr<tMatrix3d> mInvDm;                               // DmInv, used to calculate F, please read the SIGGRAPH 2012 course for more details
    tVectorXd mGravityForce, mExtForce, mIntForce, mUserForce; // internal & external force vector on each node, \in R^{3n}
    tVectorXd mXcur, mXprev, mXInit;                           // timestep previous and current
    tVectorXd mInvLumpedMassMatrixDiag;                        // row-diagnozation-lumped mass matrix
    tVectorXd mInitTetVolume;
    double mRho; // the volume density [SI] kg/m^3
    eMaterialModelType mMaterial;
    float mRayleighDamplingA, mRayleighDamplingB; // rayleigh damping for mass mat and stiffness mat
    float mFrictionCoef, mCollisionK;
    virtual void InitInvDm();
    virtual void InitPos();
    virtual void InitDiagLumpedMassMatrix();
    virtual void InitTetVolume();
    virtual void InitForceVector();
    virtual void UpdateIntForce();
    virtual tVectorXd CalcTetIntForce(size_t tet_id);
    virtual tVectorXd CalcTetIntForceBySelectionMatrix(size_t tet_id);
    virtual tVectorXd GetTetForce(size_t tet_id, const tVectorXd &total_force);
    virtual tVectorXd GetTetVerticesPos(size_t tet_id, const tVectorXd &total_pos);
    virtual void UpdateExtForce();
    double CalcEnergy();
    virtual void UpdateTriangleNormal() override;
    virtual void UpdateVertexNormalFromTriangleNormal() override;
    virtual void UpdateDeformationGradient();
    virtual void UpdateDeformationGradientForTet(int ele);
    virtual void SolveForNextPos(float dt);
    void SyncPosToVectorArray();
    int GetNumOfTets() const;
    virtual int GetNumOfFreedoms() const;
};

SIM_DECLARE_PTR(cSoftBody);

extern std::string BuildMaterialTypeStr(eMaterialModelType type);
tMatrix3d CalcGreenStrain(const tMatrix3d &F);
tMatrix3d CalcPK1(const tMatrix3d &F);
tMatrix3d CalcPK1_part1(const tMatrix3d &F);
tMatrix3d CalcPK1_part2(const tMatrix3d &F);
void GetSelectionMatrix(tMatrixXd &Sd, tMatrixXd &Sb);
extern float gMu, gLambda;