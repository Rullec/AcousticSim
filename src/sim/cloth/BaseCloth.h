#pragma once
#include "sim/BaseObject.h"

enum eClothType
{
    SEMI_IMPLICIT_CLOTH = 0,
    IMPLICIT_CLOTH,
    PBD_CLOTH,
    PD_CLOTH,
    LINCTEX_CLOTH,
    EMPTY_CLOTH, // cannot be simulated
    FEM_CLOTH,
    NUM_OF_CLOTH_TYPE
};

struct tPerturb;
SIM_DECLARE_CLASS_AND_PTR(cCollisionDetecter);
class cBaseCloth : public cBaseObject
{
public:
    inline static const std::string DAMPING_KEY = "damping",
                                    DEFAULT_TIMESTEP_KEY = "default_timestep";
    explicit cBaseCloth(eClothType cloth_type, int id_);
    virtual ~cBaseCloth();
    virtual void Init(const Json::Value &conf);
    virtual void Reset();
    virtual void SetCollisionDetecter(cCollisionDetecterPtr);
    virtual void CalcTriangleDrawBuffer(Eigen::Map<tVectorXf> &res,
                                        int &st) const;
    virtual void CalcEdgeDrawBuffer(Eigen::Map<tVectorXf> &res, int &st) const;
    static eClothType BuildClothType(std::string str);
    // virtual void Update(double dt);
    virtual void ClearForce();
    virtual void ApplyPerturb(tPerturb *pert);
    virtual void Update(float dt) override;
    virtual void UpdatePos(double dt) = 0;
    virtual void ApplyUserPerturbForceOnce(tPerturb *) override;
    virtual void SetPos(const tVectorXd &newpos);
    virtual const tVectorXd &GetPos() const;
    virtual double GetDefaultTimestep() { return mIdealDefaultTimestep; }
    virtual tVector2d GetClothShape() const;
    virtual void ClearConstraintStaticVertices();
    virtual void AddConstraintStaticVertices(const std::vector<int> &vertices);
    virtual const std::vector<int> &GetConstraintStaticVertices() const;
    virtual const tVectorXd &GetInitPos() const;
    virtual void MoveTranslation(const tVector3d &incremental_move);
    tVector CalcCOM() const;
    
protected:
    eClothType mClothType;
    cCollisionDetecterPtr mColDetecter;
    double mIdealDefaultTimestep; // default substep dt
    tVector2d mClothSizes;
    double mClothMass;            // SI key
    float mClothDensity;          // SI kg/m^2
    tVectorXd mInvMassMatrixDiag; // diag inv mass matrix
    std::string mGeometryType;
    double mDamping;                             // damping coeff
    tVectorXd mIntForce, mGravityForce, mUserForce, mCollisionForce, mDampingForce;
    
    bool mEnableClothFromObj;                    // the cloth mesh comes from obj or not
    std::string mClothObjPath;                   // if the cloth mesh comes from obj, the path
    tVectorXd mXpre, mXcur;                      // previous node position & current node position
    std::vector<int> mConstraint_StaticPointIds; // fixed constraint point
    tVectorXd mClothInitPos;                     // init position of the cloth
    virtual void InitGeometry(
        const Json::Value &conf); // discretazation from square cloth to
    virtual void InitMass(const Json::Value &conf);
    virtual void InitConstraint(const Json::Value &root);
    virtual void CalcIntForce(const tVectorXd &xcur,
                              tVectorXd &int_force) const;
    virtual void CalcExtForce(tVectorXd &ext_force) const;
    virtual void CalcDampingForce(const tVectorXd &vel,
                                  tVectorXd &damping) const;
    int GetNumOfFreedom() const;
    void CalcNodePositionVector(tVectorXd &pos) const;
};