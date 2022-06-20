#pragma once
#include "sim/BaseObject.h"
class cViscoMassSpring : public cBaseObject
{
public:
    explicit cViscoMassSpring(int id_);
    virtual ~cViscoMassSpring();
    virtual void Init(const Json::Value &conf) override;
    virtual void Reset();
    virtual void Update(float dt) override;
    virtual void UpdateImGui() override;
    virtual void ApplyUserPerturbForceOnce(tPerturb *) override;
protected:
    float mCurX, mCurV; // pos and vel
    float mMass;
    float mK, mEta; // spring stiffness and viscosity coef

    struct
    {
        float mCurX, mCurV; // pos and vel
    } mInitProp;
    float mDt;
    std::vector<tVector> mInitPos;
    virtual void UpdateMesh();
};