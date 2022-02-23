#pragma once
#include "sim/softbody/SoftBodyImplicit.h"
#include "sim/AcousticBody.h"

class cAcousticSoftBody : public cAcousticBody, public cSoftBodyImplicit
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    explicit cAcousticSoftBody(int id_);
    virtual void Init(const Json::Value &conf) override;
    virtual void Update(float dt) override;
};