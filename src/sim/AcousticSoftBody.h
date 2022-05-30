#pragma once
#include "sim/AcousticBody.h"
#include "sim/softbody/SoftBodyImplicit.h"

class cAcousticSoftBody : public cAcousticBody, public cSoftBodyImplicit
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    explicit cAcousticSoftBody(int id_);
    virtual void Init(const Json::Value &conf) override;
    virtual void Update(float dt) override;

protected:
    int mAcousticSamplingFreq;
    float mAcousticDuration;
    void SolveForMonopole();
    std::string GetWaveName() const;
};