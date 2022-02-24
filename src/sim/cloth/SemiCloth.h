#pragma once
#include "BaseCloth.h"

class cSemiCloth : public cBaseCloth
{
public:
    inline static const std::string ENABLEQBENDING_KEY = "enable_Q_bending",
                                    BENDINGSTIFFNESS_KEY = "bending_stiffness",
                                    STIFFNESS_KEY = "stiffness";
    cSemiCloth(int id_);
    virtual ~cSemiCloth();
    virtual void Init(const Json::Value &conf);
    // virtual void Update(double dt) override final;
    virtual void UpdatePos(double dt) override;

protected:
    virtual void InitGeometry(const Json::Value &conf);
    double mStiffness;
    double mEnableQBending;
    double mBendingStiffness;
    tVectorXd mXpre;
    tSparseMat mBendingHessianQ;
    void InitBendingHessian();
    tVectorXd CalcNextPositionSemiImplicit() const;
};