#pragma once
#include "BaseCloth.h"

class cImplicitCloth : public cBaseCloth
{
public:
    inline static const std::string MAX_NEWTON_ITERS_KEY = "max_newton_iters",
                                    STIFFNESS_KEY = "stiffness";
    explicit cImplicitCloth(int id_);
    virtual ~cImplicitCloth();

    virtual void Init(const Json::Value &conf);
    virtual void UpdatePos(double dt) override;

protected:
    int mMaxNewtonIters; // max newton iterations in implicit integration
    double mStiffness;   // spring stiffness
    tVectorXd CalcNextPositionImplicit();
    void CalcGxImplicit(const tVectorXd &xcur, tVectorXd &Gx,
                        tVectorXd &fint_buf, tVectorXd &fext_buf,
                        tVectorXd &fdamp_buffer) const;
    virtual void InitGeometry(const Json::Value &conf) override final;
    void CalcdGxdxImplicit(const tVectorXd &xcur, tMatrixXd &Gx) const;
    void CalcdGxdxImplicitSparse(const tVectorXd &xcur, tSparseMat &Gx) const;
    void TestdGxdxImplicit(const tVectorXd &x0, const tMatrixXd &Gx_ana);
};