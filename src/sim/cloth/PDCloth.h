#pragma once
#include "sim/cloth/BaseCloth.h"
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>

class cPDCloth : public cBaseCloth
{
public:
    inline static const std::string MAX_STEPS_OPT_KEY = "max_steps_opt",
                                    ENABLE_BENDING_KEY = "enable_bending",
                                    BENDING_STIFFNESS_KEY = "bending_stiffness",
                                    PD_CONTACT_FORCE_KP_KEY =
                                        "pd_contact_force_kp",
                                    PD_CONTACT_FORCE_KD_KEY =
                                        "pd_contact_force_kd",
                                    PD_CONTACT_FRICTION_KEY =
                                        "pd_contact_friction";
    explicit cPDCloth(int id_);
    virtual ~cPDCloth();
    virtual void Init(const Json::Value &conf) override final;
    virtual void UpdatePos(double dt) override final;

protected:
    int mMaxSteps_Opt; // iterations used in fast simulation
    bool mEnableBending;
    double mBendingStiffness; // bending stiffness

    tVectorXd CalcNextPosition() const;
    // void InitVarsOptImplicit();
    void InitVarsOptImplicitSparse();
    void InitVarsOptImplicitSparseFast();
    virtual void InitGeometry(const Json::Value &conf) override final;
    const tEigenArr<tTriplet> &GetStretchTriplet() const;
    void AddBendTriplet(tEigenArr<tTriplet> &) const;
    // tMatrixXd J, I_plus_dt2_Minv_L_inv; // vars used in fast simulation
    tSparseMat J_sparse,
        I_plus_dt2_Minv_L_sparse; // vars used in fast simulation
    Eigen::SparseLU<tSparseMat> I_plus_dt2_Minv_L_sparse_solver;
    Eigen::SparseLU<tSparseMat> I_plus_dt2_Minv_L_sparse_fast_solver;
    tSparseMat I_plus_dt2_Minv_L_sparse_fast; // vars used in fast

    double mPDContactForceKp, mPDContactForceKd, mPDContactForceFriction; // explicit normal / tangential contact force coeff
};