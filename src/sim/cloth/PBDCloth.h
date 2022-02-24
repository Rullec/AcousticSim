#pragma once
#include "sim/cloth/BaseCloth.h"

class cPBDCloth : public cBaseCloth
{
public:
    inline static const std::string PBD_CONFIG_KEY = "pbd_config",
                                    MAX_PBD_ITERS_KEY = "max_pbd_iters",
                                    STIFFNESS_PBD_KEY = "stiffness_pbd",
                                    ENABLE_PARALLEL_PBD_KEY =
                                        "enable_parallel_pbd",
                                    ENABLE_BENDING_PBD_KEY =
                                        "enable_bending_pbd",
                                    BENDING_STIFFNESS_PBD_KEY =
                                        "bending_stiffness_pbd";

    explicit cPBDCloth(int id_);
    virtual ~cPBDCloth();
    virtual void Init(const Json::Value &conf);
    virtual void UpdatePos(double dt) override;

protected:
    tVectorXd mVcur; // velocity vector
    int mItersPBD;
    double mStiffnessPBD;
    bool mEnableParallelPBD;                 // enable parallel pbd
    bool mEnableBendingPBD;                  // enable bending constraint
    double mBendingStiffnessPBD;             // bending constraint stiffness PBD
    tEigenArr<tVector> mBendingMatrixKArray; // the array of bending matrix "K"

    void UpdateSubstepPBD();
    void UpdateVelAndPosUnconstrained(const tVectorXd &fext);
    virtual void CalcExtForce(tVectorXd &ext_force) const override;
    // void ConstraintSetupPBD();
    void ConstraintProcessPBD();
    void PostProcessPBD();

    void StretchConstraintProcessPBD(double final_k);
    void BendingConstraintProcessPBD(double final_k);
    void InitBendingMatrixPBD();
};