#pragma once
#include "utils/MathUtil.h"
#include "utils/DefUtil.h"
#include "utils/SparseUtil.h"
SIM_DECLARE_CLASS_AND_PTR(tVertex);
SIM_DECLARE_CLASS_AND_PTR(tEdge);
SIM_DECLARE_CLASS_AND_PTR(tTriangle);

class cQBendingMaterial : std::enable_shared_from_this<cQBendingMaterial>
{
public:
    explicit cQBendingMaterial();
    virtual void Init(
        const std::vector<tVertexPtr> &v_array,
        const std::vector<tEdgePtr> &e_array,
        const std::vector<tTrianglePtr> &t_array, const tVector3d &bending_stiffness_warpweftbias);

    double CalcEnergy(const tVectorXd & xcur);
    tVectorXd CalcForce(const tVectorXd & xcur);
    tSparseMat GetStiffnessMatrix() const;

    // static void CheckForce();
    // static void CheckStiffnessMatrix();

protected:
    tSparseMat mStiffnessMat;
    static tSparseMat CalcKCoef(const tVector3d &v0,
                               const tVector3d &v1,
                               const tVector3d &v2,
                               const tVector3d &v3);
};
