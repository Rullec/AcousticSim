#pragma once
#include "utils/MathUtil.h"
class cQBendingMaterial : std::enable_shared_from_this<cQBendingMaterial>
{
public:
    static double CalcEnergy(
        const tVector3d &v0,
        const tVector3d &v1,
        const tVector3d &v2,
        const tVector3d &v3, double K);
    static tVectorXd CalcForce(const tVector3d &v0,
                               const tVector3d &v1,
                               const tVector3d &v2,
                               const tVector3d &v3, double K);
    static tMatrixXd CalcStiffnessMatrix(const tVector3d &v0,
                                         const tVector3d &v1,
                                         const tVector3d &v2,
                                         const tVector3d &v3, double K);
    // static void CheckForce();
    // static void CheckStiffnessMatrix();

protected:
    static tMatrixXd CalcKCoef(const tVector3d &v0,
                               const tVector3d &v1,
                               const tVector3d &v2,
                               const tVector3d &v3);
};
