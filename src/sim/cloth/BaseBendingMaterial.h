#pragma once
#include "utils/DefUtil.h"
#include "utils/EigenUtil.h"
#include "utils/MathUtil.h"
#include "utils/SparseUtil.h"
SIM_DECLARE_STRUCT_AND_PTR(tVertex);
SIM_DECLARE_STRUCT_AND_PTR(tEdge);
SIM_DECLARE_STRUCT_AND_PTR(tTriangle);
using tMatrix12f = Eigen::Matrix<float, 12, 12>;

class cBaseBendingMaterial : std::enable_shared_from_this<cBaseBendingMaterial>
{
public:
    cBaseBendingMaterial();
    virtual void Init(const std::vector<tVertexPtr> &v_array,
                      const std::vector<tEdgePtr> &e_array,
                      const std::vector<tTrianglePtr> &t_array,
                      const tVector3d &bending_stiffness_warpweftbias);
    virtual void Update() = 0;
    virtual double CalcEnergy(const tVectorXd &xcur) = 0;
    virtual tVectorXd CalcForce(const tVectorXd &xcur) = 0;
    virtual tSparseMatd GetStiffnessMatrix() const;
    virtual void CheckForce() = 0;
    virtual void CheckStiffnessMatrix() = 0;

    std::vector<tMatrix12f> GetEleStiffnessMatrixLst() const;
    std::vector<tVector4i> GetEdgeConstraintVertex() const;

protected:
    std::vector<tVertexPtr> mVertexArray;
    std::vector<tEdgePtr> mEdgeArray;
    std::vector<tTrianglePtr> mTriArray;
    double mEdgeK;
    tSparseMatd mStiffnessMat;
    std::vector<tMatrix12f> mEleKLst;
    std::vector<tVector4i> mEdgeConstraintVertexLst;
};
