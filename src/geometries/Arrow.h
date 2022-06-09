#pragma once
#include "utils/DefUtil.h"
#include "utils/EigenUtil.h"
#include <memory>
#include <vector>

SIM_DECLARE_STRUCT_AND_PTR(tVertex);
SIM_DECLARE_STRUCT_AND_PTR(tEdge);
SIM_DECLARE_STRUCT_AND_PTR(tTriangle);
// advanced primitive: arrow
class cArrow : public std::enable_shared_from_this<cArrow>
{
public:
    explicit cArrow();
    virtual void SetStEd(const tVector &st, const tVector &ed);
    virtual void GetAABB(tVector &aabb_min, tVector &aabb_max);
    virtual int GetNumOfDrawEdges() const;
    virtual int GetNumOfDrawTriangles() const;
    virtual void CalcTriangleDrawBuffer(Eigen::Map<tVectorXf> &res,
                                        int &st) const;
    virtual void CalcEdgeDrawBuffer(Eigen::Map<tVectorXf> &res, int &st) const;

protected:
    std::vector<tVertexPtr> mVertexArray;
    std::vector<tEdgePtr> mEdgeArray;
    std::vector<tTrianglePtr> mTriangleArray;
    std::vector<tVector> mVertexInitPos;
    double mInitLength;   // meter
    double mInitDiameter; // meter
    tVector mSt, mEd;
};

SIM_DECLARE_PTR(tVertex)