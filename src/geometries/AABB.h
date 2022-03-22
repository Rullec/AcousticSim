#pragma once
#include "utils/DefUtil.h"
#include "utils/MathUtil.h"

SIM_DECLARE_STRUCT_AND_PTR(tVertex);

struct tAABB
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    explicit tAABB();
    tAABB(const tAABB &old_AABB);
    int GetMaxExtent() const;
    tVector GetExtent() const;
    tVector GetMiddle() const;
    void Reset();
    void Expand(const tVector3d &);
    void Expand(const tVector &);
    void Expand(const tVertexPtr &);
    void Expand(const tAABB &);
    bool IsInvalid() const;
    bool Intersect(const tAABB &other_AABB);
    tVector mMin, mMax;
};