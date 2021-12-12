#pragma once
#include "geometries/Primitives.h"
#include "utils/DefUtil.h"
#include <vector>
#include <map>

struct tTriangle2D : std::enable_shared_from_this<tTriangle2D>
{
    tTriangle2D(int mTriangleId_, const tVector2f &p0, const tVector2f &p1,
                const tVector2f &p2);
    int mTriangleId = -1;
    tVector2f pos0;
    tVector2f pos1;
    tVector2f pos2;
};
SIM_DECLARE_STRUCT_AND_PTR(tTriangle2D);

struct tGrid
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    tGrid(const tVector2f &mAABBSt, const tVector2f &mAABBEd);
    bool IsInsideAABB(const tVector2f &cur_pos) const;
    int FindInsideTriangle(const tVector2f &cur_pos, tVector3f &bary);
    int FindNearestTriangleWhenOutside(const tVector2f &cur_pos,
                                       tVector3f &bary) const;
    void AddTriangle(int tri_id, tTriangle2DPtr cur_tri);
    float CalcDistanceFromPosToSquareGrid(const tVector2f &pos) const;
    std::vector<tTriangle2DPtr> mTrianglesArray;
    tVector2f mAABBSt, mAABBEd;
    std::map<int, int> mLocalTriId2GlobalTriId;
};
SIM_DECLARE_PTR(tGrid);

class cBoundVolume2D
{
public:
    cBoundVolume2D();
    void AddTriangle(const tVector2f &pos0, const tVector2f &pos1,
                     const tVector2f &pos2);
    void InitVolume();
    int FindIncludedTriangle(const tVector2f &pos, tVector3f &bary);
    int FindNearestTriangle(const tVector2f &pos, tVector3f &bary);
    int FindNearestTriangleAll(const tVector2f &pos, tVector3f &bary);
    void Shift(const tVector2f &shift);
protected:
    std::vector<tTriangle2DPtr> mTriangleArray;
    std::vector<tGridPtr> mGridArray;

    tVector2f mAABBMin, mAABBMax;
    void UpdateAABB(const tVector2f &pos);
};