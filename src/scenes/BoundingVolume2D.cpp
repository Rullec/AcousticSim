#include "BoundingVolume2D.h"
#include "utils/LogUtil.h"
#include <iostream>
bool tGrid::IsInsideAABB(const tVector2f &cur_pos) const
{
    return (cur_pos[0] > mAABBSt[0]) && (cur_pos[0] < mAABBEd[0]) &&
           (cur_pos[1] > mAABBSt[1]) && (cur_pos[1] < mAABBEd[1]);
}

tTriangle2D::tTriangle2D(int id_, const tVector2f &p0, const tVector2f &p1,
                         const tVector2f &p2)
{
    mTriangleId = id_;
    pos0 = p0;
    pos1 = p1;
    pos2 = p2;
}

cBoundVolume2D::cBoundVolume2D()
{
    mTriangleArrayShared.clear();
    mGridArray.clear();
    mAABBMin = tVector2f::Ones() * 1e9;
    mAABBMax = -tVector2f::Ones() * 1e9;
}

tGrid::tGrid(const tVector2f &AABBSt, const tVector2f &AABBEd)
{
    mTrianglesArray.clear();
    mAABBSt = AABBSt;
    mAABBEd = AABBEd;
}

void cBoundVolume2D::AddTriangle(const tVector2f &pos0, const tVector2f &pos1,
                                 const tVector2f &pos2)
{
    auto tri =
        std::make_shared<tTriangle2D>(mTriangleArrayShared.size(), pos0, pos1, pos2);
    mTriangleArrayShared.push_back(tri);
    UpdateAABB(pos0);
    UpdateAABB(pos1);
    UpdateAABB(pos2);
}

void tGrid::AddTriangle(int tri_id, tTriangle2DPtr cur_tri)
{
    mLocalTriId2GlobalTriId[mTrianglesArray.size()] = tri_id;
    mTrianglesArray.push_back(cur_tri);
}
void cBoundVolume2D::InitVolume()
{
    const int Gap = 5;
    // 1. create grid
    mGridArray.clear();
    double x_gap = (mAABBMax - mAABBMin)[0] / (Gap - 1);
    double y_gap = (mAABBMax - mAABBMin)[1] / (Gap - 1);
    for (int x_id = 0; x_id < Gap; x_id++)
    {
        for (int y_id = 0; y_id < Gap; y_id++)
        {
            mGridArray.push_back(std::make_shared<tGrid>(
                tVector2f(x_id * x_gap, y_id * y_gap),
                tVector2f((x_id + 1) * x_gap, (y_id + 1) * y_gap)));
        }
    }

    // 2. dispatch grid
    for (int tri_id = 0; tri_id < mTriangleArrayShared.size(); tri_id++)
    {
        auto t = mTriangleArrayShared[tri_id];
        for (auto &cur_grid : mGridArray)
        {
            if ((true == cur_grid->IsInsideAABB(t->pos0)) ||
                (true == cur_grid->IsInsideAABB(t->pos1)) ||
                (true == cur_grid->IsInsideAABB(t->pos2)))
            {
                cur_grid->AddTriangle(tri_id, t);
            }
        }
    }
}

void cBoundVolume2D::UpdateAABB(const tVector2f &pos)
{
    for (int i = 0; i < 2; i++)
    {
        if (pos[i] > mAABBMax[i])
        {
            mAABBMax[i] = pos[i];
        }
        if (pos[i] < mAABBMin[i])
        {
            mAABBMin[i] = pos[i];
        }
    }
}

int cBoundVolume2D::FindIncludedTriangle(const tVector2f &pos, tVector3f &bary)
{
    for (auto &grid : mGridArray)
    {
        if (grid->IsInsideAABB(pos))
        {
            int tri_id = grid->FindInsideTriangle(pos, bary);
            if (tri_id != -1)
            {
                return tri_id;
            }
        }
    }
    bary.noalias() = tVector3f::Ones() * std::nanf("");

    return -1;
}

/**
 * \brief           Given a vertex, find the nearest triangle
 */
int cBoundVolume2D::FindNearestTriangle(const tVector2f &pos, tVector3f &bary)
{
    // 1. find the nearest grid
    float min_dist = 1e9;
    int grid_id = -1;
    for (int i = 0; i < mGridArray.size(); i++)
    {
        auto grid = this->mGridArray[i];
        float cur_dist = grid->CalcDistanceFromPosToSquareGrid(pos);
        if (cur_dist < min_dist)
        {
            min_dist = cur_dist;
            grid_id = i;
        }
    }
    // std::cout << "[debug] find nearest grid: pos " << pos.transpose()
    //           << " grid_id = " << grid_id << " range "
    //           << mGridArray[grid_id]->mAABBSt.transpose() << " to "
    //           << mGridArray[grid_id]->mAABBEd.transpose()
    //           << ", dist = " << min_dist << std::endl;

    // 2. find the nearest triangle in this grid
    int tri_id = mGridArray[grid_id]->FindNearestTriangleWhenOutside(pos, bary);
    return tri_id;
}

int cBoundVolume2D::FindNearestTriangleAll(const tVector2f &pos,
                                           tVector3f &bary)
{
    SIM_ERROR("hasn't manage\n");
    // for (auto &grid : mGridArray)
    // {
    //     int tri_id = grid->FindNearestTriangleWhenOutside(pos, bary);
    // }
}

float sign(const tVector2f &p1, const tVector2f &p2, const tVector2f &p3)
{
    return (p1.x() - p3.x()) * (p2.y() - p3.y()) -
           (p2.x() - p3.x()) * (p1.y() - p3.y());
}

bool PointInTriangle(const tVector2f &pt, const tVector2f &v1,
                     const tVector2f &v2, const tVector2f &v3)
{
    float d1, d2, d3;
    bool has_neg, has_pos;

    d1 = sign(pt, v1, v2);
    d2 = sign(pt, v2, v3);
    d3 = sign(pt, v3, v1);

    has_neg = (d1 < 0) || (d2 < 0) || (d3 < 0);
    has_pos = (d1 > 0) || (d2 > 0) || (d3 > 0);

    return !(has_neg && has_pos);
}

/**
 * \brief
 */
int tGrid::FindInsideTriangle(const tVector2f &cur_pos, tVector3f &bary)
{
    for (auto &tri : mTrianglesArray)
    {
        if (true == PointInTriangle(cur_pos, tri->pos0, tri->pos1, tri->pos2))
        {
            bary = cMathUtil::CalcBarycentric(
                       tVector(cur_pos[0], cur_pos[1], 0, 1),
                       tVector(tri->pos0[0], tri->pos0[1], 0, 1),
                       tVector(tri->pos1[0], tri->pos1[1], 0, 1),
                       tVector(tri->pos2[0], tri->pos2[1], 0, 1))
                       .cast<float>()
                       .segment(0, 3);
            // std::cout << "get bary = " << bary.transpose() << std::endl;
            return tri->mTriangleId;
        }
    }
    bary.noalias() = std::nanf("") * tVector3f::Ones();
    return -1;
}

float CalcPointTriangleDist(const tVector2f &cur_pos, tTriangle2DPtr cur_tri)
{
    float tri_dist = 1e9;
    for (int i = 0; i < 3; i++)
    {
        tVector2f tri_st = (i == 0)
                               ? (cur_tri->pos0)
                               : ((i == 1) ? cur_tri->pos1 : (cur_tri->pos2)),
                  tri_ed = (i == 0)
                               ? (cur_tri->pos1)
                               : ((i == 1) ? cur_tri->pos2 : (cur_tri->pos0));
        tVector2f v1 = tri_ed - tri_st;
        tVector2f v2 = cur_pos - tri_st;
        float cos = v1.dot(v2) / (v1.norm() * v2.norm());
        if (v1.dot(v2) < 0)
        {
            float cur_dist = v2.norm();
            tri_dist = std::min(tri_dist, cur_dist);
        }
        else if (v2.norm() * cos > v1.norm())
        {
            float cur_dist = (cur_pos - tri_ed).norm();
            tri_dist = std::min(tri_dist, cur_dist);
        }
        else
        {
            float edge_dist =
                (v2 - (v1.normalized().dot(v2)) * v1.normalized()).norm();
            tri_dist = std::min(tri_dist, edge_dist);
        }
    }
    return tri_dist;
}
/**
 * \brief           when the cur pos is outside of any triangle, calculate the
 * pos
 */
int tGrid::FindNearestTriangleWhenOutside(const tVector2f &cur_pos,
                                          tVector3f &bary) const
{
    // 1. find the min tri id
    float min_tri_dist = 1e9;
    int min_tri_id = -1;
    for (int tri_id = 0; tri_id < this->mTrianglesArray.size(); tri_id++)
    {
        auto cur_tri = mTrianglesArray[tri_id];

        float tri_dist = CalcPointTriangleDist(cur_pos, cur_tri);
        if (tri_dist < min_tri_dist)
        {
            min_tri_dist = tri_dist;
            min_tri_id = tri_id;
        }
    }

    auto cur_tri = mTrianglesArray[min_tri_id];
    // std::cout << "[debug] for pos " << cur_pos.transpose()
    //           << " the min tri dist = " << min_tri_dist
    //           << " min tri id = " << min_tri_id
    //           << " \ntri pos = " << cur_tri->pos0.transpose() << " \n"
    //           << cur_tri->pos1.transpose() << " \n"
    //           << cur_tri->pos2.transpose() << "\n";

    // 2. begin to calculate the bary
    bary = cMathUtil::CalcBarycentric(
               tVector(cur_pos[0], cur_pos[1], 0, 1),
               tVector(cur_tri->pos0[0], cur_tri->pos0[1], 0, 1),
               tVector(cur_tri->pos1[0], cur_tri->pos1[1], 0, 1),
               tVector(cur_tri->pos2[0], cur_tri->pos2[1], 0, 1))
               .segment(0, 3)
               .cast<float>();
    // std::cout << "bary = " << bary.transpose() << std::endl;
    int g_id = mLocalTriId2GlobalTriId.find(min_tri_id)->second;
    return g_id;
}
#include "utils/LogUtil.h"
/**
 * \brief           calculate the external distance between the pos to the grid
 *      if the pos is inside the grid, return 0
 */
float tGrid::CalcDistanceFromPosToSquareGrid(const tVector2f &pos) const
{
    if (IsInsideAABB(pos) == true)
        return 0;
    else
    {
        //  x_start
        float min_dist = ((mAABBEd + mAABBSt) / 2 - pos).norm();
        return min_dist;
    }
}
void cBoundVolume2D::Shift(const tVector2f &shift)
{
    mAABBMax += shift;
    mAABBMin += shift;
    for (auto &g : mGridArray)
    {
        g->mAABBEd += shift;
        g->mAABBSt += shift;
    }
    for (auto &t : mTriangleArray)
    {
        t->pos0 += shift;
        t->pos1 += shift;
        t->pos2 += shift;
    }
}