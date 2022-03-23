#pragma once
#include "utils/DefUtil.h"
#include "utils/MathUtil.h"

SIM_DECLARE_STRUCT_AND_PTR(tVertex);

class cRenderUtil
{
public:
    static void CalcTriangleDrawBufferSingle(tVertexPtr v0, tVertexPtr v1,
                                             tVertexPtr v2,
                                             const tVector &color,
                                             Eigen::Map<tVectorXf> &buffer,
                                             int &st_pos);

    static void CalcEdgeDrawBufferSingle(tVertexPtr v0, tVertexPtr v1,
                                         const tVector &edge_normal,
                                         Eigen::Map<tVectorXf> &buffer,
                                         int &st_pos, const tVector &color);

    static void CalcEdgeDrawBufferSingle(const tVector &v0, const tVector &v1,
                                         const tVector &edge_normal,
                                         Eigen::Map<tVectorXf> &buffer,
                                         int &st_pos, const tVector &color);
    static void CalcPointDrawBufferSingle(const tVector &pos,
                                          const tVector &color,
                                          Eigen::Map<tVectorXf> &buffer,
                                          int &st_pos);
};