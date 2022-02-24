#pragma once
#include "utils/MathUtil.h"
#include "utils/DefUtil.h"

SIM_DECLARE_STRUCT_AND_PTR(tVertex);

class cRenderUtil
{
public:
    static void CalcTriangleDrawBufferSingle(tVertexPtr v0, tVertexPtr v1, tVertexPtr v2,
                                      Eigen::Map<tVectorXf> &buffer,
                                      int &st_pos);

    static void CalcEdgeDrawBufferSingle(tVertexPtr v0, tVertexPtr v1,
                                  const tVector &edge_normal,
                                  Eigen::Map<tVectorXf> &buffer, int &st_pos,
                                  const tVector &color);

    static void CalcEdgeDrawBufferSingle(const tVector &v0, const tVector &v1,
                                  const tVector &edge_normal,
                                  Eigen::Map<tVectorXf> &buffer, int &st_pos,
                                  const tVector &color);
};