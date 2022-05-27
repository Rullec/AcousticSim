#pragma once
#include "geometries/Primitives.h"

class cInnerForceCalculator
{
public:
    static std::vector<std::pair<int, int>> FindInvolvedVertices(const std::vector<tVertexPtr> &v_array,
                                     const std::vector<tEdgePtr> &e_array,
   const std::vector<tTrianglePtr> &t_array,
   const tVector2d & uv_min,
   const tVector2d & uv_max,
                                     double tar_v);
};