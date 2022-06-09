#pragma once
#include "Primitives.h"
#include "utils/DefUtil.h"

/**
 * \brief           data strucutre for tetrahedron
*/
struct tTet
{
    tTet();
    tVector4i mVertexId; // four vertices
    tVector4i mTriangleId; // four triangles
};

SIM_DECLARE_PTR(tTet);