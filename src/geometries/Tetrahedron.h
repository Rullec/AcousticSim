#pragma once
#include "Primitives.h"
#include "utils/DefUtil.h"

/**
 * \brief           data strucutre for tetrahedron
*/
struct tTet
{
    tTet();
    /**
     * v0, v1, v2, v3 index
    */
    tVector4i mVertexId; // four vertices

    /*    
    each (positively-oriented) tetrahedral element with vertices
    has triangular faces
    v0 v1 v2
    v1 v3 v2
    v2 v3 v0
    v3 v1 v0

    Here we remember the triangle id and opposite
    */
    tVector4i mTriangleId;
    bool mTriangleOpposite;
};

SIM_DECLARE_PTR(tTet);