#pragma once
#include "Primitives.h"
#include "utils/DefUtil.h"

/**
 * \brief           data strucutre for tetrahedron
*/
struct tTet
{
    tTet();
    tVector4i mTriangleFaceId;
};

SIM_DECLARE_PTR(tTet);