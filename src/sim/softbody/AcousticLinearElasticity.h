#pragma once
#include "geometries/Primitives.h"
#include "geometries/Tetrahedron.h"
#include "utils/DefUtil.h"
#include "utils/SparseUtil.h"

SIM_DECLARE_CLASS_AND_PTR(tVertex);
SIM_DECLARE_CLASS_AND_PTR(tTet);
class cAcousticLinearElasticity
{
public:
    static tSparseMatd
    CalcGlobalStiffness(const std::vector<tVertexPtr> &v_array,
                        const std::vector<tTetPtr> &tet_array,
                        double youngs_modulus, double poisson_ratio,
                        double damping_a, double damping_b);
};