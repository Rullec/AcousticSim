#pragma once
#include "utils/MathUtil.h"
class cCGSolver
{
public:
    static tVectorXd Solve(const tMatrixXd &A, const tVectorXd &b);
};