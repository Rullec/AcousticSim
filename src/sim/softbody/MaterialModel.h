#pragma once
#include "utils/MathUtil.h"
#include "sim/softbody/ThreeOrderTensor.h"
#include "sim/softbody/FourOrderTensor.h"
class cMaterialModel
{
public:
    cMaterialModel();
    virtual tMatrix3d CalcP(const tMatrix3d &F); // PK1
    virtual cFourOrderTensor CalcDPDF(const tMatrix3d &F);
    
};