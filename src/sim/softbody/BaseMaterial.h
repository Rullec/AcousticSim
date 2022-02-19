#pragma once
#include "utils/MathUtil.h"
#include "sim/softbody/ThreeOrderTensor.h"
#include "sim/softbody/FourOrderTensor.h"
namespace Json
{
    class Value;
};
class cBaseMaterial
{
public:
    explicit cBaseMaterial();
    virtual void Init(const Json::Value &conf) = 0;
    virtual tMatrix3d CalcP(const tMatrix3d &F) const = 0; // PK1
    virtual cFourOrderTensor CalcDPDF(const tMatrix3d &F) const = 0;
    virtual void CheckDPDF(const tMatrix3d &F) const;

protected:
};