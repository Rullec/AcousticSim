#pragma once
#include "sim/BaseObject.h"
#include "utils/MathUtil.h"

class cAcousticBody : public cBaseObject
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    explicit cAcousticBody(int obj_id);
    virtual ~cAcousticBody();
};