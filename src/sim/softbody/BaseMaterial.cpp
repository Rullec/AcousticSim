#include "BaseMaterial.h"
#include "utils/LogUtil.h"

cBaseMaterial::cBaseMaterial()
{
}

void cBaseMaterial::CheckDPDF(const tMatrix3d &F) const
{
    SIM_ERROR("hasn't been impled");
    exit(1);
}
