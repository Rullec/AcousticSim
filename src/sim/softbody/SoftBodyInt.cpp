#include "sim/softbody/SoftBody.h"
#include "geometries/Tetrahedron.h"
#include "utils/LogUtil.h"

void cSoftBody::UpdateIntForce()
{
    // update internal force
    int num_of_tet = this->mTetArrayShared.size();
    // internal force H = - W P(F) D_m^{-T}, iter on each tet

    mIntForce.setZero();
    tMatrix3f P = tMatrix3f::Zero();
    for (int i = 0; i < num_of_tet; i++)
    {
        auto tet = mTetArrayShared[i];
        // 1.1 get W: tet volume
        float W = mInitTetVolume[i];
        // 1.2 get deformation gradient F
        const tMatrix3f &F = mF[i];
        // 1.3 get P(F) in linear elasticity

        P.noalias() = CalcPiolaKirchoff(F);
        // 1.4 calculate force on nodes
        tMatrix3f H = -W * P * mInvDm[i].transpose();
        // std::cout << "tet " << i << " H = \n"
        //           << H << std::endl;
        tVector3f f3 = -(H.col(0) + H.col(1) + H.col(2));
        // std::cout << "f0 = " << H.col(0).transpose() << std::endl;
        // std::cout << "f1 = " << H.col(1).transpose() << std::endl;
        // std::cout << "f2 = " << H.col(2).transpose() << std::endl;
        // std::cout << "f3 = " << f3.transpose() << std::endl;
        mIntForce.segment(tet->mVertexId[0] * 3, 3) += H.col(0);
        mIntForce.segment(tet->mVertexId[1] * 3, 3) += H.col(1);
        mIntForce.segment(tet->mVertexId[2] * 3, 3) += H.col(2);
        mIntForce.segment(tet->mVertexId[3] * 3, 3) += f3;
    }
    // std::cout << "fint = " << mIntForce.transpose() << std::endl;
}

/**
 * \brief       Given deformation gradient F, calculate P(F)
*/
tMatrix3f cSoftBody::CalcPiolaKirchoff(const tMatrix3f &F)
{
    tMatrix3f P = tMatrix3f::Zero();
    tMatrix3f I = tMatrix3f::Identity();
    float mu = 1e5;
    float lambda = 0.5;

    switch (mMaterial)
    {
    // case eMaterialModelType::COROTATED:
    // {
    //     P.setZero();
    //     break;
    // }
    case eMaterialModelType::LINEAR_ELASTICITY:
    {
        /*
            P(F) = \mu * (FT + F -2I) + \lambda * tr(F - I) I
        */
        P.noalias() = mu * (F.transpose() + F - 2 * I) + lambda * (F - I).trace() * I;
        break;
    }
    // case eMaterialModelType::FIX_COROTATED:
    // {
    //     P.setZero();
    //     break;
    // }
    case eMaterialModelType::STVK:
    {
        tMatrix3f E = 0.5 * (F.transpose() * F - I);
        P = F * (2 * mu * E + lambda * E.trace() * I);
        break;
    }
    case eMaterialModelType::NEO_HOOKEAN:
    {
        /*
        P(F) = mu * F - mu * F^{-T} + lambda * log(I3) / 2 * F^{-T}
        */
        float I3 = 100;
        tMatrix3f F_inv_T = F.inverse().transpose();
        P.noalias() = mu * (F - F_inv_T) + lambda * std::log(I3) / 2 * F_inv_T;
        break;
    }
    default:
    {
        SIM_ERROR("do not support material model {}", BuildMaterialTypeStr(mMaterial));
    }
    break;
    }
    return P;
}