#include "SoftBodyImplicit.h"
#include <iostream>
#include "utils/LogUtil.h"
#include "geometries/Tetrahedron.h"

cSoftBodyImplicit::cSoftBodyImplicit(int id_) : cSoftBody(id_)
{
}

cSoftBodyImplicit::~cSoftBodyImplicit()
{
}

void cSoftBodyImplicit::Update(float dt)
{
}
void cSoftBodyImplicit::Init(const Json::Value &conf)
{
    InitDDsDxTensor();
    cSoftBody::Init(conf);

    CheckDFDx();
}

tEigenArr<tMatrix3d> cSoftBodyImplicit::CalcDFDx(int ele_id)
{
    tEigenArr<tMatrix3d> dfdx = this->mDDsDxTensor_const;
    for (int i = 0; i < dfdx.size(); i++)
    {
        dfdx[i] = (dfdx[i] * mInvDm[ele_id]).eval();
    }
    return dfdx;
}

void cSoftBodyImplicit::CalcDPDF()
{
}

void cSoftBodyImplicit::UpdateIntForce()
{
    // SIM_ERROR("for implicit scheme, we do not calculate int force explicitly.");
    // exit(1);
    CalculateStiffnessMatrix();
}

void cSoftBodyImplicit::SolveForNextPos(float dt)
{
}

/**
 * \brief           calculate sitffness matrix K
*/
void cSoftBodyImplicit::CalculateStiffnessMatrix()
{
}

void cSoftBodyImplicit::CheckDFDx()
{
    UpdateDeformationGradient();
    auto F_old_lst = this->mF;

    for (int ele_id = 0; ele_id < GetNumOfTets(); ele_id++)
    {
        // calculate dfdx analytically
        tEigenArr<tMatrix3d> dfdx_ana_lst = CalcDFDx(ele_id);
        tMatrix3d F_old = mF[ele_id];

        // calculate dfdx numerically
        auto cur_tet = mTetArrayShared[ele_id];
        tVector4i v_id = cur_tet->mVertexId;
        float eps = 1e-5;
        // calculate for a point
        for (int j = 0; j < 4; j++)
        {
            int cur_v = v_id[j];
            for (int comp = 0; comp < 3; comp += 1)
            {
                // calculate new deformation gradient
                mXcur[3 * cur_v + comp] += eps;
                SyncPosToVectorArray();
                UpdateDeformationGradientForTet(ele_id);

                tMatrix3d F_new = mF[ele_id];
                tMatrix3d dFdx_num = (F_new - F_old) / eps;
                tMatrix3d dFdx_ana = dfdx_ana_lst[3 * j + comp];
                tMatrix3d dFdx_diff = dFdx_num - dFdx_ana;
                float dFdx_diff_norm = dFdx_diff.norm();
                if (dFdx_diff_norm > 10 * eps)
                {
                    std::cout << "[error] failed, dFdx diff norm = " << dFdx_diff_norm << std::endl;
                    std::cout << "dFdx num = \n" << dFdx_num << std::endl;
                    std::cout << "dFdx ana = \n" << dFdx_ana << std::endl;
                    std::cout << "dFdx diff = \n" << dFdx_diff << std::endl;
                }

                mXcur[3 * cur_v + comp] -= eps;
                SyncPosToVectorArray();
                UpdateDeformationGradientForTet(ele_id);
            }
        }
    }

    SIM_INFO("check df/dx succ");
}
void cSoftBodyImplicit::CheckDPDF()
{
}

void cSoftBodyImplicit::InitDDsDxTensor()
{
    mDDsDxTensor_const.resize(12, tMatrix3d::Zero());
    for (int i = 0; i < 12; i++)
    {
        // for the first 3 points
        if (i < 9)
        {
            int row_id = i % 3;
            int col_id = std::floor(i / 3);
            mDDsDxTensor_const[i](row_id, col_id) = 1;
        }
        else
        {
            // for the last point
            mDDsDxTensor_const[i].row(i - 9).noalias() = tVector3d::Ones() * -1;
        }
        printf("d(Ds)/dx[%d] = ", i);
        std::cout << mDDsDxTensor_const[i] << std::endl;
    }
}