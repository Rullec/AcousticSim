#include <iostream>
#include "FEMCloth.h"
#include "utils/TimeUtil.hpp"
#include "geometries/Primitives.h"


cFEMCloth::cFEMCloth(int id_) : cBaseCloth(eClothType::FEM_CLOTH, id_)
{
    mF.clear();
    mJ.resize(0);
    mPK1.clear();
    mdFdx.clear();
}

cFEMCloth::~cFEMCloth() {}

void cFEMCloth::Init(const Json::Value &conf)
{
    cBaseCloth::Init(conf);
    mXcur.noalias() = mClothInitPos;
    mXpre.noalias() = mClothInitPos;
    InitBuffer();
}

/**
 * \brief       Update the nodal position
 */
void cFEMCloth::UpdatePos(double dt)
{
    // 1. calculate the deformation gradient
    CalculateF();
    // 1. calculate the internal force
    CalcIntForce(mXcur, mIntForce);
    // 2. update the velocity
    // 3. update the nodal position
}

/**
 * \brief       Init FEM buffer,
 */
void cFEMCloth::InitBuffer()
{
    int element_size = this->GetSingleElementFreedom();
    int num_of_triangles = GetNumOfTriangles();
    mF.resize(num_of_triangles, tMatrix32d::Zero());
    mJ.noalias() = tVectorXd::Zero(num_of_triangles);
    mPK1.resize(num_of_triangles, tMatrixXd::Zero(element_size, element_size));
    mdFdx.resize(
        num_of_triangles,
        tEigenArr<tMatrixXd>(element_size,
                             tMatrixXd::Zero(element_size, element_size)));

    InitMaterialCoords();
}

/**
 * \brief           Init matrix coords
 */
void cFEMCloth::InitMaterialCoords()
{
    // calculate material coordinates
    mVertexMateralCoords.noalias() = tMatrixXd::Zero(GetNumOfVertices(), 2);
    for (int i = 0; i < mVertexArrayShared.size(); i++)
    {
        mVertexMateralCoords.row(i).noalias() =
            mVertexArrayShared[i]->muv.cast<double>();
    }

    // calculate the DInv (used in material point)
    mDInv.resize(mTriangleArrayShared.size(), tMatrix2d::Zero());

    tMatrix2d mat1;
    for (int i = 0; i < this->mTriangleArrayShared.size(); i++)
    {
        auto tri = mTriangleArrayShared[i];
        mat1.col(0).noalias() = (mVertexMateralCoords.row(tri->mId1) -
                                 mVertexMateralCoords.row(tri->mId0))
                                    .transpose();
        mat1.col(1).noalias() = (mVertexMateralCoords.row(tri->mId2) -
                                 mVertexMateralCoords.row(tri->mId0))
                                    .transpose();

        /*
            mat0 = [b-a; c-a;] \in R 3 \times 2
            mat1 = [B-A; C-A;] \in R 2 \times 2
            mat0 = F * mat1, F \in R 3 \times 2
        */
        mDInv[i] = mat1.inverse();
    }
}
/**
 * \brief       the freedom of a triangle, 3 nodes in cartesian space, the dof =
 * 3*3 = 9
 */
int cFEMCloth::GetSingleElementFreedom() const { return 9; }

/**
 * \brief           calculate internal force in
 */
void cFEMCloth::CalcIntForce(const tVectorXd &xcur, tVectorXd &int_force) const
{
}

/**
 * \brief           calculate deformation gradient
 */
void cFEMCloth::CalculateF()
{
    // cTimeUtil::Begin("CalculateF");
    tMatrix32d mat0 = tMatrix32d::Zero();
    tMatrix2d mat1 = tMatrix2d::Zero();
    for (int i = 0; i < this->mTriangleArrayShared.size(); i++)
    {
        auto tri = mTriangleArrayShared[i];
        const tVector3d &a = mXcur.segment(3 * tri->mId0, 3);
        const tVector3d &b = mXcur.segment(3 * tri->mId1, 3);
        const tVector3d &c = mXcur.segment(3 * tri->mId2, 3);

        mat0.col(0).noalias() = b - a;
        mat0.col(1).noalias() = c - a;

        /*
            mat0 = [b-a; c-a;] \in R 3 \times 2
            mat1 = [B-A; C-A;] \in R 2 \times 2
            mat0 = F * mat1, F \in R 3 \times 2
        */
        mF[i] = mat0 * mDInv[i];
    }
    // cTimeUtil::End("CalculateF");
}