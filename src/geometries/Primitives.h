#pragma once
#include "utils/MathUtil.h"
struct tVertex
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    tVertex();
    double mMass;
    tVector mPos;
    tVector mNormal;
    tVector2f muv; // "texture" coordinate 2d, it means the plane coordinate for
                   // a vertex over a fabric, but now the texture in rendering
    tVector mColor;
};

struct tEdge
{
    tEdge();
    int mId0, mId1;
    double mRawLength; // raw length of this edge
    bool mIsBoundary;  // does this edge locate in the boundary?
    int mTriangleId0,
        mTriangleId1; // The indices of the two triangles to which this side
                      // belongs. If this edge is a boundary, the mTriangleId1
                      // is -1
    double mK_spring; // stiffness for springs
};

// struct tEdge : public tEdge
// {
//     EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
//     tEdge();
//     double mK;
// };

struct tTriangle
{
    explicit tTriangle();
    explicit tTriangle(int a, int b, int c);
    int mId0, mId1, mId2;
    tVector mNormal;
};

/**
 * \brief       an origin + a directed ray
 */
struct tRay
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    explicit tRay(const tVector &ori, const tVector &end);
    tVector mOrigin;
    tVector mDir;
};

struct tRectangle
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    tRectangle();
    tVector mVertex[4];
};