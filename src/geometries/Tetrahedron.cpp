#include "Tetrahedron.h"

tTet::tTet()
{
    this->mTriangleId.fill(-1);
    for (int i = 0; i < 4; i++)
    {
        mTriangleOpposite[i] = false;
    }
}