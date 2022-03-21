#include "geometries/ObjectBVH.h"

cObjBVH::cObjBVH()
{
    mObjId = -1;
    mVertexArray.clear();
    mEdgeArray.clear();
    mTriangleArray.clear();
    mNodes.clear();
}
void cObjBVH::Init(int obj_id, const std::vector<tVertexPtr> &v_array,
                   const std::vector<tEdgePtr> &e_array,
                   const std::vector<tTrianglePtr> &t_array)
{
    // 1. create 

}
virtual void Update();