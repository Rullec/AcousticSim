#pragma once
#include "geometries/Primitives.h"
#include <memory>
/**
 * \brief           A BVH tree for a single object
 */
struct tBVHNode
{
    tVector mAABBMin, mAABBMax;
};
class cObjBVH : std::enable_shared_from_this<cObjBVH>
{
public:
    explicit cObjBVH();
    virtual void Init(int obj_id, const std::vector<tVertexPtr> &v_array,
                      const std::vector<tEdgePtr> &e_array,
                      const std::vector<tTrianglePtr> &t_array);
    virtual void Update();

protected:
    int mObjId;
    std::vector<tVertexPtr> v_array;
    std::vector<tEdgePtr> e_array;
    std::vector<tTrianglePtr> t_array;
};