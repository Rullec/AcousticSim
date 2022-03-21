#pragma once
#include "geometries/Primitives.h"
#include <memory>
/**
 * \brief           A BVH tree for a single object
 */
struct tBVHNode : std::enable_shared_from_this<tBVHNode>
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    tBVHNode();
    tVector mAABBMin, mAABBMax;
    int mObjId, mTriangleId;
    tBVHNode *mLeft, *mRight;
};
SIM_DECLARE_PTR(tBVHNode)
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
    std::vector<tVertexPtr> mVertexArray;
    std::vector<tEdgePtr> mEdgeArray;
    std::vector<tTrianglePtr> mTriangleArray;
    std::vector<tBVHNodePtr> mNodes;
};