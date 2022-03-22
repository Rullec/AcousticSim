#pragma once
#include "geometries/AABB.h"
#include "geometries/Primitives.h"
#include <memory>

/**
 * \brief           A BVH tree for a single object
 */
SIM_DECLARE_STRUCT_AND_PTR(tBVHNode);
struct tBVHNode : std::enable_shared_from_this<tBVHNode>
{
    tBVHNode();
    int mId;
    bool mIsLeaf;
    tAABB mAABB;
    int mObjId, mTriangleId;
    tBVHNodePtr mLeft, mRight;
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
    virtual void RebuildTree();
    virtual void Print() const;
protected:
    int mObjId;
    std::vector<tVertexPtr> mVertexArray;
    std::vector<tEdgePtr> mEdgeArray;
    std::vector<tTrianglePtr> mTriangleArray;
    std::vector<tBVHNodePtr> mNodes;
    tBVHNodePtr CreateSubTree(const tAABB &node_ideal_AABB_used_for_split,
                          const std::vector<int> &vertices_array_in_this_node,
                          const std::vector<int> *local_vertex_id_sorted_xyz);
};