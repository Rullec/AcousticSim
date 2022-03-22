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
    int mTriangleId;
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
    virtual void UpdateAABB();
    virtual void RebuildTree();
    virtual void Print() const;
    virtual int GetNumOfLeaves() const;
    virtual const std::vector<tBVHNodePtr> GetLeaves() const;
    virtual const tBVHNodePtr GetRootNode() const;

    virtual std::vector<int> Intersect(tBVHNodePtr outer_node) const;

protected:
    int mObjId;
    std::vector<tVertexPtr> mVertexArray;
    std::vector<tEdgePtr> mEdgeArray;
    std::vector<tTrianglePtr> mTriangleArray;
    std::vector<tBVHNodePtr> mNodes;
    std::vector<tBVHNodePtr> mLeafNodes;
    tBVHNodePtr
    CreateSubTree(const tAABB &node_ideal_AABB_used_for_split,
                  const std::vector<int> &vertices_array_in_this_node,
                  const std::vector<int> *local_vertex_id_sorted_xyz);
};