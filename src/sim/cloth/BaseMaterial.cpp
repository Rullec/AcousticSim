#include "BaseMaterial.h"
#include "geometries/Primitives.h"
#include "utils/LogUtil.h"
#include <set>

cBaseBendingMaterial::cBaseBendingMaterial()
{
    mStiffnessMat.resize(0, 0);
    mEleKLst.clear();
    mEdgeConstraintVertexLst.clear();
}

tSparseMatd cBaseBendingMaterial::GetStiffnessMatrix() const
{
    return -1 * this->mStiffnessMat;
}

std::vector<tMatrix12f> cBaseBendingMaterial::GetEleStiffnessMatrixLst() const
{
    return mEleKLst;
}

std::vector<tVector4i> cBaseBendingMaterial::GetEdgeConstraintVertex() const
{
    return mEdgeConstraintVertexLst;
}

int SelectAnotherVertex(tTrianglePtr tri, int v0, int v1)
{
    SIM_ASSERT(tri != nullptr);
    std::set<int> vid_set = {tri->mId0, tri->mId1, tri->mId2};
    // printf("[debug] select another vertex in triangle 3 vertices (%d, %d, %d)
    // besides %d %d\n", tri->mId0, tri->mId1, tri->mId2, v0, v1);
    vid_set.erase(vid_set.find(v0));
    vid_set.erase(vid_set.find(v1));
    return *vid_set.begin();
};

void cBaseBendingMaterial::Init(const std::vector<tVertexPtr> &v_array,
                         const std::vector<tEdgePtr> &e_array,
                         const std::vector<tTrianglePtr> &t_array,
                         const tVector3d &bending_stiffness_warpweftbias)
{
    mVertexArray = v_array;
    mEdgeArray = e_array;
    mTriArray = t_array;
    mEdgeK = bending_stiffness_warpweftbias.mean();
    mEleKLst.resize(e_array.size());
    mEdgeConstraintVertexLst.clear();
    mStiffnessMat.resize(3 * v_array.size(), 3 * v_array.size());
    for (auto &e : e_array)
    {
        if (e->mIsBoundary == true)
        {
            mEdgeConstraintVertexLst.push_back(tVector4i(-1, -1, -1, -1));
        }
        else
        {
            int v0 = e->mId0, v1 = e->mId1;
            int v2 = SelectAnotherVertex(t_array[e->mTriangleId0], v0, v1),
                v3 = SelectAnotherVertex(t_array[e->mTriangleId1], v0, v1);
            int v_lst[4] = {v0, v1, v2, v3};
            mEdgeConstraintVertexLst.push_back(tVector4i(v0, v1, v2, v3));
        }
    }
}