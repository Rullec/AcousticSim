#include "InnerForceCalculator.h"
#include <iostream>
#include <set>
std::vector<std::pair<int, int>> cInnerForceCalculator::FindInvolvedVertices(
    const std::vector<tVertexPtr> &v_array,
    const std::vector<tEdgePtr> &e_array,
    const std::vector<tTrianglePtr> &t_array, const tVector2d &aabb_min,
    const tVector2d &aabb_max, double tar_v)
{
    // 1. check v range (must inside the interval)

    std::cout << "aabb min = " << aabb_min.transpose() << std ::endl;
    std::cout << "aabb max = " << aabb_max.transpose() << std ::endl;
    double v_min = aabb_min[1];
    double v_max = aabb_max[1];
    if (tar_v >= v_max || tar_v <= v_min)
    {
        printf("[warn] cur tar v %.2f, but v in [%.2f, %.2f], return\n", tar_v,
               v_min, v_max);
        return {};
    }
    std::set<std::pair<int, int>> cross_info_tri_v_set = {};
    // 2. for each edge, get the v range, and get the result.
    for (int t_id = 0; t_id < t_array.size(); t_id++)
    {
        auto t = t_array[t_id];
        int e_id_array[3] = {t->mEId0, t->mEId1, t->mEId2};

        // 2.1 if judge this edge is crossed, get the vertices who have minimium
        // v
        for (int i = 0; i < 3; i++)
        {
            int cur_e = e_id_array[i];
            int vid_min = e_array[cur_e]->mId0, vid_max = e_array[cur_e]->mId1;
            double v_min = v_array[vid_min]->muv[1];
            double v_max = v_array[vid_max]->muv[1];

            if (v_min > v_max)
            {
                SIM_SWAP(v_min, v_max);
                SIM_SWAP(vid_min, vid_max);
            }

            if (tar_v >= v_min && tar_v <= v_max)
            {
                // selected!
                auto info = std::make_pair(t_id, vid_min);
                cross_info_tri_v_set.insert(info);
            }
        }
    }
    std::vector<std::pair<int, int>> cross_info_tri_v(
        cross_info_tri_v_set.begin(), cross_info_tri_v_set.end());
    // 3. which triangle, which vertices are involved?
    return cross_info_tri_v;
}