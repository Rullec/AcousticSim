#include "TetUtil.h"
#include "geometries/Primitives.h"
#include "geometries/Tetrahedron.h"
#include "utils/FileUtil.h"
#include "utils/LogUtil.h"
#include "utils/StringUtil.h"
#include <iostream>
#include <sstream>

/**
 * \brief       Calculate the volume for given tet
 */
float cTetUtil::CalculateTetVolume(const tVector &pos0, const tVector &pos1,
                                   const tVector &pos2, const tVector &pos3)
{
    // 1/6 * (AB X AC) \cdot (AD)
    tVector AB = pos1 - pos0;
    tVector AC = pos2 - pos0;
    tVector AD = pos3 - pos0;
    return (AB.cross3(AC)).dot(AD) / 6.0;
}

/**
 * \brief               load tet mesh from ".node" files (custom format)
 * \param path          mesh file path
 * \param vertex_vec    all vertices
 * \param edge_vec      all edges
 * \param tri_vec       all triangles vector
 * \param tet_vec       all tets vector
 */

tVector3i GetSortedTriplet(int a, int b, int c)
{
    int tmp;
    if (b < a)
    {
        tmp = a;
        a = b;
        b = tmp;
        // swap ab
    }
    if (c < b)
    {
        tmp = c;
        c = b;
        b = tmp;
        // swap bc
    }
    if (b < a)
    {
        tmp = b;
        b = a;
        a = tmp;
    }
    return tVector3i(a, b, c);
}
tVector2i GetSortedPair(int a, int b)
{
    int tmp;
    if (b < a)
    {
        tmp = a;
        a = b;
        b = tmp;
        // swap ab
    }
    return tVector2i(a, b);
}
int FindTriangle(tTrianglePtrVector tri_vec, int v0, int v1, int v2)
{
    tVector3i target = GetSortedTriplet(v0, v1, v2);
    for (int i = 0; i < tri_vec.size(); i++)
    {

        tVector3i cur = GetSortedTriplet(tri_vec[i]->mId0, tri_vec[i]->mId1,
                                         tri_vec[i]->mId2);
        if ((cur - target).norm() < 1e-3)
        {
            return i;
        }
    }
    return -1;
}
int FindEdge(tEdgePtrVector edge_vec, int v0, int v1)
{

    for (int i = 0; i < edge_vec.size(); i++)
    {
        auto cur_e = edge_vec[i];
        if ((v0 == cur_e->mId0 && v1 == cur_e->mId1) ||
            (v1 == cur_e->mId0 && v0 == cur_e->mId1))
        {
            return i;
        }
    }
    return -1;
}
#include "utils/ColorUtil.h"
void cTetUtil::LoadTet(const std::string &path, tVertexPtrVector &vertex_vec,
                       tEdgePtrVector &edge_vec, tTrianglePtrVector &tri_vec,
                       tTetPtrVector &tet_vec)

{
    vertex_vec.clear();
    edge_vec.clear();
    tri_vec.clear();
    tet_vec.clear();
    SIM_DEBUG("begin to load tet from {}", path);
    // SIM_ASSERT(cFileUtil::ExistsFile(path.c_str()));
    std::string vert_file = path + ".node";
    std::string edge_file = path + ".edge";
    std::string face_file = path + ".face";
    std::string ele_file = path + ".ele";
    // 1. load vertices
    auto vertices_lines = cFileUtil::ReadFileAllLines(vert_file);
    cStringUtil::RemoveCommentLine(vertices_lines);
    cStringUtil::RemoveEmptyLine(vertices_lines);
    // 1.1 load the first line
    std::istringstream istr;
    int num_of_vertices = 0;
    for (int line_id = 0; line_id < vertices_lines.size(); line_id++)
    {
        auto cur_line = vertices_lines[line_id];
        istr.clear();
        istr.str(cur_line);
        if (line_id == 0)
        {
            istr >> num_of_vertices;
            SIM_ASSERT(vertices_lines.size() - 1 == num_of_vertices);
            std::cout << "num of vertices = " << num_of_vertices << std::endl;
        }
        else
        {
            int v_id;
            auto cur_v = std::make_shared<tVertex>();
            cur_v->mPos = tVector::Ones();
            cur_v->mColor = ColorAn;
            // cur_v->mMass = 1.0;
            istr >> v_id >> cur_v->mPos[0] >> cur_v->mPos[1] >> cur_v->mPos[2];
            // std::cout << "v " << v_id << " = " << cur_v->mPos.transpose() <<
            // std::endl;
            vertex_vec.push_back(cur_v);
        }
    }

    int num_of_edges = 0;
    auto edge_lines = cFileUtil::ReadFileAllLines(edge_file);
    cStringUtil::RemoveCommentLine(edge_lines);
    cStringUtil::RemoveEmptyLine(edge_lines);

    for (int line_id = 0; line_id < edge_lines.size(); line_id++)
    {
        auto cur_line = edge_lines[line_id];
        istr.clear();
        istr.str(cur_line);
        if (line_id == 0)
        {
            istr >> num_of_edges;
            SIM_ASSERT(edge_lines.size() - 1 == num_of_edges);
            // std::cout << "num of edges = " << num_of_edges << std::endl;
        }
        else
        {
            int e_id;
            auto cur_e = std::make_shared<tEdge>();
            istr >> e_id >> cur_e->mId0 >> cur_e->mId1;
            // printf("edge %d from %d to %d\n", e_id, cur_e->mId0,
            // cur_e->mId1);
            edge_vec.push_back(cur_e);
        }
    }

    // 3. build triangles
    int num_of_triangles = 0;
    auto tri_lines = cFileUtil::ReadFileAllLines(face_file);
    cStringUtil::RemoveCommentLine(tri_lines);
    cStringUtil::RemoveEmptyLine(tri_lines);
    for (int line_id = 0; line_id < tri_lines.size(); line_id++)
    {
        auto cur_line = tri_lines[line_id];
        istr.clear();
        istr.str(cur_line);
        if (line_id == 0)
        {
            istr >> num_of_triangles;
            SIM_ASSERT(tri_lines.size() - 1 == num_of_triangles);
            // std::cout << "num of triangles = " << num_of_triangles <<
            // std::endl;
        }
        else
        {
            int f_id;
            auto cur_t = std::make_shared<tTriangle>();
            int boundary_marker = 0; // 0: internal, -1: boundary
            istr >> f_id >> cur_t->mId0 >> cur_t->mId1 >> cur_t->mId2 >>
                boundary_marker;

            {
                int edge_id = FindEdge(edge_vec, cur_t->mId0, cur_t->mId1);
                if (edge_vec[edge_id]->mTriangleId0 == -1)
                    edge_vec[edge_id]->mTriangleId0 = f_id;
                else
                {
                    if (edge_vec[edge_id]->mTriangleId1 == -1)
                        edge_vec[edge_id]->mTriangleId1 = f_id;
                }
            }
            {
                int edge_id = FindEdge(edge_vec, cur_t->mId1, cur_t->mId2);
                if (edge_vec[edge_id]->mTriangleId0 == -1)
                    edge_vec[edge_id]->mTriangleId0 = f_id;
                else
                {
                    if (edge_vec[edge_id]->mTriangleId1 == -1)
                        edge_vec[edge_id]->mTriangleId1 = f_id;
                }
            }
            {
                int edge_id = FindEdge(edge_vec, cur_t->mId0, cur_t->mId2);
                if (edge_vec[edge_id]->mTriangleId0 == -1)
                    edge_vec[edge_id]->mTriangleId0 = f_id;
                else
                {
                    if (edge_vec[edge_id]->mTriangleId1 == -1)
                        edge_vec[edge_id]->mTriangleId1 = f_id;
                }
            }
            // printf("triangle %d contains node %d %d %d\n",
            //        f_id,
            //        cur_t->mId0,
            //        cur_t->mId1,
            //        cur_t->mId2);
            cur_t->mColor = ColorBlue;
            tri_vec.push_back(cur_t);
        }
    }

    // 4. build tets
    int num_of_tets = 0;
    auto tet_lines = cFileUtil::ReadFileAllLines(ele_file);
    cStringUtil::RemoveCommentLine(tet_lines);
    cStringUtil::RemoveEmptyLine(tet_lines);
    for (int line_id = 0; line_id < tet_lines.size(); line_id++)
    {
        auto cur_line = tet_lines[line_id];
        istr.clear();
        istr.str(cur_line);
        if (line_id == 0)
        {
            istr >> num_of_tets;
            SIM_ASSERT(tet_lines.size() - 1 == num_of_tets);
            std::cout << "num of tets = " << num_of_tets << std::endl;
        }
        else
        {
            int t_id;
            auto cur_t = std::make_shared<tTet>();
            istr >> t_id >> cur_t->mVertexId[0] >> cur_t->mVertexId[1] >>
                cur_t->mVertexId[2] >> cur_t->mVertexId[3];
            // printf("tet %d contains node %d %d %d %d\n",
            //        t_id,
            //        cur_t->mVertexId[0],
            //        cur_t->mVertexId[1],
            //        cur_t->mVertexId[2],
            //        cur_t->mVertexId[3]);
            // triangle0 : v0, v1, v2
            cur_t->mTriangleId[0] =
                FindTriangle(tri_vec, cur_t->mVertexId[0], cur_t->mVertexId[1],
                             cur_t->mVertexId[2]);
            cur_t->mTriangleId[1] =
                FindTriangle(tri_vec, cur_t->mVertexId[1], cur_t->mVertexId[3],
                             cur_t->mVertexId[2]);
            cur_t->mTriangleId[2] =
                FindTriangle(tri_vec, cur_t->mVertexId[2], cur_t->mVertexId[3],
                             cur_t->mVertexId[0]);
            cur_t->mTriangleId[3] =
                FindTriangle(tri_vec, cur_t->mVertexId[3], cur_t->mVertexId[1],
                             cur_t->mVertexId[0]);
            tet_vec.push_back(cur_t);
            // std::cout << "tet triangle = " << cur_t->mTriangleId.transpose()
            // << std::endl;
        }
    }
    AdjustTriangleIndexOrderForSurface(vertex_vec, tri_vec, tet_vec);
}
#include <set>
void cTetUtil::AdjustTriangleIndexOrderForSurface(tVertexPtrVector &vertex_vec,
                                                  tTrianglePtrVector &tri_vec,
                                                  tTetPtrVector &tet_vec)
{
    // 1. recognize boundary triangle
    std::vector<int> triangle_involved_times(tri_vec.size(), 0);
    std::vector<int> triangles_involved_tet(tri_vec.size(), -1);
    for (int i = 0; i < tet_vec.size(); i++)
    {
        auto tet = tet_vec[i];
        triangle_involved_times[tet->mTriangleId[0]] += 1;
        triangle_involved_times[tet->mTriangleId[1]] += 1;
        triangle_involved_times[tet->mTriangleId[2]] += 1;
        triangle_involved_times[tet->mTriangleId[3]] += 1;

        triangles_involved_tet[tet->mTriangleId[0]] = i;
        triangles_involved_tet[tet->mTriangleId[1]] = i;
        triangles_involved_tet[tet->mTriangleId[2]] = i;
        triangles_involved_tet[tet->mTriangleId[3]] = i;
    }

    // 2. check surface triangle
    for (int i = 0; i < tri_vec.size(); i++)
    {
        if (triangle_involved_times[i] == 1)
        {
            // surface triangle
            auto cur_tet = tet_vec[triangles_involved_tet[i]];
            auto cur_tri = tri_vec[i];
            int tri_vid0 = cur_tri->mId0;
            int tri_vid1 = cur_tri->mId1;
            int tri_vid2 = cur_tri->mId2;
            tVector tri_v0 = vertex_vec[tri_vid0]->mPos;
            tVector tri_v1 = vertex_vec[tri_vid1]->mPos;
            tVector tri_v2 = vertex_vec[tri_vid2]->mPos;
            std::set<int> triangle_vid = {};
            std::set<int> tet_vid = {
                cur_tet->mVertexId[0], cur_tet->mVertexId[1],
                cur_tet->mVertexId[2], cur_tet->mVertexId[3]};
            tet_vid.erase(cur_tri->mId0);
            tet_vid.erase(cur_tri->mId1);
            tet_vid.erase(cur_tri->mId2);
            int another_vid = *tet_vid.begin();

            // judge whether this triangle is inner pointed
            bool inner_surface =
                ((tri_v1 - tri_v0).cross3(tri_v2 - tri_v1))
                    .dot(vertex_vec[another_vid]->mPos - tri_v0) > 0;
            if (inner_surface)
            {
                // printf("[check] surface tri %d is oppo, correct.\n", i);
                SIM_SWAP(cur_tri->mId0, cur_tri->mId2);
                SIM_SWAP(cur_tri->mEId0, cur_tri->mEId1);
            }
        }
        // auto cur_tri = tri_vec[i];
    }
}