#include "Triangulator.h"
#include "utils/ColorUtil.h"
#include "utils/FileUtil.h"
#include "utils/JsonUtil.h"
#include "utils/ObjUtil.h"
#include <iostream>
#include <map>
void cTriangulator::BuildGeometry(const Json::Value &config,
                                  std::vector<tVertex *> &vertices_array,
                                  std::vector<tEdge *> &edges_array,
                                  std::vector<tTriangle *> &triangles_array)
{
    std::string geo_type =
        cJsonUtil::ParseAsString(cTriangulator::GEOMETRY_TYPE_KEY, config);
    tVector2d cloth_shape =
        cJsonUtil::ReadVectorJson(cJsonUtil::ParseAsValue("cloth_size", config))
            .segment(0, 2);
    tVector cloth_init_pos = tVector::Zero();
    tVector cloth_init_orientation = tVector::Zero();
    float cloth_uv_rotation = cJsonUtil::ParseAsFloat(
        "cloth_uv_rotation", config); // rotate the martieral coordinates.

    /*
        the default setting:
        warp: x+
        weft: y+
    */
    {
        cJsonUtil::ReadVectorJson(
            cJsonUtil::ParseAsValue("cloth_init_pos", config), cloth_init_pos);
        cJsonUtil::ReadVectorJson(
            cJsonUtil::ParseAsValue("cloth_init_orientation", config),
            cloth_init_orientation);
    }
    // std::cout << "cloth init pos = " << cloth_init_pos.transpose() <<
    // std::endl; std::cout << "cloth init orientation = " <<
    // cloth_init_orientation.transpose() << std::endl;
    tMatrix init_trans_mat =
        cMathUtil::TransformMat(cloth_init_pos, cloth_init_orientation);
    // std::cout << init_trans_mat << std::endl;
    // exit(0);
    Json::Value subdivision_json =
        cJsonUtil::ParseAsValue("subdivision", config);
    tVector2i subdivision =
        cJsonUtil::ReadVectorJson(subdivision_json).segment(0, 2).cast<int>();
    // if (geo_type == "uniform_square")
    // {
    //     SIM_ERROR("geo type uniform_square has been deprecated, because it "
    //               "doesn't support bending");
    //     exit(0);
    //     cTriangulator::BuildGeometry_UniformSquare(cloth_shape, subdivision,
    //                                                vertices_array,
    //                                                edges_array,
    //                                                triangles_array);
    // }
    // else if (geo_type == "skew_triangle")
    // {
    //     cTriangulator::BuildGeometry_SkewTriangle(cloth_shape, subdivision,
    //                                               vertices_array,
    //                                               edges_array,
    //                                               triangles_array);
    // }
    // else
    if (geo_type == "regular_triangle")
    {
        cTriangulator::BuildGeometry_UniformTriangle(
            cloth_shape, subdivision, vertices_array, edges_array,
            triangles_array, false);
        // exit(1);
    }
    else if (geo_type == "regular_triangle_perturb")
    {
        cTriangulator::BuildGeometry_UniformTriangle(
            cloth_shape, subdivision, vertices_array, edges_array,
            triangles_array, true);
    }
    else
    {
        SIM_ERROR("unsupported geo type {}", geo_type);
    }
    ValidateGeometry(vertices_array, edges_array, triangles_array);
    for (auto &v : vertices_array)
    {
        v->mPos = init_trans_mat * v->mPos;
        // v->mPos.segment(0, 3) += cloth_init_pos.segment(0, 3);
    }
    RotateMaterialCoordsAfterReset(init_trans_mat.inverse(), vertices_array,
                                   cloth_uv_rotation);
    // support vertices
    // printf(
    //     "[debug] init geometry type %s, create %d vertices, %d edges, %d
    //     triangles\n", geo_type.c_str(), vertices_array.size(),
    //     edges_array.size(), triangles_array.size());
    // exit(0);
}

void cTriangulator::BuildGeometry_UniformTriangle(
    const tVector2d &mesh_shape, const tVector2i &subdivision,
    std::vector<tVertex *> &vertices_array, std::vector<tEdge *> &edges_array,
    std::vector<tTriangle *> &triangles_array, bool add_vertices_perturb)
{
    // 1. clear all
    vertices_array.clear();
    edges_array.clear();
    triangles_array.clear();

    double height = mesh_shape.x();
    double width = mesh_shape.y();
    int num_of_height_div = subdivision.x();
    int num_of_width_div = subdivision.y();
    double unit_edge_h = height / num_of_height_div;
    double unit_edge_w = width / num_of_width_div;
    double unit_edge_skew =
        std::sqrt(unit_edge_h * unit_edge_h + unit_edge_w * unit_edge_w);
    /*
    (0, 0), HEIGHT dimension, col
    --------------------------- y+ world frame y axis, texture y axis
    |                           cartesian pos (num_of_height_div, 0)
    |
    |
    |
    |
    |
    |
    |
    | WIDTH, world frame x axis, texture x axis, row
    x+
    cartesian pos (0, num_of_width_div)
    */

    // 2. create vertices
    BuildRectVertices(height, width, num_of_height_div, num_of_width_div,
                      vertices_array, add_vertices_perturb);

    // 3. create triangles
    int num_of_width_lines = num_of_width_div + 1;
    int num_of_height_lines = num_of_height_div + 1;

    int num_of_vertices = num_of_width_lines *
                          num_of_height_lines; // the number of chessboard nodes
    int num_of_edges = num_of_width_lines * num_of_height_div +
                       num_of_height_lines * num_of_width_div +
                       num_of_height_div * num_of_width_div;
    // horizontal lines + vertical lines + skew lines
    int num_of_triangles = num_of_width_div * num_of_height_div * 2;
    edges_array.resize(num_of_edges, nullptr);

    // 1. init the triangles
    int num_edges_per_row = 2 * num_of_height_div + num_of_height_lines;
    for (int row_id = 0; row_id < num_of_width_div; row_id++)
    {
        for (int col_id = 0; col_id < num_of_height_div; col_id++)
        {
            // for even number, from upleft to downright
            int left_up_vid = row_id * num_of_height_lines + col_id;
            int right_up_vid = left_up_vid + 1;
            int left_down_vid = left_up_vid + num_of_height_lines;
            int right_down_vid = left_down_vid + 1;
            bool is_even = (row_id + col_id) % 2 == 0;
            // printf("-----[debug] for block row %d col %d---\n", row_id,
            // col_id);
            /*
            Even case
            ---------
            | \     |
            |   \   |
            |     \ |
            --------
            */
            int top_edge_id = row_id * num_edges_per_row + col_id;
            int left_edge_id =
                num_edges_per_row * row_id + num_of_height_div + col_id * 2;
            int skew_edge_id = left_edge_id + 1;
            int right_edge_id = skew_edge_id + 1;
            int bottom_edge_id = top_edge_id + num_edges_per_row;

            bool need_top_edge = false, need_left_edge = false;
            if (row_id == 0)
            {
                need_top_edge = true;
            }
            if (col_id == 0)
            {
                need_left_edge = true;
            }

            tEdge *skew_edge = new tEdge();
            if (edges_array[skew_edge_id] != nullptr)
            {
                SIM_ERROR("wrong visit in skew edge {}", skew_edge_id);
            }
            edges_array[skew_edge_id] = skew_edge;
            if (is_even)
            {
                auto tri1 =
                    new tTriangle(left_up_vid, left_down_vid, right_down_vid);
                auto tri2 =
                    new tTriangle(left_up_vid, right_down_vid, right_up_vid);

                triangles_array.push_back(tri1);
                triangles_array.push_back(tri2);

                skew_edge->mId0 = left_up_vid;
                skew_edge->mId1 = right_down_vid;
            }
            else
            {
                /*
                Odd case
                ---------
                |     / |
                |   /   |
                | /     |
                ---------
                */
                // for odd number, from upright to downleft
                auto tri1 =
                    new tTriangle(left_up_vid, left_down_vid, right_up_vid);
                auto tri2 =
                    new tTriangle(left_down_vid, right_down_vid, right_up_vid);
                triangles_array.push_back(tri1);
                triangles_array.push_back(tri2);
                skew_edge->mId0 = right_up_vid;
                skew_edge->mId1 = left_down_vid;
            }
            skew_edge->mRawLength = unit_edge_skew;
            skew_edge->mIsBoundary = false;
            skew_edge->mTriangleId0 = triangles_array.size() - 2;
            skew_edge->mTriangleId1 = triangles_array.size() - 1;
            // printf("[debug] add skew %d edge v_id %d to %d, tri id %d %d\n",
            //        skew_edge_id, skew_edge->mId0, skew_edge->mId1,
            //        skew_edge->mTriangleId0, skew_edge->mTriangleId1);
            // add bottom edge
            {
                tEdge *bottom_edge = new tEdge();
                SIM_ASSERT(edges_array[bottom_edge_id] == nullptr);
                edges_array[bottom_edge_id] = bottom_edge;
                bottom_edge->mId0 = left_down_vid;
                bottom_edge->mId1 = right_down_vid;
                bottom_edge->mRawLength = unit_edge_h;

                if (is_even)
                {
                    bottom_edge->mTriangleId0 = triangles_array.size() - 2;
                }
                else
                {
                    bottom_edge->mTriangleId0 = triangles_array.size() - 1;
                }
                bottom_edge->mTriangleId1 =
                    bottom_edge->mTriangleId0 + num_of_height_div * 2;
                if (row_id == num_of_width_div - 1)
                {
                    bottom_edge->mIsBoundary = true;
                    bottom_edge->mTriangleId1 = -1;
                }
                else
                {
                    bottom_edge->mIsBoundary = false;
                }
                // printf(
                //     "[debug] add bottom %d edge v_id %d to %d, tri id %d
                //     %d\n", bottom_edge_id, bottom_edge->mId0,
                //     bottom_edge->mId1, bottom_edge->mTriangleId0,
                //     bottom_edge->mTriangleId1);
            }

            // add right edge
            {
                tEdge *right_edge = new tEdge();
                if (edges_array[right_edge_id] != nullptr)
                {
                    SIM_ERROR("wrong visit in right edge {}", right_edge_id);
                }
                edges_array[right_edge_id] = right_edge;
                if (col_id == num_of_height_div - 1)
                {
                    right_edge->mIsBoundary = true;
                }
                else
                {
                    right_edge->mIsBoundary = false;
                }
                right_edge->mId0 = right_up_vid;
                right_edge->mId1 = right_down_vid;
                right_edge->mRawLength = unit_edge_w;
                right_edge->mTriangleId0 = triangles_array.size() - 1;
                if (col_id == num_of_height_div - 1)
                {
                    right_edge->mTriangleId1 = -1;
                }
                else
                {
                    right_edge->mTriangleId1 = right_edge->mTriangleId0 + 1;
                }
                // printf(
                //     "[debug] add right %d edge v_id %d to %d, tri id %d
                //     %d\n", right_edge_id, right_edge->mId0, right_edge->mId1,
                //     right_edge->mTriangleId0, right_edge->mTriangleId1);
            }

            // add top edge
            if (need_top_edge)
            {
                tEdge *top_edge = new tEdge();
                SIM_ASSERT(edges_array[top_edge_id] == nullptr);
                edges_array[top_edge_id] = top_edge;
                top_edge->mId0 = left_up_vid;
                top_edge->mId1 = right_up_vid;
                top_edge->mRawLength = unit_edge_h;
                top_edge->mIsBoundary = row_id == 0;
                if (is_even)
                {
                    top_edge->mTriangleId0 = triangles_array.size() - 1;
                }
                else
                {
                    top_edge->mTriangleId0 = triangles_array.size() - 2;
                }
                top_edge->mTriangleId1 = -1;
                // printf("[debug] add top %d edge v_id %d to %d, tri id %d
                // %d\n",
                //        top_edge_id, top_edge->mId0, top_edge->mId1,
                //        top_edge->mTriangleId0, top_edge->mTriangleId1);
            }

            // add left edge
            if (need_left_edge)
            {
                tEdge *left_edge = new tEdge();
                SIM_ASSERT(edges_array[left_edge_id] == nullptr);
                edges_array[left_edge_id] = left_edge;
                left_edge->mId0 = left_up_vid;
                left_edge->mId1 = left_down_vid;
                left_edge->mRawLength = unit_edge_w;
                left_edge->mIsBoundary = col_id == 0;
                left_edge->mTriangleId0 = triangles_array.size() - 2;
                left_edge->mTriangleId1 = -1;
                // printf("[debug] add left %d edge v_id %d to %d, tri id %d
                // %d\n",
                //        left_edge_id, left_edge->mId0, left_edge->mId1,
                //        left_edge->mTriangleId0, left_edge->mTriangleId1);
            }
        }
    }
}

/**
 * \brief                   create vertices as a uniform, rectangle vertices
 */
void cTriangulator::BuildRectVertices(double height, double width,
                                      int height_div, int width_div,
                                      std::vector<tVertex *> &vertices_array,
                                      bool add_vertices_perturb)
{
    vertices_array.clear();
    /*
    (0, 0), height dimension
    --------------------------- y+ world frame y axis, texture y axis
    |                           cartesian pos (height_div, 0)
    |
    |
    |
    |
    |
    |
    |
    | world frame x axis, texture x axis

    x+ width dimension
    cartesian pos (0, width_div)
    */
    double unit_edge_h = height / height_div;
    double unit_edge_w = width / width_div;

    double noise_radius_height = unit_edge_h / 5,
           noise_radius_width = unit_edge_w / 5;
    // double noise_max_radius = 0;

    tVector center_pos = tVector(width / 2, height / 2, 0, 0);
    for (int i = 0; i < width_div + 1; i++)
        for (int j = 0; j < height_div + 1; j++)
        {
            tVertex *v = new tVertex();

            // 1. first set the cartesian pos, in order to get the texture
            // coords
            v->mPos = tVector(i * unit_edge_w, j * unit_edge_h, 0, 1);
            v->mColor = ColorBlue;
            v->muv = tVector2f::Zero();

            // move the center
            v->mPos -= center_pos;
            // then add perturb
            if (add_vertices_perturb)
            {
                if (i != 0 && i != height_div && j != 0 && j != width_div)
                {
                    // std::cout << "add_vertices_perturb\n";
                    v->mPos[0] += cMathUtil::RandDouble(-noise_radius_width,
                                                        noise_radius_width);
                    v->mPos[1] += cMathUtil::RandDouble(-noise_radius_height,
                                                        noise_radius_height);
                }
            }
            // printf("[debug] add vertex (%d, %d) at pos (%.3f, %.3f, %.3f),
            // tex "
            //        "(%.2f, %.2f)\n",
            //        i, j, v->mPos[0], v->mPos[1], v->mPos[2], v->muv[0],
            //        v->muv[1]);
            vertices_array.push_back(v);
        }
}

bool ConfirmVertexInTriangles(tTriangle *tri, int vid)
{
    return (tri->mId0 == vid) || (tri->mId1 == vid) || (tri->mId2 == vid);
};
void cTriangulator::ValidateGeometry(std::vector<tVertex *> &vertices_array,
                                     std::vector<tEdge *> &edges_array,
                                     std::vector<tTriangle *> &triangles_array)
{
    // confirm the edges is really shared by triangles
    for (int i = 0; i < edges_array.size(); i++)
    {
        auto &e = edges_array[i];
        if (e->mTriangleId0 != -1)
        {
            auto tri = triangles_array[e->mTriangleId0];
            if (ConfirmVertexInTriangles(tri, e->mId0

                                         ) &&
                ConfirmVertexInTriangles(tri, e->mId1) == false)
            {
                printf("[error] validate boundary edge %d's two vertices %d "
                       "and %d doesn't located in triangle %d\n",
                       i, e->mId0, e->mId1, e->mTriangleId0);
                std::cout << "triangle vertices idx list = " << tri->mId0
                          << ", " << tri->mId1 << ", " << tri->mId2 << "\n";
                exit(0);
            }
        }
        if (e->mTriangleId1 != -1)
        {
            auto tri = triangles_array[e->mTriangleId1];
            if ((ConfirmVertexInTriangles(tri, e->mId0) &&
                 ConfirmVertexInTriangles(tri, e->mId1)) == false)
            {
                printf("[error] validate boundary edge %d's two vertices %d "
                       "and %d doesn't located in triangle %d\n",
                       i, e->mId0, e->mId1, e->mTriangleId1);
                std::cout << "triangle vertices idx list = " << tri->mId0
                          << ", " << tri->mId1 << ", " << tri->mId2 << "\n";
                exit(0);
            }
        }
    }
}

/**
 * \brief           Given the geometry info, save them to the given "path"
 *
 *          Only save a basic info
 */
void cTriangulator::SaveGeometry(std::vector<tVertex *> &vertices_array,
                                 std::vector<tEdge *> &edges_array,
                                 std::vector<tTriangle *> &triangles_array,
                                 const std::string &path)
{
    Json::Value root;
    // 1. the vertices info
    root[NUM_OF_VERTICES_KEY] = static_cast<int>(vertices_array.size());

    // 2. the edge info

    root[EDGE_ARRAY_KEY] = Json::arrayValue;
    for (auto &x : edges_array)
    {
        root[EDGE_ARRAY_KEY].append(x->mId0);
        root[EDGE_ARRAY_KEY].append(x->mId1);
    }
    // 3. the triangle info
    root[TRIANGLE_ARRAY_KEY] = Json::arrayValue;
    for (auto &x : triangles_array)
    {
        Json::Value tri = Json::arrayValue;

        tri.append(x->mId0);
        tri.append(x->mId1);
        tri.append(x->mId2);
        root[TRIANGLE_ARRAY_KEY].append(tri);
    }
    std::cout << "[debug] save geometry to " << path << std::endl;
    cJsonUtil::WriteJson(path, root, true);
}

void cTriangulator::LoadGeometry(std::vector<tVertex *> &vertices_array,
                                 std::vector<tEdge *> &edges_array,
                                 std::vector<tTriangle *> &triangles_array,
                                 const std::string &path)
{
    vertices_array.clear();
    edges_array.clear();
    triangles_array.clear();
    Json::Value root;
    cJsonUtil::LoadJson(path, root);
    int num_of_vertices = cJsonUtil::ParseAsInt(NUM_OF_VERTICES_KEY, root);
    for (int i = 0; i < num_of_vertices; i++)
    {
        vertices_array.push_back(new tVertex());
        vertices_array[vertices_array.size() - 1]->mColor = ColorBlue;
    }

    const tVectorXd &edge_info = cJsonUtil::ReadVectorJson(
        cJsonUtil::ParseAsValue(EDGE_ARRAY_KEY, root));

    for (int i = 0; i < edge_info.size() / 2; i++)
    {
        tEdge *edge = new tEdge();
        edge->mId0 = edge_info[2 * i + 0];
        edge->mId1 = edge_info[2 * i + 1];
        edges_array.push_back(edge);
    }

    const Json::Value &tri_json =
        cJsonUtil::ParseAsValue(TRIANGLE_ARRAY_KEY, root);
    for (int i = 0; i < tri_json.size(); i++)
    {
        tTriangle *tri = new tTriangle();
        tri->mId0 = tri_json[i][0].asInt();
        tri->mId1 = tri_json[i][1].asInt();
        tri->mId2 = tri_json[i][2].asInt();
        triangles_array.push_back(tri);
    }

    // 1. build edges's triangle
    {
        BuildEdgesTriangleId(edges_array, triangles_array);
    }
    printf("[debug] Load Geometry from %s done, vertices %d, edges %d, "
           "triangles %d\n",
           path.c_str(), vertices_array.size(), edges_array.size(),
           triangles_array.size());
}

void cTriangulator::RotateMaterialCoordsAfterReset(
    const tMatrix &init_mat_inv, std::vector<tVertex *> &vertices_array,
    float cloth_uv_rotation_)
{
    // degrees
    SIM_ASSERT(int(cloth_uv_rotation_) == 0 || int(cloth_uv_rotation_) == 45 ||
               int(cloth_uv_rotation_) == 90);

    // rad

    float cloth_uv_rotation = cloth_uv_rotation_ / 180.0 * M_PI;
    tMatrix2f rot_mat = cMathUtil::RotMat2D(cloth_uv_rotation).cast<float>();
    // std::cout << "angle = " << cloth_uv_rotation << ", rotmat = \n"
    //           << rot_mat << std::endl;

    // std::cout << "warp dir = " << warp_dir.transpose() << std::endl;
    // std::cout << "weft dir = " << weft_dir.transpose() << std::endl;
    // std::cout << "bias dir = " << bias_dir.transpose() << std::endl;
    // exit(1);
    bool origin_has_been_set = false;
    // for (auto &x : vertices_array)
    // std::cout << "init_mat_inv = \n" << init_mat_inv << std::endl;
    for (int i = 0; i < vertices_array.size(); i++)
    {
        auto x = vertices_array[i];
        // default origin is (0, 0)
        tVector2f cur_pos =
            (init_mat_inv * x->mPos).segment(0, 2).cast<float>();
        x->muv = rot_mat * cur_pos;
        // if (i < 10)
        //     std::cout << "i " << i << " cur pos " << cur_pos.transpose()
        //               << " uv = " << x->muv.transpose() << std::endl;
        // std::cout << "x pos = " << cur_pos.transpose()
        //           << ", uv = " << x->muv.transpose() << std::endl;
        // std::cout << x->muv.transpose() << std::endl;
    }
    // exit(1);
}

/**
 * \brief       given edges, calculate edges' triangle id
 */

void cTriangulator::BuildEdgesTriangleId(
    std::vector<tEdge *> &edges_array,
    std::vector<tTriangle *> &triangles_array)
{

    // 1. build the map from edge's vertices id to edge id
    typedef std::pair<int, int> tEdgeVerticeId;
    std::map<tEdgeVerticeId, int> edgeverticesid_to_edgeid;
    for (int i = 0; i < edges_array.size(); i++)
    {
        auto &e = edges_array[i];
        int id0 = e->mId0, id1 = e->mId1;
        if (id0 > id1)
        {
            int tmp = id1;
            id1 = id0;
            id0 = tmp;
        }

        // id0 < id1
        if ((id0 < id1) == false)
        {
            printf("edge %d id0 %d id1 %d, illegal\n", i, id0, id1);
            exit(1);
        }
        auto edge_vertice_info = tEdgeVerticeId(id0, id1);
        SIM_ASSERT(edgeverticesid_to_edgeid.end() ==
                   edgeverticesid_to_edgeid.find(edge_vertice_info));

        edgeverticesid_to_edgeid[edge_vertice_info] = i;
    }

    std::map<int, std::vector<int>> edgeid_to_triid;
    for (int i = 0; i < triangles_array.size(); i++)
    {
        auto &tri = triangles_array[i];
        int id0 = tri->mId0, id1 = tri->mId1, id2 = tri->mId2;
        std::vector<tEdgeVerticeId> mVerticeInfo;
        mVerticeInfo.push_back(tEdgeVerticeId(id0, id1));
        mVerticeInfo.push_back(tEdgeVerticeId(id1, id2));
        mVerticeInfo.push_back(tEdgeVerticeId(id2, id0));

        for (int _idx = 0; _idx < mVerticeInfo.size(); _idx++)
        {
            int vid0 = mVerticeInfo[_idx].first;
            int vid1 = mVerticeInfo[_idx].second;
            if (vid0 > vid1)
            {
                std::swap(vid0, vid1);
            }
            int edge_id = edgeverticesid_to_edgeid[tEdgeVerticeId(vid0, vid1)];
            if (edgeid_to_triid.find(edge_id) == edgeid_to_triid.end())
            {
                auto tmp = std::vector<int>();
                tmp.push_back(i);
                edgeid_to_triid[edge_id] = tmp;
            }
            else
            {
                edgeid_to_triid[edge_id].push_back(i);
            }
        }
    }

    // configure edge's info
    for (int i = 0; i < edgeid_to_triid.size(); i++)
    {
        auto &tri_id = edgeid_to_triid[i];
        SIM_ASSERT(tri_id.size() >= 1);
        SIM_ASSERT(tri_id.size() <= 2);
        // printf("edge %d belongs to triangle %d", i, tri_id[0]);
        if (tri_id.size() == 2)
        {
            int tri0 = tri_id[0];
            int tri1 = tri_id[1];
            if (tri0 > tri1)
            {
                int tmp = tri1;
                tri1 = tri0;
                tri0 = tmp;
            }
            edges_array[i]->mTriangleId0 = tri0;
            edges_array[i]->mTriangleId1 = tri1;
            edges_array[i]->mIsBoundary = false;
        }
        else
        {
            int tri0 = tri_id[0];
            edges_array[i]->mTriangleId0 = tri0;
            edges_array[i]->mTriangleId1 = -1;
            edges_array[i]->mIsBoundary = true;
        }
    }
}
void cTriangulator::RotateMaterialCoords(float cur_uv_rot_deg,
                                         float tar_uv_rot_deg,
                                         std::vector<tVertex *> &vertices_array)
{
    // std::cout << "---------------------\n";
    tMatrix2f convert_mat =
        (cMathUtil::RotMat2D(tar_uv_rot_deg / 180 * M_PI) *
         cMathUtil::RotMat2D(cur_uv_rot_deg / 180 * M_PI).inverse())
            .cast<float>();
    // printf("begin to rot from %.1f to %.1f\n", cur_uv_rot_deg,
    // tar_uv_rot_deg); std::cout << "rotmat = " << convert_mat << std::endl;

    for (int idx = 0; idx < vertices_array.size(); idx++)
    {
        bool output = false;
        if ((idx == 36) || (idx == 37) || (idx == 776) || (idx == 777))
        {
            output = true;
        }
        auto v = vertices_array[idx];

        if (output)
        {
            std::cout << "[tri] vertex " << idx
                      << " raw uv = " << v->muv.transpose() << std::endl;
        }
        v->muv = convert_mat * v->muv;
        if (output)
        {
            std::cout << "[tri] vertex " << idx
                      << " new uv = " << v->muv.transpose() << std::endl;
        }
        // std::cout << "new uv = " << v->muv.transpose() << std::endl;
    }
    // std::cout << "---------------------\n";
    // exit(1);
}
