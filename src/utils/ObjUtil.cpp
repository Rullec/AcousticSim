#include "utils/ObjUtil.h"
#include "utils/RotUtil.h"
#define TINYOBJLOADER_IMPLEMENTATION
#include "geometries/Primitives.h"
#include "geometries/Triangulator.h"
#include "tinyobjloader/tiny_obj_loader.h"
#include "utils/LogUtil.h"
#include <iostream>

void cObjUtil::LoadObj(const cObjUtil::tParams &param,
                       std::vector<tVertexPtr> &v_array,
                       std::vector<tEdgePtr> &e_array,
                       std::vector<tTrianglePtr> &t_array)
{

    std::string path = param.mPath;

    v_array.clear();
    e_array.clear();
    t_array.clear();
    // example code from https://github.com/tinyobjloader/tinyobjloader
    std::string inputfile = path;
    tinyobj::ObjReaderConfig reader_config;
    reader_config.mtl_search_path = "./"; // Path to material files

    tinyobj::ObjReader reader;

    if (!reader.ParseFromFile(inputfile, reader_config))
    {
        if (!reader.Error().empty())
        {
            std::cerr << "TinyObjReader: " << reader.Error();
        }
        exit(1);
    }

    if (!reader.Warning().empty())
    {
        std::cout << "TinyObjReader: " << reader.Warning();
    }

    auto &attrib = reader.GetAttrib();
    auto &shapes = reader.GetShapes();
    auto &materials = reader.GetMaterials();

    SIM_ASSERT(shapes.size() == 1);
    // Loop over shapes
    // for (size_t s = 0; s < shapes.size(); s++)
    auto &shape = shapes[0];

    // int num_of_vertices = 0;
    // int num_of_edges = shape.lines.indices.size();
    // v_array.resize(num_of_vertices, nullptr);
    // e_array.resize(num_of_edges, nullptr);
    // printf("[debug] %d vertices, %d edges\n", v_array.size(),
    // e_array.size()); exit(0);
    {
        // Loop over faces(polygon)
        size_t index_offset = 0;
        // SIM_ASSERT(shape.mesh.num_face_vertices.size() == 3);
        for (size_t f = 0; f < shape.mesh.num_face_vertices.size(); f++)
        {
            // for this triangle
            size_t fv = size_t(shape.mesh.num_face_vertices[f]);
            tTrianglePtr t = std::make_shared<tTriangle>();
            t->mId0 = shape.mesh.indices[index_offset + 0].vertex_index;
            t->mId1 = shape.mesh.indices[index_offset + 1].vertex_index;
            t->mId2 = shape.mesh.indices[index_offset + 2].vertex_index;
            // we only support triangles
            SIM_ASSERT(fv == 3);
            // Loop over vertices in the face.
            for (size_t v = 0; v < fv; v++)
            {
                // access to vertex
                tinyobj::index_t idx = shape.mesh.indices[index_offset + v];
                tinyobj::real_t vx =
                    attrib.vertices[3 * size_t(idx.vertex_index) + 0];
                tinyobj::real_t vy =
                    attrib.vertices[3 * size_t(idx.vertex_index) + 1];
                tinyobj::real_t vz =
                    attrib.vertices[3 * size_t(idx.vertex_index) + 2];
                int vertex_id = idx.vertex_index;
                while (vertex_id >= v_array.size())
                    v_array.push_back(nullptr);
                if (v_array[vertex_id] == nullptr)
                {
                    v_array[vertex_id] = std::make_shared<tVertex>();
                    v_array[vertex_id]->mPos = tVector(vx, vy, vz, 1);
                }
                // // Check if `normal_index` is zero or positive. negative = no
                // normal data
                if (idx.normal_index >= 0)
                {
                    tinyobj::real_t nx =
                        attrib.normals[3 * size_t(idx.normal_index) + 0];
                    tinyobj::real_t ny =
                        attrib.normals[3 * size_t(idx.normal_index) + 1];
                    tinyobj::real_t nz =
                        attrib.normals[3 * size_t(idx.normal_index) + 2];
                    v_array[vertex_id]->mNormal =
                        tVector(nx, ny, nz, 0).normalized();
                }

                // // Check if `texcoord_index` is zero or positive. negative =
                // no texcoord data
                if (idx.texcoord_index >= 0)
                {
                    tinyobj::real_t tx =
                        attrib.texcoords[2 * size_t(idx.texcoord_index) + 0];
                    tinyobj::real_t ty =
                        attrib.texcoords[2 * size_t(idx.texcoord_index) + 1];
                    v_array[vertex_id]->muv = tVector2f(tx, ty);
                }

                // // Optional: vertex colors
                // // tinyobj::real_t red   =
                // attrib.colors[3*size_t(idx.vertex_index)+0];
                // // tinyobj::real_t green =
                // attrib.colors[3*size_t(idx.vertex_index)+1];
                // // tinyobj::real_t blue  =
                // attrib.colors[3*size_t(idx.vertex_index)+2];
            }
            index_offset += fv;

            // per-face material
            // shape.mesh.material_ids[f];
            t_array.push_back(t);
        }
    }

    for (int i = 0; i < v_array.size(); i++)
    {
        if (v_array[i] == nullptr)
        {
            SIM_ERROR("vertex {} is empty, exit", i);
            exit(1);
        }
    }
    cObjUtil::BuildEdge(v_array, e_array, t_array);
}
/**
 * \brief       Given vertex array and triangle array, build the edge list
 */
#include <set>
typedef std::pair<int, int> int_pair;
void cObjUtil::BuildEdge(const std::vector<tVertexPtr> &v_array,
                         std::vector<tEdgePtr> &e_array,
                         const std::vector<tTrianglePtr> &t_array)
{
    e_array.clear();

    // 1. build duplicate edge array
    // edge_id_pair -> triangle_id_pair
    std::map<int_pair, int_pair> edge_info;
    edge_info.clear();

    // for each triangle
    for (int t_id = 0; t_id < t_array.size(); t_id++)
    {
        tTrianglePtr tri = t_array[t_id];

        // check three edges
        for (int i = 0; i < 3; i++)
        {
            // for each edge
            int id0 = (i == 0) ? (tri->mId0)
                               : ((i == 1) ? tri->mId1 : (tri->mId2)),
                id1 = (i == 0) ? (tri->mId1)
                               : ((i == 1) ? tri->mId2 : (tri->mId0));

            if (id0 > id1)
            {
                std::swap(id0, id1);
            }
            auto edge_id_pairs = int_pair(id0, id1);
            std::map<int_pair, int_pair>::iterator it =
                edge_info.find(edge_id_pairs);

            // create new edge
            if (it == edge_info.end())
            {
                edge_info[edge_id_pairs] = int_pair(t_id, -1);
            }
            else
            {
                // use old edge
                auto &triangle_pair = it->second;
                if (false ==
                    (triangle_pair.first != -1 && triangle_pair.second == -1))
                {
                    std::cout << "triangle pair first = " << triangle_pair.first
                              << std::endl;
                    std::cout
                        << "triangle pair second = " << triangle_pair.second
                        << std::endl;

                    SIM_ERROR("build edge illegal case");
                    exit(1);
                }
                triangle_pair.second = t_id;
                // it->second = t_id;

                // SIM_ASSERT(it->second.first != -1 && it->second.second ==
                // -1); it->second.second = t_id; it->second = t_id;
            }
        }
    }

    // set dataset for edges
    std::map<int_pair, int> edge_info_2_edge_id = {};
    int t_id = 0;
    for (auto t = edge_info.begin(); t != edge_info.end(); t++, t_id++)
    {
        int v0 = t->first.first, v1 = t->first.second;
        SIM_ASSERT(v0 < v1);
        int tid0 = t->second.first, tid1 = t->second.second;
        tEdgePtr edge = std::make_shared<tEdge>();
        edge->mId0 = v0;
        edge->mId1 = v1;
        edge->mRawLength = (v_array[v0]->mPos - v_array[v1]->mPos).norm();
        edge->mTriangleId0 = tid0;
        edge->mTriangleId1 = tid1;
        edge->mIsBoundary = (tid1 == -1);
        e_array.push_back(edge);
        edge_info_2_edge_id[int_pair(v0, v1)] = t_id;
        // printf("[debug] edge %d, v0 %d, v1 %d, raw length %.3f, t0 %d, t1 %d,
        // is_boud %d\n",

        //    e_array.size() - 1, v0, v1, edge->mRawLength, tid0, tid1,
        //    edge->mIsBoundary);
    }

    for (int t = 0; t < t_array.size(); t++)
    {
        int v_id[3] = {t_array[t]->mId0, t_array[t]->mId1, t_array[t]->mId2};
        t_array[t]->mEId0 = edge_info_2_edge_id[int_pair(
            SIM_MIN(v_id[0], v_id[1]), SIM_MAX(v_id[0], v_id[1]))];
        t_array[t]->mEId1 = edge_info_2_edge_id[int_pair(
            SIM_MIN(v_id[1], v_id[2]), SIM_MAX(v_id[1], v_id[2]))];
        t_array[t]->mEId2 = edge_info_2_edge_id[int_pair(
            SIM_MIN(v_id[2], v_id[0]), SIM_MAX(v_id[2], v_id[0]))];
    }
    // std::cout << "[debug] build " << e_array.size() << " edges\n";
}

/**
 * \brief           Build plane geometry data
 */
void cObjUtil::BuildPlaneGeometryData(const double scale,
                                      const tVector &plane_equation,
                                      std::vector<tVertexPtr> &vertex_array,
                                      std::vector<tEdgePtr> &edge_array,
                                      std::vector<tTrianglePtr> &triangle_array)
{
    vertex_array.clear();
    edge_array.clear();
    triangle_array.clear();
    // 1. calculate a general vertex array
    tVector cur_normal = tVector(0, 1, 0, 0);
    tEigenArr<tVector> pos_lst = {tVector(1, 0, -1, 1), tVector(-1, 0, -1, 1),
                                  tVector(-1, 0, 1, 1), tVector(1, 0, 1, 1)};
    tEigenArr<tVector3i> triangle_idx_lst = {tVector3i(0, 1, 3),
                                             tVector3i(3, 1, 2)};

    // build vertices
    for (auto &x : pos_lst)
    {
        tVertexPtr v = std::make_shared<tVertex>();
        v->mPos.noalias() = x;
        v->mPos.segment(0, 3) *= scale;
        vertex_array.push_back(v);
    }

    // build triangles
    for (auto &x : triangle_idx_lst)
    {
        tTrianglePtr tri = std::make_shared<tTriangle>();
        tri->mId0 = x[0];
        tri->mId1 = x[1];
        tri->mId2 = x[2];
        triangle_array.push_back(tri);
    }

    cObjUtil::BuildEdge(vertex_array, edge_array, triangle_array);
    cTriangulator::ValidateGeometry(vertex_array, edge_array, triangle_array);

    // rotation

    tVector normal = cMathUtil::CalcNormalFromPlane(plane_equation);
    tMatrix transform = cRotUtil::AxisAngleToRotmat(
        cMathUtil::CalcAxisAngleFromOneVectorToAnother(cur_normal, normal));

    // translation
    {
        tVector new_pt = transform * vertex_array[0]->mPos;
        tVector3d abc = plane_equation.segment(0, 3);
        double k = (-plane_equation[3] - abc.dot(new_pt.segment(0, 3))) /
                   (abc.dot(normal.segment(0, 3)));
        transform.block(0, 3, 3, 1) = k * normal.segment(0, 3);
    }

    for (auto &x : vertex_array)
    {
        // std::cout << "old pos0 = " << x->mPos.transpose() << std::endl;
        x->mPos = transform * x->mPos;
        // std::cout << "eval = " << cMathUtil::EvaluatePlane(plane_equation,
        // x->mPos) << std::endl;
    }
    // exit(0);
}