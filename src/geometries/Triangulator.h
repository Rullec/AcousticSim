#pragma once
#include "Primitives.h"
namespace Json
{
class Value;
};
class cTriangulator
{
public:
    inline static const std::string GEOMETRY_TYPE_KEY = "geometry_type";
    static void BuildGeometry(const Json::Value &config,
                              std::vector<tVertex *> &vertices_array,
                              std::vector<tEdge *> &edges_array,
                              std::vector<tTriangle *> &triangles_array);
    static void ValidateGeometry(std::vector<tVertex *> &vertices_array,
                                 std::vector<tEdge *> &edges_array,
                                 std::vector<tTriangle *> &triangles_array);

    static void SaveGeometry(std::vector<tVertex *> &vertices_array,
                             std::vector<tEdge *> &edges_array,
                             std::vector<tTriangle *> &triangles_array,
                             const std::string &path);
    static void LoadGeometry(std::vector<tVertex *> &vertices_array,
                             std::vector<tEdge *> &edges_array,
                             std::vector<tTriangle *> &triangles_array,
                             const std::string &path);
    static void
    RotateMaterialCoordsAfterReset(const tMatrix &init_mat_inv,
                                   std::vector<tVertex *> &vertices_array,
                                   float fabric_uv_rotation_deg);

    static void RotateMaterialCoords(float cur_uv_rot_deg, float tar_uv_rot_deg,
                                     std::vector<tVertex *> &vertices_array);

protected:
    // static void
    // BuildGeometry_UniformSquare(const tVector2d &mesh_shape,
    //                             const tVector2i &subdivision,
    //                             std::vector<tVertex *> &vertices_array,
    //                             std::vector<tEdge *> &edges_array,
    //                             std::vector<tTriangle *> &triangles_array);
    // static void
    // BuildGeometry_SkewTriangle(const tVector2d &mesh_shape,
    //                            const tVector2i &subdivision,
    //                            std::vector<tVertex *> &vertices_array,
    //                            std::vector<tEdge *> &edges_array,
    //                            std::vector<tTriangle *> &triangles_array);
    static void BuildGeometry_UniformTriangle(
        const tVector2d &mesh_shape, const tVector2i &subdivision,
        std::vector<tVertex *> &vertices_array,
        std::vector<tEdge *> &edges_array,
        std::vector<tTriangle *> &triangles_array, bool add_vertices_perturb);

    static void BuildRectVertices(double height, double width, int height_div,
                                  int width_div,
                                  std::vector<tVertex *> &edges_array,
                                  bool add_vertices_perturb);

    static void BuildEdgesTriangleId(std::vector<tEdge *> &edges_array,
                                     std::vector<tTriangle *> &triangles_array);
    inline static const std::string NUM_OF_VERTICES_KEY = "num_of_vertices",
                                    EDGE_ARRAY_KEY = "edge_array",
                                    TRIANGLE_ARRAY_KEY = "triangle_array";
};