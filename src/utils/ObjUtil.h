#pragma once
#include "utils/MathUtil.h"
#include "utils/DefUtil.h"
#include <string>
#include <vector>

struct tTriangle;
struct tEdge;
struct tVertex;
SIM_DECLARE_PTR(tTriangle);
SIM_DECLARE_PTR(tEdge);
SIM_DECLARE_PTR(tVertex);


/**
 * \brief           handle everything about obj
 */
class cObjUtil
{
public:
    // load obj
    static void LoadObj(const std::string &path,
                        std::vector<tVertexPtr> &mVertexArray,
                        std::vector<tEdgePtr> &mEdgeArray,
                        std::vector<tTrianglePtr> &mTriangleArray);
    static void
    BuildPlaneGeometryData(const double scale, const tVector &plane_equation,
                           std::vector<tVertexPtr> &mVertexArray,
                           std::vector<tEdgePtr> &mEdgeArray,
                           std::vector<tTrianglePtr> &mTriangleArray);

    static void BuildEdge(const std::vector<tVertexPtr> &mVertexArray,
                          std::vector<tEdgePtr> &mEdgeArray,
                          const std::vector<tTrianglePtr> &mTriangleArray);
    
    // export obj
    static bool ExportObj(std::string export_path,
                          const std::vector<tVertexPtr> &vertices_array,
                          const std::vector<tTrianglePtr> &triangles_array,
                          bool silent = false);
    static bool ExportObj(std::string export_path,
                          const tEigenArr<tVector> &vertices_pos_array,
                          const tEigenArr<tVector2f> &vertices_uv_array,
                          const tEigenArr<tVector3i> &triangle_vid_array,
                          bool silent = false);
};