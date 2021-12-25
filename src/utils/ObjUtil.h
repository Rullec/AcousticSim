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
    struct tParams
    {
        std::string mPath; // obj file path
    };

    static void LoadObj(const tParams &param,
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

};