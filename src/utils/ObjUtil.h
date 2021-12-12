#pragma once
#include "utils/MathUtil.h"
#include <string>
#include <vector>

struct tTriangle;
struct tEdge;
struct tVertex;

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
                        std::vector<tVertex *> &mVertexArray,
                        std::vector<tEdge *> &mEdgeArray,
                        std::vector<tTriangle *> &mTriangleArray);
    static void
    BuildPlaneGeometryData(const double scale, const tVector &plane_equation,
                           std::vector<tVertex *> &mVertexArray,
                           std::vector<tEdge *> &mEdgeArray,
                           std::vector<tTriangle *> &mTriangleArray);

    static void BuildEdge(const std::vector<tVertex *> &mVertexArray,
                          std::vector<tEdge *> &mEdgeArray,
                          const std::vector<tTriangle *> &mTriangleArray);

};