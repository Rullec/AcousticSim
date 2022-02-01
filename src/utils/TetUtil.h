#pragma once
#include "utils/MathUtil.h"
#include "utils/DefUtil.h"
#include <string>
#include <vector>

struct tTriangle;
struct tEdge;
struct tVertex;
struct tTet;
SIM_DECLARE_PTR(tTriangle);
SIM_DECLARE_PTR(tEdge);
SIM_DECLARE_PTR(tVertex);
SIM_DECLARE_PTR(tTet);

using tTrianglePtrVector = std::vector<tTrianglePtr>;
using tEdgePtrVector = std::vector<tEdgePtr>;
using tVertexPtrVector = std::vector<tVertexPtr>;
using tTetPtrVector = std::vector<tTetPtr>;
class cTetUtil
{
public:
    static void LoadTet(const std::string &path,
                        tVertexPtrVector &vertex_vec,
                        tEdgePtrVector &edge_vec,
                        tTrianglePtrVector &tri_vec,
                        tTetPtrVector &tet_vec);
    static float CalculateTetVolume(
        const tVector &pos0,
        const tVector &pos1,
        const tVector &pos2,
        const tVector &pos3);
};