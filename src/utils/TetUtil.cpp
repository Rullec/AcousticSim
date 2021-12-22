#include "TetUtil.h"
#include "utils/LogUtil.h"

/**
 * \brief               load tet mesh from ".node" files (custom format)
 * \param path          mesh file path
 * \param vertex_vec    all vertices
 * \param edge_vec      all edges
 * \param tri_vec       all triangles vector
 * \param tet_vec       all tets vector
*/
void cTetUtil::LoadTet(const std::string &path,
                       tVertexPtrVector &vertex_vec,
                       tEdgePtrVector &edge_vec,
                       tTrianglePtrVector &tri_vec,
                       tTetPtrVector &tet_vec)

{
    SIM_DEBUG("begin to load tet from {}", path);
    return;
}