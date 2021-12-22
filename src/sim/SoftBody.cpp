#include "SoftBody.h"
#include "utils/TetUtil.h"
#include "geometries/Tetrahedron.h"
cSoftBody::cSoftBody()
{
}

cSoftBody::~cSoftBody()
{
}

/**
 * \brief           create and allocate the softbody data
*/
void cSoftBody::Init(const Json::Value &conf)
{
    cBaseObject::Init(conf);

    mTetMeshPath = cJsonUtil::ParseAsString(conf, "tet_path");
    tVertexPtrVector vertex_vec = {};
    tEdgePtrVector edge_vec = {};
    tTrianglePtrVector tri_vec = {};
    tTetPtrVector tet_ve = {};
    cTetUtil::LoadTet(path,
                      vertex_vec,
                      edge_vec,
                      tri_vec,
                      tet_vec);
    std::cout << "load from path {} " << path << "done\n";
    exit(1);
}