#include "SoftBody.h"
#include "utils/JsonUtil.h"
#include "utils/TetUtil.h"
#include "geometries/Tetrahedron.h"
#include <iostream>
cSoftBody::cSoftBody(int id) : cBaseObject(eObjectType::SOFTBODY_TYPE, id)
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

    mTetMeshPath = cJsonUtil::ParseAsString("tet_path", conf);
    tVertexPtrVector vertex_vec = {};
    tEdgePtrVector edge_vec = {};
    tTrianglePtrVector tri_vec = {};
    tTetPtrVector tet_vec = {};
    cTetUtil::LoadTet(mTetMeshPath,
                      vertex_vec,
                      edge_vec,
                      tri_vec,
                      tet_vec);
    std::cout << "load from path " << mTetMeshPath << "done\n";
    exit(1);
}

void cSoftBody::CalcTriangleDrawBuffer(Eigen::Map<tVectorXf> &res,
                                       int &st) const
{
    
}
void cSoftBody::CalcEdgeDrawBuffer(Eigen::Map<tVectorXf> &res,
                                   int &st) const
{
}