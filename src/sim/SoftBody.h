#pragma once
#include "sim/BaseObject.h"

/**
 * \brief           softbody object (mostly based on 3D FEM)
*/

struct tTriangle;
struct tEdge;
struct tVertex;
struct tTet;
class cSoftBody : public cBaseObject
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    explicit cSoftBody(int id);
    virtual ~cSoftBody();

    virtual void Init(const Json::Value &conf);
    virtual void CalcTriangleDrawBuffer(Eigen::Map<tVectorXf> &res,
                                        int &st) const override;
    virtual void CalcEdgeDrawBuffer(Eigen::Map<tVectorXf> &res,
                                    int &st) const override;

protected:
    std::string mTetMeshPath;
};

SIM_DECLARE_PTR(cSoftBody);