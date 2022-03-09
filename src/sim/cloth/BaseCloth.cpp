#include "BaseCloth.h"
#include "geometries/Triangulator.h"
#include "utils/ColorUtil.h"
#include "utils/DefUtil.h"
#include "utils/FileUtil.h"
#include "utils/JsonUtil.h"
#include "utils/LogUtil.h"
#include "utils/RenderUtil.h"

#include <iostream>
const std::string gClothTypeStr[eClothType::NUM_OF_CLOTH_TYPE] = {
    "semi_implicit", "implicit", "pbd", "pd", "linctex", "empty", "fem", "fem_gpu"};
cBaseCloth::cBaseCloth(eClothType cloth_type, int id_)
    : cBaseObject(eObjectType::CLOTH_TYPE, id_), mClothType(cloth_type)
{

    mTriangleArrayShared.clear();
    mEdgeArrayShared.clear();
    mVertexArrayShared.clear();
    mConstraint_StaticPointIds.clear();
}

/**
 * \brief           Calculate the center of mass
 */
tVector cBaseCloth::CalcCOM() const
{
    tVector com = tVector::Zero();
    double total_mass = 0;
    for (auto &x : mVertexArrayShared)
    {
        total_mass += x->mMass;
        com += x->mMass * x->mPos;
        // std::cout << "pos = " << x->mPos.transpose() << std::endl;
    }
    com /= total_mass;
    com[3] = 1;
    return com;
}
cBaseCloth::~cBaseCloth() {}

void cBaseCloth::Init(const Json::Value &conf)
{
    cBaseObject::Init(conf);
    mGeometryType =
        cJsonUtil::ParseAsString(cTriangulator::GEOMETRY_TYPE_KEY, conf);
    mDamping = cJsonUtil::ParseAsDouble(cBaseCloth::DAMPING_KEY, conf);
    mIdealDefaultTimestep =
        cJsonUtil::ParseAsDouble(cBaseCloth::DEFAULT_TIMESTEP_KEY, conf);
    InitGeometry(conf);
    InitMass(conf);
    InitConstraint(conf);
}

void cBaseCloth::Reset()
{
    std::cout << "reset\n";
    SetPos(mClothInitPos);
    mXpre.noalias() = mClothInitPos;
    ClearForce();
}
void cBaseCloth::CalcTriangleDrawBuffer(Eigen::Map<tVectorXf> &res,
                                        int &st) const
{
    for (auto &tri : mTriangleArrayShared)
    {
        cRenderUtil::CalcTriangleDrawBufferSingle(this->mVertexArrayShared[tri->mId0],
                                                  this->mVertexArrayShared[tri->mId1],
                                                  this->mVertexArrayShared[tri->mId2], res, st);
    }
}

/**
 * \brief           calculate edge draw buffer
 */
void cBaseCloth::CalcEdgeDrawBuffer(Eigen::Map<tVectorXf> &res, int &st) const
{
    float alpha_channel = mVertexArrayShared[0]->mColor[3];

    tVector WARP_COLOR = ColorShoJoHi, WEFT_COLOR = ColorShiQing;
    tVector BLACK_COLOR = ColorAn;
    tVector normal = tVector::Zero();
    tVector cur_color = tVector::Zero();
    for (int idx = 0; idx < mEdgeArrayShared.size(); idx++)
    {
        auto e = mEdgeArrayShared[idx];
        // 1. get the averge normal direction for this edge
        normal = mTriangleArrayShared[e->mTriangleId0]->mNormal;
        if (e->mTriangleId1 != -1)
        {
            normal += mTriangleArrayShared[e->mTriangleId1]->mNormal;
            normal /= 2;
        }
        tVector2f uv_move =
            mVertexArrayShared[e->mId1]->muv - mVertexArrayShared[e->mId0]->muv;
        float thre = 1e-3;
        if (std::abs(uv_move[1]) < thre)
        {
            cur_color = WEFT_COLOR;
        }
        else if (std::abs(uv_move[0]) < thre)
        {
            cur_color = WARP_COLOR;
        }
        else if (std::fabs(std::fabs(uv_move[0]) - std::fabs(uv_move[1])) <
                 thre)
        {
            cur_color = BLACK_COLOR;
        }
        else
        {
            // std::cout << "v0 tex = " <<
            // mVertexArrayShared[e->mId0]->muv.transpose()
            //           << std::endl;
            // std::cout << "v1 tex = " <<
            // mVertexArrayShared[e->mId1]->muv.transpose()
            //           << std::endl;
            // std::cout << "uv move = " << uv_move.transpose() << std::endl;
            // std::cout << "edge id = " << idx << std::endl;
            // std::cout << "vertex 0 id = " << e->mId0 << std::endl;
            // std::cout << "vertex 1 id = " << e->mId1 << std::endl;
            // std::cout << "[warn] failed to determine the edge color\n";
        }
        // 2. calculate the edge draw position
        cRenderUtil::CalcEdgeDrawBufferSingle(mVertexArrayShared[e->mId0], mVertexArrayShared[e->mId1],
                                              normal, res, st, cur_color);
    }
}

void cBaseCloth::ClearForce()
{
    int dof = GetNumOfFreedom();
    mIntForce.noalias() = tVectorXd::Zero(dof);
    mUserForce.noalias() = tVectorXd::Zero(dof);
    mDampingForce.noalias() = tVectorXd::Zero(dof);
    mCollisionForce.noalias() = tVectorXd::Zero(dof);
}
#include "sim/Perturb.h"

void cBaseCloth::ApplyPerturb(tPerturb *pert)
{
    if (pert == nullptr)
        return;
    tVector force = pert->GetPerturbForce();
    mUserForce.segment(mTriangleArrayShared[pert->mAffectedTriId]->mId0 * 3, 3) +=
        force.segment(0, 3) / 3;
    mUserForce.segment(mTriangleArrayShared[pert->mAffectedTriId]->mId1 * 3, 3) +=
        force.segment(0, 3) / 3;
    mUserForce.segment(mTriangleArrayShared[pert->mAffectedTriId]->mId2 * 3, 3) +=
        force.segment(0, 3) / 3;
}
/**
 * \brief           add damping forces
 */
void cBaseCloth::CalcDampingForce(const tVectorXd &vel,
                                  tVectorXd &damping) const
{
    damping.noalias() = -vel * mDamping;
}

void cBaseCloth::SetPos(const tVectorXd &newpos)
{
    mXcur = newpos;
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < mVertexArrayShared.size(); i++)
    {
        mVertexArrayShared[i]->mPos.segment(0, 3).noalias() = mXcur.segment(i * 3, 3);
    }
}

// extern tVector gGravity;
const tVectorXd &cBaseCloth::GetPos() const { return this->mXcur; }
void cBaseCloth::CalcExtForce(tVectorXd &ext_force) const
{
    // 1. apply gravity
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < mVertexArrayShared.size(); i++)
    {
        ext_force.segment(3 * i, 3) +=
            mGravity * mVertexArrayShared[i]->mMass;
    }

    // std::cout << "add ext noise\n";
    // ext_force.segment(3 * (mVertexArrayShared.size() - 1), 3) += tVector3d(0, 0,
    // 10);

    //  2. add perturb force
    // if (mPerturb != nullptr)
    // {
    //     tVector perturb_force = mPerturb->GetPerturbForce();
    //     // printf(
    //     //     "[debug] perturb vid %d %d %d, ",
    //     //     mPerturb->mAffectedVerticesId[0],
    //     //     mPerturb->mAffectedVerticesId[1], mPerturb        // std::cout
    //     <<
    //     //     "perturb force = " << perturb_force.transpose()
    //     //           << std::endl;->mAffectedVerticesId[2]);

    //     ext_force.segment(mPerturb->mAffectedVerticesId[0] * 3, 3) +=
    //         perturb_force.segment(0, 3) / 3;
    //     ext_force.segment(mPerturb->mAffectedVerticesId[1] * 3, 3) +=
    //         perturb_force.segment(0, 3) / 3;
    //     ext_force.segment(mPerturb->mAffectedVerticesId[2] * 3, 3) +=
    //         perturb_force.segment(0, 3) / 3;
    //     // 2. give the ray to the perturb, calculate force on each vertices
    //     // 3. apply the force
    // }
}

std::vector<int>
FindVerticesIdFromUV(const tEigenArr<tVector2f> &uv_list,
                     const std::vector<tVertexPtr> &vertices_array)
{
    std::vector<int> vertices_result(uv_list.size(), -1);
    std::vector<double> SelectedConstraintStaticPointApproxDist(uv_list.size(),
                                                                std::nan(""));
    std::vector<tVector2f> result_uv(uv_list.size());
    for (int v_id = 0; v_id < vertices_array.size(); v_id++)
    {
        const tVector2f &v_uv = vertices_array[v_id]->muv;
        // std::cout << v_uv.transpose() << std::endl;
        for (int j = 0; j < uv_list.size(); j++)
        {
            double dist = (v_uv - uv_list[j]).norm();
            if (std::isnan(SelectedConstraintStaticPointApproxDist[j]) ||
                dist < SelectedConstraintStaticPointApproxDist[j])
            {
                result_uv[j] = v_uv;
                vertices_result[j] = v_id;
                SelectedConstraintStaticPointApproxDist[j] = dist;
            }
        }
    }
    // for (int i = 0; i < uv_list.size(); i++)
    // {

    //     std::cout << i << " [debug] find uv " << result_uv[i].transpose()
    //               << " for given " << uv_list[i].transpose() << std::endl;
    //     std::cout << "gap dist = " <<
    //     SelectedConstraintStaticPointApproxDist[i]
    //               << std::endl;
    //     std::cout << "vertex id = " << vertices_result[i] << std::endl;
    //     std::cout << "real id's uv = "
    //               << vertices_array[vertices_result[i]]->muv.transpose()
    //               << std::endl;
    // }
    return vertices_result;
}
void cBaseCloth::InitConstraint(const Json::Value &root)
{
    if (root.isMember("constraints") == false)
        return;
    auto cons = cJsonUtil::ParseAsValue("constraints", root);
    if (cons.isMember("static_point") == true)
    {
        // 1. read all 2d constraint for constraint_static point
        auto constraint_static_cons =
            cJsonUtil::ParseAsValue("static_point", cons);
        int num_of_constraint_static_pts = constraint_static_cons.size();
        tEigenArr<tVector2f> constraint_static_tex_coords(
            num_of_constraint_static_pts);
        for (int i = 0; i < num_of_constraint_static_pts; i++)
        {
            SIM_ASSERT(constraint_static_cons[i].size() == 2);
            constraint_static_tex_coords[i] =
                tVector2f(constraint_static_cons[i][0].asDouble(),
                          constraint_static_cons[i][1].asDouble());
        }

        // 2. iterate over all vertices to find which point should be finally
        // constraint_static
        mConstraint_StaticPointIds =
            FindVerticesIdFromUV(constraint_static_tex_coords, mVertexArrayShared);
        // output
        for (int i = 0; i < num_of_constraint_static_pts; i++)
        {
            printf("[debug] constraint_static uv (%.3f %.3f) selected v_id %d uv (%.3f, %.3f)\n",
                   constraint_static_tex_coords[i][0],
                   constraint_static_tex_coords[i][1],
                   mConstraint_StaticPointIds[i],
                   mVertexArrayShared[mConstraint_StaticPointIds[i]]->muv[0],
                   mVertexArrayShared[mConstraint_StaticPointIds[i]]->muv[1]);
        }
    }
    // for (auto &i : mConstraint_StaticPointIds)
    // {
    //     mInvMassMatrixDiag.segment(i * 3, 3).setZero();
    //     // printf("[debug] constraint_static point id %d at ", i);
    //     // exit(0);
    //     // next_pos.segment(i * 3, 3) = mXcur.segment(i * 3, 3);
    //     // std::cout << mXcur.segment(i * 3, 3).transpose() << std::endl;
    // }
}
#include "utils/ObjUtil.h"
void MoveObjPos(std::vector<tVertexPtr> &mVertexArrayShared)
{
    // 1. detect the biggest y
    float highest_y = -1e3;

    for (auto &x : mVertexArrayShared)
    {
        if (x->mPos[1] > highest_y)
        {
            highest_y = x->mPos[1];
        }
    }
    // 2. move this to given height
    double given_height = 0.27;
    double shift_value = given_height - highest_y;
    for (auto &x : mVertexArrayShared)
    {
        x->mPos[1] += shift_value;
    }
    double circle_height_threshold = 0.01; // 1cm
    tVector circle_origin = tVector::Zero();
    int counter = 0;
    for (auto &x : mVertexArrayShared)
    {
        if ((given_height - x->mPos[1]) < circle_height_threshold)
        {
            circle_origin = (circle_origin * counter + x->mPos) / (counter + 1);
            counter += 1;
        }
    }
    std::cout << "circle_origin = " << circle_origin.transpose() << std::endl;
    std::cout << "counter = " << counter << std::endl;

    for (auto &x : mVertexArrayShared)
    {
        x->mPos[0] -= circle_origin[0];
        x->mPos[2] -= circle_origin[2];
    }
    // 3. judge the circle around the most height pos

    // 4. judge the origin
    // 5. add the shift value
}

void cBaseCloth::InitGeometry(const Json::Value &conf)
{
    // 1. build the geometry
    bool enable_cloth_from_obj =
        cJsonUtil::ParseAsBool("enable_cloth_from_obj", conf);
    // std::cout << enable_cloth_from_obj << std::endl;
    // exit(1);
    if (enable_cloth_from_obj == true)
    {
        mClothObjPath = cJsonUtil::ParseAsString("cloth_obj_path", conf);
        bool enabel_relocate_obj =
            cJsonUtil::ParseAsBool("enable_relocate_obj", conf);
        std::cout << "enable cloth from obj from " << mClothObjPath
                  << std::endl;
        tVector3d cloth_obj_translation =
            cJsonUtil::ReadVectorJson(
                cJsonUtil::ParseAsValue("cloth_obj_translation", conf))
                .segment(0, 3);
        std::cout << "obj translation = " << cloth_obj_translation.transpose()
                  << std::endl;
        if (false == cFileUtil::ExistsFile(mClothObjPath))
        {
            SIM_ERROR("cloth obj path {} doesn't exist", mClothObjPath);
            exit(1);
        }
        cObjUtil::tParams params;
        params.mPath = mClothObjPath;
        cObjUtil::LoadObj(params, mVertexArrayShared, mEdgeArrayShared, mTriangleArrayShared);
        for (auto &x : mVertexArrayShared)
        {
            // x->mPos /= 1e3;
            x->mPos /= 1;
            x->mPos[1] += 0.2;
            // x->muv /= (10 / 3.0);
            x->muv *= (10 / 3.0);
            x->mNormal *= -1;
            x->mColor = ColorBlue;
            x->mPos += cMathUtil::Expand(cloth_obj_translation, 0);
        }
        for (auto &x : mEdgeArrayShared)
        {
            auto v0 = mVertexArrayShared[x->mId0];
            auto v1 = mVertexArrayShared[x->mId1];
            tVector cartesian_dist = v0->mPos - v1->mPos;
            tVector2f uv_dist = v0->muv - v1->muv;
            // std::cout << "cartesian dist = " << cartesian_dist.norm()
            //           << " uv_dist = " << uv_dist.norm() << std::endl;
        }
        // exit(1);
        if (enabel_relocate_obj == true)
            MoveObjPos(mVertexArrayShared);
        tVectorXd res = cJsonUtil::ReadVectorJson(
            cJsonUtil::ParseAsValue("cloth_obj_size", conf));
        SIM_ASSERT(res.size() == 2);
        SIM_DEBUG("Load obj as cloth from {}", cloth_obj_translation);
        mClothSizes = res.segment(0, 2);
        std::cout << "load obj done, num of ver = " << mVertexArrayShared.size()
                  << " num of edge " << mEdgeArrayShared.size()
                  << " num of triangles " << mTriangleArrayShared.size() << std::endl;
        tVector aabbmin, aabbmax;
        this->CalcAABB(aabbmin, aabbmax);
        std::cout << "aabbmin = " << aabbmin.transpose()
                  << " aabbmax = " << aabbmax.transpose() << std::endl;
    }
    else
    {
        // build from custom
        tVectorXd res = cJsonUtil::ReadVectorJson(
            cJsonUtil::ParseAsValue("cloth_size", conf));
        SIM_ASSERT(res.size() == 2);
        mClothSizes = res.segment(0, 2);

        cTriangulator::BuildGeometry(conf, mVertexArrayShared, mEdgeArrayShared,
                                     mTriangleArrayShared);
        for (int e_id = 0; e_id < mEdgeArrayShared.size(); e_id++)
        {
            auto &e = mEdgeArrayShared[e_id];
            if ((e->mId0 >= e->mId1) == true)
            {
                printf("[error] edge id %d, vid0 %d, vid1 %d\n", e_id, e->mId0,
                       e->mId1);
                exit(1);
            }
        }

        // mVertexArrayShared.clear();
        // mEdgeArrayShared.clear();
        // mTriangleArrayShared.clear();
        // {
        //     auto v0 = std::make_shared<tVertex>();
        //     v0->mColor = ColorBlue;
        //     v0->mPos = tVector(0, 0.1, 0, 1);
        //     v0->muv = tVector2f(0, 0);

        //     auto v1 = std::make_shared<tVertex>();
        //     v1->mColor = ColorBlue;
        //     v1->mPos = tVector(0.1, 0.1, 0, 1);
        //     v1->muv = tVector2f(0.1, 0);

        //     auto v2 = std::make_shared<tVertex>();
        //     v2->mColor = ColorBlue;
        //     v2->mPos = tVector(0.05, 0.2, 0, 1);
        //     v2->muv = tVector2f(0.05, 0.1);
        //     mVertexArrayShared = {v0, v1, v2};

        //     auto e0 = std::make_shared<tEdge>();
        //     e0->mId0 = 0;
        //     e0->mId1 = 1;
        //     e0->mTriangleId0 = 0;
        //     auto e1 = std::make_shared<tEdge>();
        //     e1->mId0 = 1;
        //     e1->mId1 = 2;
        //     e1->mTriangleId0 = 0;

        //     auto e2 = std::make_shared<tEdge>();
        //     e2->mId0 = 2;
        //     e2->mId1 = 0;
        //     e2->mTriangleId0 = 0;

        //     mEdgeArrayShared = {e0, e1, e2};
        //     auto t = std::make_shared<tTriangle>();
        //     t->mId0 = 0;
        //     t->mId1 = 1;
        //     t->mId2 = 2;

        //     mTriangleArrayShared = {t};
        // }
    }

    CalcNodePositionVector(mClothInitPos);

    // init the inv mass vector
    mMassMatrixDiag.noalias() = tVectorXd::Zero(GetNumOfFreedom());
    for (int i = 0; i < mVertexArrayShared.size(); i++)
    {
        mMassMatrixDiag.segment(i * 3, 3).fill(mVertexArrayShared[i]->mMass);
    }

    // update the normal information
    UpdateTriangleNormal();
    UpdateVertexNormalFromTriangleNormal();
    mXcur.noalias() = mClothInitPos;
    mXpre.noalias() = mClothInitPos;
}

void cBaseCloth::InitMass(const Json::Value &conf)
{
    // 1. get the mass
    mClothMass = cJsonUtil::ParseAsDouble("cloth_mass", conf);
    mClothDensity = cJsonUtil::ParseAsDouble("cloth_density", conf);

    // 2. get the density, confirm the correspondence
    if (std::fabs(mClothSizes.x() * mClothSizes.y() * mClothDensity -
                  mClothMass) > 1e-6)
    {
        std::cout << "cloth area = " << mClothSizes.x() * mClothSizes.y()
                  << std::endl;
        std::cout << "cloth density = " << mClothDensity << std::endl;
        std::cout << "cloth mass = " << mClothMass << std::endl;
        SIM_ERROR("mass {} != density * area {} mismatch", mClothMass,
                  mClothSizes.x() * mClothSizes.y() * mClothDensity);
        exit(1);
    }
    // 3. averge the mass to the vertices
    for (auto &x : mVertexArrayShared)
    {
        x->mMass = mClothMass / mVertexArrayShared.size();
    }
}
void cBaseCloth::CalcNodePositionVector(tVectorXd &pos) const
{
    if (pos.size() != GetNumOfFreedom())
    {
        pos.noalias() = tVectorXd::Zero(GetNumOfFreedom());
    }
    for (int i = 0; i < mVertexArrayShared.size(); i++)
    {
        pos.segment(i * 3, 3) = mVertexArrayShared[i]->mPos.segment(0, 3);
    }
}

/**
 * \brief       internal force
 */
void cBaseCloth::CalcIntForce(const tVectorXd &xcur, tVectorXd &int_force) const
{
    // std::vector<std::atomic<double>> int_force_atomic(int_force.size());
    // for (int i = 0; i < int_force.size(); i++)
    //     int_force_atomic[i] = 0;
    // double res = 1;
    // std::vector<double> int_force_atomic(int_force.size());

    // std::cout << "input fint = " << int_force.transpose() << std::endl;
    int id0, id1;
    double dist;
#ifdef USE_OPENMP
#pragma omp parallel for private(id0, id1, dist)
#endif
    for (int i = 0; i < mEdgeArrayShared.size(); i++)
    {
        const auto &spr = mEdgeArrayShared[i];
        // 1. calcualte internal force for each spring
        id0 = spr->mId0;
        id1 = spr->mId1;
        tVector3d pos0 = xcur.segment(id0 * 3, 3);
        tVector3d pos1 = xcur.segment(id1 * 3, 3);
        dist = (pos0 - pos1).norm();
        tVector3d force0 = spr->mK_spring * (spr->mRawLength - dist) *
                           (pos0 - pos1).segment(0, 3) / dist;
        // tVector3d force1 = -force0;
        // const tVectorXd &inf_force_0 = int_force.segment(3 * id0, 3);
        // const tVectorXd &inf_force_1 = int_force.segment(3 * id1, 3);
        //         std::cout << "spring " << i << " force = " <<
        //         force0.transpose() << ", dist " << dist << ", v0 " << id0 <<
        //         " v1 " << id1 << std::endl;
        // std::cout << "spring " << i << ", v0 = " << id0 << " v1 = " << id1 <<
        // std::endl;
        // 2. add force
        {
#ifdef USE_OPENMP
#pragma omp atomic
#endif
            int_force[3 * id0 + 0] += force0[0];
#ifdef USE_OPENMP
#pragma omp atomic
#endif
            int_force[3 * id0 + 1] += force0[1];
#ifdef USE_OPENMP
#pragma omp atomic
#endif
            int_force[3 * id0 + 2] += force0[2];
#ifdef USE_OPENMP
#pragma omp atomic
#endif
            int_force[3 * id1 + 0] += -force0[0];
#ifdef USE_OPENMP
#pragma omp atomic
#endif
            int_force[3 * id1 + 1] += -force0[1];
#ifdef USE_OPENMP
#pragma omp atomic
#endif
            int_force[3 * id1 + 2] += -force0[2];
        }
    }
    // std::cout << "output fint = " << int_force.transpose() << std::endl;
    // exit(0);
}

int cBaseCloth::GetNumOfFreedom() const { return 3 * GetNumOfVertices(); }

eClothType cBaseCloth::BuildClothType(std::string str)
{
    for (int i = 0; i < eClothType::NUM_OF_CLOTH_TYPE; i++)
    {
        if (gClothTypeStr[i] == str)
        {
            return static_cast<eClothType>(i);
        }
    }
    SIM_ERROR("unsupported cloth type {}", str);
    return eClothType::NUM_OF_CLOTH_TYPE;
}

void cBaseCloth::SetCollisionDetecter(cCollisionDetecterPtr ptr)
{
    mColDetecter = ptr;
}

/**
 * \brief           Get cloth shape: [height, width] = (Y_axis, X_axis)
 * For more details please check the triangulator
 */
tVector2d cBaseCloth::GetClothShape() const { return this->mClothSizes; }

void cBaseCloth::ClearConstraintStaticVertices()
{
    this->mConstraint_StaticPointIds.clear();
}
void cBaseCloth::AddConstraintStaticVertices(const std::vector<int> &vertices)
{
    for (auto &x : vertices)
    {
        mConstraint_StaticPointIds.push_back(x);
    }
}

const std::vector<int> &cBaseCloth::GetConstraintStaticVertices() const
{
    return mConstraint_StaticPointIds;
}

const tVectorXd &cBaseCloth::GetInitPos() const { return this->mClothInitPos; }

void cBaseCloth::MoveTranslation(const tVector3d &incremental_move)
{
    for (int i = 0; i < mVertexArrayShared.size(); i++)
    {
        mXcur.segment(3 * i, 3) += incremental_move;
    }
    SetPos(mXcur);
}

void cBaseCloth::ApplyUserPerturbForceOnce(tPerturb *)
{
}

void cBaseCloth::Update(float dt)
{
    this->UpdatePos(dt);
}
