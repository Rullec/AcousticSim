#ifdef _WIN32
#include "LinctexCloth.h"
#include "SePhysicalProperties.h"
#include "SePiece.h"
#include "SeScene.h"
#include "SeSceneOptions.h"
#include "SeSimParameters.h"
#include "SeSimulationProperties.h"
#include "geometries/Triangulator.h"
#include "sim/cloth/ClothProperty.h"
#include "utils/FileUtil.h"
#include "utils/JsonUtil.h"
#include <SeFeatureVertices.h>
#include <algorithm>
#include <iostream>

void tTargetDragPoint::Init(const Json::Value &conf)
{
    mUV = cJsonUtil::ReadVectorJson(cJsonUtil::ParseAsValue("uv", conf))
              .segment(0, 2)
              .cast<float>();
    std::cout << "[target_drag] uv = " << mUV.transpose() << std::endl;
    mTargetCartesian =
        cJsonUtil::ReadVectorJson(
            cJsonUtil::ParseAsValue("target_cartesian_pos", conf))
            .segment(0, 3)
            .cast<float>();
    std::cout << "[target_drag] target cartesian = "
              << mTargetCartesian.transpose() << std::endl;
}

SE_USING_NAMESPACE

std::string strip(std::string old_str)
{

    while (old_str[0] == ' ' || old_str[0] == '\n' || old_str[0] == '\t')
    {
        old_str = old_str.substr(1, old_str.size() - 1);
    }
    while (old_str[old_str.size() - 1] == ' ' ||
           old_str[old_str.size() - 1] == '\n' ||
           old_str[old_str.size() - 1] == '\t')
    {
        old_str = old_str.substr(0, old_str.size() - 1);
    }
    return old_str;
}

cLinctexCloth ::cLinctexCloth(int id_)
    : cBaseCloth(eClothType::LINCTEX_CLOTH, id_)
{
    SIM_DEBUG(
        "linctex build {}, build name {}, stretch scale {}, bending scale {}",
        LINCTEX_BUILD_VERSION, LINCTEX_BUILD_NAME, STRETCH_SCALE,
        BENDING_SCALE);
}
cLinctexCloth ::~cLinctexCloth() {}
void cLinctexCloth ::Init(const Json::Value &conf)
{
    mEnableDumpGeometryInfo =
        cJsonUtil::ParseAsBool("enable_dump_geometry_info", conf);
    if (mEnableDumpGeometryInfo)
    {
        mDumpGeometryInfoPath =
            cJsonUtil::ParseAsString("dump_triangle_info_path", conf);
    }
    mThickness = cJsonUtil::ParseAsFloat("cloth_thickness", conf);
    // if (mThickness > 1e-6)
    // {
    //     // SIM_ERROR("invalid thickness {}", mThickness);
    //     // exit(1);
    // }
    cBaseCloth::Init(conf);

    // 1. load the physical property
    mClothProp = std::make_shared<tPhyProperty>();
    mClothProp->Init(conf);
    AddPiece();
    UpdateLinctexConstraintStaticVertices(mConstraint_StaticPointIds);
    InitLinctexConstraintTargetDragVertices();

    InitClothFeatureVector();
}
void cLinctexCloth::UpdatePos(double dt)
{
    auto &pos = mSeCloth->FetchPositions();
    for (int i = 0; i < mVertexArrayShared.size(); i++)
    {
        mVertexArrayShared[i]->mPos.noalias() =
            tVector(pos[i][0], pos[i][1], pos[i][2], 1);
    }

    UpdateClothFeatureVector();
    UpdateTriangleNormal();
    UpdateVertexNormalFromTriangleNormal();
    UpdateLinctexConstraintTargetDragVertices();
    // CalcMaxCurvature();
    SyncSimPropertyToInternal();
    SyncDensityToInternal();
    // std::cout << mSeCloth->GetPhysicalProperties()->GetBendingBias()
    //           << std::endl;
    // std::cout << "stretch = "
    //           << mSeCloth->GetPhysicalProperties()->GetStretchWarp()
    //           << std::endl;
}

void cLinctexCloth::InitConstraint(const Json::Value &root)
{
    // we must rotate the fixed point's uv as well...
    if (root.isMember("constraints") == false)
        return;
    Json::Value new_root = root;
    Json::Value cons_conf = cJsonUtil::ParseAsValue("constraints", new_root);
    if (cons_conf.isMember("static_point"))
    {
        float angle = cJsonUtil::ParseAsDouble("cloth_uv_rotation", new_root) /
                      180 * M_PI;
        // angle = 0;
        // std::cout << "raw_angle = " << angle << std::endl;
        tMatrix2f rotmat = cRotUtil::RotMat2D(angle).cast<float>();
        tMatrix2f rotmat_inv = rotmat.inverse();
        // std::cout << "rotmat = \n" << rotmat << std::endl;
        // std::cout << "rotmat inv = \n" << rotmat_inv << std::endl;

        Json::Value fixed_cons_conf =
            cJsonUtil::ParseAsValue("static_point", cons_conf);

        int num_of_fixed_pts = fixed_cons_conf.size();

        for (int i = 0; i < num_of_fixed_pts; i++)
        {
            SIM_ASSERT(fixed_cons_conf[i].size() == 2);
            tVector2f raw_fixed = tVector2f(fixed_cons_conf[i][0].asDouble(),
                                            fixed_cons_conf[i][1].asDouble());
            // tVector2f new_fixed = rotmat_inv * raw_fixed;
            tVector2f new_fixed = rotmat * raw_fixed;
            // std::cout << "raw fixed = " << raw_fixed.transpose() << " "
            //           << "new fixed = " << new_fixed.transpose() <<
            //           std::endl;
            fixed_cons_conf[i][0] = new_fixed[0];
            fixed_cons_conf[i][1] = new_fixed[1];
        }
        cons_conf["static_point"] = fixed_cons_conf;
        new_root["constraints"] = cons_conf;
    }
    // 2. create t target drag

    if (true == cons_conf.isMember("target_drag_point"))
    {
        for (auto &v : cJsonUtil::ParseAsValue("target_drag_point", cons_conf))
        {
            auto pt = std::make_shared<tTargetDragPoint>();
            pt->Init(v);
            mTargetDragPointsArray.push_back(pt);
        }
    }
    // std::cout << "target drag points array = " <<
    // mTargetDragPointsArray.size()
    //           << std::endl;
    // exit(1);

    cBaseCloth::InitConstraint(new_root);
}

/**
 * \brief               Init the geometry and set the init positions
 */
void cLinctexCloth::InitGeometry(const Json::Value &conf)
{
    cBaseCloth::InitGeometry(conf);
    std::string init_geo =
        cJsonUtil::ParseAsString("cloth_init_nodal_position", conf);
    bool is_illegal = false;
    if (cFileUtil::ExistsFile(init_geo) == true)
    {
        Json::Value value;
        cJsonUtil::LoadJson(init_geo, value);
        tVectorXd vec =
            cJsonUtil::ReadVectorJson(cJsonUtil::ParseAsValue("input", value));
        // std::cout << "vec size = " << vec.size();
        if (vec.size() == GetNumOfFreedom())
        {
            mClothInitPos = vec;
            mXcur = vec;
            SetPos(mXcur);
            std::cout << "[debug] set init vec from " << init_geo << std::endl;
        }
        else
        {
            SIM_WARN("nodal {} freedom {} != {}", init_geo, vec.size(),
                     GetNumOfFreedom());
            is_illegal = true;
        }
    }
    else
    {
        SIM_WARN("nodal {} doesn't exist", init_geo);
        is_illegal = true;
    }

    // 1. check whehter it's illegal
    if (is_illegal)
    {
        std::cout << "[warn] init geometry file " << init_geo
                  << " is illegal, ignore\n";
    }

    if (mEnableDumpGeometryInfo == true)
    {
        cTriangulator::SaveGeometry(this->mVertexArray, this->mEdgeArray,
                                    this->mTriangleArray,
                                    this->mDumpGeometryInfoPath);
    }
}

void cLinctexCloth::AddPiece()
{
    /*
SePiecePtr SePiece::Create(const std::vector<Int3> & triangles,
                    const std::vector<Float3> & positions,
                    const std::vector<Float2> & materialCoords,
                    SePhysicalPropertiesPtr pPhysicalProperties)
    */

    std::vector<Int3> indices(0);
    std::vector<Float3> pos3D(0);
    std::vector<Float2> pos2D(0);
    for (int i = 0; i < mVertexArrayShared.size(); i++)
    {
        auto v = mVertexArrayShared[i];
        pos3D.push_back(Float3(v->mPos[0], v->mPos[1], v->mPos[2]));
        pos2D.push_back(Float2(v->muv[0], v->muv[1]));
    }

    for (int i = 0; i < mTriangleArrayShared.size(); i++)
    {
        auto tri = mTriangleArrayShared[i];
        indices.push_back(Int3(tri->mId0, tri->mId1, tri->mId2));
        // triangle_array_se[i] = ;
        // printf("[debug] triangle %d: vertices: %d, %d, %d\n", i, tri->mId0,
        // tri->mId1, tri->mId2);
    }

    auto phyProp = SePhysicalProperties::Create();
    phyProp->SetMassDensity(mClothDensity);
    // std::cout << "linctex mass density = " << phyProp->GetMassDensity() <<
    // std::endl; std::cout << "cur mass density = " << mClothDensity <<
    // std::endl; SIM_ASSERT((std::fabs(phyProp->GetMassDensity() -
    // mClothDensity) < 1e-6) ==
    //            true);

    // std::cout << "mass density = " << phyProp->GetMassDensity() << std::endl;
    mSeCloth = SePiece::Create(indices, pos3D, pos2D, phyProp);
    SetSimProperty(mClothProp);

    SeSimulationPropertiesPtr cloth_sim_prop = mSeCloth->GetSimulationProperties();
    // 碰撞厚度
    cloth_sim_prop->SetCollisionThickness(mThickness);

    // 面料厚度
    mSeCloth->GetPhysicalProperties()->SetThickness(mThickness);
    // std::cout << "[debug] se col thickness = "
    //           << cloth_sim_prop->GetCollisionThickness() << std::endl;
}

/**
 *
 * \brief           let M = a K^2 + b K,
 *
 * dMdk = 2 a K + b \ge 0
 * a >= - b / (2 K)
 */
float clamp_bending_2nd(float term_1st, float term_2nd)
{
    if (term_1st < 0)
    {
        SIM_ERROR("term 1st {} <0, illegal", term_1st);
        exit(1);
    }
    else
    {
        float max_K = 200; // d\theta / dl
        float fenmu = -2 * max_K;
        float min_term_2nd = term_1st / fenmu;
        if (term_2nd < min_term_2nd)
        {
            // std::cout << "fenmu = " << fenmu << std::endl;
            // std::cout << "for linear = " << term_1st
            //           << " the threshold = " << min_term_2nd << std::endl;
            SIM_WARN("nonlinear 2nd term {} < min {} when 1st term {}, cliped",
                     term_2nd, min_term_2nd, term_1st);
            std::cout << "-b / (2a) = " << -term_1st / (2 * term_2nd)
                      << std::endl;
        }

        term_2nd = std::max(term_2nd, min_term_2nd);
        return term_2nd;
    }
}

/**
 * \brief           Given the external sim property, calculate the result
 */
void cLinctexCloth::SetSimProperty(const tPhyPropertyPtr &prop)
{
    mClothProp = prop;
    auto phyProp = mSeCloth->GetPhysicalProperties();

    // SIM_ERROR("cannot be set directly");
    // exit(0);
    auto want_to_set_simulation_stretch_warp_value =
        tPhyProperty::ConvertStretchCoefFromGUIToSim(
            mClothProp->GetStretchWarp());

    // set stretch values
    phyProp->SetStretchWarp(want_to_set_simulation_stretch_warp_value);

    phyProp->SetStretchWeft(tPhyProperty::ConvertStretchCoefFromGUIToSim(
        mClothProp->GetStretchWeft()));
    phyProp->SetShearing(tPhyProperty::ConvertStretchCoefFromGUIToSim(
        mClothProp->GetStretchBias()));
    // set linear bending values
    float bending_warp_1st = tPhyProperty::ConvertLinearBendingCoefFromGUIToSim(
        mClothProp->GetLinearBendingWarp());
    // std::cout << "[debug] linear bending warp gui = "
    //           << mClothProp->GetLinearBendingWarp()
    //           << " sim = " << bending_warp_1st << std::endl;
    phyProp->SetBendingWarp(bending_warp_1st);
    float bending_weft_1st = tPhyProperty::ConvertLinearBendingCoefFromGUIToSim(
        mClothProp->GetLinearBendingWeft());
    phyProp->SetBendingWeft(bending_weft_1st);
    float bending_bias_1st = tPhyProperty::ConvertLinearBendingCoefFromGUIToSim(
        mClothProp->GetLinearBendingBias());
    phyProp->SetBendingBias(bending_bias_1st);

    // set nonlinear bending values
#ifdef SE_HAVE_NONLINEAR_PROPERTIES
    float bending_warp_2nd =
        tPhyProperty::ConvertNonLinearBendingCoefFromGUIToSim(
            mClothProp->GetNonLinearBendingWarp());
    // bending_warp_2nd = clamp_bending_2nd(bending_warp_1st, bending_warp_2nd);
    // std::cout << "[debug] nonlinear bending warp gui = "
    //           << mClothProp->GetNonLinearBendingWarp()
    //           << " sim = " << bending_warp_2nd << std::endl;

    mClothProp->SetNonLinearBendingWarp(
        tPhyProperty::ConvertNonLinearBendingCoefFromSimToGUI(
            bending_warp_2nd));

    phyProp->SetBending2ndWarp(bending_warp_2nd);
    float bending_weft_2nd =
        tPhyProperty::ConvertNonLinearBendingCoefFromGUIToSim(
            mClothProp->GetNonLinearBendingWeft());

    // bending_weft_2nd = clamp_bending_2nd(bending_weft_1st, bending_weft_2nd);
    mClothProp->SetNonLinearBendingWeft(
        tPhyProperty::ConvertNonLinearBendingCoefFromSimToGUI(
            bending_weft_2nd));

    phyProp->SetBending2ndWeft(bending_weft_2nd);

    float bending_bias_2nd =
        tPhyProperty::ConvertNonLinearBendingCoefFromGUIToSim(
            mClothProp->GetNonLinearBendingBias());
    // bending_bias_2nd = clamp_bending_2nd(bending_bias_1st, bending_bias_2nd);
    mClothProp->SetNonLinearBendingBias(
        tPhyProperty::ConvertNonLinearBendingCoefFromSimToGUI(
            bending_bias_2nd));
    phyProp->SetBending2ndBias(bending_bias_2nd);
#endif
}

tPhyPropertyPtr cLinctexCloth::GetSimProperty() const { return mClothProp; }

/**
 * \brief           Get the feature vector of this cloth
 *
 *  Current it's all nodal position of current time
 */
const tVectorXd &cLinctexCloth::GetClothFeatureVector() const { return mXcur; }

int cLinctexCloth::GetClothFeatureSize() const { return mXcur.size(); }

/**
 * \brief               Init the feature vector
 */
void cLinctexCloth::InitClothFeatureVector()
{
    mXcur.noalias() = tVectorXd::Zero(mVertexArrayShared.size() * 3);
    UpdateClothFeatureVector();
    // std::cout << "mXcur size = " << mXcur.size() << std::endl;
}

/**
 * \brief               Calculate the feature vector of the cloth
 */
void cLinctexCloth::UpdateClothFeatureVector()
{
    for (int i = 0; i < mVertexArrayShared.size(); i++)
    {
        mXcur.segment(3 * i, 3).noalias() = mVertexArrayShared[i]->mPos.segment(0, 3);
    }
}

/**
 * \brief               Update nodal position from a vector
 */
void cLinctexCloth::SetPos(const tVectorXd &xcur)
{
    cBaseCloth::SetPos(xcur);
    if (mSeCloth)
    {
        std::vector<Float3> pos(0);
        for (int i = 0; i < mVertexArrayShared.size(); i++)
        {
            pos.push_back(
                Float3(xcur[3 * i + 0], xcur[3 * i + 1], xcur[3 * i + 2]));
        }
        mSeCloth->SetPositions(pos);
    }
}
bool cLinctexCloth::IsNonlinearBendingEnabled() const
{
#ifdef SE_HAVE_NONLINEAR_PROPERTIES
    return true;
#else
    return false;
#endif
}
void cLinctexCloth::Reset()
{
    cBaseCloth::Reset();
    UpdateClothFeatureVector();
}

std::shared_ptr<StyleEngine::SePiece> cLinctexCloth::GetPiece() const
{
    return this->mSeCloth;
}

#include <cmath>
#include <math.h>

/**
 * \brief           Apply random gaussian noise, point-wise
 */
void cLinctexCloth::ApplyNoise(bool enable_y_random_rotation,
                               double &rotation_angle, bool enable_y_random_pos,
                               const double random_ypos_std)
{
    rotation_angle = 0;
    if (enable_y_random_rotation == true)
    {
        rotation_angle = cMathUtil::RandDouble(0, 2 * M_PI);
        // std::cout << "rotation angle = " << rotation_angle << std::endl;
    }

    tMatrix mat = cRotUtil::EulerAnglesToRotMat(
        tVector(0, rotation_angle, 0, 0), eRotationOrder::XYZ);
    for (int i = 0; i < mVertexArrayShared.size(); i++)
    {
        // 1. apply rotation
        tVector cur_pos = tVector::Ones();
        cur_pos.segment(0, 3).noalias() = mXcur.segment(i * 3, 3);
        mXcur.segment(i * 3, 3).noalias() = (mat * cur_pos).segment(0, 3);

        // 2. apply noise translation
        if (enable_y_random_pos == true)
        {
            mXcur[i * 3 + 1] += cMathUtil::RandDoubleNorm(0, random_ypos_std);
        }

        // if (enable_y_random_pos)
        // {
        // }
        // exit(0);
    }

    // 3. final update
    SetPos(mXcur);
}

/**
 * \brief          Apply manual-hold noise
 *          When the hands are holding clothes, there will be a folding axis
 * going through the center of cloth
 *          uv-based method
 */
#include "geometries/ObjExport.h"
void cLinctexCloth::ApplyManualHoldNoise(const tVector2d &random_vector_xoz_,
                                         float one_side_gap_m,
                                         float bending_angle_rad)
{
    {
        tVector aabb_min, aabb_max;
        CalcAABB(aabb_min, aabb_max);
        double y_gap = aabb_max[1] - aabb_min[1];
        if (y_gap > 0.05) // y> 5cm
        {
            SIM_WARN("manual hold noise: Y gap for AABB is {} >= 5cm, fail to "
                     "apply the noise",
                     y_gap);
            return;
        }
    }

    // 1. random an axis on XOZ plane
    tVector2d principal_axis_2d = random_vector_xoz_.normalized();
    // principal_axis_2d[0] = std::fabs(principal_axis_2d[0]);
    tVector3d principal_axis_3d =
        tVector3d(principal_axis_2d[0], 0, principal_axis_2d[1]);
    if (principal_axis_2d.hasNaN() == true)
    {
        principal_axis_2d = tVector2d::Random().normalized();
        SIM_WARN("random set principal_axis_2d to {}",
                 principal_axis_2d.transpose());
    }

    int num_of_vertices = mXcur.size() / 3;

    // 2. get the highest Y axis
    tVector min_aabb, max_aabb;
    CalcAABB(min_aabb, max_aabb);
    double height_highest = max_aabb[1];
    tVector mid_point = (min_aabb + max_aabb) / 2;

    // 2.1 calculate the plane
    tVector2d principal_plane_normal_2d =
        cRotUtil::RotMat2D(M_PI / 2) * principal_axis_2d;
    tVector3d principal_plane_normal_3d = tVector3d(
        principal_plane_normal_2d[0], 0, principal_plane_normal_2d[1]);

    // 2.2 construct the plane
    tVector plane_equation = tVector::Zero();
    plane_equation.segment(0, 3) = principal_plane_normal_3d;
    plane_equation[3] = -principal_plane_normal_3d.dot(mid_point.segment(0, 3));

    SIM_ASSERT(cMathUtil::CalcPlanePointDist(plane_equation,
                                             mid_point.segment(0, 3)) < 1e-10);

    // 3. set the left - right gap = 6cm  = 0.06 m
    // 3.1 clamp gap
    one_side_gap_m = std::fabs(one_side_gap_m);
    double one_side_gap_min = 0.06; // m
    double one_side_gap_max = 0.1;  // m
    double one_side_gap =
        cMathUtil::Clamp(one_side_gap_m, one_side_gap_min, one_side_gap_max);
    // 3.2 clamp value
    bending_angle_rad = std::fabs(bending_angle_rad);
    bending_angle_rad =
        cMathUtil::Clamp(bending_angle_rad, M_PI / 4.0, M_PI / 2.0);
    // printf("[log] apply manual-fold noise on axis %.2f %.2f, angle %.2f, gap
    // "
    //        "%.2f m\n",
    //        principal_axis_2d[0], principal_axis_2d[1], bending_angle_rad,
    //        one_side_gap);
    tMatrix bending_mat_negative = cRotUtil::AxisAngleToRotmat(
        -cMathUtil::Expand(principal_axis_3d, 0) * bending_angle_rad);
    tMatrix bending_mat_positive = cRotUtil::AxisAngleToRotmat(
        cMathUtil::Expand(principal_axis_3d, 0) * bending_angle_rad);

    // 4. fold the vertices
    for (int i = 0; i < num_of_vertices; i++)
    {
        tVector3d pos = mXcur.segment(3 * i, 3);
        // 4.1 calculate distance
        double dist_to_plane =
            cMathUtil::CalcPlanePointDist(plane_equation, pos);
        // give plane
        if (dist_to_plane < one_side_gap)
        {
            pos[1] = height_highest;
        }
        else
        {
            // std::cout << ""
            // give bending
            tVector mid_to_pos = cMathUtil::Expand(pos, 1) - mid_point;
            mid_to_pos[3] = 0;
            tVector mid_to_pos_clamped_at_edge =
                one_side_gap / dist_to_plane * mid_to_pos;
            tVector left_vec = mid_to_pos - mid_to_pos_clamped_at_edge;
            left_vec[3] = 0;
            tVector new_left_vec = bending_mat_positive * left_vec;
            if (new_left_vec[1] > 0)
            {
                new_left_vec = bending_mat_negative * left_vec;
            }

            // new left vec should perpendicular to principal normal
            float should_be_zero =
                new_left_vec.segment(0, 3).dot(principal_plane_normal_3d);
            // if (should_be_zero > 1e-6)
            // {
            //     std::cout << "should be zero = " << should_be_zero <<
            //     std::endl;
            // }

            tVector new_pos =
                new_left_vec + mid_to_pos_clamped_at_edge + mid_point;
            // double dist_to_edge = (dist_to_plane - one_side_gap);
            // double new_highest = height_highest - ;
            pos = new_pos.segment(0, 3);
            // exit(1);
        }

        mXcur.segment(3 * i, 3) = pos;
    }
    SetPos(mXcur);

    // cObjExporter::ExportObj("after_manual_noise.obj", mVertexArray,
    //                         mTriangleArray);
}

/**
 * \brief                   Apply multiple folds noise
 * \param num_of_folds      Given number of folds
 */
void cLinctexCloth::ApplyMultiFoldsNoise(int num_of_folds, double max_amp)
{
    SIM_ASSERT(num_of_folds >= 2 && num_of_folds <= 10);
    // 1. calculate the fold cycle (theta)
    double theta = 2 * M_PI / num_of_folds;

    // 2. calculate the fold direction, random the amptitude
    double st_bias = 0;

    // lying on the XOZ plane
    tEigenArr<tVector2d> fold_directions_array(0);
    std::vector<double> fold_st_angle_array(0);
    std::vector<double> fold_amp_array(0);
    for (int i = 0; i < num_of_folds; i++)
    {
        // 2 * 1 vector
        double angle = theta * i + st_bias;
        angle = cMathUtil::NormalizeAngle(angle);
        // std::cout << "angle " << i << " = " << angle << std::endl;
        // double amp = cMathUtil::RandDouble(0, 0.1); // up to 10 cm amp
        // double amp = 0.1;
        double amp = cMathUtil::RandDouble(0, max_amp);
        // double amp = cMathUtil::RandDoubleNorm(0.1, 0.1);
        tVector fold_dir =
            cRotUtil::AxisAngleToRotmat(tVector(0, 1, 0, 0) * angle) *
            tVector(1, 0, 0, 0);
        // printf("[debug] angle %d = %.3f, dir = ", i, angle);
        fold_st_angle_array.push_back(angle);
        fold_amp_array.push_back(amp);
        // project to XOZ plane
        fold_directions_array.push_back(tVector2d(fold_dir[0], fold_dir[2]));
        // printf("[info] angle %d = %.3f, amp = %.3f\n", i, angle, amp);
        // std::cout << fold_directions_array[fold_directions_array.size() -
        // 1].transpose() << std::endl;
    }

    /*
        3. calculate the noise for each point
            3.1 for each point, calculate the pole coordinate
            3.2 confirm two fold direction
            3.3 calculate the height field for these 2 fold. (clamped cos
       function) 3.4 averaging these 2 values, apply this height
    */
    auto calc_angle_distance = [](double first, double second)
    {
        first = cMathUtil::NormalizeAngle(first);
        second = cMathUtil::NormalizeAngle(second);
        double dist = std::fabs(first - second);
        if (dist > M_PI)
            return 2 * M_PI - dist;
        else
            return dist;
    };

    std::vector<int> times(num_of_folds, 0);
    tVector com = CalcCOM();
    // int idx = 0;
    for (int v_id = 0; v_id < mVertexArrayShared.size(); v_id++)
    {
        auto &v = mVertexArrayShared[v_id];
        // project the nodal vector to XOZ plane, calculate the "angle"
        tVector node_vec = v->mPos - com;
        // double theta = 2 * M_PI / mVertexArrayShared.size() * (idx++);
        // tVector node_vec =
        //     tVector(
        //         std::cos(theta),
        //         0,
        //         std::sin(theta), 0);
        node_vec[1] = 0;
        node_vec.normalize();
        // std::cout << tVector(1, 0, 0, 0).cross3(tVector(0, 0, -1, 0)) <<
        // std::endl; exit(0);
        tVector res = cMathUtil::CalcAxisAngleFromOneVectorToAnother(
            tVector(1, 0, 0, 0), node_vec);
        if (res[1] < 0)
        {
            res = tVector(0, 2 * M_PI + res[1], 0, 0);
        }
        // {
        //     tVector residual = tVector(0, 1, 0, 0) * res.norm() - res;
        //     SIM_ASSERT(residual.norm() < 1e-6);
        //     if (residual.norm() > 1e-6)
        //     {
        //         std::cout << "node_vec = " << node_vec.transpose() <<
        //         std::endl; std::cout << "res = " << res.transpose() <<
        //         std::endl; exit(1);
        //     }
        // }
        double cur_angle = cMathUtil::NormalizeAngle(res[1]);
        // double cur_angle = 3;
        // double cur_angle = 1.24081;
        // std::cout << "for point " << node_vec.transpose() << " its angle = "
        // << cur_angle << std::endl;
        int interval0 = -1, interval1 = -1;
        for (int i = 0; i < num_of_folds; i++)
        {
            // for fold 1
            int fold0_id = i, fold1_id = (i + 1) % num_of_folds;
            double angle0 = fold_st_angle_array[fold0_id],
                   angle1 = fold_st_angle_array[fold1_id];

            // if the value is on the boundary, include them
            if (std::fabs(angle0 - cur_angle) < 1e-3 ||
                std::fabs(angle1 - cur_angle) < 1e-3)
            {
                interval0 = fold0_id;
                interval1 = fold1_id;
                break;
            }
            else
            {
                double interval = calc_angle_distance(angle0, angle1);
                if (calc_angle_distance(angle0, cur_angle) < interval &&
                    calc_angle_distance(angle1, cur_angle) < interval)
                {
                    interval0 = fold0_id;
                    interval1 = fold1_id;
                    break;
                }
            }
        }

        if (interval0 == -1)
            interval0 = 0;
        if (interval1 == -1)
            interval1 = 0;
        // {
        //     // std::cout << "[error] for angle " << cur_angle
        //     //           << " failed to judge the interval. cur interval
        //     are:";
        //     // for (auto &x : fold_st_angle_array)
        //     //     std::cout << x << " ";
        //     // std::cout << std::endl;
        //     interval0 = 0;
        //     interval1 = 0;
        //     // exit(0);
        // }
        // else
        {
            times[interval0] += 1;

            double amp0 = fold_amp_array[interval0],
                   amp1 = fold_amp_array[interval1];
            double int_angle0 = fold_st_angle_array[interval0],
                   int_angle1 = fold_st_angle_array[interval1];
            double angle_with0 = calc_angle_distance(cur_angle, int_angle0);
            double angle_with1 = calc_angle_distance(cur_angle, int_angle1);
            double bias = 0;
            if (angle_with0 < theta / 2)
            {
                bias += std::cos(angle_with0 / (theta / 2) * M_PI) * amp0;
            }
            else
            {
                bias += -amp0 * std::pow(angle_with1 / theta, 2);
            }
            if (angle_with1 < theta / 2)
            {
                bias += std::cos(angle_with1 / (theta / 2) * M_PI) * amp1;
            }
            else
            {
                bias += -amp1 * std::pow(angle_with0 / theta, 2);
            }

            // remove stretch
            double raw_length = (v->mPos - com).segment(0, 3).norm();
            double cloth_width = (mClothSizes[0] + mClothSizes[1]) / 2;
            double max_length = cloth_width / 2 * std::sqrt(2);
            bias *= std::pow(raw_length / max_length, 3);
            mXcur[3 * v_id + 1] += bias;
            // tVector3d ref = mXcur.segment(3 * (v_id - 1), 3);
            // tVector3d cur = mXcur.segment(3 * (v_id), 3);

            // mXcur.segment(3 * v_id, 3) = (cur - ref).normalized() *
            // (v->mPos.segment(0, 3) - ref).norm() + ref.segment(0, 3);
            // v->mPos[1] ;
            // printf("[info] angle %.4f is in [%.3f, %.3f]\n", cur_angle,
            // fold_st_angle_array[interval0],
            //        fold_st_angle_array[interval1]);
        }
    }
    SetPos(mXcur);
    // for (auto &x : times)
    // {
    //     std::cout << x << std::endl;
    // }
    // exit(0);
}

void cLinctexCloth::RotateMaterialCoordsAfterReset(const tMatrix &init_mat_inv,
                                                   int angle)
{
    cTriangulator::RotateMaterialCoordsAfterReset(init_mat_inv, mVertexArray,
                                                  angle);
    auto props = mSeCloth->GetSimulationProperties();
    props->SetWeftWarpAngle(float(angle) / 180 * M_PI);
}

void cLinctexCloth::RotateMaterialCoords(int raw_angle, int new_angle)
{
    cTriangulator::RotateMaterialCoords(raw_angle, new_angle, mVertexArrayShared);
    auto props = mSeCloth->GetSimulationProperties();
    props->SetWeftWarpAngle(float(new_angle) / 180 * M_PI);
}

void cLinctexCloth::SyncSimPropertyToInternal()
{
    if (false == NeedToSyncProperty())
    {
        // we don't need to sync
        return;
    }
    else
    {
        // SIM_WARN(
        //     "physical prop is different from the internal value in linctex "
        //     "dll:\nformal bending warp {} != internal bending warp {}\nformal
        //     " "bending weft {} != internal bending weft {}\nformal bending
        //     bias {} "
        //     "!= internal bending bias {}\nformal stretch warp {} != internal
        //     " "stretch warp {}\nformal stretch weft {} != internal stretch
        //     weft {}\n" "formal stretch bias {} != internal stretch bias {} ",
        //     formal_bending_warp, internal_bending_warp, formal_bending_weft,
        //     internal_bending_weft, formal_bending_bias,
        //     internal_bending_bias,

        //     formal_stretch_warp, internal_stretch_warp, formal_stretch_weft,
        //     internal_stretch_weft, formal_stretch_bias,
        //     internal_stretch_bias);
        // std::cout << "set sim prop = "
        //           << mClothProp->BuildFullFeatureVector().transpose()
        //           << std::endl;
        SetSimProperty(mClothProp);
    }
}

void cLinctexCloth::ClearConstraintStaticVertices()
{
    cBaseCloth::ClearConstraintStaticVertices();
    mSeCloth->Remove(mLinctexConstraintStaticVertices);
}
void cLinctexCloth::AddConstraintStaticVertices(
    const std::vector<int> &vertices)
{
    cBaseCloth::AddConstraintStaticVertices(vertices);
    this->UpdateLinctexConstraintStaticVertices(mConstraint_StaticPointIds);
}

void cLinctexCloth::UpdateLinctexConstraintStaticVertices(
    const std::vector<int> &vertices)
{
    mLinctexConstraintStaticVertices = mSeCloth->AddFixedVertices(vertices);
}

void cLinctexCloth::SetVerticesPos(const std::vector<int> &vertices,
                                   const tEigenArr<tVector3d> &pos_array)
{
    SIM_ASSERT(vertices.size() == pos_array.size());
    std::vector<StyleEngine::Float3> se_array(0);
    for (int i = 0; i < vertices.size(); i++)
    {
        mXcur.segment(3 * i, 3) = pos_array[i];
        se_array.push_back(StyleEngine::Float3(pos_array[i][0], pos_array[i][1],
                                               pos_array[i][2]));
    }
    mSeCloth->SetPositionsByIndices(vertices, se_array);
    cBaseCloth::SetPos(mXcur);
}

std::vector<int>
FindTrianglesIdFromUV(const tEigenArr<tVector2f> &uv_list,
                      const std::vector<tVertexPtr > &vertices_array,
                      const std::vector<tTriangle *> &tri_array)
{
    std::vector<int> tri_id_result(uv_list.size(), -1);
    std::vector<double> SelectedConstraintStaticPointApproxDist(uv_list.size(),
                                                                std::nan(""));
    std::vector<tVector2f> result_uv(uv_list.size());
    for (int t_id = 0; t_id < tri_array.size(); t_id++)
    {
        tVector2f v_uv = (vertices_array[tri_array[t_id]->mId0]->muv +
                          vertices_array[tri_array[t_id]->mId1]->muv +
                          vertices_array[tri_array[t_id]->mId2]->muv) /
                         3;
        // std::cout << v_uv.transpose() << std::endl;
        for (int j = 0; j < uv_list.size(); j++)
        {
            double dist = (v_uv - uv_list[j]).norm();
            if (std::isnan(SelectedConstraintStaticPointApproxDist[j]) ||
                dist < SelectedConstraintStaticPointApproxDist[j])
            {
                result_uv[j] = v_uv;
                tri_id_result[j] = t_id;
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
    //     std::cout << "triangle id = " << tri_id_result[i] << std::endl;
    //     std::cout << "real id's uv = "
    //               << vertices_array[tri_array[tri_id_result[i]]->mId0]
    //                      ->muv.transpose()
    //               << std::endl;
    // }
    return tri_id_result;
}

tVector3d GetTriangleMidPoint(tTriangle *tri,
                              const std::vector<tVertexPtr > &v_array)
{
    return (v_array[tri->mId0]->mPos + v_array[tri->mId1]->mPos +
            v_array[tri->mId2]->mPos)
               .segment(0, 3) /
           3;
}
/**
 * \brief           set the constraint target drag vertices
 */
tEigenArr<tVector3f> init_pos_array;
void cLinctexCloth::InitLinctexConstraintTargetDragVertices()
{
    tVector2f min_uv = tVector2f::Ones() * 1e10;
    tVector2f max_uv = -tVector2f::Ones() * 1e10;
    for (auto &v : mVertexArrayShared)
    {
        if (v->muv[0] > max_uv[0])
            max_uv[0] = v->muv[0];
        if (v->muv[1] > max_uv[1])
            max_uv[1] = v->muv[1];

        if (v->muv[0] < min_uv[0])
            min_uv[0] = v->muv[0];
        if (v->muv[1] < min_uv[1])
            min_uv[1] = v->muv[1];
    }
    // std::cout << min_uv.transpose() << std::endl;
    // std::cout << max_uv.transpose() << std::endl;
    init_pos_array.clear();
    int i = 0;
    for (auto &v : mTargetDragPointsArray)
    {

        // std::cout << "--------begin to add target drag pt---------\n";
        StyleEngine::TriangleBaryCoord bary;
        bary.coord = Vec3f(0.5, 0.5, 0.5);
        // 1. begin to find the vertices
        bary.index =
            FindTrianglesIdFromUV({v->mUV}, mVertexArray, mTriangleArray)[0];
        mTargetDragTriIdArray.push_back(bary.index);
        // std::cout << "set bary index = " << bary.index << std::endl;
        // std::cout << "v 10098 uv = " << mVertexArrayShared[10098]->muv.transpose()
        //           << std::endl;
        // std::cout << "v 10200 uv = " << mVertexArrayShared[10200]->muv.transpose()
        //           << std::endl;
        SIM_ASSERT((bary.index >= 0) && (bary.index < mTriangleArrayShared.size()));
        tVector3f mid_pt =
            GetTriangleMidPoint(mTriangleArrayShared[bary.index], mVertexArrayShared)
                .cast<float>();
        init_pos_array.push_back(mid_pt);
        std::shared_ptr<StyleEngine::SeDraggedPoints> drag_pt =
            GetPiece()->AddDraggedPoint(
                bary, StyleEngine::Float3(mid_pt[0], mid_pt[1], mid_pt[2]));
        // std::cout << "init = " << mid_pt.transpose() << std::endl;
        mLinctexTargetDragPointsArray.push_back(drag_pt);
        i++;
    }
}
int global_fid = 0;
/**
 * \brief           update dragged vertices
 */
void cLinctexCloth::UpdateLinctexConstraintTargetDragVertices()
{
    // for(auto &v : mLinctexTargetDragPointsArray)
    for (int i = 0; i < mLinctexTargetDragPointsArray.size(); i++)
    {

        tVector3f cur_pos = init_pos_array[i];
        tVector3f tar_pos = mTargetDragPointsArray[i]->mTargetCartesian;
        float total_f = 3e2;

        tVector3f new_tar =
            global_fid / 3e2 * tar_pos + (3e2 - global_fid) / 3e2 * cur_pos;
        // tVector3f vec = tar_pos - cur_pos;
        // for (int i = 0; i < 3; i++)
        // {
        //     float new_val = std::fabs(vec[i]) > 1e-3 ? 1e-3 :
        //     std::fabs(vec[i]); vec[i] = new_val * (vec[i] > 0 ? 1 : -1);
        // }
        // tVector3f tar_pos_now = cur_pos + vec;
        // std::cout << "-----\n";
        // std::cout << "cur pos = " << cur_pos.transpose() << std::endl;
        // std::cout << "tar pos = " << tar_pos.transpose() << std::endl;
        // // std::cout << "move vec = " << vec.transpose() << std::endl;
        // std::cout << "new tar pos = " << new_tar.transpose() << std::endl;
        mLinctexTargetDragPointsArray[i]->SetPositions(
            {StyleEngine::Float3(new_tar[0], new_tar[1], new_tar[2])});
    }
    if (global_fid < 3e2)
        global_fid += 1;
}
/**
 * \brief           judge whether do we need to sync the property into the
 * internal simulation property
 */
bool cLinctexCloth::NeedToSyncProperty()
{
    auto judge_float_equal = [](float a, float b) -> bool
    { return std::fabs(a - b) < 1e-6; };
    // 1. get current sim property from world
    auto phyProp = mSeCloth->GetPhysicalProperties();
    float internal_linear_bending_warp =
              tPhyProperty::ConvertLinearBendingCoefFromSimToGUI(
                  phyProp->GetBendingWarp()),
          internal_linear_bending_weft =
              tPhyProperty::ConvertLinearBendingCoefFromSimToGUI(
                  phyProp->GetBendingWeft()),
          internal_linear_bending_bias =
              tPhyProperty::ConvertLinearBendingCoefFromSimToGUI(
                  phyProp->GetBendingBias());

    float internal_stretch_warp = tPhyProperty::ConvertStretchCoefFromSimToGUI(
              phyProp->GetStretchWarp()),
          internal_stretch_weft = tPhyProperty::ConvertStretchCoefFromSimToGUI(
              phyProp->GetStretchWeft()),
          internal_stretch_bias = tPhyProperty::ConvertStretchCoefFromSimToGUI(
              phyProp->GetShearing());

    float formal_linear_bending_warp = mClothProp->GetLinearBendingWarp(),
          formal_linear_bending_weft = mClothProp->GetLinearBendingWeft(),
          formal_linear_bending_bias = mClothProp->GetLinearBendingBias();

    float formal_stretch_warp = mClothProp->GetStretchWarp(),
          formal_stretch_weft = mClothProp->GetStretchWeft(),
          formal_stretch_bias = mClothProp->GetStretchBias();

    bool need_update = false;

    // if stretch is not synced, set need_update
    if (false ==
            judge_float_equal(internal_stretch_warp, formal_stretch_warp) ||
        false ==
            judge_float_equal(internal_stretch_weft, formal_stretch_weft) ||
        false == judge_float_equal(internal_stretch_bias, formal_stretch_bias))
    {
        need_update = true;
    }

    // if linear is not synced, set need_update
    if (false == judge_float_equal(internal_linear_bending_warp,
                                   formal_linear_bending_warp) ||
        false == judge_float_equal(internal_linear_bending_weft,
                                   formal_linear_bending_weft) ||
        false == judge_float_equal(internal_linear_bending_bias,
                                   formal_linear_bending_bias))
    {
        need_update = true;
    }

#ifdef SE_HAVE_NONLINEAR_PROPERTIES

    // 1. get formal nonlinear property
    float formal_nonlinear_bending_bias = mClothProp->GetNonLinearBendingBias();
    float formal_nonlinear_bending_weft = mClothProp->GetNonLinearBendingWeft();
    float formal_nonlinear_bending_warp = mClothProp->GetNonLinearBendingWarp();

    // 2. get internal nonlinear property
    float internal_nonlinear_bending_bias =
        tPhyProperty::ConvertNonLinearBendingCoefFromSimToGUI(
            phyProp->GetBending2ndBias());
    float internal_nonlinear_bending_warp =
        tPhyProperty::ConvertNonLinearBendingCoefFromSimToGUI(
            phyProp->GetBending2ndWarp());
    float internal_nonlinear_bending_weft =
        tPhyProperty::ConvertNonLinearBendingCoefFromSimToGUI(
            phyProp->GetBending2ndWeft());

    if (

        false == judge_float_equal(formal_nonlinear_bending_bias,
                                   internal_nonlinear_bending_bias) ||
        false == judge_float_equal(formal_nonlinear_bending_warp,
                                   internal_nonlinear_bending_warp) ||
        false == judge_float_equal(formal_nonlinear_bending_weft,
                                   internal_nonlinear_bending_weft)

    )
    {
        // printf(
        //     "[nonlinear_judge_failed] formal nl bias %f, internal nl bias
        //     %f\n", formal_nonlinear_bending_bias,
        //     internal_nonlinear_bending_bias);
        // printf(
        //     "[nonlinear_judge_failed] formal nl warp %f, internal nl warp
        //     %f\n", formal_nonlinear_bending_warp,
        //     formal_nonlinear_bending_warp);
        // printf(
        //     "[nonlinear_judge_failed] formal nl weft %f, internal nl weft
        //     %f\n", formal_nonlinear_bending_weft,
        //     formal_nonlinear_bending_weft);
        need_update = true;
    }
#endif

    return need_update;
}

float &cLinctexCloth::GetClothDensityRef() { return this->mClothDensity; }

/**
 * \brief           do we need to update the density
 */
bool cLinctexCloth::NeedToSyncDensity()
{
    float internal_density =
        GetPiece()->GetPhysicalProperties()->GetMassDensity();
    return std::fabs(internal_density - this->mClothDensity) > 1e-6;
}

/**
 * \brief
 */
void cLinctexCloth::SyncDensityToInternal()
{
    if (NeedToSyncDensity() == true)
    {
        auto prop = GetPiece()->GetPhysicalProperties();
        prop->SetMassDensity(mClothDensity);
    }
}

/**
 * \brief           calculate max curvature
 */

void cLinctexCloth::CalcMaxCurvature() const
{
    float max_K = -1;
    float max_K_color = 200;
    for (auto &cur_e : mEdgeArray)
    {
        if (cur_e->mTriangleId0 == -1 || cur_e->mTriangleId1 == -1)
            continue;
        auto tri0 = mTriangleArrayShared[cur_e->mTriangleId0];
        auto tri1 = mTriangleArrayShared[cur_e->mTriangleId1];
        // 1. calc normal0 and normal1
        tVector normal0 =
            (mVertexArrayShared[tri0->mId1]->mPos - mVertexArrayShared[tri0->mId0]->mPos)
                .cross3(mVertexArrayShared[tri0->mId2]->mPos -
                        mVertexArrayShared[tri0->mId1]->mPos)
                .normalized();

        tVector normal1 =
            (mVertexArrayShared[tri1->mId1]->mPos - mVertexArrayShared[tri1->mId0]->mPos)
                .cross3(mVertexArrayShared[tri1->mId2]->mPos -
                        mVertexArrayShared[tri1->mId1]->mPos)
                .normalized();
        float theta_rad = std::acos(normal0.dot(normal1));
        float edge_length =
            (mVertexArrayShared[cur_e->mId0]->mPos - mVertexArrayShared[cur_e->mId1]->mPos)
                .norm();
        float area0 = cMathUtil::CalcTriangleArea(
            mVertexArrayShared[tri0->mId0]->mPos, mVertexArrayShared[tri0->mId1]->mPos,
            mVertexArrayShared[tri0->mId2]->mPos);
        float area1 = cMathUtil::CalcTriangleArea(
            mVertexArrayShared[tri1->mId0]->mPos, mVertexArrayShared[tri1->mId1]->mPos,
            mVertexArrayShared[tri1->mId2]->mPos);
        float height0 = 2 * area0 / edge_length;
        float height1 = 2 * area1 / edge_length;
        float curvature = 2 * theta_rad / (height0 + height1);
        max_K = std::max(curvature, max_K);

        float hot = curvature / max_K_color;
        tVector cur_color =
            tVector(1, 0, 0, 1) * hot + tVector(0, 1, 0, 1) * (1 - hot);

        cur_color[3] = GetVertexColorAlpha();
        {
            mVertexArrayShared[cur_e->mId0]->mColor = cur_color;
            mVertexArrayShared[cur_e->mId1]->mColor = cur_color;
        }
    }
    std::cout << "max K = " << max_K << std::endl;
}
#endif
