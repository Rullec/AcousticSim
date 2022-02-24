#ifdef _WIN32
#pragma once
#include "sim/cloth/BaseCloth.h"
#include "utils/DefUtil.h"

namespace StyleEngine
{
class SePiece;
class SeScene;
class SeDraggedPoints;
class SeFixedVertices;
}; // namespace StyleEngine

SIM_DECLARE_CLASS_AND_PTR(tPhyProperty);

struct tTargetDragPoint
{
    void Init(const Json::Value &conf);
    tVector2f mUV = tVector2f::Zero();
    tVector3f mTargetCartesian = tVector3f::Zero();
};

SIM_DECLARE_PTR(tTargetDragPoint);

class cLinctexCloth : public cBaseCloth
{
public:
    explicit cLinctexCloth(int id_);
    virtual ~cLinctexCloth();
    virtual void Init(const Json::Value &conf);
    virtual void Reset() override final;
    virtual void UpdatePos(double dt) override;
    virtual void SyncSimPropertyToInternal();
    virtual void SyncDensityToInternal();
    void SetSimProperty(const tPhyPropertyPtr &prop);
    tPhyPropertyPtr GetSimProperty() const;
    std::shared_ptr<StyleEngine::SePiece> GetPiece() const;
    const tVectorXd &GetClothFeatureVector() const;
    // tVector CalcCOM() const;
    int GetClothFeatureSize() const;

    // apply the noise
    void ApplyNoise(bool enable_y_random_rotation, double &rotation_angle,
                    bool enable_y_random_pos, const double random_ypos_std);
    void ApplyManualHoldNoise(const tVector2d &random_vector_xoz,
                              float one_side_gap_m, float bending_angle_rad);
    void ApplyMultiFoldsNoise(int num_of_folds, double max_amp);
    void RotateMaterialCoordsAfterReset(const tMatrix &init_mat_inv, int angle);
    void RotateMaterialCoords(int raw_angle, int new_angle);
    virtual void ClearConstraintStaticVertices() override;
    virtual void
    AddConstraintStaticVertices(const std::vector<int> &vertices) override;
    virtual void SetVerticesPos(const std::vector<int> &vertices,
                                const tEigenArr<tVector3d> &pos_array);
    virtual void SetPos(const tVectorXd &xcur) override final;
    bool IsNonlinearBendingEnabled() const;
    float &GetClothDensityRef();

protected:
    std::shared_ptr<StyleEngine::SePiece> mSeCloth;
    tPhyPropertyPtr mClothProp;   // cloth property
    bool mEnableDumpGeometryInfo; // if true, we save the geometry information
                                  // after the initialization
    std::string mDumpGeometryInfoPath; // save path for initial geometry
    // tVectorXd mClothFeature;
    std::shared_ptr<StyleEngine::SeFixedVertices>
        mLinctexConstraintStaticVertices;

    std::vector<tTargetDragPointPtr> mTargetDragPointsArray = {};
    std::vector<std::shared_ptr<StyleEngine::SeDraggedPoints>>
        mLinctexTargetDragPointsArray = {};
    std::vector<int> mTargetDragTriIdArray = {};
    float mThickness;
    virtual void InitConstraint(const Json::Value &root) override final;
    virtual void InitGeometry(const Json::Value &conf);
    void AddPiece(); // add the simulation data into the se engine

    void InitClothFeatureVector();
    void UpdateClothFeatureVector();
    void
    UpdateLinctexConstraintStaticVertices(const std::vector<int> &vertices);
    void InitLinctexConstraintTargetDragVertices();
    void UpdateLinctexConstraintTargetDragVertices();
    bool NeedToSyncProperty();
    bool NeedToSyncDensity();
    void CalcMaxCurvature() const;
};
#endif