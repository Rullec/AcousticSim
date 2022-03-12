#pragma once
#include "utils/DefUtil.h"
#include "utils/MathUtil.h"
#include "utils/SparseUtil.h"
template <typename data_type> using tMatrix32 = Eigen::Matrix<data_type, 3, 2>;
typedef tMatrix32<double> tMatrix32d;
typedef tMatrix32<float> tMatrix32f;
typedef Eigen::Matrix<double, 9, 1> tVector9d;
typedef Eigen::Matrix<float, 9, 1> tVector9f;
typedef Eigen::Matrix<double, 9, 9> tMatrix9d;
typedef Eigen::Matrix<double, 9, 2> tMatrix92d;
typedef Eigen::Matrix<float, 9, 2> tMatrix92f;
typedef Eigen::Matrix<float, 9, 9> tMatrix9f;
SIM_DECLARE_CLASS_AND_PTR(cBaseObject);

class cBaraffMaterial
{
public:
    explicit cBaraffMaterial();
    virtual void Init(cBaseObjectPtr object, const tVector3d &Kwarpweftbias);
    virtual void SetK(const tVector3d &Kwarpweftbias);
    virtual void
    Update(bool calc_energy = true, bool calc_fint = true,
           bool calc_stiffnessmatrix = true); // update deformation gra

    virtual double GetEnergy(int tri_id) const;
    virtual tVector9d GetForce(int tri_id) const;
    virtual tMatrix9d GetStiffMatrix(int tri_id) const;
    virtual double CalcTotalEnergy() const;
    virtual tVectorXd CalcTotalForce() const;
    virtual tSparseMatd CalcTotalStiffnessMatrix();
    virtual void CheckForce();
    virtual void CheckStiffnessMatrix();

    void GetFLst(std::vector<tMatrix32f> &vec);
    void GetnLst(std::vector<tMatrix32f> &vec);
    void GetgLst(std::vector<tMatrix92f> &vec);
    void GetCLst(std::vector<tVector2f> &vec);

    void GetFprimeLst(std::vector<tMatrix32f> &vec);
    void GetnprimeLst(std::vector<tMatrix32f> &vec);
    void GetgprimeLst(std::vector<tMatrix92f> &vec);
    void GetCprimeLst(std::vector<tVector2f> &vec);
    void GetEleKLst(std::vector<tMatrix9f> &mat);
    void GetEleFintLst(std::vector<tVector9f> &mat);

protected:
    cBaseObjectPtr mObject;
    int mNumOfVertices, mNumOfTriangles;
    tMatrix32d mS; // selection matrix from x to \delta x
    double mBu = 1.0, mBv = 1.0, mBshear = 1.0;
    tVector3d mKwarpweftshear; // warp and weft stiffness
    tEigenArr<tMatrix32d> mNLst, mFLst, mNprimeLst,
        mFprimeLst; // F = X * S * D_m^{-1}, N = S * D_m^{-1}
    tEigenArr<tMatrix32d> mnLst, mnprimeLst; // Fi / |Fi|, F colwise normalized
    tEigenArr<tMatrix92d> mgLst, mgprimeLst; // gi = Ni \otimes ni
    std::vector<double> mELst;               // energy
    tEigenArr<tVector2d> mCLst, mCprimeLst;
    tEigenArr<tVector9d> mIntForceLst;
    tEigenArr<tMatrix9d> mKLst;
    tSparseMatd global_K_buf;
    tEigenArr<tTriplet> total_triplets_buf = {};
    virtual void Allocate();
    virtual void InitN();
    virtual void UpdateFnC();

    // virtual tVector2d CalcCondition(const tMatrix32d &F) const;
    virtual void CalcFAndN(const tMatrix3d &pos, const tMatrix32d &uv_coords,
                           tMatrix32d &F, tMatrix32d &N) const;

    static tVector2f CalcC(tMatrix32d &F);
    static tMatrix92d Calcg(const tMatrix32d &N, const tMatrix32d &n);
    virtual tMatrix32d Calcn(const tMatrix3d &pos, tMatrix32d &N) const;
    static void CalcPi(const tMatrix32d &n, tMatrix3d &P0, tMatrix3d &P1);

    virtual void CalcFAndN_shearing(const tMatrix3d &pos,
                                    const tMatrix32d &uv_coords, tMatrix32d &F,
                                    tMatrix32d &N) const;

    virtual void CheckStretchForce();
    virtual void CheckStretchStiffnessMatrix();
    virtual double CalcStretchEnergy(const tMatrix3d &pos,
                                     const tMatrix32d &uv_coords) const;
    virtual tVector9d CalcStretchForce(const tMatrix3d &pos,
                                       const tMatrix32d &uv_coords) const;
    virtual tMatrix9d
    CalcStretchStiffMatrix(const tMatrix3d &pos,
                           const tMatrix32d &rest_texture_coords) const;

    // shearing
    virtual void CheckShearingForce();
    virtual void CheckShearingStiffnessMatrix();
    virtual double CalcShearingEnergy(const tMatrix3d &pos,
                                      const tMatrix32d &uv_coords) const;
    virtual tVector9d CalcShearingForce(const tMatrix3d &pos,
                                        const tMatrix32d &uv_coords) const;
    virtual tMatrix9d
    CalcShearingStiffMatrix(const tMatrix3d &pos,
                            const tMatrix32d &rest_texture_coords) const;
};