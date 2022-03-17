#include "BaraffMaterial.h"
#include "geometries/Primitives.h"
#include "sim/BaseObject.h"
#include "utils/RotUtil.h"

cBaraffMaterial::cBaraffMaterial()
{
    mObject = nullptr;
    mKwarpweftbias.setZero();
    mCoefFu_warp_weft.clear();
    mCoefFv_warp_weft.clear();
    mCoefFu_diag.clear();
    mCoefFv_diag.clear();
    mNumOfVertices = 0;
    mNumOfTriangles = 0;
}

void CalcCoef(const tMatrix2d &mat, tVector3d &coef_u, tVector3d &coef_v)
{
    float a = mat(0, 0);
    float b = mat(0, 1);
    float c = mat(1, 0);
    float d = mat(1, 1);

    coef_u[0] = -a - c;
    coef_u[1] = a;
    coef_u[2] = c;
    coef_v[0] = -b - d;
    coef_v[1] = b;
    coef_v[2] = d;
}
void cBaraffMaterial::Init(cBaseObjectPtr object,
                           const tVector3d &Kwarpweftbias)
{
    mObject = object;
    mKwarpweftbias = Kwarpweftbias;

    mNumOfVertices = object->GetNumOfVertices();
    mNumOfTriangles = object->GetNumOfTriangles();

    // begin to init vector
    mCoefFu_warp_weft.resize(mNumOfTriangles, tVector3d::Zero());
    mCoefFv_warp_weft.resize(mNumOfTriangles, tVector3d::Zero());
    mCoefFu_diag.resize(mNumOfTriangles, tVector3d::Zero());
    mCoefFv_diag.resize(mNumOfTriangles, tVector3d::Zero());

    // calculate for each triangle
    auto tri_array = mObject->GetTriangleArray();
    auto v_array = mObject->GetVertexArray();
    tMatrix2d Dminv;
    tMatrix2d R45 = cRotUtil::RotMat2D(M_PI / 4);
    mTriangleAreaArray.resize(mNumOfTriangles);
    for (int i = 0; i < mNumOfTriangles; i++)
    {
        /*
        Dm^{-1} * [R45.T] = [
            a, b
            c, d
        ]
        */
        auto tri = tri_array[i];
        auto v0 = v_array[tri->mId0], v1 = v_array[tri->mId1],
             v2 = v_array[tri->mId2];
        Dminv.col(0) = (v1->muv - v0->muv).cast<double>();
        Dminv.col(1) = (v2->muv - v0->muv).cast<double>();
        Dminv = Dminv.inverse().eval();

        CalcCoef(Dminv, mCoefFu_warp_weft[i], mCoefFv_warp_weft[i]);
        Dminv = Dminv * R45.transpose();
        CalcCoef(Dminv, mCoefFu_diag[i], mCoefFv_diag[i]);
        mTriangleAreaArray[i] =
            cMathUtil::CalcTriangleArea(v0->mPos, v1->mPos, v2->mPos);
    }
}
void cBaraffMaterial::SetK(const tVector3d &Kwarpweftbias)
{
    this->mKwarpweftbias = Kwarpweftbias;
}

tVector3d CalcFuFv(const tVector3d &v0, const tVector3d &v1,
                   const tVector3d &v2, const tVector3d &coef)
{
    tVector3d F = tVector3d::Zero();
    F = v0 * coef[0] + v1 * coef[1] + v2 * coef[2];
    return F;
};
void CalcFuAndFv(const tVector3d &v0, const tVector3d &v1, const tVector3d &v2,
                 const tVector3d &Fu_coef, const tVector3d &Fv_coef,
                 tVector3d &Fu, tVector3d &Fv)
{
    Fu = Fu_coef[0] * v0 + Fu_coef[1] * v1 + Fu_coef[2] * v2;
    Fv = Fv_coef[0] * v0 + Fv_coef[1] * v1 + Fv_coef[2] * v2;
}
/*
E = S / 2 * (Ku * Cu^2 + Kv * Cv^2)
*/
double CalculateEnergy(double tri_area, double Ku, double Kv, double Cu,
                       double Cv)
{
    return tri_area / 2 * (Ku * Cu * Cu + Kv * Cv * Cv);
}

/*
Fint = -S * (Ku * Cu * ciu * Fu_normed + Kv * Cv * civ * Fv_normed)
    = -S * Ku * Cu * Fu_normed * ciu - S * Kv * Cv * Fv_normed * civ
    = u_base * ciu + v_base * civ


Fint_mat = [f0; f1; f2]
*/
tMatrix3d CalcFint(double tri_area, double Ku, double Kv, double Cu, double Cv,
                   const tVector3d &coef_u, const tVector3d &coef_v,
                   const tVector3d &Fu, const tVector3d &Fv)
{
    tVector3d u_base = -tri_area * Ku * Cu * Fu.normalized();
    tVector3d v_base = -tri_area * Kv * Cv * Fv.normalized();
    tMatrix3d fint = tMatrix3d::Zero();
    for (int i = 0; i < 3; i++)
    {
        fint.col(i) = u_base * coef_u[i] + v_base * coef_v[i];
    }
    return fint;
}

/*

H = -S * (Ku * ciu * cju * Dju + Kv * civ * cjv * Djv)

Du = 1/|Fu| * Fu_normed * Fu_normed^T
if |Fu| >=1, DJu += (1 - 1.0 / |Fu|) * I
*/
#include "utils/LogUtil.h"
tEigenArr<tMatrix3d> CalcHessian(double tri_area, double Ku, double Kv,
                                 const tVector3d &coef_u,
                                 const tVector3d &coef_v, const tVector3d &Fu,
                                 const tVector3d &Fv)
{
    // 1. calculate Du and Dv
    double Fu_norm = Fu.norm();
    double Fv_norm = Fu.norm();
    SIM_ASSERT(Fu_norm > 1e-6);
    SIM_ASSERT(Fv_norm > 1e-6);

    tVector3d Fu_normalized = Fu.normalized();
    tVector3d Fv_normalized = Fv.normalized();

    tMatrix3d Du = 1.0 / Fu_norm * Fu_normalized * Fu_normalized.transpose();
    tMatrix3d Dv = 1.0 / Fv_norm * Fv_normalized * Fv_normalized.transpose();

    if (Fu_norm > 1)
        Du += (1 - 1.0 / Fu_norm) * tMatrix3d::Identity();
    if (Fv_norm > 1)
        Dv += (1 - 1.0 / Fv_norm) * tMatrix3d::Identity();

    // 2. calculate array
    // H = -S * (Ku * ciu * cju * Du + Kv * civ * cjv * Dv)
    // row major storage in a list
    tEigenArr<tMatrix3d> array = {};
    for (int i = 0; i < 3; i++)
    {
        double ciu = coef_u[i], civ = coef_v[i];
        for (int j = 0; j < 3; j++)
        {
            double cju = coef_u[j], cjv = coef_v[j];
            array.push_back(-tri_area *
                            (Ku * ciu * cju * Du + Kv * civ * cjv * Dv));
        }
    }
    return array;
}
void cBaraffMaterial::Update()
{
    mTriEnergy.reserve(mNumOfTriangles);
    mTriEnergy.clear();
    mTriFint.reserve(mNumOfTriangles);
    mTriFint.clear();
    mTriHessian.reserve(mNumOfTriangles);
    mTriHessian.clear();

    auto v_array = mObject->GetVertexArray();
    auto tri_array = mObject->GetTriangleArray();
    OMP_PARALLEL_FOR
    for (int tri_id = 0; tri_id < mNumOfTriangles; tri_id++)
    {
        auto cur_tri = tri_array[tri_id];
        tVector3d v0 = v_array[cur_tri->mId0]->mPos.segment(0, 3);
        tVector3d v1 = v_array[cur_tri->mId1]->mPos.segment(0, 3);
        tVector3d v2 = v_array[cur_tri->mId2]->mPos.segment(0, 3);
        double tri_area = mTriangleAreaArray[tri_id];
        /* 1. update deformation gradient
                Fu, Fv for warp weft
                Fu, Fv for diag antidiag
        */
        // 2. update condition
        // 3. update energy
        // 4. update fint
        // 5. update stiffnessmatrix
        double E_total = 0;
        tMatrix3d Fint_total = tMatrix3d::Zero();
        tEigenArr<tMatrix3d> Hessian_total(
            9, tMatrix3d::Zero()); // dfdx = - d^2E /dx
        // warp weft direction
        {
            double Ku = mKwarpweftbias[0];
            double Kv = mKwarpweftbias[1];
            tVector3d coef_Fu = mCoefFu_warp_weft[tri_id],
                      coef_Fv = mCoefFv_warp_weft[tri_id];
            tVector3d Fu, Fv;
            CalcFuAndFv(v0, v1, v2, coef_Fu, coef_Fv, Fu, Fv);

            double Cu = Fu.norm() - 1;
            double Cv = Fv.norm() - 1;
            E_total += CalculateEnergy(tri_area, Ku, Kv, Cu, Cv);
            Fint_total +=
                CalcFint(tri_area, Ku, Kv, Cu, Cv, coef_Fu, coef_Fv, Fu, Fv);

            auto Hessian_stretch =
                CalcHessian(tri_area, Ku, Kv, coef_Fu, coef_Fv, Fu, Fv);
            for (int i = 0; i < 9; i++)
            {
                Hessian_total[i] += Hessian_stretch[i];
            }
        }
        // diag antidiag direction
        {
            double Ku = mKwarpweftbias[2];
            double Kv = mKwarpweftbias[2];
            tVector3d coef_Fu = mCoefFu_diag[tri_id],
                      coef_Fv = mCoefFv_diag[tri_id];
            tVector3d Fu, Fv;
            CalcFuAndFv(v0, v1, v2, coef_Fu, coef_Fv, Fu, Fv);

            double Cu = Fu.norm() - 1;
            double Cv = Fv.norm() - 1;
            E_total += CalculateEnergy(tri_area, Ku, Kv, Cu, Cv);
            Fint_total +=
                CalcFint(tri_area, Ku, Kv, Cu, Cv, coef_Fu, coef_Fv, Fu, Fv);

            auto Hessian_stretch =
                CalcHessian(tri_area, Ku, Kv, coef_Fu, coef_Fv, Fu, Fv);
            for (int i = 0; i < 9; i++)
            {
                Hessian_total[i] += Hessian_stretch[i];
            }
        }
#pragma omp critical
        mTriEnergy.push_back(E_total);
#pragma omp critical
        mTriFint.push_back(Fint_total);
#pragma omp critical
        mTriHessian.insert(mTriHessian.end(), Hessian_total.begin(),
                           Hessian_total.end());
    }
}

double cBaraffMaterial::CalcTotalEnergy() const
{
    double E_total = 0;
    for (auto &e : mTriEnergy)
    {
        E_total += e;
    }
    return E_total;
}
tVectorXd cBaraffMaterial::CalcTotalForce() const
{
    tVectorXd fint = tVectorXd::Zero(mNumOfVertices * 3);
    auto &tri_array = mObject->GetTriangleArray();
    for (int tri_id = 0; tri_id < mNumOfTriangles; tri_id++)
    {
        int v0 = tri_array[tri_id]->mId0;
        int v1 = tri_array[tri_id]->mId1;
        int v2 = tri_array[tri_id]->mId2;
        const tMatrix3d &Fint_mat = mTriFint[tri_id];
        fint.segment(3 * v0, 3) += Fint_mat.col(0);
        fint.segment(3 * v1, 3) += Fint_mat.col(1);
        fint.segment(3 * v2, 3) += Fint_mat.col(2);
    }
    return fint;
}
void Dispatch(tSparseMatd &sparse_mat, const tMatrix3d &mat, int v0_global,
              int v1_global)
{
    int dof = sparse_mat.outerSize();
    SIM_ASSERT(3 * v0_global < dof);
    SIM_ASSERT(3 * v1_global < dof);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
        {
            sparse_mat.coeffRef(3 * v0_global + i, 3 * v1_global + j) +=
                mat(i, j);
        }
}
tEigenArr<tTriplet> total_triplets_buf = {};
tSparseMatd cBaraffMaterial::CalcTotalStiffnessMatrix()
{
#define num_divide 15
    auto &tri_lst = mObject->GetTriangleArray();
    mHessian.resize(3 * mNumOfVertices, 3 * mNumOfVertices);
    mHessian.setZero();
    total_triplets_buf.clear();
    std::vector<tTriplet> sub_triples_set[num_divide];
    for (int i = 0; i < num_divide; i++)
    {
        sub_triples_set[i].reserve(3 * 3 * 9);
    }
#pragma omp parallel for num_threads(num_divide)
    for (int t = 0; t < mNumOfTriangles; t++)
    {
        int thread_num = omp_get_thread_num();
        auto &sub_triples = sub_triples_set[thread_num];
        sub_triples.clear();
        auto cur_tri = tri_lst[t];
        int v_id[3] = {cur_tri->mId0, cur_tri->mId1, cur_tri->mId2};

        for (size_t i = 0; i < 3; i++)
        {
            size_t global_vi_idx = v_id[i];
            for (size_t j = 0; j < 3; j++)
            {
                size_t global_vj_idx = v_id[j];

                tMatrix3d ele_K_part = mTriHessian[t * 9 + i * 3 + j];
                size_t i3 = 3 * global_vi_idx, j3 = 3 * global_vj_idx;
                sub_triples.emplace_back(i3, j3, ele_K_part(0, 0));
                sub_triples.emplace_back(i3, j3 + 1, ele_K_part(0, 1));
                sub_triples.emplace_back(i3, j3 + 2, ele_K_part(0, 2));

                sub_triples.emplace_back(i3 + 1, j3, ele_K_part(1, 0));
                sub_triples.emplace_back(i3 + 1, j3 + 1, ele_K_part(1, 1));
                sub_triples.emplace_back(i3 + 1, j3 + 2, ele_K_part(1, 2));

                sub_triples.emplace_back(i3 + 2, j3, ele_K_part(2, 0));
                sub_triples.emplace_back(i3 + 2, j3 + 1, ele_K_part(2, 1));
                sub_triples.emplace_back(i3 + 2, j3 + 2, ele_K_part(2, 2));
            }
        }
#pragma omp critical
        total_triplets_buf.insert(total_triplets_buf.end(), sub_triples.begin(),
                                  sub_triples.end());
    }
    mHessian.setFromTriplets(total_triplets_buf.begin(),
                             total_triplets_buf.end());
    return mHessian;
}