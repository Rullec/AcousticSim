#include "DihedralBending.h"
#include "geometries/Primitives.h"
#include <iostream>

double CalcTriAreaFromEdgeLenght(double a, double b, double c)
{
    double p = (a + b + c) / 2;
    double S = std::sqrt(p * (p - a) * (p - b) * (p - c));
    return S;
}
tVector12d CalcDBetaDx(const tVector3d &v0, const tVector3d &v1,
                       const tVector3d &v2, const tVector3d &v3);
void VerifyDBetaDx();
void VerifyDBetaDnormal();
cDihedralMaterial::cDihedralMaterial()
{
    // VerifyDBetaDx();
    // VerifyDBetaDnormal();
    // exit(1);
}

void cDihedralMaterial::Init(const std::vector<tVertexPtr> &v_array,
                             const std::vector<tEdgePtr> &e_array,
                             const std::vector<tTrianglePtr> &t_array,
                             const tVector3d &bending_stiffness_warpweftbias)
{
    cBaseMaterial::Init(v_array, e_array, t_array,
                        bending_stiffness_warpweftbias);

    int numOfVertices = v_array.size();
    int numOfEdges = e_array.size();
    int numOfTris = t_array.size();

    // 1. raw edge length
    mRawEdgeLengthArray.resize(numOfEdges);
    mEleForceArray.resize(numOfEdges);
    mGlobalForce.resize(3 * numOfVertices);
    for (int i = 0; i < numOfEdges; i++)
    {
        auto edge = e_array[i];
        int v0 = edge->mId0;
        int v1 = edge->mId1;
        mRawEdgeLengthArray[i] = (v_array[v0]->mPos - v_array[v1]->mPos).norm();
        std::cout << "raw edge length " << i << " = " << mRawEdgeLengthArray[i]
                  << std::endl;
    }

    // 2. raw triangle area
    mRawTriangleAreaArray.resize(numOfTris);
    for (int i = 0; i < numOfTris; i++)
    {
        double e0_length = mRawEdgeLengthArray[t_array[i]->mEId0];
        double e1_length = mRawEdgeLengthArray[t_array[i]->mEId1];
        double e2_length = mRawEdgeLengthArray[t_array[i]->mEId2];
        mRawTriangleAreaArray[i] =
            CalcTriAreaFromEdgeLenght(e0_length, e1_length, e2_length);
        std::cout << "raw tri area " << i << " = " << mRawTriangleAreaArray[i]
                  << std::endl;
    }

    // 3. raw height (for inner triangles)
    mRawHeightArray.resize(numOfEdges);
    tVector2f height = tVector2f::Zero();
    for (int i = 0; i < numOfEdges; i++)
    {
        auto e = mEdgeArray[i];
        if (e->mIsBoundary == false)
        {
            // 3.1 get triangle id (left and right)
            int t0 = e->mTriangleId0;
            int t1 = e->mTriangleId1;
            // 3.2 get edge length and triangle area
            double a0 = mRawTriangleAreaArray[t0];
            double a1 = mRawTriangleAreaArray[t1];
            double e_length = mRawEdgeLengthArray[i];
            // 3.3 get height. S = 0.5 * e * h; h = 2 * S / e
            height = 2 * tVector2f(a0, a1) / e_length;
            std::cout << "edge height " << i << " =  " << height.transpose()
                      << std::endl;
        }
        else
        {
            height.setZero();
        }
        mRawHeightArray[i] = height;
    }
}

double cDihedralMaterial::CalcEnergy(const tVectorXd &xcur) { return 0; }

tVectorXd cDihedralMaterial::CalcForce(const tVectorXd &xcur)
{
    std::cout << "[f] force size = " << mGlobalForce.size() << std::endl;
    return this->mGlobalForce;
}

void cDihedralMaterial::CheckForce() {}

void cDihedralMaterial::CheckStiffnessMatrix() {}

/*
    beta = arctan2(y, x) - \pi
    n0 = (x1 - x0) \times (x2 - x1)
    n1 = (x0 - x1) \times (x3 - x0)

    \bar n = n / |n|

    I'_vec = (I3 - bar_v * bar_v^T) / |v|

    d(beta)/dn0 = - sign(sin \beta) / sqrt(1 - x*x) * I'_{n0}^T * \bar n_1
    d(beta)/dn1 = - sign(sin \beta) / sqrt(1 - x*x) * I'_{n1}^T * \bar n_0
*/
tMatrix3d CalcIprime(const tVector3d &vec)
{
    tMatrix3d res = tMatrix3d::Identity();
    return (res - vec.normalized() * vec.normalized().transpose()) / vec.norm();
}
double CalcBeta(const tVector3d &n0_bar, const tVector3d &n1_bar,
                const tVector3d &e_bar)
{
    /*
        x = n0_bar . n1_bar
        y = e_bar . (n0_bar \times n1_bar)
        beta = arctan2(y, x) - pi
    */
    double x = n0_bar.dot(n1_bar);
    double y = e_bar.dot(n0_bar.cross(n1_bar));
    double beta = std::atan2(y, x);
    return beta;
}

void CalcNormalAndEdge(const tVector3d &v0, const tVector3d &v1,
                       const tVector3d &v2, const tVector3d &v3, tVector3d &n0,
                       tVector3d &n1, tVector3d &e)
{
    n0 = (v1 - v0).cross(v2 - v1);
    n1 = (v0 - v1).cross(v3 - v0);
    e = (v1 - v0);
}

double CalcBeta(const tVector3d &v0, const tVector3d &v1, const tVector3d &v2,
                const tVector3d &v3)
{
    tVector3d n0, n1, e;
    CalcNormalAndEdge(v0, v1, v2, v3, n0, n1, e);
    return CalcBeta(n0.normalized(), n1.normalized(), e.normalized());
}

void CalcDBetaDnormal_and_De(const tVector3d &n0, const tVector3d &n1,
                             const tVector3d &e, tVector3d &dbeta_dn0,
                             tVector3d &dbeta_dn1, tVector3d &dbeta_de)
{
    // 1. prepare
    tVector3d n0_bar = n0.normalized();
    tVector3d n1_bar = n1.normalized();
    tVector3d e_bar = e.normalized();
    double x = n0_bar.dot(n1_bar);
    double y = e_bar.dot(n0_bar.cross(n1_bar));
    double beta = CalcBeta(n0_bar, n1_bar, e_bar);
    // std::cout << "beta = " << beta << std::endl;
    double sign_sinb = std::sin(beta) > 0 ? 1.0 : -1.0;
    double sign_cosb = std::cos(beta) > 0 ? 1.0 : -1.0;
    // 2. I'
    tMatrix3d Iprime0 = CalcIprime(n0), Iprime1 = CalcIprime(n1);
    double deno = std::sqrt(1 - std::min(x * x, 1.0));
    if (deno < 1e-10)
    {
        deno = 1e-10;
    }
    double factor = -sign_sinb * 1.0 / deno;
    dbeta_dn0 = factor * Iprime0.transpose() * n1.normalized();
    dbeta_dn1 = factor * Iprime1.transpose() * n0.normalized();
    dbeta_de = sign_cosb / std::sqrt(1 - y * y) * CalcIprime(e).transpose() *
               (n0_bar.cross(n1_bar));
}

tVector12d CalcDBetaDx(const tVector3d &v0, const tVector3d &v1,
                       const tVector3d &v2, const tVector3d &v3)
{
    // 1. calculate beta
    tVector3d n0, n1, e;
    CalcNormalAndEdge(v0, v1, v2, v3, n0, n1, e);
    tVector3d n0_bar = n0.normalized(), n1_bar = n1.normalized(),
              e_bar = e.normalized();
    double beta = CalcBeta(n0_bar, n1_bar, e_bar);

    // 2. calculate dbeta/dn0, dbeta/dn1
    tVector3d dbeta_dn0, dbeta_dn1, dbeta_de;
    CalcDBetaDnormal_and_De(n0, n1, e, dbeta_dn0, dbeta_dn1, dbeta_de);
    // std::cout << "dbeta_dn0 = " << dbeta_dn0.transpose() << std::endl;
    // std::cout << "dbeta_dn1 = " << dbeta_dn1.transpose() << std::endl;
    /*
        3. calcualte
            dbeta / dv0 =
                (v1 - v2) \times dbdn0
                +
                (v1 - v3) \times dbdn1

            dbeta / dv1 =
                (v2 - v0) \times dbdn0
                +
                (v0 - v3) \times dbdn1

            dbeta / dv2 =
                (v0 - v1) \times dbdn0

            dbeta / dv3 =
                (v1 - v0) \times dbdn1
    */
    // std::cout << "dbeta_de = " << dbeta_de.transpose() << std::endl;
    tVector3d dbeta_dv0 =
        (v1 - v2).cross(dbeta_dn0) + (v3 - v1).cross(dbeta_dn1);
    tVector3d dbeta_dv1 =
        (v2 - v0).cross(dbeta_dn0) + (v0 - v3).cross(dbeta_dn1);
    tVector3d dbeta_dv2 = (v0 - v1).cross(dbeta_dn0);
    tVector3d dbeta_dv3 = (v1 - v0).cross(dbeta_dn1);
    tVector12d res(12);
    res.segment(0, 3) = dbeta_dv0;
    res.segment(3, 3) = dbeta_dv1;
    res.segment(6, 3) = dbeta_dv2;
    res.segment(9, 3) = dbeta_dv3;
    return res;
}
void VerifyDBetaDnormal()
{
    for (int i = 0; i < 10; i++)
    {
        tVector3d n0 = tVector3d::Random(), n1 = tVector3d::Random();
        tVector3d e = n0.cross(n1);

        // tVector3d n0 = tVector3d(-0.11283083, -0.00337882, -0.03312034);
        // tVector3d n1 = tVector3d(0.17173063, 0.06385464, -0.03608133);
        // tVector3d e = tVector3d(-0.10729597, 0.4681182, 0.31776865);

        // 1. calculate old beta
        double old_beta =
            CalcBeta(n0.normalized(), n1.normalized(), e.normalized());

        // 2. calc ana deriv
        tVector3d dbeta_dn0_ana = tVector3d::Zero();
        tVector3d dbeta_dn1_ana = tVector3d::Zero();
        tVector3d _dbeta_de = tVector3d::Zero();
        CalcDBetaDnormal_and_De(n0, n1, e, dbeta_dn0_ana, dbeta_dn1_ana,
                                _dbeta_de);

        // 3. calc num deriv
        double eps = 1e-5;
        tVector3d dbeta_dn0_num = tVector3d::Zero();
        tVector3d dbeta_dn1_num = tVector3d::Zero();
        for (int i = 0; i < 3; i++)
        {
            n0[i] += eps;
            double new_beta =
                CalcBeta(n0.normalized(), n1.normalized(), e.normalized());
            dbeta_dn0_num[i] = (new_beta - old_beta) / eps;
            n0[i] -= eps;
        }
        for (int i = 0; i < 3; i++)
        {
            n1[i] += eps;
            double new_beta =
                CalcBeta(n0.normalized(), n1.normalized(), e.normalized());
            dbeta_dn1_num[i] = (new_beta - old_beta) / eps;
            n1[i] -= eps;
        }
        tVector3d diff_n0 = (dbeta_dn0_ana - dbeta_dn0_num);
        tVector3d diff_n1 = (dbeta_dn1_ana - dbeta_dn1_num);
        double diff_n0_norm = diff_n0.norm();
        double diff_n1_norm = diff_n1.norm();
        if (diff_n0_norm > 0.1 || diff_n1_norm > 0.1)
        {
            std::cout << "dbeta_dn0_ana = " << dbeta_dn0_ana.transpose()
                      << std::endl;
            std::cout << "dbeta_dn0_num = " << dbeta_dn0_num.transpose()
                      << std::endl;
            std::cout << "dbeta_dn0_diff = " << diff_n0.transpose()
                      << " norm = " << diff_n0_norm << std::endl;
            std::cout << "dbeta_dn1_ana = " << dbeta_dn1_ana.transpose()
                      << std::endl;
            std::cout << "dbeta_dn1_num = " << dbeta_dn1_num.transpose()
                      << std::endl;
            std::cout << "dbeta_dn1_diff = " << diff_n1.transpose()
                      << " norm = " << diff_n1_norm << std::endl;
            exit(1);
        }
    }
    std::cout << "VerifyDBetaDnormal succ\n";
}
void VerifyDBetaDx()
{
    tVectorXd pos = tVectorXd::Random(12);
    // pos.segment(0, 3) = tVector3d(0.5488135, 0.71518937, 0.60276338);
    // pos.segment(3, 3) = tVector3d(0.54488318, 0.4236548, 0.64589411);
    // pos.segment(6, 3) = tVector3d(0.43758721, 0.891773, 0.96366276);
    // pos.segment(9, 3) = tVector3d(0.38344152, 0.79172504, 0.52889492);

    // pos.segment(0, 3) = tVector3d(0.54488318, 0.4236548, 0.64589411);
    // pos.segment(3, 3) = tVector3d(0.43758721, 0.891773, 0.96366276);
    // pos.segment(6, 3) = tVector3d(0.5488135, 0.71518937, 0.60276338);
    // pos.segment(9, 3) = tVector3d(0.38344152, 0.79172504, 0.52889492);

    // pos[0] =
    // 1. get old value
    tVectorXd v0 = pos.segment(0, 3), v1 = pos.segment(3, 3),
              v2 = pos.segment(6, 3), v3 = pos.segment(9, 3);
    double old_beta = CalcBeta(v0, v1, v2, v3);
    std::cout << "old_beta = " << old_beta << std::endl;
    // 2. get ana deriv
    tVectorXd dbeta_dx_ana = CalcDBetaDx(v0, v1, v2, v3);
    tVectorXd dbeta_dx_num = tVectorXd::Zero(12);

    double eps = 1e-5;
    for (int i = 0; i < 12; i++)
    {
        pos[i] += eps;
        v0 = pos.segment(0, 3);
        v1 = pos.segment(3, 3);
        v2 = pos.segment(6, 3);
        v3 = pos.segment(9, 3);
        double new_beta = CalcBeta(v0, v1, v2, v3);
        dbeta_dx_num[i] = (new_beta - old_beta) / eps;
        pos[i] -= eps;
    }
    // 3. get num deriv
    tVectorXd diff = dbeta_dx_ana - dbeta_dx_num;
    std::cout << "dbeta_dx_ana = " << dbeta_dx_ana.transpose() << std::endl;
    std::cout << "dbeta_dx_num = " << dbeta_dx_num.transpose() << std::endl;
    std::cout << "diff = " << diff.transpose() << std::endl;
}
void cDihedralMaterial::Update()
{
    printf("dihedral materia update now!\n");

    int numOfEdges = this->mEdgeArray.size();

    for (int i = 0; i < numOfEdges; i++)
    {
        auto cur_e = mEdgeArray[i];
        if (cur_e->mIsBoundary)
            continue;
        if (i == 4)
        {
            i = 4;
        }
        // printf("-------------[dih] edge %d ------------\n", i);
        // 1. calculate force for each element, storage into the vector

        /*
        all variables are static, except dbeta (Gauss-Newton method, quasi
        Newton)

        shape   = e^2 / (2 * (A1 + A2))
                = e / (h1 + h2)
        dt0 = -1.0 / h1 * n0;
        dt1 = cos(alpha_2) / h1 * n0 + cos(alpha_2_tile) / h1 * n1
        dt2 = cos(alpha_1) / h2 * n0 + cos(alpha_1_tile) / h2 * n1

        dbeta = [dt0, dt1, dt2, dt3] \in R^{12}

        f = -modulus * shape * beta  * dbeta / 2
        H = -modulus * shape * dbeta * dbeta.T
        */
        tVector4i v_id_lst = mEdgeConstraintVertexLst[i];
        tVector3d v0 = mVertexArray[v_id_lst[0]]->mPos.segment(0, 3);
        tVector3d v1 = mVertexArray[v_id_lst[1]]->mPos.segment(0, 3);
        tVector3d v2 = mVertexArray[v_id_lst[2]]->mPos.segment(0, 3);
        tVector3d v3 = mVertexArray[v_id_lst[3]]->mPos.segment(0, 3);

        double beta = CalcBeta(v0, v1, v2, v3);
        tVectorXd dAdx = CalcDBetaDx(v0, v1, v2, v3);
        // 2. calcualte shape
        double e = mRawEdgeLengthArray[i];
        tVector2f height = mRawHeightArray[i];
        double factor = -mEdgeK * e / (height[0] + height[1]);
        tVectorXd force = factor * beta * dAdx / 2;
        tMatrix12d hessian = factor * dAdx * dAdx.transpose();
        if (hessian.hasNaN() == true)
        {
            std::cout << "edge " << i << " hessian nan = \n"
                      << hessian << std::endl;
            exit(1);
        }
        // std::cout << "force = " << force.transpose() << std::endl;
        // std::cout << "hessian = " << hessian.transpose() << std::endl;
        mEleKLst[i] = -hessian.cast<float>();

        // assign the force and hessian into the global stiffness matrix
        mEleForceArray[i] = force;
        // exit(1);
    }

    // 2. calculate elementwise K, form the triplet, then set from
    AssignForceAndMatrix();
    // exit(1);
}

void cDihedralMaterial::AssignForceAndMatrix()
{
    int numOfEdges = mEdgeArray.size();
    mGlobalForce.setZero();
    mStiffnessMat.setZero();
    std::vector<tTriplet> triplets = {};
    for (int i = 0; i < numOfEdges; i++)
    {
        // 1. get edge vertices
        auto cur_e = mEdgeArray[i];
        if (cur_e->mIsBoundary)
            continue;
        tVector4i vid = mEdgeConstraintVertexLst[i];
        const tVector12d &ele_f = mEleForceArray[i];
        const tMatrix12d &ele_H = mEleKLst[i].cast<double>();
        // local id
        for (int v_lid = 0; v_lid < 4; v_lid++)
        {
            // global id
            int v_gid = vid[v_lid];

            // 2. dispatch force vector
            mGlobalForce.segment(3 * v_gid, 3) +=
                ele_f.segment(3 * v_lid, 3).cast<double>();

            // 3. dispatch sparse hessian (3 * 3 block)
            for (int v_lid2 = 0; v_lid2 < 4; v_lid2++)
            {
                int v_gid2 = vid[v_lid2];
                tMatrix3d H_block = ele_H.block(3 * v_lid, 3 * v_lid2, 3, 3);
                for (int a = 0; a < 3; a++)
                    for (int b = 0; b < 3; b++)
                    {
                        triplets.push_back(tTriplet(
                            3 * v_gid + a, 3 * v_gid2 + b, H_block(a, b)));
                    }
            }
        }
    }
    this->mStiffnessMat.setFromTriplets(triplets.begin(), triplets.end());
}