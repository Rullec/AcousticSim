#include "BezierCurve.h"
#include "utils/DefUtil.h"
#include "utils/LogUtil.h"
#include <iostream>
tVectorXd GetLinspace(float low, float high, int num_of_points)
{
    SIM_ASSERT(num_of_points > 1);
    float step = (high - low) / (num_of_points - 1);
    tVectorXd res = tVectorXd::Zero(num_of_points);
    for (int i = 0; i < num_of_points; i++)
    {
        res[i] = i * step;
    }
    return res;
}

extern void CalcEdgeDrawBufferSingle(const tVector &v0, const tVector &v1,
                                     const tVector &edge_normal,
                                     Eigen::Map<tVectorXf> &buffer, int &st_pos,
                                     const tVector &color);

cBezierCurve::cBezierCurve(int num_of_div, const tVector2d &A,
                           const tVector2d &B, const tVector2d &C,
                           const tVector2d &D)
    : A(A), B(B), C(C), D(D), mNumOfDiv(num_of_div)
{
    // init the edge buffer
    InitPointlist(mPointList);
    // printf("point shape %d, %d\n", mPointList.rows(), mPointList.cols());
    // mPointList.row(0) = GetLinspace(0, 1, mNumOfDiv);
    // mPointList.row(1) = GetLinspace(0, 1, mNumOfDiv);

    int st = 0;
    int num_of_edges = GetNumOfDrawEdges();
    int num_of_vertices = num_of_edges * 2;
    int buffer_size = num_of_vertices * RENDERING_SIZE_PER_VERTICE;
    mDrawBuffer.noalias() = tVectorXf::Zero(buffer_size);
    Eigen::Map<tVectorXf> map(mDrawBuffer.data(), mDrawBuffer.size());
    tVector BLACK_COLOR = tVector(0, 0, 0, 1);
    for (int i = 0; i < mNumOfDiv - 1; i++)
    {
        tVector st_pos = tVector(0, mPointList(0, i), mPointList(1, i), 0);
        tVector ed_pos =
            tVector(0, mPointList(0, i + 1), mPointList(1, i + 1), 0);
        CalcEdgeDrawBufferSingle(st_pos, ed_pos, tVector::Zero(), map, st,
                                 BLACK_COLOR);
    }
}

int cBezierCurve::GetNumOfDrawEdges() const { return mNumOfDiv - 1; }

const tVectorXf &cBezierCurve::GetDrawBuffer() { return mDrawBuffer; }

int factorial(int n)
{
    if (n == 0)
        return 1;
    int val = 1;
    for (int i = 1; i <= n; i++)
        val *= i;
    return val;
}

int Com(int n, int k)
{
    SIM_ASSERT((n >= 1) && ((k >= 0) && (k <= n)));
    return int(factorial(n) / (factorial(k) * factorial(n - k)));
}

void cBezierCurve::InitPointlist(tMatrixXd &point_lst)
{
    int order = 3;
    tVectorXd u = GetLinspace(0, 1, mNumOfDiv);
    // std::cout << "u = " << u.size() << std::endl;
    tVectorXd one_minus_u = tVectorXd::Ones(mNumOfDiv) - u;
    // std::cout << "1 - u = " << one_minus_u.size() << std::endl;
    point_lst.noalias() = tMatrixXd::Zero(2, mNumOfDiv);
    tMatrixXd ctrl_point_lst = tMatrixXd::Zero(2, 4);
    ctrl_point_lst.col(0) = A;
    ctrl_point_lst.col(1) = B;
    ctrl_point_lst.col(2) = C;
    ctrl_point_lst.col(3) = D;
    for (int i = 0; i <= order; i++)
    {
        // (N \times 1)
        tVectorXd res =
            (u.array().pow(i) * one_minus_u.array().pow(order - i)).transpose();
        tVectorXd ctrl_pt = ctrl_point_lst.col(i);
        tMatrixXd new_point_incre = Com(order, i) * ctrl_pt * res.transpose();
        point_lst.noalias() += new_point_incre;
    }
}

double cBezierCurve::GetTotalLength() const
{
    double total_length = 0;
    for (int i = 0; i < this->mPointList.cols() - 1; i++)
    {
        total_length += (mPointList.col(i + 1) - mPointList.col(i)).norm();
    }
    return total_length;
}
