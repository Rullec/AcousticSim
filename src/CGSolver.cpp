#include "CGSolver.h"
int max_iters = 3000;
float threshold = 1e-5;
int recalc_gap = 1;
tVectorXd cCGSolver::Solve(const tMatrixXd &A, const tVectorXd &b)
{
    tVectorXd x = tVectorXd::Random(b.size());
    tVectorXd r_cur = b - A * x;
    tVectorXd d_cur = r_cur;
    tVectorXd r_next;
    for (int i = 0; i < max_iters; i++)
    {
        float r_cur_max = r_cur.cwiseAbs().maxCoeff();
        printf("iter %d res max %.5f\n", i, r_cur_max);
        if (r_cur_max < threshold)
            break;
        float alpha_i = r_cur.dot(r_cur) / (d_cur.dot(A * d_cur));

        x.noalias() = x + alpha_i * d_cur;

        // recalculate r
        if (i != 0 && i % recalc_gap == 0)
            r_next = b - A * x;
        else
            r_next.noalias() = r_cur - alpha_i * A * d_cur;

        float beta_iplus1 = r_next.dot(r_next) / r_cur.dot(r_cur);

        d_cur = r_next + beta_iplus1 * d_cur;
        r_cur = r_next;
    }
    return x;
}