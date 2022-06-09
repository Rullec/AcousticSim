#include "sim/softbody/AcousticLinearElasticity.h"
#include "utils/LogUtil.h"
// calculate M matrix: its determinant is 6 * volume
tMatrixXd CalcM_6V(const tEigenArr<tVector> &v_array)
{
    tMatrix M = tMatrix::Zero();
    M.col(0).array() = 1;
    for (int i = 0; i < 4; i++)
    {
        M.block(i, 1, 1, 3) = v_array[i].segment(0, 3).transpose();
    };
    return M;
}

tMatrixXd CalcB(const tEigenArr<tVector> &v_array)
{
    // 1. calculate coef
    tMatrixXd M = CalcM_6V(v_array);
    double Mdet_6V = M.determinant();
    SIM_ASSERT(Mdet_6V > 0);

    tMatrix C = Mdet_6V * M.inverse().transpose();
    tVector beta = C.col(1), gamma = C.col(2), delta = C.col(3);

    tMatrixXd B = tMatrixXd::Zero(6, 12);
    for (int i = 0; i < 4; i++)
    {
        B(0, 3 * i) = beta[i];
        B(1, 3 * i + 1) = gamma[i];
        B(2, 3 * i + 2) = delta[i];
    }

    for (int i = 0; i < 4; i++)
    {
        B(3, 3 * i) = gamma[i];
        B(3, 3 * i + 1) = beta[i];

        B(4, 3 * i + 1) = delta[i];
        B(4, 3 * i + 2) = gamma[i];

        B(5, 3 * i + 0) = delta[i];
        B(5, 3 * i + 2) = beta[i];
    }

    B /= Mdet_6V;
    return B;
}

tMatrixXd CalcD(double youngs, double poisson_ratio)
{
    tMatrixXd D = tMatrixXd::Zero(6, 6);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
        {
            if (i == j)
            {
                D(i, j) = 1 - poisson_ratio;
            }
            else
            {
                D(i, j) = poisson_ratio;
            }
        }
    for (int j = 3; j < 6; j++)
    {
        D(j, j) = (1 - 2 * poisson_ratio) / 2;
    }
    D *= youngs / (1 + poisson_ratio) / (1 - 2 * poisson_ratio);
    return D;
}

tSparseMatd cAcousticLinearElasticity::CalcGlobalStiffness(
    const std::vector<tVertexPtr> &v_array,
    const std::vector<tTetPtr> &tet_array, double youngs_modulus,
    double poisson_ratio, double damping_a, double damping_b)
{
    // 1. calculate D
    tMatrixXd D = CalcD(youngs_modulus, poisson_ratio);
    int num_of_tet = tet_array.size();
    tMatrixXd B;
    std::vector<tTriplet> tri_array = {};
    for (int i = 0; i < num_of_tet; i++)
    {
        // 2. for each tet
        auto cur_tet = tet_array[i];
        tEigenArr<tVector> vertex_pos = {};
        for (int j = 0; j < 4; j++)
            vertex_pos.push_back(v_array[cur_tet->mVertexId[j]]->mPos);

        // 3. calculate B
        B.noalias() = CalcB(vertex_pos);

        double sign_volume = CalcM_6V(vertex_pos).determinant() / 6;
        SIM_ASSERT(sign_volume > 0);

        // 4. calculate K_ele
        tMatrixXd K_ele = sign_volume * B.transpose() * D * B;

        // 5. assign to K_global
        for (int a = 0; a < 4; a++)
        {
            int va = cur_tet->mVertexId[a];
            for (int b = 0; b < 4; b++)
            {
                int vb = cur_tet->mVertexId[b];
                for (int j = 0; j < 3; j++)
                    for (int k = 0; k < 3; k++)
                    {
                        int global_row_id = 3 * va + j;
                        int global_col_id = 3 * vb + k;

                        tri_array.push_back(
                            tTriplet(global_row_id, global_col_id,
                                     K_ele(3 * a + j, 3 * b + k)));
                    }
            }
        }
    }
    int num_of_dof = v_array.size() * 3;
    tSparseMatd K_global(num_of_dof, num_of_dof);
    K_global.setFromTriplets(tri_array.begin(), tri_array.end());
    return K_global;
}