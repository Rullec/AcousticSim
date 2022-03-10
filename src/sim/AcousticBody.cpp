#include "AcousticBody.h"
#include <iostream>

cAcousticBody::cAcousticBody()
{
}
#include <fstream>
tVectorXd sum_wave;
void cAcousticBody::SolveVibration(
    const tVectorXd &MassDiag,
    const tSparseMatd &StiffMat,
    const tVector2f &rayleigh_damping,
    const tVectorXd &xcur_vec,
    const tVectorXd &xprev_vec)
{
    tVectorXd eigenvalues;
    tMatrixXd eigenvecs;
    {
        Eigen::GeneralizedEigenSolver<tMatrixXd> ges;

        tMatrixXd dense_M = tMatrixXd::Zero(MassDiag.size(), MassDiag.size());
        dense_M.diagonal() = MassDiag;
        tMatrixXd dense_K = StiffMat.toDense();

        ges.compute(-dense_K, dense_M);
        eigenvalues = ges.eigenvalues().real();
        eigenvecs = ges.eigenvectors().real().transpose();
        std::cout << "eigenvalues = " << eigenvalues.transpose() << std::endl;
    }

    // solve decoupled linear system
    // auto output = "vib.txt";
    // std::ofstream fout(output);
    double dt = 2e-3;
    size_t steps = 1 / dt;
    sum_wave = tVectorXd::Zero(steps);
    int num_of_dof = MassDiag.size();
    {
        double damp_a = rayleigh_damping[0];
        double damp_b = rayleigh_damping[1];
        tVectorXd m_vec = tVectorXd::Ones(num_of_dof);

        tVectorXd c_vec = damp_a * tVectorXd::Ones(num_of_dof) + damp_b * eigenvalues;
        tVectorXd k_vec = eigenvalues;

        double force_amp = 100;
        tVectorXd fprev_vec = tVectorXd::Ones(num_of_dof) * force_amp;
        tVectorXd UTf0_vec = eigenvecs.transpose() * fprev_vec;
        for (size_t i = 0; i < num_of_dof; i++)
        {
            double m = m_vec[i],
                   c = c_vec[i],
                   k = k_vec[i],
                   fprev = UTf0_vec[i];
            if (k < 1e-6)
                continue;

            /*  
                    w = sqrt(k / m)
                    xi = c / (2 * m * w)
                    wd = w * sqrt(1 - xi^2)
                    eps = e^{- xi * w * h}
                    theta = w_d * h
                    gamma = arcsin(xi)
                */
            double w = std::sqrt(k / m);
            double xi = c / (2 * m * w);
            double wd = w * std::sqrt(1 - xi * xi);
            std::cout << "w = " << w << std::endl;
            std::cout << "wd = " << wd << std::endl;
            double eps = std::exp(-xi * w * dt);
            double theta = wd * dt;
            double gamma = std::asin(xi);
            double xpprev = 0;
            double xprev = 0;
            std::vector<double> x_lst = {xprev};
            for (size_t t = 0; t < steps - 1; t++)
            {
                /*
                xcur = 2 * eps * cos(theta) * xprev
                        - eps^2 * xpprev
                        + 2 * fprev * 
                            (
                                eps * cos(theta + gamma) 
                                - eps^2 * cos(2 * theta + gamma)
                            )
                            /
                            (3 * m * w * wd)
                */
                double xcur = 2 * eps * std::cos(theta) * xprev - eps * eps * xpprev + 2 * fprev * (eps * std::cos(theta + gamma) - eps * eps * std::cos(2 * theta + gamma)) / (3 * m * w * wd);

                // update xprev and xcur, fprev
                xpprev = xprev;
                xprev = xcur;
                // if (t == steps / 2)
                // {
                //     fprev = UTf0_vec[i];
                // }
                // else
                // {
                fprev = 0;
                // }
                x_lst.push_back(xcur);
                if (std::isnan(xcur))
                {
                    std::cout << "m " << m << std::endl;
                    std::cout << "w " << w << std::endl;
                    std::cout << "wd " << wd << std::endl;
                    std::cout << "xi " << xi << std::endl;
                    std::cout << "c " << c << std::endl;
                    exit(1);
                }
            }
            tVectorXd x_vec = tVectorXd(x_lst.size());
            for (size_t j = 0; j < x_vec.size(); j++)
                x_vec[j] = x_lst[j];
            // fout << "dof " << i << " x vec = " << x_vec.transpose() << std::endl;
            sum_wave += x_vec;
        }
    }
    // fout << "sum wave now is = " << sum_wave.transpose() << std::endl;
    // fout.close();
    // std::cout << "write vibrations to " << output << std::endl;
}