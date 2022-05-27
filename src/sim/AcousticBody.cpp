#include "AcousticBody.h"
#include "sim/AudioWave.h"
#include <iostream>

cAcousticBody::cAcousticBody() {}
#include <fstream>
tVectorXd CalcAccel(const tVectorXd &x, float dt)
{
    tVectorXd accel(x.size() - 2);
    for (int i = 1; i < x.size() - 1; i++)
    {
        accel[i - 1] = (x[i + 1] + x[i - 1] - 2 * x[i]) / (dt * dt);
    }
    return accel;
}
tVectorXd SolveModeITimeSeries(double m, double c, double k, double fprev,
                               double dt, double steps)
{

    if (k < 1e-6)
        return tVectorXd::Zero(steps);

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
    // std::cout << "mode freq = " << w << std::endl;
    // std::cout << "wd = " << wd << std::endl;

    double qprev = 0;

    tVectorXd q_vec = tVectorXd::Zero(steps);
    if (std::isnan(wd) == false)
    {

        double eps = std::exp(-xi * w * dt);
        double theta = wd * dt;
        double gamma = std::asin(xi);
        double qpprev = 0;

        for (size_t t = 0; t < steps - 1; t++)
        {
            /*
            xcur = 2 * eps * cos(theta) * qprev
                    - eps^2 * qpprev
                    + 2 * fprev *
                        (
                            eps * cos(theta + gamma)
                            - eps^2 * cos(2 * theta + gamma)
                        )
                        /
                        (3 * m * w * wd)
            */
            double qcur = 2 * eps * std::cos(theta) * qprev -
                          eps * eps * qpprev +
                          2 * fprev *
                              (eps * std::cos(theta + gamma) -
                               eps * eps * std::cos(2 * theta + gamma)) /
                              (3 * m * w * wd);

            // update qprev and qcur, fprev
            qpprev = qprev;
            qprev = qcur;

            // the force is applied at first
            fprev = 0;

            q_vec[t + 1] = qcur;
            if (std::isnan(qcur))
            {
                std::cout << "m " << m << std::endl;
                std::cout << "w " << w << std::endl;
                std::cout << "wd " << wd << std::endl;
                std::cout << "xi " << xi << std::endl;
                std::cout << "c " << c << std::endl;
                exit(1);
            }
        }
    }
    else
    {
        printf("[warn] w %.1f wd %.1f nan, return\n", w, wd);
    }

    return q_vec;
}

tDiscretedWavePtr cAcousticBody::SolveVibration(
    const tEigenArr<tVector> &normal_array, const tVectorXd &MassDiag,
    const tSparseMatd &StiffMat, const tVector2f &rayleigh_damping,
    const tVectorXd &xcur_vec, const tVectorXd &qprev_vec, int sampling_freq,
    float duration)
{
    std::cout << "mass = " << MassDiag.transpose() << std::endl;

    tVectorXd eigenvalues;
    tMatrixXd eigenvecs;
    {
        Eigen::GeneralizedEigenSolver<tMatrixXd> ges;

        tMatrixXd dense_M = tMatrixXd::Zero(MassDiag.size(), MassDiag.size());
        dense_M.diagonal() = MassDiag;
        tMatrixXd dense_K = StiffMat.toDense();
        std::cout << "stiffness diag = " << dense_K.diagonal().transpose()
                  << std::endl;

        ges.compute(-dense_K, dense_M);
        eigenvalues = ges.eigenvalues().real();
        eigenvecs = ges.eigenvectors().real().transpose();
        // std::cout << "eigenvalues = " << eigenvalues.transpose() <<
        // std::endl;
    }

    // solve decoupled linear system
    auto mode_output = "mode_vib.txt";
    auto v_output = "v_vib.txt";
    std::ofstream fout_mode(mode_output);
    std::ofstream fout_v(v_output);

    float dt = 1.0 / (sampling_freq * 1.0);
    tDiscretedWavePtr wave = std::make_shared<tDiscretedWave>(dt, duration);

    size_t steps = wave->GetNumOfData();
    int num_of_vertices = MassDiag.size() / 3;
    int num_of_modes = MassDiag.size();
    tMatrixXd vertex_displacement_time_series(num_of_vertices * 3, steps);
    tMatrixXd mode_displacement_time_series(num_of_modes, steps);

    {
        double damp_a = rayleigh_damping[0];
        double damp_b = rayleigh_damping[1];
        tVectorXd m_vec = tVectorXd::Ones(num_of_modes);

        tVectorXd c_vec =
            damp_a * tVectorXd::Ones(num_of_modes) + damp_b * eigenvalues;
        tVectorXd k_vec = eigenvalues;

        double force_amp = 100.0;

        // apply force on the x axis of the first vertex
        tVectorXd fprev_vec = tVectorXd::Zero(num_of_modes);
        fprev_vec[0] = force_amp;

        tVectorXd UTf0_vec = eigenvecs.transpose() * fprev_vec;
        std::cout << "modal k = " << k_vec.transpose() << std::endl;
        tVectorXd freq =
            (k_vec.cwiseProduct(m_vec.cwiseInverse())).cwiseSqrt() / (2 * M_PI);
        std::cout << "freq = " << freq.transpose() << std::endl;
        for (size_t i = 0; i < num_of_modes; i++)
        {
            double m = m_vec[i], c = c_vec[i], k = k_vec[i],
                   fprev = UTf0_vec[i];
            tVectorXd mode_i_time_series =
                SolveModeITimeSeries(m, c, k, fprev, dt, steps);
            mode_displacement_time_series.row(i) = mode_i_time_series;
        }
    }

    // convert mode vibration to vertex vibration (displacement)
    vertex_displacement_time_series = eigenvecs * mode_displacement_time_series;

    fout_mode << sampling_freq << " Hz\n";
    fout_mode << mode_displacement_time_series << std::endl;
    fout_v << sampling_freq << " Hz\n";
    fout_v << vertex_displacement_time_series << std::endl;

    fout_mode.close();
    fout_v.close();

    // set data
    // wave->SetData(vertex_displacement_time_series.row(0).cast<float>());
    float air_density = 1.225; // kg / m^3
    int num_of_v = normal_array.size();
    int num_of_dt = vertex_displacement_time_series.cols();
    tEigenArr<tMatrixXd> sound_pressure_array = {};
    for (int i = 0; i < num_of_v; i++)
    {
        // for vertex i
        tMatrixXd vertex_disp =
            vertex_displacement_time_series.block(3 * i, 0, 3, num_of_dt);
        tVector3d normal = normal_array[i].segment(0, 3);
        // std::cout << "v " << i << " normal = " << normal.transpose() <<
        // std::endl;
        tMatrixXd sound_pressure_alltime = tMatrixXd::Zero(3, num_of_dt - 2);
        for (int t = 1; t < num_of_dt - 1; t++)
        {
            // its accel is
            tVector3d accel = (2 * vertex_disp.col(t) - vertex_disp.col(t - 1) -
                               vertex_disp.col(t + 1)) /
                              (dt * dt);
            tVector3d accel_along_normal = accel.dot(normal) * normal;
            tVector3d sound_pressure = accel_along_normal * air_density;
            // std::cout << "sound pressure = " << sound_pressure.transpose() <<
            // std::endl;
            sound_pressure_alltime.col(t - 1) = sound_pressure;
        }
        sound_pressure_array.push_back(sound_pressure_alltime);
        break;
    }
    // std::cout << "wave = " << sound_pressure_array[0].row(0) << std::endl;
    wave->SetData(sound_pressure_array[0].row(0).cast<float>());
    return wave;
}