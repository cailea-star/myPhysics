/**
 * @file    mcmc_test.cpp
 * @author  cailea
 * @date    2026-05-27
 * @brief   Random-walk Metropolis test for a standard Gaussian target.
 */

#include <Eigen/Core>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include "bayesian_mcmc.hpp"

/**
 * @brief  Evaluate the standard Gaussian log density.
 * @math   ℓ(x) = -xᵀx / 2
 * @output Log density at x_F1D_x.
 */
double gaussian_logp(const Eigen::VectorXd& x_F1D_x) {
    return -0.5 * x_F1D_x.squaredNorm();
}

/**
 * @brief  Test random-walk Metropolis sampling of a standard Gaussian target.
 * @math   x_h ∼ N(0, I)
 * @output Prints labeled samples and statistics; asserts key sampling accuracy.
 */
int main() {
    // Function usage.
    Eigen::VectorXd min_F1D_x(2);
    Eigen::VectorXd max_F1D_x(2);
    min_F1D_x << -4.0, -4.0;
    max_F1D_x << 4.0, 4.0;
    int Nwarm_I = 1000;
    int Nsample_I = 100000;
    int thinWarm_I = 1001;
    int thin_I = 1;
    double tolAvg_F = 5.0e-2;
    double tolCov_F = 5.0e-2;
    MCMCSampler sampler(gaussian_logp, min_F1D_x, max_F1D_x);
    sampler.run(Nwarm_I, thinWarm_I);
    sampler.run(Nsample_I, thin_I);
    double acceptance_F = 0.0;
    for (int h_I = 1; h_I < Nsample_I; ++h_I) {acceptance_F += (sampler.samples_F2D_h_x[h_I] - sampler.samples_F2D_h_x[h_I - 1]).squaredNorm() > 0.0;}
    acceptance_F /= Nsample_I - 1;
    bool inside_B = true;
    for (const Eigen::VectorXd& sample_F1D_x : sampler.samples_F2D_h_x) {
        if ((sample_F1D_x.array() < min_F1D_x.array()).any() || (sample_F1D_x.array() > max_F1D_x.array()).any()) {
            inside_B = false;
            break;
        }
    }

    // Labeled input and output.
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "[Input] minimum bounds = " << min_F1D_x.transpose() << "\n";
    std::cout << "[Input] maximum bounds = " << max_F1D_x.transpose() << "\n";
    std::cout << "[Input] warmup steps = " << Nwarm_I << ", sampling steps = " << Nsample_I << ", intervals = " << thinWarm_I << ", " << thin_I << "\n";
    std::cout << "[Reference] Gaussian posterior mean = " << Eigen::VectorXd::Zero(2).transpose() << "\n";
    std::cout << "[Reference] Gaussian posterior covariance =\n";
    std::cout << Eigen::MatrixXd::Identity(2, 2) << "\n";
    std::cout << "[Computed] samples (h, x) =\n";
    for (int h_I = 0; h_I < 5; ++h_I) {std::cout << h_I << " " << sampler.samples_F2D_h_x[h_I].transpose() << "\n";}
    std::cout << "...\n";
    for (int h_I = Nsample_I - 5; h_I < Nsample_I; ++h_I) {std::cout << h_I << " " << sampler.samples_F2D_h_x[h_I].transpose() << "\n";}
    std::cout << "[Computed] MCMC sample count = " << sampler.samples_F2D_h_x.size() << "\n";
    std::cout << "[Computed] MCMC acceptance rate = " << acceptance_F << "\n";
    std::cout << "[Computed] Gaussian posterior mean = " << sampler.avg_F1D_x.transpose() << "\n";
    std::cout << "[Computed] Gaussian posterior covariance =\n";
    std::cout << sampler.cov_F2D_x_x << "\n";

    // Acceptance asserts.
    assert(sampler.samples_F2D_h_x.size() == static_cast<std::size_t>(Nsample_I));
    assert(inside_B);
    assert(acceptance_F > 0.1 && acceptance_F < 0.4);
    assert(sampler.avg_F1D_x.cwiseAbs().maxCoeff() < tolAvg_F);
    assert((sampler.cov_F2D_x_x - Eigen::MatrixXd::Identity(2, 2)).cwiseAbs().maxCoeff() < tolCov_F);
    return 0;
}
