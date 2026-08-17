#include <Eigen/Core>
#include <cmath>
#include <iomanip>
#include <iostream>
#include "bayesian_mcmc.hpp"

double gaussian_logp(const Eigen::VectorXd& x_F1D_x) {
    return -0.5 * x_F1D_x.squaredNorm();
}

int main() {
    Eigen::VectorXd min_F1D_x(2);
    Eigen::VectorXd max_F1D_x(2);
    min_F1D_x << -10.0, -10.0;
    max_F1D_x << 10.0, 10.0;

    MCMCSampler sampler(gaussian_logp, min_F1D_x, max_F1D_x);
    sampler.run(1000, 1001); // 预热1000步, 不记录样本
    sampler.run(5000, 1);   // 采样5000步, 每步记录样本

    std::cout << std::scientific << std::setprecision(6);
    std::cout << "MCMC sample count = " << sampler.samples_F2D_h_x.size() << "\n";
    std::cout << "Gaussian posterior mean = " << sampler.avg_F1D_x.transpose() << "\n";
    std::cout << "Gaussian posterior covariance =\n" << sampler.cov_F2D_x_x << "\n";
}


