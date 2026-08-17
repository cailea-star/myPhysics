#include <Eigen/Core>
#include <ios>
#include <cmath>
#include <iomanip>
#include <ostream>
#include <iostream>
#include <Eigen/Dense>
#include <gsl/gsl_sf_legendre.h>
#include "spherical_fourier.hpp"
#include "spherical_harmonics.hpp"

/**
 * Gaussian f(r) = exp(-r^2/2), its l=0 spherical FT is (2pi)^(3/2) * exp(-k^2/2)
 */
double gaussian_r(double r_F) {
    return std::exp(-r_F * r_F / 2.0);
}

double gaussian_k(double k_F) {
    return std::pow(2.0 * M_PI, 1.5) * std::exp(-k_F * k_F / 2.0);
}

/**
 * Y_lm = P_l^|m|(cosθ) * cos(mφ), P_l^m是归一化的Associated Legendre多项式, 适合测试球面函数积分的正交性;
 */
double ReYlm(int l, int m, double theta, double phi) {
    double y = gsl_sf_legendre_sphPlm(l, std::abs(m), std::cos(theta));
    return  y * std::cos(m * phi);
}

int main() {
    int l_I = 0;
    int n_I = 1001;
    double k_F = 1.0, r_F = 1.0;
    double kmin_F = 0.0, kmax_F = 20.0;
    double rmin_F = 0.0, rmax_F = 20.0;
    Eigen::VectorXd k_F1D_k = Eigen::VectorXd::LinSpaced(n_I, kmin_F, kmax_F);
    Eigen::VectorXd r_F1D_r = Eigen::VectorXd::LinSpaced(n_I, rmin_F, rmax_F);
    Eigen::VectorXd f_F1D_r(n_I);
    Eigen::VectorXd F_F1D_k(n_I);
    for (int i = 0; i < n_I; ++i) {
        f_F1D_r(i) = gaussian_r(r_F1D_r(i));
        F_F1D_k(i) = gaussian_k(k_F1D_k(i));
    }
    Eigen::VectorXd F_simpson_F1D_k(n_I);
    Eigen::VectorXd f_simpson_F1D_r(n_I);
    Eigen::VectorXd F_qag_F1D_k(n_I);
    Eigen::VectorXd f_qag_F1D_r(n_I);
    for (int i = 0; i < n_I; ++i) {
        F_simpson_F1D_k(i) = spherical_FT_r2k(f_F1D_r, r_F1D_r, k_F1D_k(i), l_I);
        f_simpson_F1D_r(i) = spherical_FT_k2r(F_F1D_k, k_F1D_k, r_F1D_r(i), l_I);
        F_qag_F1D_k(i) = spherical_FT_r2k(gaussian_r, k_F1D_k(i), l_I, rmin_F, rmax_F);
        f_qag_F1D_r(i) = spherical_FT_k2r(gaussian_k, r_F1D_r(i), l_I, kmin_F, kmax_F);
    }

    std::cout << std::scientific << std::setprecision(4);

    std::cout << "\n Test spherical_FT_r2k (Simpson vs qag):\n";
    std::cout << std::setw(20) << "k_F" 
    << std::setw(20) << "FT_exact" << std::setw(20) << "err(FT_qag)" << std::setw(20) << "err(FT_simpson)" << "\n";
    for (int i = 0; i < n_I; i += 100) {
        double ft_exact_F = gaussian_k(k_F1D_k(i));
        std::cout << std::setw(20) << k_F1D_k(i) 
        << std::setw(20) << ft_exact_F 
        << std::setw(20) << F_qag_F1D_k(i) - ft_exact_F 
        << std::setw(20) << F_simpson_F1D_k(i) - ft_exact_F << "\n";
    }

    std::cout << "\n Test spherical_FT_k2r (Simpson vs qag):\n";
    std::cout << std::setw(20) << "r_F" 
    << std::setw(20) << "IFT_exact" << std::setw(20) << "err(IFT_qag)" << std::setw(20) << "err(IFT_simpson)" << "\n";
    for (int i = 0; i < n_I; i += 100) {
        double ift_exact_F = gaussian_r(r_F1D_r(i));
        std::cout << std::setw(20) << r_F1D_r(i) 
        << std::setw(20) << ift_exact_F 
        << std::setw(20) << f_qag_F1D_r(i) - ift_exact_F 
        << std::setw(20) << f_simpson_F1D_r(i) - ift_exact_F << "\n";
    }

    std::cout << std::scientific << std::setprecision(3) << std::right;

    std::cout << "\n Test spherical_Omega2lm:\n";
    for (int l1 = 4; l1 <=6; l1++) for (int m1 = 0; m1 <= l1/2; m1++) {
        auto ReYlm_l1m1 = [&](double theta,double phi){return ReYlm(l1,m1,theta,phi);};
        for (int l2 = 4; l2 <=6; l2++) for (int m2 = 0; m2 <= l2/2; m2++) {
            int nTheta = (l1 + l2) / 2 + 1;
            int nPhi = 2 * (m1 + m2) + 1;
            double result = spherical_Omega2lm(ReYlm_l1m1, l2, m2, nTheta, nPhi);
            std::cout << std::setw(12) << result;
        }
        std::cout << "\n";
    }

    std::cout << "\n Test spherical_lm2Omega:\n";
    for (int l1 = 4; l1 <=6; l1++) for (int m1 = 0; m1 <= l1/2; m1++) {
        auto ReYlm_l1m1 = [&](double theta,double phi){return ReYlm(l1,m1,theta,phi);};
        for (int l2 = 4; l2 <=6; l2++) for (int m2 = 0; m2 <= l2/2; m2++) {
            int nTheta = (l1 + l2) / 2 + 1;
            int nPhi = 2 * (m1 + m2) + 1;
            double result = spherical_Omega2lm(ReYlm_l1m1, l2, m2, nTheta, nPhi);
            auto ReYlm_result = [&](double theta,double phi){return spherical_lm2Omega(result, l2, m2, theta, phi);};
            result = spherical_Omega2lm(ReYlm_result, l1, m1, nTheta, nPhi);
            std::cout << std::setw(12) << result;
        }
        std::cout << "\n";
    }

    return 0;
}
