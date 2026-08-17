/**
 * @file    spherical_fourier_test.cpp
 * @author  cailea
 * @date    2026-08-17
 * @brief   Tests for spherical Fourier transforms.
 */

#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <Eigen/Core>
#include "spherical_fourier.hpp"

/**
 * @brief  Evaluate the radial Gaussian reference function.
 * @math   f(r) = exp(-r²/2)
 * @output f(r).
 */
double gaussian_r(double r_F) {
    return std::exp(-r_F * r_F / 2.0);
}

/**
 * @brief  Evaluate the momentum-space Gaussian reference function.
 * @math   F(k) = (2π)^{3/2} exp(-k²/2)
 * @output F(k).
 */
double gaussian_k(double k_F) {
    return std::pow(2.0 * M_PI, 1.5) * std::exp(-k_F * k_F / 2.0);
}

/**
 * @brief  Test spherical Fourier transforms using GSL QAG and composite Simpson quadrature.
 * @math   f(r) = exp(-r²/2) ↔ F(k) = (2π)^{3/2} exp(-k²/2)
 * @output Labeled transform comparisons with acceptance assertions.
 * @note   Uses the l = 0 Gaussian transform pair on uniform grids.
 */
int main() {
    int l_I = 0;
    int Ngrid_I = 1001;
    double kmin_F = 0.0;
    double kmax_F = 20.0;
    double rmin_F = 0.0;
    double rmax_F = 20.0;
    Eigen::VectorXd k_F1D_k = Eigen::VectorXd::LinSpaced(Ngrid_I, kmin_F, kmax_F);
    Eigen::VectorXd r_F1D_r = Eigen::VectorXd::LinSpaced(Ngrid_I, rmin_F, rmax_F);
    Eigen::VectorXd f_F1D_r(Ngrid_I);
    Eigen::VectorXd F_F1D_k(Ngrid_I);
    for (int grid_I = 0; grid_I < Ngrid_I; ++grid_I) {
        f_F1D_r(grid_I) = gaussian_r(r_F1D_r(grid_I));
        F_F1D_k(grid_I) = gaussian_k(k_F1D_k(grid_I));
    }

    Eigen::VectorXd F_simpson_F1D_k(Ngrid_I);
    Eigen::VectorXd F_qag_F1D_k(Ngrid_I);
    for (int k_I = 0; k_I < Ngrid_I; ++k_I) {
        F_simpson_F1D_k(k_I) = spherical_FT_r2k(f_F1D_r, r_F1D_r, k_F1D_k(k_I), l_I);
        F_qag_F1D_k(k_I) = spherical_FT_r2k(gaussian_r, k_F1D_k(k_I), l_I, rmin_F, rmax_F);
    }

    Eigen::VectorXd f_simpson_F1D_r(Ngrid_I);
    Eigen::VectorXd f_qag_F1D_r(Ngrid_I);
    for (int r_I = 0; r_I < Ngrid_I; ++r_I) {
        f_simpson_F1D_r(r_I) = spherical_FT_k2r(F_F1D_k, k_F1D_k, r_F1D_r(r_I), l_I);
        f_qag_F1D_r(r_I) = spherical_FT_k2r(gaussian_k, r_F1D_r(r_I), l_I, kmin_F, kmax_F);
    }

    std::cout << std::scientific << std::setprecision(4);

    std::cout << "\nInput: f_l(r) = exp(-r^2/2), l = 0, r in [0, 20], N = 1001\n";
    std::cout << "Reference/Computed: F_l(k), QAG and Simpson\n";
    std::cout << std::setw(20) << "k_F" << std::setw(20) << "FT_reference" << std::setw(20) << "FT_qag" << std::setw(20) << "FT_simpson" << std::setw(20) << "err_qag" << std::setw(20) << "err_simpson" << "\n";
    for (int k_I = 0; k_I < Ngrid_I; ++k_I) {
        if (Ngrid_I > 10 && k_I == 5) {std::cout << std::setw(20) << "..." << "\n";}
        if (Ngrid_I > 10 && k_I >= 5 && k_I < Ngrid_I - 5) {continue;}
        std::cout << std::setw(20) << k_F1D_k(k_I) << std::setw(20) << F_F1D_k(k_I) << std::setw(20) << F_qag_F1D_k(k_I) << std::setw(20) << F_simpson_F1D_k(k_I) << std::setw(20) << F_qag_F1D_k(k_I) - F_F1D_k(k_I) << std::setw(20) << F_simpson_F1D_k(k_I) - F_F1D_k(k_I) << "\n";
    }

    std::cout << "\nInput: F_l(k) = (2pi)^(3/2) exp(-k^2/2), l = 0, k in [0, 20], N = 1001\n";
    std::cout << "Reference/Computed: f_l(r), QAG and Simpson\n";
    std::cout << std::setw(20) << "r_F" << std::setw(20) << "IFT_reference" << std::setw(20) << "IFT_qag" << std::setw(20) << "IFT_simpson" << std::setw(20) << "err_qag" << std::setw(20) << "err_simpson" << "\n";
    for (int r_I = 0; r_I < Ngrid_I; ++r_I) {
        if (Ngrid_I > 10 && r_I == 5) {std::cout << std::setw(20) << "..." << "\n";}
        if (Ngrid_I > 10 && r_I >= 5 && r_I < Ngrid_I - 5) {continue;}
        std::cout << std::setw(20) << r_F1D_r(r_I) << std::setw(20) << f_F1D_r(r_I) << std::setw(20) << f_qag_F1D_r(r_I) << std::setw(20) << f_simpson_F1D_r(r_I) << std::setw(20) << f_qag_F1D_r(r_I) - f_F1D_r(r_I) << std::setw(20) << f_simpson_F1D_r(r_I) - f_F1D_r(r_I) << "\n";
    }

    assert((F_qag_F1D_k - F_F1D_k).cwiseAbs().maxCoeff() < 1.0e-8);
    assert((F_simpson_F1D_k - F_F1D_k).cwiseAbs().maxCoeff() < 1.0e-8);
    assert((f_qag_F1D_r - f_F1D_r).cwiseAbs().maxCoeff() < 1.0e-8);
    assert((f_simpson_F1D_r - f_F1D_r).cwiseAbs().maxCoeff() < 1.0e-8);
    return 0;
}
