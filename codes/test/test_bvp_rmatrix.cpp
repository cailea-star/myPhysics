/**
 * @file    test_bvp_rmatrix.cpp
 * @author  cailea
 * @date    2026-08-18
 * @brief   Validate the B-spline R-matrix against a constant two-channel solution.
 */

#include <cassert>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <Eigen/Dense>
#include "bvp_rmatrix.hpp"

/**
 * @brief  Evaluate the analytic constant-region scattering amplitude.
 * @math   S = e⁻²ⁱᵏᴸ(λ + ik)/(λ - ik)
 * @output Analytic single-mode scattering amplitude.
 * @note   Assumes x_min = 0 and ψ′(0) = -ikψ(0).
 */
std::complex<double> constant_region_S(double k_F, double q_F, double L_F) {
    const std::complex<double> I_C(0.0, 1.0);
    const std::complex<double> rho_C = (q_F - k_F) / (q_F + k_F);
    const std::complex<double> exp_p_C = std::exp(I_C * q_F * L_F);
    const std::complex<double> exp_m_C = std::exp(-I_C * q_F * L_F);
    const std::complex<double> lambda_C = I_C * q_F * (rho_C * exp_p_C - exp_m_C) / (rho_C * exp_p_C + exp_m_C);
    return std::exp(-2.0 * I_C * k_F * L_F) * (lambda_C + I_C * k_F) / (lambda_C - I_C * k_F);
}

/**
 * @brief  Validate the B-spline R-matrix using a constant two-channel reference.
 * @math   S = U diag(S_1, S_2) U^T
 * @output Labeled reference and computed scattering matrices.
 */
int main() {
    const int Nx_I = 81;
    const int Nquad_I = 5;
    const double xmin_F = 0.0;
    const double xmax_F = 8.0;
    const double k_F = 1.1;

    // B_i(x): cubic B-spline basis on [x_min, x_max].
    Eigen::VectorXd x_F1D_x = Eigen::VectorXd::LinSpaced(Nx_I, xmin_F, xmax_F);
    BSplineBasisFunction b_basis_func(x_F1D_x);
    BSplineBasis b_basis(b_basis_func, Nquad_I, false);

    // q_mode = {q_1, q_2}: internal constant-region wave numbers.
    Eigen::VectorXcd k_C1D_ch = Eigen::VectorXcd::Constant(2, k_F);
    Eigen::Vector2d q_F1D_mode;
    q_F1D_mode << 0.7, 1.5;

    // U: channel-to-mode rotation.
    const double theta_F = 0.35;
    Eigen::Matrix2d U_F2D_ch_mode;
    U_F2D_ch_mode << std::cos(theta_F), -std::sin(theta_F), std::sin(theta_F), std::cos(theta_F);

    // F = -U diag(q_1², q_2²) Uᵀ.
    Eigen::Matrix2cd q2_C2D_mode_mode = q_F1D_mode.array().square().matrix().asDiagonal().toDenseMatrix().cast<std::complex<double>>();
    Eigen::MatrixXcd F0_C2D_ch_ch = -U_F2D_ch_mode.cast<std::complex<double>>() * q2_C2D_mode_mode * U_F2D_ch_mode.transpose().cast<std::complex<double>>();
    auto F_Func = [F0_C2D_ch_ch](double) {return F0_C2D_ch_ch;};

    // Σ_xmin = -ik, Σ_xmax = 0.
    Eigen::VectorXcd kzero_C1D_ch = Eigen::VectorXcd::Zero(k_C1D_ch.size());
    RMatrix r_matrix(b_basis, F_Func, k_C1D_ch, kzero_C1D_ch);

    // O = eⁱᵏˣ, O′ = ikO, I = e⁻ⁱᵏˣ, I′ = -ikI at x_max.
    const std::complex<double> I_C(0.0, 1.0);
    Eigen::VectorXcd O_C1D_ch(2);
    Eigen::VectorXcd dO_C1D_ch(2);
    Eigen::VectorXcd I_C1D_ch(2);
    Eigen::VectorXcd dI_C1D_ch(2);
    for (int ch_I = 0; ch_I < k_C1D_ch.size(); ++ch_I) {
        O_C1D_ch(ch_I) = std::exp(I_C * k_C1D_ch(ch_I) * xmax_F);
        dO_C1D_ch(ch_I) = I_C * k_C1D_ch(ch_I) * O_C1D_ch(ch_I);
        I_C1D_ch(ch_I) = std::exp(-I_C * k_C1D_ch(ch_I) * xmax_F);
        dI_C1D_ch(ch_I) = -I_C * k_C1D_ch(ch_I) * I_C1D_ch(ch_I);
    }
    Eigen::MatrixXcd Scalc_C2D_ch_ch = r_matrix.calc_SMatrix(O_C1D_ch, dO_C1D_ch, I_C1D_ch, dI_C1D_ch);

    // S_ref = U diag(S_1, S_2) Uᵀ.
    Eigen::Vector2cd S_C1D_mode;
    S_C1D_mode << constant_region_S(k_F, q_F1D_mode(0), xmax_F), constant_region_S(k_F, q_F1D_mode(1), xmax_F);
    Eigen::MatrixXcd Sref_C2D_ch_ch = U_F2D_ch_mode.cast<std::complex<double>>() * S_C1D_mode.asDiagonal() * U_F2D_ch_mode.transpose().cast<std::complex<double>>();
    const double Serror_F = (Scalc_C2D_ch_ch - Sref_C2D_ch_ch).cwiseAbs().maxCoeff();

    // (x, k, q, θ) → S_ref, S.
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "[Input] x = [" << xmin_F << ", " << xmax_F << "], k = " << k_F << ", q = " << q_F1D_mode.transpose() << ", theta = " << theta_F << "\n";
    std::cout << "[Reference] S =\n" << Sref_C2D_ch_ch << "\n";
    std::cout << "[Computed] S =\n" << Scalc_C2D_ch_ch << "\n";
    std::cout << "[Computed] max |S - S_ref| = " << Serror_F << "\n";

    // ||S - S_ref||_max < ε.
    assert(Serror_F < 1.0e-5);
    return 0;
}
