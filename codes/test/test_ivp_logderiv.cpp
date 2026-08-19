/**
 * @file    test_ivp_logderiv.cpp
 * @author  cailea
 * @date    2026-08-19
 * @brief   Test QR-stabilized matrix log-derivative propagation.
 */

#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>

#include <Eigen/Dense>

#include "ivp_logderiv.hpp"

using doubleC = std::complex<double>;

/**
 * @brief  Test QR-stabilized RK4 propagation against a coupled analytic solution.
 * @math   F=U diag(κ²)U†, Y=U diag(κ tanh(κx))U†
 * @output Process exit status.
 */
int main() {
    // Standard coupled complex input.
    int Nx_I = 4001;
    double xmin_F = 0.0;
    double xmax_F = 4.0;
    double theta_F = 0.4;
    double phase_F = 0.3;
    double costheta_F = std::cos(theta_F);
    double sintheta_F = std::sin(theta_F);
    doubleC phase_C = std::exp(doubleC(0.0, phase_F));
    double tol_F = 1.0e-5;
    Eigen::VectorXd x_F1D_x = Eigen::VectorXd::LinSpaced(Nx_I, xmin_F, xmax_F);
    Eigen::Vector2d k_F1D_mode;
    k_F1D_mode << 1.0, 12.0;
    Eigen::Matrix2cd U_C2D_ch_mode;
    U_C2D_ch_mode << costheta_F, sintheta_F * phase_C, -sintheta_F * std::conj(phase_C), costheta_F;
    Eigen::Matrix2cd Fdiag_C2D_mode_mode = k_F1D_mode.array().square().matrix().cast<doubleC>().asDiagonal();
    Eigen::Matrix2cd F_C2D_ch_ch;
    F_C2D_ch_ch.noalias() = U_C2D_ch_mode * Fdiag_C2D_mode_mode * U_C2D_ch_mode.adjoint();
    Eigen::Matrix2cd Y0_C2D_ch_ch = Eigen::Matrix2cd::Zero();
    Real2TMatFunc<doubleC> F_Func = [&](double, Eigen::Ref<Eigen::MatrixXcd> Fout_C2D_ch_ch) {Fout_C2D_ch_ch = F_C2D_ch_ch;};

    // Y(x) = U diag(κ tanh(κx)) U†.
    Eigen::Tensor<doubleC, 3, Eigen::ColMajor> Y_C3D_ch_ch_x = ivp_logderiv_rk4<doubleC>(F_Func, Y0_C2D_ch_ch, x_F1D_x);
    Eigen::Matrix2cd Ydiag_C2D_mode_mode = Eigen::Matrix2cd::Zero();
    Eigen::Matrix2cd Yref_C2D_ch_ch;
    double maxerror_F = 0.0;
    for (int x_I = 0; x_I < Nx_I; ++x_I) {
        Ydiag_C2D_mode_mode.diagonal() = (k_F1D_mode.array() * (k_F1D_mode.array() * x_F1D_x(x_I)).tanh()).matrix().cast<doubleC>();
        Yref_C2D_ch_ch.noalias() = U_C2D_ch_mode * Ydiag_C2D_mode_mode * U_C2D_ch_mode.adjoint();
        Eigen::Map<const Eigen::Matrix2cd> Yx_C2D_ch_ch(Y_C3D_ch_ch_x.data() + 4 * x_I);
        maxerror_F = std::max(maxerror_F, (Yx_C2D_ch_ch - Yref_C2D_ch_ch).cwiseAbs().maxCoeff());
    }
    Eigen::Map<const Eigen::Matrix2cd> Y0calc_C2D_ch_ch(Y_C3D_ch_ch_x.data());
    Eigen::Map<const Eigen::Matrix2cd> Yfinal_C2D_ch_ch(Y_C3D_ch_ch_x.data() + 4 * (Nx_I - 1));
    double initialerror_F = (Y0calc_C2D_ch_ch - Y0_C2D_ch_ch).cwiseAbs().maxCoeff();

    // Inputs → reference → computed → errors.
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "[Input] x = [" << xmin_F << ", " << xmax_F << "], Nx = " << Nx_I << ", kappa = " << k_F1D_mode.transpose() << "\n";
    std::cout << "[Input] F =\n" << F_C2D_ch_ch << "\n";
    std::cout << "[Input] Y0 =\n" << Y0_C2D_ch_ch << "\n";
    std::cout << "[Reference] Y(xmax) =\n" << Yref_C2D_ch_ch << "\n";
    std::cout << "[Computed] Y(xmax) =\n" << Yfinal_C2D_ch_ch << "\n";
    std::cout << "[Error] initial max error = " << initialerror_F << "\n";
    std::cout << "[Error] global max error = " << maxerror_F << "\n";

    // Acceptance.
    assert(initialerror_F < tol_F);
    assert(maxerror_F < tol_F);
    return 0;
}
