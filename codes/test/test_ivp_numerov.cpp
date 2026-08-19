/**
 * @file    test_ivp_numerov.cpp
 * @author  cailea
 * @date    2026-08-17
 * @brief   Test scalar, vector, and matrix Numerov propagation.
 */

#include <Eigen/Core>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include "ivp_numerov.hpp"

/**
 * @brief  Return the scalar coefficient for oscillatory Numerov propagation.
 * @math   F(x)=-4
 * @output Scalar coefficient F(x).
 */
double F_scal(double) {
    return -4.0;
}

/**
 * @brief  Fill the diagonal coefficient matrix for two-channel Numerov propagation.
 * @math   F(x)=diag(-4,4)
 * @output Two-channel coefficient matrix F(x).
 */
void F_mat(double, Eigen::Ref<Eigen::MatrixXd> F_F2D_ch_ch) {
    F_F2D_ch_ch << -4.0, 0.0, 0.0, 4.0;
}

/**
 * @brief  Test scalar, vector, and matrix Numerov propagation against analytic solutions.
 * @math   y''(x)=F(x)y(x)
 * @output Process exit status.
 * @note   Forward Numerov propagation is unstable for exponentially decaying modes.
 */
int main() {
    int Nx_I = 100;
    double xmin_F = 0.0;
    double xmax_F = 10.0;
    double tol_F = 1.0e-4;
    Eigen::VectorXd x_F1D_x = Eigen::VectorXd::LinSpaced(Nx_I, xmin_F, xmax_F);
    double ycos0_F = std::cos(2.0 * x_F1D_x(0));
    double ycos1_F = std::cos(2.0 * x_F1D_x(1));
    double ysin0_F = std::sin(2.0 * x_F1D_x(0));
    double ysin1_F = std::sin(2.0 * x_F1D_x(1));
    Eigen::VectorXd ycos_F1D_x = (x_F1D_x * 2.0).array().cos();
    Eigen::VectorXd ysin_F1D_x = (x_F1D_x * 2.0).array().sin();
    double yexpn0_F = std::exp(-2.0 * x_F1D_x(0));
    double yexpn1_F = std::exp(-2.0 * x_F1D_x(1));
    double yexpp0_F = std::exp(2.0 * x_F1D_x(0));
    double yexpp1_F = std::exp(2.0 * x_F1D_x(1));
    Eigen::VectorXd yexpp_F1D_x = (x_F1D_x * 2.0).array().exp();
    Eigen::VectorXd yexpn_F1D_x = (-x_F1D_x * 2.0).array().exp();
    std::cout << std::scientific << std::setprecision(4) << std::right;

    // Scalar Numerov propagation.
    Eigen::VectorXd y_F1D_x = ivp_numerov<double>(F_scal, ycos0_F, ycos1_F, x_F1D_x);
    Eigen::VectorXd err_F1D_x = y_F1D_x - ycos_F1D_x;
    std::cout << "\n[INPUT][scalar] F(x)=-4, x=[" << xmin_F << "," << xmax_F << "], y0=" << ycos0_F << ", y1=" << ycos1_F << "\n";
    std::cout << "[REFERENCE/COMPUTED][scalar]" << std::setw(15) << "x" << std::setw(15) << "y_ref" << std::setw(15) << "y_computed" << std::setw(15) << "error" << "\n";
    for (int x_I = 0; x_I < 5; ++x_I) {std::cout << std::setw(15) << x_F1D_x(x_I) << std::setw(15) << ycos_F1D_x(x_I) << std::setw(15) << y_F1D_x(x_I) << std::setw(15) << err_F1D_x(x_I) << "\n";}
    std::cout << "...\n";
    for (int x_I = Nx_I - 5; x_I < Nx_I; ++x_I) {std::cout << std::setw(15) << x_F1D_x(x_I) << std::setw(15) << ycos_F1D_x(x_I) << std::setw(15) << y_F1D_x(x_I) << std::setw(15) << err_F1D_x(x_I) << "\n";}
    assert(err_F1D_x.cwiseAbs().maxCoeff() < tol_F);

    // Matrix Numerov propagation with Nsol = 1.
    Eigen::MatrixXd y0one_F2D_ch_sol(2, 1);
    Eigen::MatrixXd y1one_F2D_ch_sol(2, 1);
    y0one_F2D_ch_sol << ycos0_F, yexpp0_F;
    y1one_F2D_ch_sol << ycos1_F, yexpp1_F;
    auto yone_F3D_ch_sol_x = ivp_numerov_mat<double>(F_mat, y0one_F2D_ch_sol, y1one_F2D_ch_sol, x_F1D_x);
    Eigen::Map<const Eigen::MatrixXd> yone_F2D_ch_x(yone_F3D_ch_sol_x.data(), 2, Nx_I);
    Eigen::MatrixXd errone_F2D_ch_x = yone_F2D_ch_x;
    errone_F2D_ch_x.row(0) = yone_F2D_ch_x.row(0) - ycos_F1D_x.transpose();
    errone_F2D_ch_x.row(1) = yone_F2D_ch_x.row(1) - yexpp_F1D_x.transpose();
    std::cout << "\n[INPUT][matrix Nsol=1] y0=[" << y0one_F2D_ch_sol.transpose() << "], y1=[" << y1one_F2D_ch_sol.transpose() << "]\n";
    std::cout << "[REFERENCE/COMPUTED][matrix Nsol=1]" << std::setw(15) << "x" << std::setw(15) << "y0_ref" << std::setw(15) << "y1_ref" << std::setw(15) << "y0_computed" << std::setw(15) << "y1_computed" << std::setw(15) << "error0" << std::setw(15) << "error1" << "\n";
    for (int x_I = 0; x_I < 5; ++x_I) {std::cout << std::setw(15) << x_F1D_x(x_I) << std::setw(15) << ycos_F1D_x(x_I) << std::setw(15) << yexpp_F1D_x(x_I) << std::setw(15) << yone_F2D_ch_x(0, x_I) << std::setw(15) << yone_F2D_ch_x(1, x_I) << std::setw(15) << errone_F2D_ch_x(0, x_I) << std::setw(15) << errone_F2D_ch_x(1, x_I) << "\n";}
    std::cout << "...\n";
    for (int x_I = Nx_I - 5; x_I < Nx_I; ++x_I) {std::cout << std::setw(15) << x_F1D_x(x_I) << std::setw(15) << ycos_F1D_x(x_I) << std::setw(15) << yexpp_F1D_x(x_I) << std::setw(15) << yone_F2D_ch_x(0, x_I) << std::setw(15) << yone_F2D_ch_x(1, x_I) << std::setw(15) << errone_F2D_ch_x(0, x_I) << std::setw(15) << errone_F2D_ch_x(1, x_I) << "\n";}
    assert(errone_F2D_ch_x.row(0).cwiseAbs().maxCoeff() < tol_F && (errone_F2D_ch_x.row(1).array() / yexpp_F1D_x.transpose().array()).abs().maxCoeff() < tol_F);

    // Matrix Numerov propagation.
    Eigen::MatrixXd y0_F2D_ch_sol(2, 2);
    Eigen::MatrixXd y1_F2D_ch_sol(2, 2);
    y0_F2D_ch_sol << ycos0_F, ysin0_F, yexpn0_F, yexpp0_F;
    y1_F2D_ch_sol << ycos1_F, ysin1_F, yexpn1_F, yexpp1_F;
    auto ymat_F3D_ch_sol_x = ivp_numerov_mat<double>(F_mat, y0_F2D_ch_sol, y1_F2D_ch_sol, x_F1D_x);
    auto errmat_F3D_ch_sol_x = ymat_F3D_ch_sol_x;
    for (int x_I = 0; x_I < Nx_I; ++x_I) {
        errmat_F3D_ch_sol_x(0, 0, x_I) = ymat_F3D_ch_sol_x(0, 0, x_I) - ycos_F1D_x(x_I);
        errmat_F3D_ch_sol_x(0, 1, x_I) = ymat_F3D_ch_sol_x(0, 1, x_I) - ysin_F1D_x(x_I);
        errmat_F3D_ch_sol_x(1, 0, x_I) = ymat_F3D_ch_sol_x(1, 0, x_I) - yexpn_F1D_x(x_I);
        errmat_F3D_ch_sol_x(1, 1, x_I) = ymat_F3D_ch_sol_x(1, 1, x_I) - yexpp_F1D_x(x_I);
    }
    std::cout << "\n[INPUT][matrix y0]\n" << y0_F2D_ch_sol << "\n";
    std::cout << "[INPUT][matrix y1]\n" << y1_F2D_ch_sol << "\n";
    std::cout << "[REFERENCE][matrix x=0]\n" << y0_F2D_ch_sol << "\n";
    std::cout << "[COMPUTED][matrix x=0]\n";
    std::cout << ymat_F3D_ch_sol_x(0, 0, 0) << " " << ymat_F3D_ch_sol_x(0, 1, 0) << "\n";
    std::cout << ymat_F3D_ch_sol_x(1, 0, 0) << " " << ymat_F3D_ch_sol_x(1, 1, 0) << "\n";
    std::cout << "[REFERENCE][matrix x=Nx-1]\n";
    std::cout << ycos_F1D_x(Nx_I - 1) << " " << ysin_F1D_x(Nx_I - 1) << "\n";
    std::cout << yexpn_F1D_x(Nx_I - 1) << " " << yexpp_F1D_x(Nx_I - 1) << "\n";
    std::cout << "[COMPUTED][matrix x=Nx-1]\n";
    std::cout << ymat_F3D_ch_sol_x(0, 0, Nx_I - 1) << " " << ymat_F3D_ch_sol_x(0, 1, Nx_I - 1) << "\n";
    std::cout << ymat_F3D_ch_sol_x(1, 0, Nx_I - 1) << " " << ymat_F3D_ch_sol_x(1, 1, Nx_I - 1) << "\n";
    std::cout << "[ERROR][matrix x=Nx-1]\n";
    std::cout << errmat_F3D_ch_sol_x(0, 0, Nx_I - 1) << " " << errmat_F3D_ch_sol_x(0, 1, Nx_I - 1) << "\n";
    std::cout << errmat_F3D_ch_sol_x(1, 0, Nx_I - 1) << " " << errmat_F3D_ch_sol_x(1, 1, Nx_I - 1) << "\n";
    assert(std::abs(errmat_F3D_ch_sol_x(1, 0, Nx_I - 1)) < tol_F);
    assert(std::abs(errmat_F3D_ch_sol_x(0, 0, Nx_I - 1)) < tol_F && std::abs(errmat_F3D_ch_sol_x(0, 1, Nx_I - 1)) < tol_F && std::abs(errmat_F3D_ch_sol_x(1, 1, Nx_I - 1) / yexpp_F1D_x(Nx_I - 1)) < tol_F);
}
