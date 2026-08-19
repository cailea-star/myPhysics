/**
 * @file    test_ivp_rungekutta.cpp
 * @author  cailea
 * @date    2026-08-17
 * @brief   Test scalar and vector fourth-order Runge-Kutta propagation.
 */

#include <Eigen/Core>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include "ivp_rungekutta.hpp"

/**
 * @brief  Evaluate the scalar initial-value equation used by RK4.
 * @math   y'=cos(x)
 * @output Scalar derivative y'.
 */
double sinfFunc(double x_F, double y_F) {
    return std::cos(x_F);
}

/**
 * @brief  Evaluate the two-channel initial-value equation used by vector RK4.
 * @math   (y₁',y₂')=(y₂,y₁)
 * @output Writes the two-channel derivative vector.
 */
void expfFunc(double x_F, const Eigen::Ref<const Eigen::VectorXd>& y_F1D_ch, Eigen::Ref<Eigen::VectorXd> dydx_F1D_ch) {
    dydx_F1D_ch(0) = y_F1D_ch(1);
    dydx_F1D_ch(1) = y_F1D_ch(0);
}

/**
 * @brief  Test scalar and vector fourth-order Runge-Kutta propagation.
 * @math   y_{n+1}=RK4[f,x_n,y_n,h]
 * @output Process exit status.
 */
int main() {
    int Nx_I = 100;
    double xmin_F = 0.0;
    double xmax_F = 10.0;
    double tol_F = 1.0e-4;
    Eigen::VectorXd x_F1D_x = Eigen::VectorXd::LinSpaced(Nx_I, xmin_F, xmax_F);
    double ysin0_F = std::sin(x_F1D_x(0));
    double yexpp0_F = std::exp(x_F1D_x(0));
    double yexpn0_F = std::exp(-x_F1D_x(0));
    Eigen::VectorXd yexppvec0_F1D_ch(2);
    yexppvec0_F1D_ch(0) = yexpp0_F;
    yexppvec0_F1D_ch(1) = yexpp0_F;
    Eigen::VectorXd yexpnvec0_F1D_ch(2);
    yexpnvec0_F1D_ch(0) = yexpn0_F;
    yexpnvec0_F1D_ch(1) = -yexpn0_F;
    Eigen::VectorXd sin_F1D_x = x_F1D_x.array().sin();
    Eigen::VectorXd expp_F1D_x = x_F1D_x.array().exp();
    Eigen::VectorXd expn_F1D_x = (-x_F1D_x).array().exp();
    std::cout << std::scientific << std::setprecision(4) << std::right;

    // Scalar RK4 propagation.
    Eigen::VectorXd ysin_F1D_x = ivp_rk4<double>(sinfFunc, ysin0_F, x_F1D_x);
    Eigen::VectorXd err_sin_F1D_x = ysin_F1D_x - sin_F1D_x;
    std::cout << "\n[INPUT][scalar RK4] x=[" << xmin_F << "," << xmax_F << "], y0=" << ysin0_F << "\n";
    std::cout << "[REFERENCE/COMPUTED][scalar RK4]\n";
    std::cout << std::setw(15) << "x" << std::setw(15) << "y_ref" << std::setw(15) << "y_computed" << std::setw(15) << "error" << "\n";
    for (int x_I = 0; x_I < 5; ++x_I) {std::cout << std::setw(15) << x_F1D_x(x_I) << std::setw(15) << sin_F1D_x(x_I) << std::setw(15) << ysin_F1D_x(x_I) << std::setw(15) << err_sin_F1D_x(x_I) << "\n";}
    std::cout << "...\n";
    for (int x_I = Nx_I - 5; x_I < Nx_I; ++x_I) {std::cout << std::setw(15) << x_F1D_x(x_I) << std::setw(15) << sin_F1D_x(x_I) << std::setw(15) << ysin_F1D_x(x_I) << std::setw(15) << err_sin_F1D_x(x_I) << "\n";}
    assert(err_sin_F1D_x.cwiseAbs().maxCoeff() < tol_F);

    // Vector RK4 propagation.
    Eigen::MatrixXd yexpp_F2D_ch_x = ivp_rk4_vec<double>(expfFunc, x_F1D_x, yexppvec0_F1D_ch);
    Eigen::MatrixXd yexpn_F2D_ch_x = ivp_rk4_vec<double>(expfFunc, x_F1D_x, yexpnvec0_F1D_ch);
    Eigen::VectorXd err_expp_F1D_x = yexpp_F2D_ch_x.row(0) - expp_F1D_x.transpose();
    Eigen::VectorXd err_expn_F1D_x = yexpn_F2D_ch_x.row(0) - expn_F1D_x.transpose();
    std::cout << "\n[INPUT][vector RK4] y0+=[" << yexppvec0_F1D_ch.transpose() << "], y0-=[" << yexpnvec0_F1D_ch.transpose() << "]\n";
    std::cout << "[REFERENCE/COMPUTED][vector RK4]\n";
    std::cout << std::setw(15) << "x" << std::setw(15) << "yp_ref" << std::setw(15) << "yp_computed" << std::setw(15) << "error+" << std::setw(15) << "yn_ref" << std::setw(15) << "yn_computed" << std::setw(15) << "error-" << "\n";
    for (int x_I = 0; x_I < 5; ++x_I) {std::cout << std::setw(15) << x_F1D_x(x_I) << std::setw(15) << expp_F1D_x(x_I) << std::setw(15) << yexpp_F2D_ch_x(0, x_I) << std::setw(15) << err_expp_F1D_x(x_I) << std::setw(15) << expn_F1D_x(x_I) << std::setw(15) << yexpn_F2D_ch_x(0, x_I) << std::setw(15) << err_expn_F1D_x(x_I) << "\n";}
    std::cout << "...\n";
    for (int x_I = Nx_I - 5; x_I < Nx_I; ++x_I) {std::cout << std::setw(15) << x_F1D_x(x_I) << std::setw(15) << expp_F1D_x(x_I) << std::setw(15) << yexpp_F2D_ch_x(0, x_I) << std::setw(15) << err_expp_F1D_x(x_I) << std::setw(15) << expn_F1D_x(x_I) << std::setw(15) << yexpn_F2D_ch_x(0, x_I) << std::setw(15) << err_expn_F1D_x(x_I) << "\n";}
    assert((err_expp_F1D_x.array() / expp_F1D_x.array()).abs().maxCoeff() < tol_F && err_expn_F1D_x.cwiseAbs().maxCoeff() < tol_F);
    return 0;
}
