/**
 * @file    test_integration.cpp
 * @author  cailea
 * @date    2026-08-17
 * @brief   Tests for one-dimensional integration utilities.
 */

#include "integration.hpp"

#include <Eigen/Dense>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>

double V0_F = 1.0;
double a0_F = 2.0;
double R0_F = 6.0;

/**
 * @brief  Evaluate the Woods-Saxon potential.
 * @math   V(x) = V_0 / [1 + exp((x - R_0) / a_0)]
 * @output Potential value at x.
 */
double WS_func(double x_F) {
    double exp_F = std::exp((x_F - R0_F) / a0_F);
    return V0_F / (1.0 + exp_F);
}

/**
 * @brief  Evaluate the derivative of the Woods-Saxon potential.
 * @math   V'(x) = -V_0 exp((x - R_0) / a_0) / {a_0 [1 + exp((x - R_0) / a_0)]^2}
 * @output Potential derivative at x.
 */
double dWS_func(double x_F) {
    double exp_F = std::exp((x_F - R0_F) / a0_F);
    return -V0_F * std::pow(1.0 + exp_F, -2) * exp_F / a0_F;
}

/**
 * @brief  Test Gauss-Kronrod, Simpson, and trapezoidal quadrature.
 * @math   ∫V'(x) dx = V(x_max) - V(x_min)
 * @output Labeled integrals and acceptance assertions.
 */
int main() {
    // (V', [x_min,x_max]) → {I_QAG,I_S,I_T}.
    double xmin_F = 0.0;
    double xmax_F = 100.0; // x_max ≈ ∞.
    int Nx_I = 1001;
    Eigen::VectorXd x_F1D_x = Eigen::VectorXd::LinSpaced(Nx_I, xmin_F, xmax_F);
    Eigen::VectorXd f_F1D_x(Nx_I);
    for (int x_I = 0; x_I < x_F1D_x.size(); ++x_I) {
        f_F1D_x(x_I) = dWS_func(x_F1D_x(x_I));
    }

    double integralExact_F = WS_func(xmax_F) - WS_func(xmin_F); // ∫V'(x) dx = V(x_max) - V(x_min).
    double integralQag_F = integrate_qag(dWS_func, xmin_F, xmax_F);
    double integralSimpson_F = integrate_simpson(f_F1D_x, x_F1D_x);
    double integralTrapezoidal_F = integrate_trapezoidal(f_F1D_x, x_F1D_x);
    double tol_F = 1.0e-4;

    // (I_ref,I_QAG,I_S,I_T) → stdout.
    std::cout << std::scientific << std::setprecision(4) << std::left;
    std::cout << "[Input] x bounds = [" << xmin_F << ", " << xmax_F << "], Nx = " << Nx_I << "\n";
    std::cout << "[Reference] Woods-Saxon derivative integral = " << integralExact_F << "\n";
    std::cout << "[Computed] QAG = " << integralQag_F << ", error = " << integralQag_F - integralExact_F << "\n";
    std::cout << "[Computed] Simpson = " << integralSimpson_F << ", error = " << integralSimpson_F - integralExact_F << "\n";
    std::cout << "[Computed] trapezoidal = " << integralTrapezoidal_F << ", error = " << integralTrapezoidal_F - integralExact_F << "\n";

    // {|I_QAG-I_ref|,|I_S-I_ref|,|I_T-I_ref|} < tolerances.
    assert(std::abs(integralQag_F - integralExact_F) < tol_F);
    assert(std::abs(integralSimpson_F - integralExact_F) < tol_F);
    assert(std::abs(integralTrapezoidal_F - integralExact_F) < tol_F);
    return 0;
}
