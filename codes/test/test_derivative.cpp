/**
 * @file    test_derivative.cpp
 * @author  cailea
 * @date    2026-08-17
 * @brief   Tests for numerical derivative utilities.
 */

#include <cassert>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <Eigen/Dense>
#include "derivative.hpp"

double V0_F = 1.0;
double a0_F = 2.0;
double R0_F = 6.0;

/**
 * @brief  Evaluate the Woods-Saxon potential directly.
 * @math   V(x) = V_0 / [1 + exp((x - R_0)/a_0)]
 * @output V(x).
 */
double WS_func(double x_F) {
    double exp_F = std::exp((x_F - R0_F) / a0_F);
    return V0_F / (1 + exp_F);
}

/**
 * @brief  Evaluate the analytic first derivative of the Woods-Saxon potential.
 * @math   V'(x) = -V_0 exp((x - R_0)/a_0) / [a_0(1 + exp((x - R_0)/a_0))²]
 * @output V'(x).
 */
double dWS_func(double x_F) {
    double exp_F = std::exp((x_F - R0_F) / a0_F);
    return V0_F * (-1) * std::pow(1 + exp_F, -2) * exp_F * (1 / a0_F);
}

/**
 * @brief  Evaluate the analytic second derivative of the Woods-Saxon potential.
 * @math   V''(x) = 2V_0 e^{2z}/[a_0²(1+e^z)³] - V_0 e^z/[a_0²(1+e^z)²], z = (x-R_0)/a_0
 * @output V''(x).
 */
double ddWS_func(double x_F) {
    // z=(x-R_0)/a_0 → (e^z,V''_1,V''_2).
    double exp_F = std::exp((x_F - R0_F) / a0_F);
    double term1_F = V0_F * (-2) * (-1) * std::pow(1 + exp_F, -3) * std::pow(exp_F, 2) * std::pow(1 / a0_F, 2);
    double term2_F = V0_F * (-1) * std::pow(1 + exp_F, -2) * exp_F * std::pow(1 / a0_F, 2);
    return term1_F + term2_F;
}

/**
 * @brief  Test fourth-order callable and uniform-grid finite differences.
 * @math   f'(x) ≈ D_h f, f''(x) ≈ D_h² f
 * @output Labeled derivative errors with acceptance assertions.
 * @note   Grid endpoints use the lower-order stencils defined by derivative.hpp.
 */
int main() {
    // (grid,V) → ({V_i},{V'_i},{V''_i}).
    int Nx_I = 100;
    double xmin_F = 1.0e-10;
    double xmax_F = 15.0;
    double tol1_F = 1.0e-3;
    double tol2_F = 5.0e-5;
    Eigen::VectorXd x_F1D_x = Eigen::VectorXd::LinSpaced(Nx_I, xmin_F, xmax_F);
    Eigen::VectorXd y_F1D_x(Nx_I);
    for (int x_I = 0; x_I < Nx_I; ++x_I) {y_F1D_x(x_I) = WS_func(x_F1D_x(x_I));}
    Eigen::VectorXd dydx_F1D_x = derivative1<double>(y_F1D_x, x_F1D_x);
    Eigen::VectorXd d2ydx2_F1D_x = derivative2<double>(y_F1D_x, x_F1D_x);

    // (V',V'',errors) → stdout; errors < tolerances.
    std::cout << std::scientific << std::setprecision(2);
    std::cout << "\nInput: Woods-Saxon potential, x in [1e-10, 15], N = 100\n";
    std::cout << "Reference/Computed: callable and uniform-grid derivatives\n";
    std::cout << std::setw(15) << "x" << std::setw(15) << "WS(x)" << std::setw(15) << "d1_ref" << std::setw(15) << "d1_func" << std::setw(15) << "d1_vec" << std::setw(15) << "err1_func" << std::setw(15) << "err1_vec" << std::setw(15) << "d2_ref" << std::setw(15) << "d2_func" << std::setw(15) << "d2_vec" << std::setw(15) << "err2_func" << std::setw(15) << "err2_vec" << "\n";
    for (int x_I = 0; x_I < Nx_I; ++x_I) {
        double x_F = x_F1D_x(x_I);
        double d1Func_F = derivative1(WS_func, x_F);
        double d1Vec_F = dydx_F1D_x(x_I);
        double d2Func_F = derivative2(WS_func, x_F);
        double d2Vec_F = d2ydx2_F1D_x(x_I);
        double err1Func_F = d1Func_F - dWS_func(x_F);
        double err1Vec_F = d1Vec_F - dWS_func(x_F);
        double err2Func_F = d2Func_F - ddWS_func(x_F);
        double err2Vec_F = d2Vec_F - ddWS_func(x_F);
        if (Nx_I > 10 && x_I == 5) {std::cout << std::setw(15) << "..." << "\n";}
        if (!(Nx_I > 10 && x_I >= 5 && x_I < Nx_I - 5)) {
            std::cout << std::setw(15) << x_F << std::setw(15) << WS_func(x_F) << std::setw(15) << dWS_func(x_F) << std::setw(15) << d1Func_F << std::setw(15) << d1Vec_F << std::setw(15) << err1Func_F << std::setw(15) << err1Vec_F << std::setw(15) << ddWS_func(x_F) << std::setw(15) << d2Func_F << std::setw(15) << d2Vec_F << std::setw(15) << err2Func_F << std::setw(15) << err2Vec_F << "\n";
        }
        assert(std::abs(err1Func_F) < tol1_F);
        assert(std::abs(err1Vec_F) < tol1_F);
        assert(std::abs(err2Func_F) < tol2_F);
        assert(std::abs(err2Vec_F) < tol2_F);
    }

    // f(x)=(1+i)x² → f''(x)=2(1+i).
    int Nx_check_I = 7;
    Eigen::VectorXd x_check_F1D_x = Eigen::VectorXd::LinSpaced(Nx_check_I, -1.0, 1.0);
    Eigen::VectorXcd f_check_C1D_x(Nx_check_I);
    std::complex<double> coeff_C(1.0, 1.0);
    for (int x_I = 0; x_I < Nx_check_I; ++x_I) {f_check_C1D_x(x_I) = coeff_C * x_check_F1D_x(x_I) * x_check_F1D_x(x_I);}
    Eigen::VectorXcd d2f_check_C1D_x = derivative2(f_check_C1D_x, x_check_F1D_x);
    std::complex<double> d2f_expected_C = 2.0 * coeff_C;

    // (x_i, f''_i, 2(1+i)) → stdout.
    std::cout << "\nInput: f(x) = (1+i)x^2, x in [-1, 1], N = 7\n";
    std::cout << "Reference/Computed: complex uniform-grid second derivative\n";
    std::cout << std::setw(15) << "x" << std::setw(25) << "d2f_reference" << std::setw(25) << "d2f_computed" << "\n";
    for (int x_I = 0; x_I < Nx_check_I; ++x_I) {std::cout << std::setw(15) << x_check_F1D_x(x_I) << std::setw(25) << d2f_expected_C << std::setw(25) << d2f_check_C1D_x(x_I) << "\n";}

    // ||D_h²f - 2(1+i)||∞ < 10⁻¹².
    for (int x_I = 0; x_I < Nx_check_I; ++x_I) {assert(std::abs(d2f_check_C1D_x(x_I) - d2f_expected_C) < 1.0e-12 && std::abs(derivative2(f_check_C1D_x, x_check_F1D_x, x_I) - d2f_expected_C) < 1.0e-12);}
    return 0;
}
