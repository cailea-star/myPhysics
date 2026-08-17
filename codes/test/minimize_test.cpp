/**
 * @file    minimize_test.cpp
 * @author  cailea
 * @date    2026-01-07
 * @brief   Manual tests for scalar minimization algorithms.
 */

#include <cmath>
#include <iomanip>
#include <iostream>
#include "minimize.hpp"

/**
 * @brief  Evaluate a quadratic test function.
 * @math   f(x) = (x - 2)^2 + 1
 * @output Returns f(x), whose minimum is f(2) = 1.
 */
double x2_1_func(double x_F) {
    return std::pow(x_F - 2.0, 2) + 1.0;
}

/**
 * @brief  Evaluate a piecewise quadratic test function.
 * @math   f(x) = (x - 2)^2 + 1 for x < 4; f(x) = (x - 6)^2 otherwise
 * @output Returns f(x), with local minima at x = 2 and x = 6.
 */
double x2_1_x6_0_func(double x_F) {
    double fleft_F = std::pow(x_F - 2.0, 2) + 1.0;
    double fright_F = std::pow(x_F - 6.0, 2);
    return x_F < 4.0 ? fleft_F : fright_F;
}

/**
 * @brief  Evaluate a monotone exponential test function.
 * @math   f(x) = exp(x)
 * @output Returns f(x), whose minimum on [0,+infinity) is f(0) = 1.
 */
double exp_func(double x_F) {
    return std::exp(x_F);
}

/**
 * @brief  Run and print scalar minimization tests.
 * @math   x^* = arg min_{x in [0,20]} f(x)
 * @output Prints golden-section and Brent minimization results.
 */
int main() {
    double xlo_F = 0.0;
    double xup_F = 20.0;
    double tol_F = 1e-6;

    std::cout << std::scientific << std::setprecision(6) << std::left;

    // Golden-section search.
    double xmin1_F = minimize_golden(x2_1_func, xlo_F, xup_F, tol_F);
    double xmin2_F = minimize_golden(x2_1_x6_0_func, xlo_F, xup_F, tol_F);
    std::cout << "Golden Section Search Results:\n";
    std::cout << std::setw(50) << "Minimum of (x-2)^2 + 1:" << "x = " << xmin1_F << ", f(x) = " << x2_1_func(xmin1_F) << "\n";
    std::cout << std::setw(50) << "Minimum of [(x-2)^2 + 1] or [(x-6)^2 + 0]: " << "x = " << xmin2_F << ", f(x) = " << x2_1_x6_0_func(xmin2_F) << "\n";

    // Brent minimization.
    xmin1_F = minimize_brent(x2_1_func, xlo_F, xup_F, tol_F);
    xmin2_F = minimize_brent(x2_1_x6_0_func, xlo_F, xup_F, tol_F);
    std::cout << "\nBrent's Method Search Results:\n";
    std::cout << std::setw(50) << "Minimum of (x-2)^2 + 1: " << "x = " << xmin1_F << ", f(x) = " << x2_1_func(xmin1_F) << "\n";
    std::cout << std::setw(50) << "Minimum of [(x-2)^2 + 1] or [(x-6)^2 + 0]: " << "x = " << xmin2_F << ", f(x) = " << x2_1_x6_0_func(xmin2_F) << "\n";

    // Boundary minimum: Brent initialization cannot find an interior bracket.
    std::cout << "\nTesting exp(x) at boundary (x=0):\n";
    double xmin3_F = minimize_golden(exp_func, xlo_F, xup_F, tol_F);
    std::cout << std::setw(50) << "Minimum of exp(x) with golden: " << "x = " << xmin3_F << ", f(x) = " << exp_func(xmin3_F) << "\n";
    xmin3_F = minimize_brent(exp_func, xlo_F, xup_F, tol_F);
    std::cout << std::setw(50) << "Minimum of exp(x) with brent: " << "x = " << xmin3_F << ", f(x) = " << exp_func(xmin3_F) << "\n";
    return 0;
}
