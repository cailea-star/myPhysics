/**
 * @file    root_test.cpp
 * @author  cailea
 * @date    2026-01-29
 * @brief   Manual tests for scalar root-finding methods.
 */

#include <cmath>
#include <iomanip>
#include <iostream>
#include "root.hpp"

/**
 * @brief  Evaluate a polynomial with roots at 1 and 30.
 * @math   f(x) = (x - 1)(x - 30)
 * @output Polynomial value at x.
 */
double polynomial_func(double x_F) {
    return (x_F - 1.0) * (x_F - 30.0);
}

/**
 * @brief  Evaluate an exponential function with a root at 2.
 * @math   f(x) = 2^x - 4
 * @output Exponential function value at x.
 */
double exponential_func(double x_F) {
    return std::pow(2.0, x_F) - 4.0;
}

/**
 * @brief  Evaluate a quadratic function with a root at 0.
 * @math   f(x) = x^2
 * @output Quadratic function value at x.
 */
double square_func(double x_F) {
    return x_F * x_F;
}

/**
 * @brief  Run and print scalar root-finding tests.
 * @math   f(x^*) = 0
 * @output Prints roots and residuals for bisection, Newton, and Brent methods.
 */
int main() {
    double xlo_F = -1.0;
    double xup_F = 20.0;
    double tol_F = 1e-6;

    std::cout << std::scientific << std::setprecision(6) << std::left;

    // Bisection method.
    double xroot1_F = root_bisection(polynomial_func, xlo_F, xup_F, tol_F);
    double xroot2_F = root_bisection(exponential_func, xlo_F, xup_F, tol_F);
    std::cout << "Bisection Method Results:\n";
    std::cout << std::setw(30) << "Root of (x-1)(x-30): " << "x = " << xroot1_F << ", f(x) = " << polynomial_func(xroot1_F) << "\n";
    std::cout << std::setw(30) << "Root of 2^x-4: " << "x = " << xroot2_F << ", f(x) = " << exponential_func(xroot2_F) << "\n";

    // Safeguarded Newton method with decreasing and increasing functions.
    std::cout << "\nNewton's Method Results:\n";
    xroot1_F = root_newton(polynomial_func, xlo_F, xup_F, tol_F);
    xroot2_F = root_newton(exponential_func, xlo_F, xup_F, tol_F);
    std::cout << std::setw(30) << "Root of (x-1)(x-30): " << "x = " << xroot1_F << ", f(x) = " << polynomial_func(xroot1_F) << "\n";
    std::cout << std::setw(30) << "Root of 2^x-4: " << "x = " << xroot2_F << ", f(x) = " << exponential_func(xroot2_F) << "\n";

    // Brent's method.
    std::cout << "\nBrent's Method Results:\n";
    xroot1_F = root_brent(polynomial_func, xlo_F, xup_F, tol_F);
    xroot2_F = root_brent(exponential_func, xlo_F, xup_F, tol_F);
    std::cout << std::setw(30) << "Root of (x-1)(x-30): " << "x = " << xroot1_F << ", f(x) = " << polynomial_func(xroot1_F) << "\n";
    std::cout << std::setw(30) << "Root of 2^x-4: " << "x = " << xroot2_F << ", f(x) = " << exponential_func(xroot2_F) << "\n";

    // Recoverable Brent failure without a sign change.
    std::cout << "\nTesting an unbracketed even-multiplicity root (x^2=0):\n";
    double xroot3_F = root_brent(square_func, xlo_F, xup_F, tol_F);
    std::cout << std::setw(30) << "Root of x^2: " << "x = " << xroot3_F << ", f(x) = " << square_func(xroot3_F) << "\n";
    return 0;
}
