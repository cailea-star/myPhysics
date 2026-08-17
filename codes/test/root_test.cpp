/**
 * @file    root_test.cpp
 * @author  cailea
 * @date    2026-01-29
 * @brief   Test scalar root-finding algorithms.
 */

#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include "root.hpp"

/**
 * @brief  Evaluate a polynomial with roots at 1 and 30.
 * @math   f(x)=(x-1)(x-30)
 * @output Function value f(x).
 */
double polynomial_func(double x_F) {
    return (x_F - 1.0) * (x_F - 30.0);
}

/**
 * @brief  Evaluate an exponential function with a root at 2.
 * @math   f(x)=2^x-4
 * @output Function value f(x).
 */
double exponential_func(double x_F) {
    return std::pow(2.0, x_F) - 4.0;
}

/**
 * @brief  Evaluate a quadratic function with a root at 0.
 * @math   f(x)=x^2
 * @output Function value f(x).
 */
double square_func(double x_F) {
    return x_F * x_F;
}

/**
 * @brief  Test bisection, safeguarded Newton, and Brent root finding.
 * @math   f(x^*)=0
 * @output Process exit status.
 */
int main() {
    double xlo_F = -1.0;
    double xup_F = 20.0;
    double tol_F = 1e-6;
    std::cout << std::scientific << std::setprecision(6) << std::left;

    // Bisection method.
    double xroot1_F = root_bisection(polynomial_func, xlo_F, xup_F, tol_F);
    double xroot2_F = root_bisection(exponential_func, xlo_F, xup_F, tol_F);
    std::cout << "\n[INPUT][bisection] x=[" << xlo_F << "," << xup_F << "], tol=" << tol_F << "\n";
    std::cout << "[REFERENCE][bisection] x1=1, f1=0; x2=2, f2=0\n";
    std::cout << "[COMPUTED][bisection] " << std::setw(22) << "root of (x-1)(x-30)" << "x=" << xroot1_F << ", f(x)=" << polynomial_func(xroot1_F) << "\n";
    std::cout << "[COMPUTED][bisection] " << std::setw(22) << "root of 2^x-4" << "x=" << xroot2_F << ", f(x)=" << exponential_func(xroot2_F) << "\n";
    assert(std::abs(polynomial_func(xroot1_F)) < 10.0 * tol_F && std::abs(exponential_func(xroot2_F)) < 10.0 * tol_F);

    // Safeguarded Newton method.
    xroot1_F = root_newton(polynomial_func, xlo_F, xup_F, tol_F);
    xroot2_F = root_newton(exponential_func, xlo_F, xup_F, tol_F);
    std::cout << "\n[INPUT][Newton] x=[" << xlo_F << "," << xup_F << "], tol=" << tol_F << "\n";
    std::cout << "[REFERENCE][Newton] x1=1, f1=0; x2=2, f2=0\n";
    std::cout << "[COMPUTED][Newton] " << std::setw(22) << "root of (x-1)(x-30)" << "x=" << xroot1_F << ", f(x)=" << polynomial_func(xroot1_F) << "\n";
    std::cout << "[COMPUTED][Newton] " << std::setw(22) << "root of 2^x-4" << "x=" << xroot2_F << ", f(x)=" << exponential_func(xroot2_F) << "\n";
    assert(std::abs(polynomial_func(xroot1_F)) < tol_F && std::abs(exponential_func(xroot2_F)) < tol_F);

    // Brent method.
    xroot1_F = root_brent(polynomial_func, xlo_F, xup_F, tol_F);
    xroot2_F = root_brent(exponential_func, xlo_F, xup_F, tol_F);
    std::cout << "\n[INPUT][Brent] x=[" << xlo_F << "," << xup_F << "], tol=" << tol_F << "\n";
    std::cout << "[REFERENCE][Brent] x1=1, f1=0; x2=2, f2=0\n";
    std::cout << "[COMPUTED][Brent] " << std::setw(22) << "root of (x-1)(x-30)" << "x=" << xroot1_F << ", f(x)=" << polynomial_func(xroot1_F) << "\n";
    std::cout << "[COMPUTED][Brent] " << std::setw(22) << "root of 2^x-4" << "x=" << xroot2_F << ", f(x)=" << exponential_func(xroot2_F) << "\n";
    assert(std::abs(polynomial_func(xroot1_F)) < tol_F && std::abs(exponential_func(xroot2_F)) < tol_F);

    // Recoverable Brent failure without a sign change.
    double xroot3_F = root_brent(square_func, xlo_F, xup_F, tol_F);
    std::cout << "\n[INPUT][Brent unbracketed] f(x)=x^2, x=[" << xlo_F << "," << xup_F << "], tol=" << tol_F << "\n";
    std::cout << "[REFERENCE][Brent unbracketed] x=nan, f(x)=nan\n";
    std::cout << "[COMPUTED][Brent unbracketed] x=" << xroot3_F << ", f(x)=" << square_func(xroot3_F) << "\n";
    assert(std::isnan(xroot3_F));
    return 0;
}
