/**
 * @file    test_minimize.cpp
 * @author  cailea
 * @date    2026-01-07
 * @brief   Test scalar minimization algorithms.
 */

#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include "minimize.hpp"

/**
 * @brief  Evaluate a quadratic test function.
 * @math   f(x)=(x-2)^2+1
 * @output Function value f(x).
 */
double x2_1_func(double x_F) {
    return std::pow(x_F - 2.0, 2) + 1.0;
}

/**
 * @brief  Evaluate a piecewise quadratic test function.
 * @math   f(x) = {(x-2)²+1, x<4; (x-6)², x≥4}
 * @output Function value f(x).
 */
double x2_1_x6_0_func(double x_F) {
    double fleft_F = std::pow(x_F - 2.0, 2) + 1.0;
    double fright_F = std::pow(x_F - 6.0, 2);
    return x_F < 4.0 ? fleft_F : fright_F;
}

/**
 * @brief  Evaluate a monotone exponential test function.
 * @math   f(x)=exp(x)
 * @output Function value f(x).
 */
double exp_func(double x_F) {
    return std::exp(x_F);
}

/**
 * @brief  Test golden-section and Brent scalar minimization.
 * @math   x^*=arg min_{x∈[0,20]}f(x)
 * @output Labeled minimizers with acceptance assertions.
 */
int main() {
    // ([x_lo,x_up],ε) → minimization inputs.
    double xlo_F = 0.0;
    double xup_F = 20.0;
    double tol_F = 1e-6;
    std::cout << std::scientific << std::setprecision(6) << std::left;

    // (argmin f₁,argmin f₂) = (2,6): golden-section search.
    double xmin1_F = minimize_golden(x2_1_func, xlo_F, xup_F, tol_F);
    double xmin2_F = minimize_golden(x2_1_x6_0_func, xlo_F, xup_F, tol_F);
    std::cout << "\n[INPUT][golden] x=[" << xlo_F << "," << xup_F << "], tol=" << tol_F << "\n";
    std::cout << "[REFERENCE][golden] x1=2, f1=1; x2=6, f2=0\n";
    std::cout << "[COMPUTED][golden] " << std::setw(34) << "minimum of (x-2)^2+1" << "x=" << xmin1_F << ", f(x)=" << x2_1_func(xmin1_F) << "\n";
    std::cout << "[COMPUTED][golden] " << std::setw(34) << "minimum of piecewise quadratic" << "x=" << xmin2_F << ", f(x)=" << x2_1_x6_0_func(xmin2_F) << "\n";
    assert(std::abs(xmin1_F - 2.0) < tol_F && std::abs(xmin2_F - 6.0) < tol_F);

    // (argmin f₁,argmin f₂) = (2,6): Brent minimization.
    xmin1_F = minimize_brent(x2_1_func, xlo_F, xup_F, tol_F);
    xmin2_F = minimize_brent(x2_1_x6_0_func, xlo_F, xup_F, tol_F);
    std::cout << "\n[INPUT][Brent] x=[" << xlo_F << "," << xup_F << "], tol=" << tol_F << "\n";
    std::cout << "[REFERENCE][Brent] x1=2, f1=1; x2=6, f2=0\n";
    std::cout << "[COMPUTED][Brent] " << std::setw(34) << "minimum of (x-2)^2+1" << "x=" << xmin1_F << ", f(x)=" << x2_1_func(xmin1_F) << "\n";
    std::cout << "[COMPUTED][Brent] " << std::setw(34) << "minimum of piecewise quadratic" << "x=" << xmin2_F << ", f(x)=" << x2_1_x6_0_func(xmin2_F) << "\n";
    assert(std::abs(xmin1_F - 2.0) < tol_F && std::abs(xmin2_F - 6.0) < tol_F);

    // f′(x)>0 ⇒ x*=x_min; golden→x_min, Brent→NaN.
    double xmin3_F = minimize_golden(exp_func, xlo_F, xup_F, tol_F);
    std::cout << "\n[INPUT][golden boundary] f(x)=exp(x), x=[" << xlo_F << "," << xup_F << "], tol=" << tol_F << "\n";
    std::cout << "[REFERENCE][golden boundary] x=0, f(x)=1\n";
    std::cout << "[COMPUTED][golden boundary] x=" << xmin3_F << ", f(x)=" << exp_func(xmin3_F) << "\n";
    assert(std::abs(xmin3_F - xlo_F) < tol_F);
    xmin3_F = minimize_brent(exp_func, xlo_F, xup_F, tol_F);
    std::cout << "\n[INPUT][Brent boundary] f(x)=exp(x), x=[" << xlo_F << "," << xup_F << "], tol=" << tol_F << "\n";
    std::cout << "[REFERENCE][Brent boundary] x=nan, f(x)=nan\n";
    std::cout << "[COMPUTED][Brent boundary] x=" << xmin3_F << ", f(x)=" << exp_func(xmin3_F) << "\n";
    assert(std::isnan(xmin3_F));
    return 0;
}
