/**
 * @file    test_root_broyden.cpp
 * @author  cailea
 * @date    2026-05-14
 * @brief   Acceptance test for the modified Broyden fixed-point iterator.
 */

#include <Eigen/Core>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include "root_broyden.hpp"

/**
 * @brief  Return the target fixed point.
 * @math   x^* = (1.5, -0.75)^T
 * @output Target fixed-point vector.
 */
Eigen::Vector2d target_x() {
    Eigen::Vector2d xstar_F1D_i;
    xstar_F1D_i << 1.5, -0.75;
    return xstar_F1D_i;
}

/**
 * @brief  Evaluate a contractive nonlinear fixed-point map.
 * @math   G(x^*) = x^*
 * @output Updates G_F1D_i with G(x).
 * @note   Requires x_F1D_i.size() = 2.
 */
void fixed_point_func(const Eigen::VectorXd& x_F1D_i, Eigen::VectorXd& G_F1D_i) {
    Eigen::Vector2d xstar_F1D_i = target_x();
    Eigen::Vector2d y_F1D_i = x_F1D_i - xstar_F1D_i;
    G_F1D_i.resize(2);
    G_F1D_i(0) = xstar_F1D_i(0) + 0.30 * y_F1D_i(0) + 0.08 * std::sin(y_F1D_i(1));
    G_F1D_i(1) = xstar_F1D_i(1) + 0.10 * std::sin(y_F1D_i(0)) + 0.25 * y_F1D_i(1);
}

/**
 * @brief  Test the modified Broyden fixed-point iteration.
 * @math   r = G(x) - x
 * @output Prints labeled inputs and outputs and asserts convergence accuracy.
 */
int main() {
    int Nh_I = 2;
    int Nitermax_I = 30;
    double alpha_F = 1.0;
    double tol_F = 1.0e-10;
    Eigen::Vector2d x0_F1D_i;
    x0_F1D_i << 4.0, -3.0;
    Eigen::VectorXd G0_F1D_i(2);
    fixed_point_func(x0_F1D_i, G0_F1D_i);
    BroydenIterator solver(Nh_I, fixed_point_func, alpha_F, x0_F1D_i, G0_F1D_i);
    double res_F = solver.rnext_F1D_i.array().abs().maxCoeff();
    int iter_I = 1;
    while (res_F > tol_F && iter_I < Nitermax_I) {
        res_F = solver.iterate(fixed_point_func, alpha_F);
        ++iter_I;
    }
    Eigen::Vector2d xstar_F1D_i = target_x();
    Eigen::VectorXd err_F1D_i = solver.xnext_F1D_i - xstar_F1D_i;
    std::cout << std::scientific << std::setprecision(6) << std::left;
    std::cout << "Input:\n";
    std::cout << std::setw(18) << "history length" << Nh_I << "\n";
    std::cout << std::setw(18) << "max iterations" << Nitermax_I << "\n";
    std::cout << std::setw(18) << "alpha" << alpha_F << "\n";
    std::cout << std::setw(18) << "tolerance" << tol_F << "\n";
    std::cout << std::setw(18) << "x0 input[0]" << x0_F1D_i(0) << "\n";
    std::cout << std::setw(18) << "x0 input[1]" << x0_F1D_i(1) << "\n";
    std::cout << std::setw(18) << "G0 input[0]" << G0_F1D_i(0) << "\n";
    std::cout << std::setw(18) << "G0 input[1]" << G0_F1D_i(1) << "\n";
    std::cout << "Reference:\n";
    std::cout << std::setw(18) << "x0 exact" << xstar_F1D_i(0) << "\n";
    std::cout << std::setw(18) << "x1 exact" << xstar_F1D_i(1) << "\n";
    std::cout << "Computed:\n";
    std::cout << "Broyden Fixed Point Results:\n";
    std::cout << std::setw(18) << "iterations" << iter_I << "\n";
    std::cout << std::setw(18) << "residual max" << res_F << "\n";
    std::cout << std::setw(18) << "x0" << solver.xnext_F1D_i(0) << "\n";
    std::cout << std::setw(18) << "x1" << solver.xnext_F1D_i(1) << "\n";
    std::cout << std::setw(18) << "error norm" << err_F1D_i.norm() << "\n";
    assert(res_F <= tol_F);
    assert(err_F1D_i.lpNorm<Eigen::Infinity>() <= tol_F);
    return 0;
}
