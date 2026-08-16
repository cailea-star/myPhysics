/**
 * @file    root_broyden_test.cpp
 * @author  cailea
 * @date    2026-05-14
 * @brief   Manual fixed-point test for the modified Broyden iterator.
 */

#include <Eigen/Core>
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
 */
void fixed_point_func(const Eigen::VectorXd& x_F1D_i, Eigen::VectorXd& G_F1D_i) {
    Eigen::Vector2d xstar_F1D_i = target_x();
    Eigen::Vector2d y_F1D_i = x_F1D_i - xstar_F1D_i;
    G_F1D_i.resize(2);
    G_F1D_i(0) = xstar_F1D_i(0) + 0.30 * y_F1D_i(0) + 0.08 * std::sin(y_F1D_i(1));
    G_F1D_i(1) = xstar_F1D_i(1) + 0.10 * std::sin(y_F1D_i(0)) + 0.25 * y_F1D_i(1);
}

/**
 * @brief  Run and print the modified Broyden fixed-point test.
 * @math   r = G(x) - x
 * @output Prints the iteration count, residual, solution, and error.
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
    std::cout << "Broyden Fixed Point Results:\n";
    std::cout << std::setw(18) << "iterations" << iter_I << "\n";
    std::cout << std::setw(18) << "residual max" << res_F << "\n";
    std::cout << std::setw(18) << "x0" << solver.xnext_F1D_i(0) << "\n";
    std::cout << std::setw(18) << "x1" << solver.xnext_F1D_i(1) << "\n";
    std::cout << std::setw(18) << "x0 exact" << xstar_F1D_i(0) << "\n";
    std::cout << std::setw(18) << "x1 exact" << xstar_F1D_i(1) << "\n";
    std::cout << std::setw(18) << "error norm" << err_F1D_i.norm() << "\n";
    return 0;
}
