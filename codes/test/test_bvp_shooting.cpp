/**
 * @file    test_bvp_shooting.cpp
 * @author  cailea
 * @date    2026-08-18
 * @brief   Acceptance test for Numerov and RK4 shooting methods.
 */

#include <cassert>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <Eigen/Dense>
#include "integration.hpp"
#include "bvp_shooting.hpp"
#include "spherical_radial_boundary.hpp"

/**
 * @brief  Compare Numerov and RK4 shooting for a spherical Woods-Saxon well.
 * @math   h_μ u''(r) = [V(r) - E]u(r)
 * @output Labeled reference and computed eigenstate diagnostics.
 */
int main() {
    int Nx_I = 1201;
    int l_I = 0;
    int node_I = 2;
    int xmatch_I = Nx_I / 2;
    double xmin_F = 0.05;
    double xmax_F = 20.0;
    double hmass_F = 197.329 * 197.329 / (2.0 * 938.0);
    double Ze2_F = 20.0 * (197.329 / 137.0);
    double tol_F = 1.0e-6;
    double V0_F = 80.0;
    double a0_F = 0.65;
    double R0_F = 1.25 * std::pow(40.0, 1.0 / 3.0);
    Eigen::VectorXd x_F1D_x = Eigen::VectorXd::LinSpaced(Nx_I, xmin_F, xmax_F);

    // V = V_WS + V_C + h_μ l(l + 1) / r².
    ShootingPotentialFunc V_Func = [=](double x_F) {
        double Vws_F = -V0_F / (1.0 + std::exp((x_F - R0_F) / a0_F));
        double ratio_F = x_F / R0_F;
        double Vc_F = x_F < R0_F ? Ze2_F * (3.0 - ratio_F * ratio_F) / (2.0 * R0_F) : Ze2_F / x_F;
        double l_F = static_cast<double>(l_I);
        return Vws_F + Vc_F + hmass_F * l_F * (l_F + 1.0) / (x_F * x_F);
    };

    // B_in = B_reg, B_out = (exp(-κx), -κ exp(-κx)).
    ShootingBoundaryFunc Bin_Func = [=](double x_F, double E_F) {return spherical_radial_boundary_regular(x_F, l_I, hmass_F, Ze2_F, E_F);};
    ShootingBoundaryFunc Bout_Func = [](double x_F, double E_F) {
        (void) E_F;
        double k_F = 0.6;
        doubleC y_C = std::exp(-k_F * x_F);
        return ShootingBoundary{y_C, -k_F * y_C};
    };

    // (V, B_in, B_out, x) → problem.
    ShootingProblem problem(V_Func, Bin_Func, Bout_Func, x_F1D_x, node_I, hmass_F, xmatch_I, tol_F);

    // E_ref = {E_Numerov, E_RK4}.
    double ErefNumerov_F = -8.929844e-1;
    double ErefRk4_F = -8.938847e-1;

    // problem → (E, R, y, N).
    ShootingSolution solutionNumerov = shooting_numerov(problem);
    ShootingSolution solutionRk4 = shooting_rk4(problem);
    double normNumerov_F = integrate_trapezoidal(solutionNumerov.y_C1D_x.cwiseAbs2(), x_F1D_x);
    double normRk4_F = integrate_trapezoidal(solutionRk4.y_C1D_x.cwiseAbs2(), x_F1D_x);

    // Input → reference → computed.
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "[Input] n = " << node_I << ", l = " << l_I << ", x = [" << xmin_F << ", " << xmax_F << "]\n";
    std::cout << "[Input] x = ";
    for (int x_I = 0; x_I < 5; ++x_I) {std::cout << x_F1D_x(x_I) << ", ";}
    std::cout << "..., ";
    for (int x_I = Nx_I - 5; x_I < Nx_I; ++x_I) {std::cout << x_F1D_x(x_I) << (x_I < Nx_I - 1 ? ", " : "\n");}
    std::cout << "[Input] potential (x, V):\n";
    std::cout << "  x  V\n";
    for (int x_I = 0; x_I < 5; ++x_I) {std::cout << "  " << x_F1D_x(x_I) << "  " << problem.V_F1D_x(x_I) << "\n";}
    std::cout << "  ...\n";
    for (int x_I = Nx_I - 5; x_I < Nx_I; ++x_I) {std::cout << "  " << x_F1D_x(x_I) << "  " << problem.V_F1D_x(x_I) << "\n";}
    std::cout << "[Reference] Numerov E = " << ErefNumerov_F << "\n";
    std::cout << "[Reference] RK4 E = " << ErefRk4_F << "\n";
    std::cout << "[Computed] Numerov E = " << solutionNumerov.E_F << ", |R| = " << std::abs(solutionNumerov.R_C) << ", norm = " << normNumerov_F << ", nodes = " << solutionNumerov.node_I << "\n";
    std::cout << "[Computed] RK4 E = " << solutionRk4.E_F << ", |R| = " << std::abs(solutionRk4.R_C) << ", norm = " << normRk4_F << ", nodes = " << solutionRk4.node_I << "\n";
    std::cout << "[Computed] Numerov y = ";
    for (int x_I = 0; x_I < 5; ++x_I) {std::cout << solutionNumerov.y_C1D_x(x_I) << ", ";}
    std::cout << "..., ";
    for (int x_I = Nx_I - 5; x_I < Nx_I; ++x_I) {std::cout << solutionNumerov.y_C1D_x(x_I) << (x_I < Nx_I - 1 ? ", " : "\n");}
    std::cout << "[Computed] RK4 y = ";
    for (int x_I = 0; x_I < 5; ++x_I) {std::cout << solutionRk4.y_C1D_x(x_I) << ", ";}
    std::cout << "..., ";
    for (int x_I = Nx_I - 5; x_I < Nx_I; ++x_I) {std::cout << solutionRk4.y_C1D_x(x_I) << (x_I < Nx_I - 1 ? ", " : "\n");}

    // |E - E_ref| < ε_E, |R| < ε_R, |Norm - 1| < ε_N, N = n.
    assert(std::abs(solutionNumerov.E_F - ErefNumerov_F) < 1.0e-5 && std::abs(solutionRk4.E_F - ErefRk4_F) < 1.0e-5);
    assert(std::abs(solutionNumerov.R_C) < tol_F && std::abs(solutionRk4.R_C) < tol_F);
    assert(std::abs(normNumerov_F - 1.0) < 1.0e-8 && std::abs(normRk4_F - 1.0) < 1.0e-8);
    assert(solutionNumerov.node_I == node_I && solutionRk4.node_I == node_I);
    return 0;
}
