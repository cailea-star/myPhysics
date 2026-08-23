/**
 * @file    test_spherical_radial_boundary.cpp
 * @author  cailea
 * @date    2026-08-18
 * @brief   Validate regular and Coulomb-Hankel spherical radial boundaries.
 */

#include <cassert>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include "spherical_radial_boundary.hpp"

/**
 * @brief  Validate spherical radial boundaries using analytic free-wave solutions.
 * @math   B_reg = (r, 1), H_0^(±) = exp(±iρ)
 * @output Labeled reference and computed boundary pairs.
 */
int main() {
    // (r,l,h_μ,Z e²,E,ε) → boundary inputs.
    double r_F = 2.0;
    int l_I = 0;
    double hmass_F = 1.0;
    double Ze2_F = 0.0;
    double E_F = 4.0;
    double tol_F = 1.0e-11;

    // k = √(E/h_μ), ρ = kr.
    double k_F = std::sqrt(E_F / hmass_F);
    double rho_F = k_F * r_F;

    // B_reg = (r,1).
    SphericalRadialBoundary Breg_ref_C1D_ydy = {doubleC(r_F, 0.0), doubleC(1.0, 0.0)};
    SphericalRadialBoundary Breg_C1D_ydy = spherical_radial_boundary_regular(r_F, l_I, hmass_F, Ze2_F, E_F);

    // H_0^(±) = exp(±iρ), dH_0^(±)/dr = ±ikH_0^(±).
    SphericalRadialBoundary Bplus_ref_C1D_ydy = {std::exp(doubleC(0.0, rho_F)), doubleC(0.0, k_F) * std::exp(doubleC(0.0, rho_F))};
    SphericalRadialBoundary Bminus_ref_C1D_ydy = {std::exp(doubleC(0.0, -rho_F)), doubleC(0.0, -k_F) * std::exp(doubleC(0.0, -rho_F))};
    SphericalRadialBoundary Bplus_C1D_ydy = spherical_radial_boundary_coulomb_hplus(r_F, l_I, hmass_F, Ze2_F, E_F);
    SphericalRadialBoundary Bminus_C1D_ydy = spherical_radial_boundary_coulomb_hminus(r_F, l_I, hmass_F, Ze2_F, E_F);

    // (r,l,h_μ,Ze²,E) → (B_ref,B).
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "[Input] r = " << r_F << ", l = " << l_I << ", h_mu = " << hmass_F << ", Ze2 = " << Ze2_F << ", E = " << E_F << "\n";
    std::cout << "[Reference] B_reg = " << Breg_ref_C1D_ydy.first << "  " << Breg_ref_C1D_ydy.second << "\n";
    std::cout << "[Reference] B_+ = " << Bplus_ref_C1D_ydy.first << "  " << Bplus_ref_C1D_ydy.second << "\n";
    std::cout << "[Reference] B_- = " << Bminus_ref_C1D_ydy.first << "  " << Bminus_ref_C1D_ydy.second << "\n";
    std::cout << "[Computed] B_reg = " << Breg_C1D_ydy.first << "  " << Breg_C1D_ydy.second << "\n";
    std::cout << "[Computed] B_+ = " << Bplus_C1D_ydy.first << "  " << Bplus_C1D_ydy.second << "\n";
    std::cout << "[Computed] B_- = " << Bminus_C1D_ydy.first << "  " << Bminus_C1D_ydy.second << "\n";

    // max |B - B_ref| < ε.
    assert(std::abs(Breg_C1D_ydy.first - Breg_ref_C1D_ydy.first) < tol_F && std::abs(Breg_C1D_ydy.second - Breg_ref_C1D_ydy.second) < tol_F);
    assert(std::abs(Bplus_C1D_ydy.first - Bplus_ref_C1D_ydy.first) < tol_F && std::abs(Bplus_C1D_ydy.second - Bplus_ref_C1D_ydy.second) < tol_F);
    assert(std::abs(Bminus_C1D_ydy.first - Bminus_ref_C1D_ydy.first) < tol_F && std::abs(Bminus_C1D_ydy.second - Bminus_ref_C1D_ydy.second) < tol_F);
    return 0;
}
