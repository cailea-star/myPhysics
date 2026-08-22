/**
 * @file    bifold_fermik.hpp
 * @author  cailea
 * @date    2026-08-22
 * @brief   Local Fermi momentum and reduced spherical Bessel factor.
 */

#pragma once

#include <cassert>
#include <cmath>
#include <functional>
#include <numbers>

#include "derivative.hpp"

using RealToRealFunc = std::function<double(double)>;
using RealRealRealToRealFunc = std::function<double(double, double, double)>;

/**
 * @brief  Compute an angular-density local Fermi momentum using fourth-order five-point central differences.
 * @math   k_F² = (3π²ρ/2)^(2/3) + 5C_s|∇ρ|²/(3ρ²) + 5∇²ρ/(36ρ)
 * @output Local Fermi momentum.
 * @note   Uses fixed coordinate steps and clamps r to 4Δr before differentiating.
 */
inline double calc_fermik(const RealRealRealToRealFunc& rho_Func, double r_F, double theta_F, double phi_F, double Cs_F = 1.0 / 36.0) {
    assert(static_cast<bool>(rho_Func));
    assert(std::isfinite(r_F));
    assert(r_F >= 0.0);
    assert(std::isfinite(theta_F));
    assert(theta_F >= 0.0 && theta_F <= std::numbers::pi);
    assert(std::isfinite(phi_F));
    assert(std::isfinite(Cs_F));
    assert(Cs_F >= 0.0);

    // Δr = 10⁻³ fm; Δθ = 10⁻³π rad; Δφ = 2 × 10⁻³π rad; r → max(r, 4Δr).
    double dr_F = 1.0e-3; // Radial step in fm.
    double dtheta_F = 1.0e-3 * std::numbers::pi; // Polar-angle step in rad.
    double dphi_F = 2.0e-3 * std::numbers::pi; // Azimuthal-angle step in rad.
    double denominator_epsilon_F = 1.0e-10; // Denominator regularizer.
    r_F = std::fmax(r_F, 4.0 * dr_F);

    // ρ_r(r′) = ρ(r′, θ, φ); ρ_θ(θ′) = ρ(r, θ′, φ); ρ_φ(φ′) = ρ(r, θ, φ′).
    auto rho_r_Func = [&](double r_F) { return rho_Func(r_F, theta_F, phi_F); };
    auto rho_theta_Func = [&](double theta_F) { return rho_Func(r_F, theta_F, phi_F); };
    auto rho_phi_Func = [&](double phi_F) { return rho_Func(r_F, theta_F, phi_F); };
    double rho_F = rho_Func(r_F, theta_F, phi_F); // Local density.
    assert(std::isfinite(rho_F));
    assert(rho_F >= 0.0);

    // ∂ρ, ∂²ρ ← D_h[ρ].
    double drho_dr_F = derivative1(rho_r_Func, r_F, dr_F);
    double d2rho_dr2_F = derivative2(rho_r_Func, r_F, dr_F);
    double drho_dtheta_F = derivative1(rho_theta_Func, theta_F, dtheta_F);
    double d2rho_dtheta2_F = derivative2(rho_theta_Func, theta_F, dtheta_F);
    double drho_dphi_F = derivative1(rho_phi_Func, phi_F, dphi_F);
    double d2rho_dphi2_F = derivative2(rho_phi_Func, phi_F, dphi_F);

    // (∂_rρ, ∂_θρ, ∂_φρ) → (|∇ρ|², ∇²ρ).
    double grad_rho2_F = drho_dr_F * drho_dr_F + drho_dtheta_F * drho_dtheta_F / (r_F * r_F + denominator_epsilon_F) + drho_dphi_F * drho_dphi_F / (r_F * r_F * std::sin(theta_F) * std::sin(theta_F) + denominator_epsilon_F);
    double laplacian_rho_F = d2rho_dr2_F + 2.0 * drho_dr_F / (r_F + denominator_epsilon_F) + (d2rho_dtheta2_F + drho_dtheta_F / (std::tan(theta_F) + denominator_epsilon_F)) / (r_F * r_F + denominator_epsilon_F) + d2rho_dphi2_F / (r_F * r_F * std::sin(theta_F) * std::sin(theta_F) + denominator_epsilon_F);

    // k_F² = term₁ + term₂ + term₃ → max(k_F², 0).
    double term1_F = std::pow(1.5 * std::numbers::pi * std::numbers::pi * rho_F, 2.0 / 3.0);
    double term2_F = 5.0 * grad_rho2_F * Cs_F / (3.0 * rho_F * rho_F + denominator_epsilon_F);
    double term3_F = 5.0 * laplacian_rho_F / (36.0 * rho_F + denominator_epsilon_F);
    double kF2_F = std::fmax(term1_F + term2_F + term3_F, 0.0); // Nonnegative squared Fermi momentum.
    return std::sqrt(kF2_F);
}

/**
 * @brief  Compute a radial-density local Fermi momentum using fourth-order five-point central differences.
 * @math   k_F² = (3π²ρ/2)^(2/3) + 5C_s(∂ρ/∂r)²/(3ρ²) + 5[∂²ρ/∂r² + 2(∂ρ/∂r)/r]/(36ρ)
 * @output Local Fermi momentum.
 * @note   Uses a fixed radial step and clamps r to 4Δr before differentiating.
 */
inline double calc_fermik(const RealToRealFunc& rho_Func, double r_F, double Cs_F = 1.0 / 36.0) {
    assert(static_cast<bool>(rho_Func));
    assert(std::isfinite(r_F));
    assert(r_F >= 0.0);
    assert(std::isfinite(Cs_F));
    assert(Cs_F >= 0.0);

    // Δr = 10⁻³ fm; r → max(r, 4Δr); ρ → ρ(r).
    double dr_F = 1.0e-3; // Radial step in fm.
    double denominator_epsilon_F = 1.0e-10; // Denominator regularizer.
    r_F = std::fmax(r_F, 4.0 * dr_F);
    auto rho_r_Func = [&](double r_F) { return rho_Func(r_F); };
    double rho_F = rho_Func(r_F); // Local density.
    assert(std::isfinite(rho_F));
    assert(rho_F >= 0.0);

    // (∂_rρ, ∂²_rρ) → (|∇ρ|², ∇²ρ).
    double drho_dr_F = derivative1(rho_r_Func, r_F, dr_F);
    double d2rho_dr2_F = derivative2(rho_r_Func, r_F, dr_F);
    double grad_rho2_F = drho_dr_F * drho_dr_F;
    double laplacian_rho_F = d2rho_dr2_F + 2.0 * drho_dr_F / (r_F + denominator_epsilon_F);

    // k_F² = term₁ + term₂ + term₃ → max(k_F², 0).
    double term1_F = std::pow(1.5 * std::numbers::pi * std::numbers::pi * rho_F, 2.0 / 3.0);
    double term2_F = 5.0 * grad_rho2_F * Cs_F / (3.0 * rho_F * rho_F + denominator_epsilon_F);
    double term3_F = 5.0 * laplacian_rho_F / (36.0 * rho_F + denominator_epsilon_F);
    double kF2_F = std::fmax(term1_F + term2_F + term3_F, 0.0); // Nonnegative squared Fermi momentum.
    return std::sqrt(kF2_F);
}

/**
 * @brief  Evaluate the reduced spherical Bessel factor using a fourth-order small-argument expansion.
 * @math   ĵ₁(x) = 3j₁(x)/x = 3[sin(x) - x cos(x)]/x³; ĵ₁(0) = 1
 * @output Reduced spherical Bessel factor.
 * @note   Uses 1 - x²/10 + x⁴/280 for |x| <= 10⁻².
 */
inline double calc_hatj1(double x_F) {
    assert(std::isfinite(x_F));

    // ĵ₁ = (1 - m)ĵ₁^series + mĵ₁^exact; m = [|x| > 10⁻²].
    double x_series_max_F = 1.0e-2; // Series-domain upper bound for |x|.
    double hat_j1_series_F = 1.0 - std::pow(x_F, 2) / 10.0 + std::pow(x_F, 4) / 280.0;
    double use_exact_F = static_cast<double>(std::abs(x_F) > x_series_max_F); // Exact-expression mask.
    double x3_safe_F = std::pow(x_F, 3) + 1.0 - use_exact_F; // Nonzero exact-form denominator.
    double hat_j1_exact_F = 3.0 * (std::sin(x_F) - x_F * std::cos(x_F)) / x3_safe_F;
    return hat_j1_series_F + use_exact_F * (hat_j1_exact_F - hat_j1_series_F);
}
