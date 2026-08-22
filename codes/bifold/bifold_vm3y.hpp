/**
 * @file    bifold_vm3y.hpp
 * @author  cailea
 * @date    2026-08-22
 * @brief   M3Y-Reid and M3Y-Paris interaction kernels.
 */

#pragma once

#include <cassert>
#include <cmath>
#include <functional>
#include <numbers>
#include <utility>

/**
 * @brief  Evaluate the coordinate-space Coulomb kernel analytically.
 * @math   v_C(r) = Ze² / r
 * @output Coordinate-space Coulomb interaction.
 * @note   Ze² is the supplied charge coupling; requires r > 0.
 */
inline double v_coulomb_r(double r_F, double Ze2_F) {
    // Evaluate outside the Coulomb singularity.
    assert(r_F > 0.0);
    return Ze2_F / r_F;
}

/**
 * @brief  Evaluate the momentum-space Coulomb kernel analytically.
 * @math   ṽ_C(k) = 4πZe² / k²
 * @output Momentum-space Coulomb interaction.
 * @note   Ze² is the supplied charge coupling; requires k > 0.
 */
inline double v_coulomb_k(double k_F, double Ze2_F) {
    // Evaluate outside the zero-momentum singularity.
    assert(k_F > 0.0);
    return 4.0 * std::numbers::pi * Ze2_F / (k_F * k_F);
}

/**
 * @brief  Evaluate a coordinate-space Yukawa kernel analytically.
 * @math   v(r) = V₀ exp(-μr) / (μr)
 * @output Coordinate-space Yukawa interaction.
 * @note   Requires r > 0 and μ > 0.
 */
inline double v_yukawa_r(double r_F, double V0_F, double mu_F) {
    // Evaluate outside the integrable origin singularity.
    assert(r_F > 0.0);
    assert(mu_F > 0.0);
    return V0_F * std::exp(-mu_F * r_F) / (mu_F * r_F);
}

/**
 * @brief  Evaluate the analytic Fourier transform of a Yukawa kernel.
 * @math   ṽ(k) = 4πV₀ / [μ(k² + μ²)]
 * @output Momentum-space Yukawa interaction.
 * @note   Requires k ≥ 0 and μ > 0.
 */
inline double v_yukawa_k(double k_F, double V0_F, double mu_F) {
    // Transform the central kernel with the spherical Fourier convention.
    assert(k_F >= 0.0);
    assert(mu_F > 0.0);
    return 4.0 * std::numbers::pi * V0_F / (mu_F * (k_F * k_F + mu_F * mu_F));
}

/**
 * @brief  Evaluate a coordinate-space Gaussian kernel analytically.
 * @math   v(r) = V₀ exp[-(r/r₀)²]
 * @output Coordinate-space Gaussian interaction.
 * @note   Requires r ≥ 0 and r₀ > 0.
 */
inline double v_gaussian_r(double r_F, double V0_F, double r0_F) {
    // Evaluate the radial Gaussian profile.
    assert(r_F >= 0.0);
    assert(r0_F > 0.0);
    return V0_F * std::exp(-std::pow(r_F / r0_F, 2));
}

/**
 * @brief  Evaluate the analytic Fourier transform of a Gaussian kernel.
 * @math   ṽ(k) = π^(3/2) V₀ r₀³ exp[-(kr₀/2)²]
 * @output Momentum-space Gaussian interaction.
 * @note   Requires k ≥ 0 and r₀ > 0.
 */
inline double v_gaussian_k(double k_F, double V0_F, double r0_F) {
    // Transform the central Gaussian with the spherical Fourier convention.
    assert(k_F >= 0.0);
    assert(r0_F > 0.0);
    return std::pow(std::numbers::pi * r0_F * r0_F, 1.5) * V0_F * std::exp(-std::pow(k_F * r0_F / 2.0, 2));
}

/**
 * @brief  Evaluate the legacy cutoff proxy for a coordinate-space delta kernel.
 * @math   v_δ(r) = V₀ for r < 10⁻⁶ fm; otherwise v_δ(r) = 0
 * @output Cutoff representation of the coordinate-space delta interaction.
 * @note   This pointwise proxy is not a normalized Dirac distribution.
 */
inline double v_delta_r(double r_F, double V0_F) {
    // Apply the reference cutoff through a branch-free mask.
    assert(r_F >= 0.0);
    return V0_F * static_cast<double>(r_F < 1.0e-6);
}

/**
 * @brief  Evaluate the legacy momentum-space delta-kernel strength.
 * @math   ṽ_δ = 4πV₀
 * @output Momentum-space delta-interaction strength.
 */
inline double v_delta_k(double V0_F) {
    // Preserve the reference-code normalization.
    return 4.0 * std::numbers::pi * V0_F;
}

/**
 * @brief  Evaluate the M3Y-Reid direct interaction in momentum space analytically.
 * @math   ṽ_D(k) = 4π[7999/(4(k²+4²)) - 2134.25/(2.5(k²+2.5²))]
 * @output Momentum-space M3Y-Reid direct interaction.
 */
inline double v_reid_d_k(double k_F) {
    // Sum Yukawa components from shorter to longer range.
    double term1_F = v_yukawa_k(k_F, 7999.0, 4.0);
    double term2_F = v_yukawa_k(k_F, -2134.25, 2.5);
    return term1_F + term2_F;
}

/**
 * @brief  Evaluate the finite-range M3Y-Reid exchange interaction analytically.
 * @math   v_E(r) = 4631.38e^(-4r)/(4r) - 1787.13e^(-2.5r)/(2.5r) - 7.8474e^(-0.7072r)/(0.7072r)
 * @output Coordinate-space M3Y-Reid exchange interaction.
 * @note   Requires r > 0.
 */
inline double v_reid_ex_s(double r_F) {
    // Sum Yukawa components from shorter to longer range.
    double term1_F = v_yukawa_r(r_F, 4631.38, 4.0);
    double term2_F = v_yukawa_r(r_F, -1787.13, 2.5);
    double term3_F = v_yukawa_r(r_F, -7.8474, 0.7072);
    return term1_F + term2_F + term3_F;
}

/**
 * @brief  Evaluate the M3Y-Paris direct interaction in momentum space analytically.
 * @math   ṽ_D(k) = 4π[11061.625/(4(k²+4²)) - 2537.5/(2.5(k²+2.5²))]
 * @output Momentum-space M3Y-Paris direct interaction.
 */
inline double v_paris_d_k(double k_F) {
    // Sum Yukawa components from shorter to longer range.
    double term1_F = v_yukawa_k(k_F, 11061.625, 4.0);
    double term2_F = v_yukawa_k(k_F, -2537.5, 2.5);
    return term1_F + term2_F;
}

/**
 * @brief  Evaluate the finite-range M3Y-Paris exchange interaction analytically.
 * @math   v_E(r) = -1524.25e^(-4r)/(4r) - 518.75e^(-2.5r)/(2.5r) - 7.8474e^(-0.7072r)/(0.7072r)
 * @output Coordinate-space M3Y-Paris exchange interaction.
 * @note   Requires r > 0.
 */
inline double v_paris_ex_s(double r_F) {
    // Sum Yukawa components from shorter to longer range.
    double term1_F = v_yukawa_r(r_F, -1524.25, 4.0);
    double term2_F = v_yukawa_r(r_F, -518.75, 2.5);
    double term3_F = v_yukawa_r(r_F, -7.8474, 0.7072);
    return term1_F + term2_F + term3_F;
}

using Real2RealFunc = std::function<double(double)>;
using RealReal2RealFunc = std::function<double(double, double)>;
using RealInt2RealFunc = std::function<double(double, int)>;

/**
 * @brief Bundle the kernels and energy factor for one M3Y parameterization.
 */
class VnnFunctions {
public:
    RealReal2RealFunc vcoul_k_Func;
    Real2RealFunc vnnd_k_Func;
    Real2RealFunc vnne_r_Func;
    RealInt2RealFunc vnne_zr_Func;
    RealInt2RealFunc gE_Func;

    /**
     * @brief  Construct an empty interaction bundle.
     * @math   𝒱 = ∅
     * @output Empty callable bundle.
     */
    VnnFunctions() = default;

    /**
     * @brief  Construct an interaction bundle from callable kernels.
     * @math   𝒱 = {v_C, v_D, v_E, v_E^ZR, g(E)}
     * @output Initialized callable bundle.
     */
    VnnFunctions(RealReal2RealFunc vcoul_k_Func_, Real2RealFunc vnnd_k_Func_, Real2RealFunc vnne_r_Func_, RealInt2RealFunc vnne_zr_Func_, RealInt2RealFunc gE_Func_)
        : vcoul_k_Func(std::move(vcoul_k_Func_)), vnnd_k_Func(std::move(vnnd_k_Func_)), vnne_r_Func(std::move(vnne_r_Func_)), vnne_zr_Func(std::move(vnne_zr_Func_)), gE_Func(std::move(gE_Func_)) {}
};

/**
 * @brief  Evaluate the linear M3Y-Reid beam-energy factor.
 * @math   g(E) = 1 - 0.002 E_lab/A_p
 * @output M3Y-Reid energy factor.
 * @note   Requires E_lab ≥ 0 and A_p > 0.
 */
inline double reid_gE(double Elab_F, int Ap_I) {
    // Apply the Reid energy correction per projectile nucleon.
    assert(Elab_F >= 0.0);
    assert(Ap_I > 0);
    return 1.0 - 2.0e-3 * Elab_F / Ap_I;
}

/**
 * @brief  Evaluate the linear M3Y-Paris beam-energy factor.
 * @math   g(E) = 1 - 0.003 E_lab/A_p
 * @output M3Y-Paris energy factor.
 * @note   Requires E_lab ≥ 0 and A_p > 0.
 */
inline double paris_gE(double Elab_F, int Ap_I) {
    // Apply the Paris energy correction per projectile nucleon.
    assert(Elab_F >= 0.0);
    assert(Ap_I > 0);
    return 1.0 - 3.0e-3 * Elab_F / Ap_I;
}

/**
 * @brief  Evaluate the legacy M3Y-Reid zero-range exchange strength.
 * @math   J₀₀(E) = -276[1 - 0.005 E_lab/A_p]
 * @output M3Y-Reid zero-range exchange strength.
 * @note   Requires E_lab ≥ 0 and A_p > 0.
 */
inline double v_reid_ex_zr(double Elab_F, int Ap_I) {
    // Apply the Reid zero-range energy correction.
    assert(Elab_F >= 0.0);
    assert(Ap_I > 0);
    return -276.0 * (1.0 - 5.0e-3 * Elab_F / Ap_I);
}

/**
 * @brief  Evaluate the legacy M3Y-Paris zero-range exchange strength.
 * @math   J₀₀(E) = -590[1 - 0.002 E_lab/A_p]
 * @output M3Y-Paris zero-range exchange strength.
 * @note   Requires E_lab ≥ 0 and A_p > 0.
 */
inline double v_paris_ex_zr(double Elab_F, int Ap_I) {
    // Apply the Paris zero-range energy correction.
    assert(Elab_F >= 0.0);
    assert(Ap_I > 0);
    return -590.0 * (1.0 - 2.0e-3 * Elab_F / Ap_I);
}

// Bind each parameterization to its interaction kernels.
inline const VnnFunctions reid_functions(v_coulomb_k, v_reid_d_k, v_reid_ex_s, v_reid_ex_zr, reid_gE);
inline const VnnFunctions paris_functions(v_coulomb_k, v_paris_d_k, v_paris_ex_s, v_paris_ex_zr, paris_gE);
