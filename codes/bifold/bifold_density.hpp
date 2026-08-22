/**
 * @file    bifold_density.hpp
 * @author  cailea
 * @date    2026-08-22
 * @brief   Projectile and target density profiles for double-folding calculations.
 */

#pragma once

#include <cassert>
#include <cmath>
#include <functional>
#include <numbers>
#include <utility>

#include <gsl/gsl_sf_legendre.h>

#include "integration.hpp"

using Real2RealFunc = std::function<double(double)>;
using RealReal2RealFunc = std::function<double(double, double)>;

/**
 * @brief Store a normalized spherical projectile density and its nuclear data.
 */
class DensityProjec {
public:
    int Ap_I = 0; // Projectile mass number.
    int Zp_I = 0; // Projectile proton number.
    double rmax_F = 0.0; // Smallest 0.1-grid cutoff enclosing 99.9% of the density.
    Real2RealFunc rho_Func; // 4π ∫₀^∞ r_p² ρ_p(r_p) dr_p = 1.

    /**
     * @brief  Construct an empty projectile-density bundle.
     * @math   𝓡_p = ∅
     * @output Empty projectile-density bundle.
     */
    DensityProjec() = default;

    /**
     * @brief  Construct a projectile-density bundle using a trapezoidal radial-cutoff search.
     * @math   4π ∫₀^∞ r_p² ρ_p(r_p) dr_p = 1
     * @output Initialized projectile-density bundle.
     * @note   Requires A_p > 0, 0 <= Z_p <= A_p, and a normalized density reaching 0.999.
     */
    DensityProjec(int Ap_I_, int Zp_I_, Real2RealFunc rho_Func_)
        : Ap_I(Ap_I_), Zp_I(Zp_I_), rmax_F(0.0), rho_Func(std::move(rho_Func_)) {
        assert(Ap_I > 0);
        assert(Zp_I >= 0);
        assert(Zp_I <= Ap_I);
        assert(static_cast<bool>(rho_Func));

        // 4π ∫₀^r_max r² ρ_p(r) dr → 0.999.
        double dr_F = 0.1; // Radial search step.
        double integral_F = 0.0; // Enclosed density normalization.
        auto radial_integrand_Func = [&](double r_F) {
            return rho_Func(r_F) * r_F * r_F * 4.0 * std::numbers::pi;
        };
        while (integral_F < 0.999) {
            integral_F += 0.5 * (radial_integrand_Func(rmax_F) + radial_integrand_Func(rmax_F + dr_F)) * dr_F;
            rmax_F += dr_F;
        }
    }

    /**
     * @brief  Evaluate the projectile matter density.
     * @math   ρ_{p,m}(r_p) = A_p ρ_p(r_p)
     * @output Projectile matter density.
     * @note   Requires r_p >= 0.
     */
    double rho_matter(double r_F) const {
        assert(r_F >= 0.0);
        assert(static_cast<bool>(rho_Func));
        return rho_Func(r_F) * Ap_I;
    }

    /**
     * @brief  Evaluate the projectile charge density.
     * @math   ρ_{p,c}(r_p) = Z_p ρ_p(r_p)
     * @output Projectile charge density.
     * @note   Requires r_p >= 0.
     */
    double rho_charge(double r_F) const {
        assert(r_F >= 0.0);
        assert(static_cast<bool>(rho_Func));
        return rho_Func(r_F) * Zp_I;
    }
};

/**
 * @brief Store a normalized angular target density and its nuclear data.
 */
class DensityTarget {
public:
    int At_I = 0; // Target mass number.
    int Zt_I = 0; // Target proton number.
    double rmax_F = 0.0; // Smallest 0.1-grid cutoff enclosing 99.9% of the density.
    RealReal2RealFunc rho_Func; // 2π ∫₀^∞ r_t² dr_t ∫₀^π ρ_t(r_t, θ_t) sin(θ_t) dθ_t = 1.

    /**
     * @brief  Construct an empty target-density bundle.
     * @math   𝓡_t = ∅
     * @output Empty target-density bundle.
     */
    DensityTarget() = default;

    /**
     * @brief  Construct a target-density bundle using spherical and trapezoidal radial-cutoff quadrature.
     * @math   2π ∫₀^∞ r_t² dr_t ∫₀^π ρ_t(r_t, θ_t) sin(θ_t) dθ_t = 1
     * @output Initialized target-density bundle.
     * @note   Requires A_t > 0, 0 <= Z_t <= A_t, and a normalized density reaching 0.999.
     */
    DensityTarget(int At_I_, int Zt_I_, RealReal2RealFunc rho_Func_)
        : At_I(At_I_), Zt_I(Zt_I_), rmax_F(0.0), rho_Func(std::move(rho_Func_)) {
        assert(At_I > 0);
        assert(Zt_I >= 0);
        assert(Zt_I <= At_I);
        assert(static_cast<bool>(rho_Func));

        // 2π ∫₀^r_max r² dr ∫₀^π ρ_t(r, θ) sin(θ) dθ → 0.999.
        double dr_F = 0.1; // Radial search step.
        double integral_F = 0.0; // Enclosed density normalization.
        auto radial_integrand_Func = [&](double r_F) {
            auto angular_integrand_Func = [&](double theta_F, double) {
                return rho_Func(r_F, theta_F);
            };
            return integrate_spherical(angular_integrand_Func, 7, 1) * r_F * r_F;
        };
        while (integral_F < 0.999) {
            integral_F += 0.5 * (radial_integrand_Func(rmax_F) + radial_integrand_Func(rmax_F + dr_F)) * dr_F;
            rmax_F += dr_F;
        }
    }

    /**
     * @brief  Evaluate the target matter density.
     * @math   ρ_{t,m}(r_t, θ_t) = A_t ρ_t(r_t, θ_t)
     * @output Target matter density.
     * @note   Requires r_t >= 0.
     */
    double rho_matter(double r_F, double theta_F) const {
        assert(r_F >= 0.0);
        assert(static_cast<bool>(rho_Func));
        return rho_Func(r_F, theta_F) * At_I;
    }

    /**
     * @brief  Evaluate the target charge density.
     * @math   ρ_{t,c}(r_t, θ_t) = Z_t ρ_t(r_t, θ_t)
     * @output Target charge density.
     * @note   Requires r_t >= 0.
     */
    double rho_charge(double r_F, double theta_F) const {
        assert(r_F >= 0.0);
        assert(static_cast<bool>(rho_Func));
        return rho_Func(r_F, theta_F) * Zt_I;
    }
};

/**
 * @brief Represent a normalized spherical Gaussian density.
 */
class DensityGaussian {
protected:
    double norm_F = 1.0; // Numerical density normalization.

public:
    double a_F = 0.0; // Gaussian width.

    /**
     * @brief  Normalize a Gaussian density with adaptive Gauss-Kronrod QAG quadrature.
     * @math   ρ(r) ∝ exp[-(r/a)²]
     * @output Normalized Gaussian density profile.
     * @note   Requires a > 0 and uses the legacy radial interval [0, 100].
     */
    explicit DensityGaussian(double a_F_)
        : norm_F(1.0), a_F(a_F_) {
        assert(a_F > 0.0);

        // 𝒩 = 4π ∫₀¹⁰⁰ r² exp[-(r/a)²] dr.
        auto radial_integrand_Func = [&](double r_F) {
            return density(r_F) * r_F * r_F;
        };
        norm_F = integrate_qag(radial_integrand_Func, 0.0, 100.0) * 4.0 * std::numbers::pi;
    }

    /**
     * @brief  Evaluate the normalized Gaussian density.
     * @math   ρ(r) = exp[-(r/a)²] / 𝒩
     * @output Normalized Gaussian density.
     * @note   Requires r >= 0.
     */
    double density(double r_F) const {
        assert(r_F >= 0.0);
        return std::exp(-(r_F * r_F) / (a_F * a_F)) / norm_F;
    }
};

/**
 * @brief Represent a normalized axially deformed Fermi density.
 */
class DensityFermi {
protected:
    double norm_F = 1.0; // Numerical density normalization.

public:
    double R0_F = 0.0; // Reference nuclear radius.
    double a0_F = 0.0; // Surface diffuseness.
    double beta2_F = 0.0; // Quadrupole deformation.
    double beta4_F = 0.0; // Hexadecapole deformation.

    /**
     * @brief  Normalize a deformed Fermi density with Gauss-Legendre, trapezoidal, and adaptive Gauss-Kronrod quadrature.
     * @math   ρ(r, θ) ∝ {1 + exp[(r - R(θ))/a₀]}⁻¹
     * @output Normalized Fermi density profile.
     * @note   Uses the legacy radial interval [0, 100] and the default spherical grid.
     */
    DensityFermi(double R0_F_, double a0_F_, double beta2_F_ = 0.0, double beta4_F_ = 0.0)
        : norm_F(1.0), R0_F(R0_F_), a0_F(a0_F_), beta2_F(beta2_F_), beta4_F(beta4_F_) {
        assert(R0_F > 0.0);
        assert(a0_F > 0.0);

        // 𝒩 = 2π ∫₀¹⁰⁰ r² dr ∫₀^π ρ(r, θ) sin(θ) dθ.
        auto radial_integrand_Func = [&](double r_F) {
            auto angular_integrand_Func = [&](double theta_F, double) {
                return density(r_F, theta_F);
            };
            return integrate_spherical(angular_integrand_Func, 7, 1) * r_F * r_F;
        };
        norm_F = integrate_qag(radial_integrand_Func, 0.0, 100.0);
    }

    /**
     * @brief  Evaluate the deformed Fermi density with GSL normalized spherical associated Legendre functions.
     * @math   R(θ) = R₀[1 + β₂Y₂₀(θ) + β₄Y₄₀(θ)]; ρ(r, θ) = {1 + exp[(r - R(θ))/a₀]}⁻¹ / 𝒩
     * @output Normalized Fermi density.
     * @note   Requires r >= 0 and 0 <= θ <= π.
     */
    double density(double r_F, double theta_F) const {
        assert(r_F >= 0.0);
        assert(theta_F >= 0.0 && theta_F <= std::numbers::pi);

        // R(θ) → ρ(r, θ).
        double Y20_F = gsl_sf_legendre_sphPlm(2, 0, std::cos(theta_F));
        double Y40_F = gsl_sf_legendre_sphPlm(4, 0, std::cos(theta_F));
        double R_F = R0_F * (1.0 + beta2_F * Y20_F + beta4_F * Y40_F);
        return 1.0 / ((1.0 + std::exp((r_F - R_F) / a0_F)) * norm_F);
    }
};
