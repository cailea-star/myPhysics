/**
 * @file    spherical_fourier.hpp
 * @author  cailea
 * @date    2026-04-09
 * @brief   Spherical Fourier transform utilities.
 */

#pragma once

#include <cassert>
#include <cmath>
#include <functional>
#include <Eigen/Core>
#include <gsl/gsl_sf_bessel.h>
#include "integration.hpp"

using RealToRealFunc = std::function<double(double)>;

/**
 * @brief  Radial-to-momentum transform by GSL QAG (61-point Gauss-Kronrod).
 * @math   F_l(k) = 4π ∫_{r_min}^{r_max} r² f_l(r) j_l(kr) dr
 * @output F_l(k).
 * @note   Finite [r_min, r_max] approximates [0, ∞).
 */
inline double spherical_FT_r2k(const RealToRealFunc& fl_Func, double k_F, int l_I = 0, double rmin_F = 1e-6, double rmax_F = 1e2) {
    // f_l(r) → r²f_l(r)j_l(kr) → F_l(k).
    assert(l_I >= 0);
    assert(rmin_F < rmax_F);
    auto integrand_Func = [&](double r_F) {
        return r_F * r_F * fl_Func(r_F) * gsl_sf_bessel_jl(l_I, k_F * r_F);
    };
    double integral_F = integrate_qag(integrand_Func, rmin_F, rmax_F);
    return 4.0 * M_PI * integral_F;
}

/**
 * @brief  Radial-to-momentum transform by composite Simpson quadrature.
 * @math   F_l(k) ≈ 4π Q_S[r² f_l(r) j_l(kr)]
 * @output F_l(k).
 * @note   Equal-length data; odd, panelwise-uniform radial grid.
 */
inline double spherical_FT_r2k(const Eigen::Ref<const Eigen::VectorXd>& fl_F1D_r, const Eigen::Ref<const Eigen::VectorXd>& r_F1D_r, double k_F, int l_I = 0) {
    // {r_i,f_i} → Q_S[r²f_l(r)j_l(kr)].
    int Nr_I = r_F1D_r.size();
    int Ninterval_I = Nr_I - 1;
    assert(l_I >= 0);
    assert(fl_F1D_r.size() == Nr_I);
    assert(Nr_I >= 3 && Nr_I % 2 == 1);

    double integrand2_F = fl_F1D_r(0) * r_F1D_r(0) * r_F1D_r(0) * gsl_sf_bessel_jl(l_I, k_F * r_F1D_r(0));
    double integral_F = 0.0;

    // [r_{2p},r_{2p+2}] → (h/3)(g_0+4g_1+g_2).
    for (int panel_I = 0; panel_I <= Ninterval_I / 2 - 1; ++panel_I) {
        int r0_I = 2 * panel_I;
        int r1_I = 2 * panel_I + 1;
        int r2_I = 2 * panel_I + 2;
        double r0_F = r_F1D_r(r0_I);
        double r1_F = r_F1D_r(r1_I);
        double r2_F = r_F1D_r(r2_I);
        assert(std::abs((r2_F - r1_F) / (r1_F - r0_F) - 1.0) < 1.0e-12);
        double integrand0_F = integrand2_F;
        double integrand1_F = fl_F1D_r(r1_I) * r1_F * r1_F * gsl_sf_bessel_jl(l_I, k_F * r1_F);
        integrand2_F = fl_F1D_r(r2_I) * r2_F * r2_F * gsl_sf_bessel_jl(l_I, k_F * r2_F);
        double dr_F = (r2_F - r0_F) / 2.0;
        integral_F += (dr_F / 3.0) * (integrand0_F + 4.0 * integrand1_F + integrand2_F);
    }
    return 4.0 * M_PI * integral_F;
}

/**
 * @brief  Momentum-to-radial transform by GSL QAG (61-point Gauss-Kronrod).
 * @math   f_l(r) = (2π²)⁻¹ ∫_{k_min}^{k_max} k² F_l(k) j_l(kr) dk
 * @output f_l(r).
 * @note   Finite [k_min, k_max] approximates [0, ∞).
 */
inline double spherical_FT_k2r(const RealToRealFunc& Fl_Func, double r_F, int l_I = 0, double kmin_F = 1e-6, double kmax_F = 1e2) {
    // F_l(k) → k²F_l(k)j_l(kr) → f_l(r).
    assert(l_I >= 0);
    assert(kmin_F < kmax_F);
    auto integrand_Func = [&](double k_F) {
        return k_F * k_F * Fl_Func(k_F) * gsl_sf_bessel_jl(l_I, k_F * r_F);
    };
    double integral_F = integrate_qag(integrand_Func, kmin_F, kmax_F);
    return integral_F / (2.0 * M_PI * M_PI);
}

/**
 * @brief  Momentum-to-radial transform by composite Simpson quadrature.
 * @math   f_l(r) ≈ (2π²)⁻¹ Q_S[k² F_l(k) j_l(kr)]
 * @output f_l(r).
 * @note   Equal-length data; odd, panelwise-uniform momentum grid.
 */
inline double spherical_FT_k2r(const Eigen::Ref<const Eigen::VectorXd>& Fl_F1D_k, const Eigen::Ref<const Eigen::VectorXd>& k_F1D_k, double r_F, int l_I = 0) {
    // {k_i,F_i} → Q_S[k²F_l(k)j_l(kr)].
    int Nk_I = k_F1D_k.size();
    int Ninterval_I = Nk_I - 1;
    assert(l_I >= 0);
    assert(Fl_F1D_k.size() == Nk_I);
    assert(Nk_I >= 3 && Nk_I % 2 == 1);

    double integrand2_F = Fl_F1D_k(0) * k_F1D_k(0) * k_F1D_k(0) * gsl_sf_bessel_jl(l_I, k_F1D_k(0) * r_F);
    double integral_F = 0.0;

    // [k_{2p},k_{2p+2}] → (h/3)(g_0+4g_1+g_2).
    for (int panel_I = 0; panel_I <= Ninterval_I / 2 - 1; ++panel_I) {
        int k0_I = 2 * panel_I;
        int k1_I = 2 * panel_I + 1;
        int k2_I = 2 * panel_I + 2;
        double k0_F = k_F1D_k(k0_I);
        double k1_F = k_F1D_k(k1_I);
        double k2_F = k_F1D_k(k2_I);
        assert(std::abs((k2_F - k1_F) / (k1_F - k0_F) - 1.0) < 1.0e-12);
        double integrand0_F = integrand2_F;
        double integrand1_F = Fl_F1D_k(k1_I) * k1_F * k1_F * gsl_sf_bessel_jl(l_I, k1_F * r_F);
        integrand2_F = Fl_F1D_k(k2_I) * k2_F * k2_F * gsl_sf_bessel_jl(l_I, k2_F * r_F);
        double dk_F = (k2_F - k0_F) / 2.0;
        integral_F += (dk_F / 3.0) * (integrand0_F + 4.0 * integrand1_F + integrand2_F);
    }
    return integral_F / (2.0 * M_PI * M_PI);
}
