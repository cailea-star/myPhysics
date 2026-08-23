/**
 * @file    spherical_harmonics.hpp
 * @author  cailea
 * @date    2026-04-09
 * @brief   Real spherical-harmonic expansion utilities.
 */

#pragma once

#include <cassert>
#include <cmath>
#include <functional>
#include <gsl/gsl_sf_legendre.h>
#include "integration.hpp"

using RealRealToRealFunc = std::function<double(double, double)>;

/**
 * @brief  Project onto Re Y_lm by Gauss-Legendre/trapezoidal quadrature.
 * @math   Re f_lm = ∫ f(Ω) Re Y_lm(Ω) dΩ
 * @output Re f_lm.
 * @note   Cosine sector; Nφ = 1 assumes axial symmetry.
 */
inline double spherical_Omega2lm(const RealRealToRealFunc& f_Func, int l_I, int m_I = 0, int Ntheta_I = 7, int Nphi_I = 1) {
    // f(Ω) → f(Ω) ReY_lm → Re f_lm.
    int mabs_I = std::abs(m_I);
    assert(l_I >= 0 && mabs_I <= l_I);
    assert(Nphi_I > 1 || mabs_I == 0);
    auto fReYlm_Func = [&](double theta_F, double phi_F) {
        double ReYlm_F = gsl_sf_legendre_sphPlm(l_I, mabs_I, std::cos(theta_F)) * std::cos(mabs_I * phi_F);
        return f_Func(theta_F, phi_F) * ReYlm_F;
    };
    double Reflm_F = integrate_spherical(fReYlm_Func, Ntheta_I, Nphi_I);
    return Reflm_F;
}

/**
 * @brief  Reconstruct one real spherical-harmonic mode directly.
 * @math   f_lm(Ω) = (2 - δ_m0) Re f_lm Re Y_lm(Ω)
 * @output Mode contribution f_lm(Ω).
 * @note   Im f_lm = 0.
 */
inline double spherical_lm2Omega(double Reflm_F, int l_I, int m_I, double theta_F, double phi_F) {
    // Re f_lm → (2-δ_m0) Re f_lm ReY_lm.
    int mabs_I = std::abs(m_I);
    assert(l_I >= 0 && mabs_I <= l_I);
    double ReYlm_F = gsl_sf_legendre_sphPlm(l_I, mabs_I, std::cos(theta_F)) * std::cos(mabs_I * phi_F);
    double flmOmega_F = Reflm_F * ReYlm_F;
    double mFactor_F = 2.0 - static_cast<double>(mabs_I == 0);
    flmOmega_F *= mFactor_F;
    return flmOmega_F;
}
