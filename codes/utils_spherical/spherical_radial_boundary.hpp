/**
 * @file    spherical_radial_boundary.hpp
 * @author  cailea
 * @date    2026-05-26
 * @brief   Spherical radial boundary-condition utilities.
 */

#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <functional>
#include <stdexcept>
#include <utility>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_coulomb.h>

using doubleC = std::complex<double>;
using SphericalRadialBoundary = std::pair<doubleC, doubleC>;
using SphericalRadialBoundaryFunc = std::function<SphericalRadialBoundary(double r_F, int l_I, double hmass_F, double Ze2_F, double E_F)>;

/**
 * @brief  Generate regular boundary by analytic power law.
 * @math   u_l(r) = A r^{l+1}, u_l'(r) = A(l+1)r^l
 * @output Boundary pair (u_l, u_l').
 * @note   hmass_F, Ze2_F, and E_F are callback-only.
 */
inline SphericalRadialBoundary spherical_radial_boundary_regular(double r_F, int l_I, double hmass_F, double Ze2_F, double E_F) {
    // (r,l) → (u_l,u_l').
    assert(std::isfinite(r_F) && r_F > 0.0);
    assert(l_I >= 0);
    (void) hmass_F;
    (void) Ze2_F;
    (void) E_F;
    double norm_F = 1.0 / std::sqrt(std::tgamma(l_I + 1.0));
    double u_F = norm_F * std::pow(r_F, l_I + 1);
    double du_F = norm_F * (l_I + 1.0) * std::pow(r_F, l_I);
    return {doubleC(u_F, 0.0), doubleC(du_F, 0.0)};
}

/**
 * @brief  Generate a Coulomb-Hankel boundary using GSL Coulomb functions.
 * @math   H_l^(±) = G_l ± iF_l, dH_l^(±)/dr = k(dG_l/dρ ± i dF_l/dρ)
 * @output Commonly scaled boundary pair (H_l^(±), dH_l^(±)/dr).
 * @note   Requires E_F > 0; sign_I = ±1 selects H_l^(±).
 * @note   Removes common GSL scale to prevent overflow.
 */
inline SphericalRadialBoundary spherical_radial_boundary_coulomb_hankel(double r_F, int l_I, double hmass_F, double Ze2_F, double E_F, int sign_I) {
    // (r,l,h_μ,Ze²,E,sign) → (k,η,ρ).
    assert(std::isfinite(r_F) && r_F > 0.0);
    assert(l_I >= 0);
    assert(std::isfinite(hmass_F) && hmass_F > 0.0);
    assert(std::isfinite(E_F) && E_F > 0.0);
    assert(sign_I == 1 || sign_I == -1);
    double k_F = std::sqrt(E_F / hmass_F);
    double eta_F = Ze2_F / (2.0 * hmass_F * k_F);
    double rho_F = k_F * r_F;
    double sign_F = static_cast<double>(sign_I);

    // (η,ρ,l) → (F_l,F_l',G_l,G_l').
    gsl_sf_result F_GSL;
    gsl_sf_result dFdrho_GSL;
    gsl_sf_result G_GSL;
    gsl_sf_result dGdrho_GSL;
    double expF_F = 0.0;
    double expG_F = 0.0;

    int status_I = gsl_sf_coulomb_wave_FG_e(eta_F, rho_F, static_cast<double>(l_I), 0, &F_GSL, &dFdrho_GSL, &G_GSL, &dGdrho_GSL, &expF_F, &expG_F);
    if (status_I != GSL_SUCCESS && status_I != GSL_EOVRFLW) {throw std::runtime_error("GSL Coulomb-Hankel evaluation failed");}

    // (F_l,F_l',G_l,G_l') → e^{-s}(F_l,F_l',G_l,G_l') → (H_l^{±},∂_rH_l^{±}).
    double expScale_F = std::max(expF_F, expG_F);
    double F_F = F_GSL.val * std::exp(expF_F - expScale_F);
    double dFdrho_F = dFdrho_GSL.val * std::exp(expF_F - expScale_F);
    double G_F = G_GSL.val * std::exp(expG_F - expScale_F);
    double dGdrho_F = dGdrho_GSL.val * std::exp(expG_F - expScale_F);
    doubleC i_C(0.0, 1.0);
    doubleC u_C = G_F + sign_F * i_C * F_F;
    doubleC du_C = k_F * (dGdrho_F + sign_F * i_C * dFdrho_F);
    return {u_C, du_C};
}

/**
 * @brief  Generate the outgoing Coulomb-Hankel boundary through the GSL kernel.
 * @math   H_l^(+) = G_l + iF_l
 * @output Commonly scaled outgoing boundary pair (H_l^(+), dH_l^(+)/dr).
 * @note   Outgoing assumes the exp(-iEt/ħ) time convention.
 */
inline SphericalRadialBoundary spherical_radial_boundary_coulomb_hplus(double r_F, int l_I, double hmass_F, double Ze2_F, double E_F) {
    return spherical_radial_boundary_coulomb_hankel(r_F, l_I, hmass_F, Ze2_F, E_F, 1);
}

/**
 * @brief  Generate the incoming Coulomb-Hankel boundary through the GSL kernel.
 * @math   H_l^(-) = G_l - iF_l
 * @output Commonly scaled incoming boundary pair (H_l^(-), dH_l^(-)/dr).
 * @note   Incoming assumes the exp(-iEt/ħ) time convention.
 */
inline SphericalRadialBoundary spherical_radial_boundary_coulomb_hminus(double r_F, int l_I, double hmass_F, double Ze2_F, double E_F) {
    return spherical_radial_boundary_coulomb_hankel(r_F, l_I, hmass_F, Ze2_F, E_F, -1);
}
