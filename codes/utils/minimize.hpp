/**
 * @file    minimize.hpp
 * @author  cailea
 * @date    2026-01-07
 * @brief   Scalar minimization algorithms.
 */

#pragma once

#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>

using Real2RealFunc = std::function<double(double)>;

// ==================== Golden-section search ====================

/**
 * @brief  Find a scalar minimizer using golden-section search.
 * @math   x^* = arg min_{x in [x_lo,x_up]} f(x)
 * @output Returns the estimated minimizer x^*.
 * @note   Assumes that f(x) is unimodal on [x_lo,x_up].
 */
inline double minimize_golden(const Real2RealFunc& f_Func, double xlo_F, double xup_F, double tol_F = 1e-6) {
    const double tau_F = (std::sqrt(5.0) - 1.0) / 2.0;
    double xlo_tmp_F = xup_F - tau_F * (xup_F - xlo_F);
    double xup_tmp_F = xlo_F + tau_F * (xup_F - xlo_F);
    double flo_tmp_F = f_Func(xlo_tmp_F);
    double fup_tmp_F = f_Func(xup_tmp_F);

    for (int h_I = 0; h_I < 100; ++h_I) {
        if (flo_tmp_F < fup_tmp_F) {
            xup_F = xup_tmp_F;
            xup_tmp_F = xlo_tmp_F;
            fup_tmp_F = flo_tmp_F;
            xlo_tmp_F = xup_F - tau_F * (xup_F - xlo_F);
            flo_tmp_F = f_Func(xlo_tmp_F);
        } else {
            xlo_F = xlo_tmp_F;
            xlo_tmp_F = xup_tmp_F;
            flo_tmp_F = fup_tmp_F;
            xup_tmp_F = xlo_F + tau_F * (xup_F - xlo_F);
            fup_tmp_F = f_Func(xup_tmp_F);
        }
        if (std::abs(xup_F - xlo_F) < tol_F) {break;}
    }
    return (xlo_F + xup_F) / 2.0;
}

// ==================== Brent minimization ====================

/**
 * @brief  Find a bracketed scalar minimizer using the GSL Brent method.
 * @math   x_lo < x_mid < x_up, f(x_mid) < min{f(x_lo),f(x_up)}
 * @output Returns the estimated minimizer x^*.
 * @note   The initial three-point bracket is supplied explicitly.
 */
inline double minimize_brent(const Real2RealFunc& f_Func, double xlo_F, double xmid_F, double xup_F, double tol_F) {
    double flo_F = f_Func(xlo_F);
    double fmid_F = f_Func(xmid_F);
    double fup_F = f_Func(xup_F);
    assert(xlo_F < xmid_F && xmid_F < xup_F);
    assert(std::isfinite(flo_F) && std::isfinite(fmid_F) && std::isfinite(fup_F) && flo_F > fmid_F && fup_F > fmid_F);

    gsl_function f_GSL;
    f_GSL.function = [](double x_F, void* params_V) { return (*static_cast<Real2RealFunc*>(params_V))(x_F); };
    f_GSL.params = const_cast<Real2RealFunc*>(&f_Func);
    gsl_min_fminimizer* minimizer_GSL = gsl_min_fminimizer_alloc(gsl_min_fminimizer_brent);
    gsl_min_fminimizer_set(minimizer_GSL, &f_GSL, xmid_F, xlo_F, xup_F);

    int status_I = GSL_CONTINUE;
    for (int h_I = 0; status_I == GSL_CONTINUE && h_I < 100; ++h_I) {
        gsl_min_fminimizer_iterate(minimizer_GSL);
        xlo_F = gsl_min_fminimizer_x_lower(minimizer_GSL);
        xup_F = gsl_min_fminimizer_x_upper(minimizer_GSL);
        status_I = gsl_min_test_interval(xlo_F, xup_F, tol_F, 0.0);
    }

    double xmin_F = gsl_min_fminimizer_x_minimum(minimizer_GSL);
    gsl_min_fminimizer_free(minimizer_GSL);
    return xmin_F;
}

/**
 * @brief  Find a scalar minimizer using an automatically initialized GSL Brent method.
 * @math   x_mid = golden(f,[x_lo,x_up]), f(x_mid) <= min{f(x_lo),f(x_up)}
 * @output Returns the estimated minimizer x^* or NaN if initialization fails.
 * @note   Golden-section search supplies the initial interior point.
 */
inline double minimize_brent(const Real2RealFunc& f_Func, double xlo_F, double xup_F, double tol_F = 1e-6) {
    double xmid_F = minimize_golden(f_Func, xlo_F, xup_F, tol_F);
    double flo_F = f_Func(xlo_F);
    double fmid_F = f_Func(xmid_F);
    double fup_F = f_Func(xup_F);

    if (!std::isfinite(flo_F) || !std::isfinite(fmid_F) || !std::isfinite(fup_F) || flo_F <= fmid_F || fup_F <= fmid_F) {
        std::cerr << "[ERROR]: In [minimize_brent] failed to find x let f(xlo_F) > f(x) < f(xup_F)." << std::endl;
        return NAN;
    }

    return minimize_brent(f_Func, xlo_F, xmid_F, xup_F, tol_F);
}
