/**
 * @file    root.hpp
 * @author  cailea
 * @date    2026-01-29
 * @brief   Scalar root-finding algorithms.
 */

#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>

using Real2RealFunc = std::function<double(double)>;

/**
 * @brief  Find a scalar root by bisection.
 * @math   f(x^*) = 0, f(xlo)f(xup) ≤ 0
 * @output Root approximation in [xlo, xup].
 */
inline double root_bisection(const Real2RealFunc& f_Func, double xlo_F, double xup_F, double tol_F = 1e-6) {
    // [x_lo,x_up] → (f_lo,f_up).
    double flo_F = f_Func(xlo_F);
    double fup_F = f_Func(xup_F);
    if (std::abs(flo_F) <= tol_F) {return xlo_F;}
    if (std::abs(fup_F) <= tol_F) {return xup_F;}

    assert(flo_F * fup_F <= 0.0);

    // [xlo_h, xup_h] → [xlo_{h+1}, xup_{h+1}].
    for (int h_I = 0; h_I < 100 && std::abs(xup_F - xlo_F) > tol_F; ++h_I) {
        double xmid_F = (xlo_F + xup_F) / 2.0;
        double fmid_F = f_Func(xmid_F);
        if (std::abs(fmid_F) < tol_F) {return xmid_F;}
        if (flo_F * fmid_F < 0.0) {
            xup_F = xmid_F;
            fup_F = fmid_F;
        } else {
            xlo_F = xmid_F;
            flo_F = fmid_F;
        }
    }
    return (xlo_F + xup_F) / 2.0;
}

/**
 * @brief  Find root by safeguarded finite-difference Newton iteration.
 * @math   x_{h+1}=x_h-f_h/f′_h, f′_h≈[f(x_h+ε_x)-f_h]/ε_x
 * @output Root or NaN for non-finite f.
 * @note   Invalid Newton steps use the bracket midpoint.
 */
inline double root_newton(const Real2RealFunc& f_Func, double xlo_F, double xup_F, double ftol_F = 1e-6, double dxmax_F = 5.0, int Nitermax_I = 50) {
    // [x_lo,x_up] → (ε_x,f_lo,f_up,x_0).
    assert(xlo_F < xup_F);

    assert(std::isfinite(dxmax_F) && dxmax_F > 0.0);
    double epsilonx_F = 0.01 * dxmax_F;

    double flo_F = f_Func(xlo_F);
    double fup_F = f_Func(xup_F);
    assert(std::isfinite(flo_F) && std::isfinite(fup_F));
    if (std::abs(flo_F) <= ftol_F) {return xlo_F;}
    if (std::abs(fup_F) <= ftol_F) {return xup_F;}
    assert(flo_F * fup_F <= 0.0);

    double x_F = 0.5 * (xlo_F + xup_F);

    // x_h → f_h → bracket_h → x_{h+1}.
    for (int h_I = 0; h_I <= Nitermax_I; ++h_I) {
        double fx_F = f_Func(x_F);

        if (!std::isfinite(fx_F)) {
            std::cerr << "[ERROR]: In [root_newton] Function value is not finite." << std::endl;
            return NAN;
        }
        if (std::abs(fx_F) <= ftol_F) {return x_F;}
        if (h_I == Nitermax_I) {break;}

        // f(x_lo)f(x_up) ≤ 0.
        if (flo_F * fx_F <= 0.0) {
            xup_F = x_F;
            fup_F = fx_F;
        } else {
            xlo_F = x_F;
            flo_F = fx_F;
        }

        // f'(x_h) ≈ [f(x_h + ε_x) - f(x_h)]/ε_x, x_{h+1} = x_h - f(x_h)/f'(x_h).
        double dfdx_F = (f_Func(x_F + epsilonx_F) - fx_F) / epsilonx_F;
        double xnext_F = 0.5 * (xlo_F + xup_F);
        if (std::isfinite(dfdx_F) && std::abs(dfdx_F) >= 1.0e-20) {
            double dx_F = - fx_F / dfdx_F;
            if (std::isfinite(dx_F)) {xnext_F = x_F + std::clamp(dx_F, -dxmax_F, dxmax_F);}
        }

        if (xnext_F < xlo_F || xnext_F > xup_F) {xnext_F = 0.5 * (xlo_F + xup_F);}

        // |x_{h+1}-x_h| ≤ 10⁻¹² max(1,|x_h|) → stop.
        if (std::abs(xnext_F - x_F) <= 1.0e-12 * std::max(1.0, std::abs(x_F))) {
            std::cerr << "[WARNING]: In [root_newton] Convergence not reached." << std::endl;
            break;
        }

        x_F = xnext_F; // x_h ← x_{h+1}.
    }
    return x_F;
}

/**
 * @brief  Find a scalar root with the GSL Brent solver.
 * @math   f(x^*) = 0, f(xlo)f(xup) ≤ 0
 * @output Root approximation in [xlo, xup].
 */
inline double root_brent(const Real2RealFunc& f_Func, double xlo_F, double xup_F, double tol_F = 1e-6) {
    // f(x_lo)f(x_up) ≤ 0.
    if (f_Func(xlo_F) * f_Func(xup_F) > 0.0) {
        std::cerr << "[ERROR]: In [root_brent] Function values at the endpoints must have opposite signs." << std::endl;
        return NAN;
    }

    // [x_lo,x_up] → GSL Brent state.
    gsl_function gsl_Func;
    gsl_Func.function = [](double x_F, void* f_Ptr) {
        return (*static_cast<Real2RealFunc*>(f_Ptr))(x_F);
    };
    gsl_Func.params = const_cast<Real2RealFunc*>(&f_Func);
    gsl_root_fsolver* gsl_Solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
    gsl_root_fsolver_set(gsl_Solver, &gsl_Func, xlo_F, xup_F);

    // [xlo_h, xup_h] → [xlo_{h+1}, xup_{h+1}].
    int status_I = GSL_CONTINUE;
    for (int h_I = 0; status_I == GSL_CONTINUE && h_I < 100; ++h_I) {
        gsl_root_fsolver_iterate(gsl_Solver);
        xlo_F = gsl_root_fsolver_x_lower(gsl_Solver);
        xup_F = gsl_root_fsolver_x_upper(gsl_Solver);
        status_I = gsl_root_test_interval(xlo_F, xup_F, tol_F, 0.0);
    }

    double xroot_F = gsl_root_fsolver_root(gsl_Solver);
    gsl_root_fsolver_free(gsl_Solver);
    return xroot_F;
}
