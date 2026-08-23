/**
 * @file    integration.hpp
 * @author  cailea
 * @date    2026-04-09
 * @brief   Numerical integration utilities.
 */

#pragma once

#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <functional>
#include <gsl/gsl_integration.h>

using RealToRealFunc = std::function<double(double)>;
using RealRealToRealFunc = std::function<double(double, double)>;

/**
 * @brief  Integrate a scalar function with adaptive Gauss-Kronrod quadrature.
 * @math   I = ∫_{x_min}^{x_max} f(x) dx
 * @output Numerical integral over [x_min, x_max].
 */
inline double integrate_qag(const RealToRealFunc& f_Func, double xmin_F, double xmax_F) {
    // QAG₆₁ → (N_sub,ε_abs,ε_rel).
    const std::size_t Nsubmax_I = 4096;
    const double epsabs_F = 1.0e-8;
    const double epsrel_F = 1.0e-8;
    const int rule_I = GSL_INTEG_GAUSS61;

    // f → gsl_function.
    gsl_function gsl_Func;
    gsl_Func.function = [](double x_F, void* f_Func_ptr) {
        return (*static_cast<const RealToRealFunc*>(f_Func_ptr))(x_F);
    };
    gsl_Func.params = const_cast<void*>(static_cast<const void*>(&f_Func));

    // (f,[x_min,x_max]) → (I,ε_I).
    gsl_integration_workspace* workspace = gsl_integration_workspace_alloc(Nsubmax_I);
    double integral_F;
    double error_F;
    gsl_integration_qag(&gsl_Func, xmin_F, xmax_F, epsabs_F, epsrel_F, Nsubmax_I, rule_I, workspace, &integral_F, &error_F);
    gsl_integration_workspace_free(workspace);
    return integral_F;
}

/**
 * @brief  Integrate discrete samples with the composite Simpson rule.
 * @math   Q_S = Σ_k h_k [f_{2k} + 4f_{2k+1} + f_{2k+2}] / 3
 * @output Composite Simpson integral.
 * @note   |f|=|x|=2p+1≥3; Δx_{2p}=Δx_{2p+1}.
 */
inline double integrate_simpson(const Eigen::Ref<const Eigen::VectorXd>& f_F1D_x, const Eigen::Ref<const Eigen::VectorXd>& x_F1D_x) {
    // ({f_i},{x_i}) → (N_x,N_interval).
    int Nx_I = f_F1D_x.size();
    int Ninterval_I = Nx_I - 1;
    assert(f_F1D_x.size() == x_F1D_x.size());
    assert((Nx_I >= 3 && Nx_I % 2 == 1) && "Simpson's rule requires an odd number of points (>= 3)");

    // [x_{2p},x_{2p+2}] → h_p(f_{2p}+4f_{2p+1}+f_{2p+2})/3.
    double integral_F = 0.0;
    for (int panel_I = 0; panel_I <= Ninterval_I / 2 - 1; ++panel_I) {
        int x0_I = 2 * panel_I;
        int x1_I = 2 * panel_I + 1;
        int x2_I = 2 * panel_I + 2;
        assert(std::abs((x_F1D_x(x2_I) - x_F1D_x(x1_I)) / (x_F1D_x(x1_I) - x_F1D_x(x0_I)) - 1.0) < 1.0e-12);
        double h_F = (x_F1D_x(x2_I) - x_F1D_x(x0_I)) / 2.0;
        integral_F += (h_F / 3.0) * (f_F1D_x(x0_I) + 4.0 * f_F1D_x(x1_I) + f_F1D_x(x2_I));
    }
    return integral_F;
}

/**
 * @brief  Integrate discrete samples with the composite trapezoidal rule.
 * @math   Q_T = Σ_i (x_{i+1} - x_i)(f_i + f_{i+1}) / 2
 * @output Composite trapezoidal integral.
 * @note   |f|=|x|≥2.
 */
inline double integrate_trapezoidal(const Eigen::Ref<const Eigen::VectorXd>& f_F1D_x, const Eigen::Ref<const Eigen::VectorXd>& x_F1D_x) {
    // ({f_i},{x_i}) → (N_x,N_interval).
    int Nx_I = f_F1D_x.size();
    int Ninterval_I = Nx_I - 1;
    assert(f_F1D_x.size() == x_F1D_x.size());
    assert(Nx_I >= 2);

    // [x_i,x_{i+1}] → Q_{T,i}.
    double integral_F = 0.0;
    for (int x_I = 0; x_I <= Ninterval_I - 1; ++x_I) {
        double h_F = x_F1D_x(x_I + 1) - x_F1D_x(x_I);
        integral_F += (h_F / 2.0) * (f_F1D_x(x_I) + f_F1D_x(x_I + 1));
    }
    return integral_F;
}

/**
 * @brief  Integrate sphere by Gauss-Legendre and trapezoidal quadrature.
 * @math   I_Ω = ∫_0^{2π} ∫_0^π f(θ, φ) sin(θ) dθ dφ
 * @output Numerical integral over the unit sphere.
 * @note   N_θ,N_φ>0; N_φ=1 requires ∂_φf=0.
 */
inline double integrate_spherical(const RealRealToRealFunc& f_Func, int Ntheta_I = 7, int Nphi_I = 25) {
    // (N_θ,N_φ) → ({φ_j},{f_j}).
    assert(Ntheta_I > 0 && Nphi_I > 0);
    Eigen::VectorXd phi_F1D_phi = Eigen::VectorXd::LinSpaced(Nphi_I, 0.0, 2.0 * M_PI);
    Eigen::VectorXd fphi_F1D_phi = Eigen::VectorXd::Zero(Nphi_I);

    // [-1,1] → {μ_i,w_i}.
    gsl_integration_fixed_workspace* workspace = gsl_integration_fixed_alloc(gsl_integration_fixed_legendre, Ntheta_I, -1.0, 1.0, 0.0, 0.0);
    const double* mu_F1D_mu = gsl_integration_fixed_nodes(workspace);
    const double* w_F1D_mu = gsl_integration_fixed_weights(workspace);

    // θ_i=arccos μ_i; I_Ω=Σ_iw_iQ_φ[f(θ_i,φ)].
    double integral_F = 0.0;
    for (int theta_I = 0; theta_I <= Ntheta_I - 1; ++theta_I) {
        double theta_F = std::acos(mu_F1D_mu[theta_I]);
        double integralPhi_F;
        if (Nphi_I == 1) {
            integralPhi_F = f_Func(theta_F, 0.0) * 2.0 * M_PI;
        } else {
            for (int phi_I = 0; phi_I <= Nphi_I - 1; ++phi_I) {fphi_F1D_phi(phi_I) = f_Func(theta_F, phi_F1D_phi(phi_I));}
            integralPhi_F = integrate_trapezoidal(fphi_F1D_phi, phi_F1D_phi);
        }
        integral_F += w_F1D_mu[theta_I] * integralPhi_F;
    }

    gsl_integration_fixed_free(workspace);
    return integral_F;
}
