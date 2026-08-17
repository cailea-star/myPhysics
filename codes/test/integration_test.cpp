/**
 * @file    integration_test.cpp
 * @author  cailea
 * @date    2026-08-17
 * @brief   Manual checks for numerical integration utilities.
 */

#include "integration.hpp"

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <string>
#include <gsl/gsl_sf_legendre.h>

double V0_F = 1.0;
double a0_F = 2.0;
double R0_F = 6.0;

/**
 * @brief  Evaluate the Woods-Saxon potential.
 * @math   V(x) = V_0 / [1 + exp((x - R_0) / a_0)]
 * @output Potential value at x.
 */
double WS_func(double x_F) {
    double exp_F = std::exp((x_F - R0_F) / a0_F);
    return V0_F / (1.0 + exp_F);
}

/**
 * @brief  Evaluate the derivative of the Woods-Saxon potential.
 * @math   V'(x) = -V_0 exp((x - R_0) / a_0) / {a_0 [1 + exp((x - R_0) / a_0)]^2}
 * @output Potential derivative at x.
 */
double dWS_func(double x_F) {
    double exp_F = std::exp((x_F - R0_F) / a0_F);
    return -V0_F * std::pow(1.0 + exp_F, -2) * exp_F / a0_F;
}

/**
 * @brief  Evaluate a real spherical-harmonic factor.
 * @math   Y_lm(θ, φ) = P_l^{|m|}(cos θ) cos(mφ)
 * @output Real spherical-harmonic factor at (θ, φ).
 */
double Ylm(int l_I, int m_I, double theta_F, double phi_F) {
    double Y_F = gsl_sf_legendre_sphPlm(l_I, std::abs(m_I), std::cos(theta_F));
    return Y_F * std::cos(m_I * phi_F);
}

/**
 * @brief  Evaluate a product of two real spherical-harmonic factors.
 * @math   F(θ, φ) = Y_lm(θ, φ) Y_l'm'(θ, φ)
 * @output Spherical-harmonic product at (θ, φ).
 */
double Ylm_product(int l1_I, int m1_I, int l2_I, int m2_I, double theta_F, double phi_F) {
    return Ylm(l1_I, m1_I, theta_F, phi_F) * Ylm(l2_I, m2_I, theta_F, phi_F);
}

/**
 * @brief  Print representative integration results.
 * @math   ∫V'(x) dx = V(x_max) - V(x_min), ∫Y_lm Y_l'm' dΩ
 * @output Printed one-dimensional errors and spherical overlap matrix.
 */
int main() {
    std::cout << std::scientific << std::setprecision(4) << std::left;
    double xmin_F = 0.0;
    double xmax_F = 100.0; // x_max ≈ ∞.
    int Nx_I = 1001;
    Eigen::VectorXd x_F1D_x = Eigen::VectorXd::LinSpaced(Nx_I, xmin_F, xmax_F);
    Eigen::VectorXd f_F1D_x(Nx_I);
    for (int x_I = 0; x_I < x_F1D_x.size(); ++x_I) {
        f_F1D_x(x_I) = dWS_func(x_F1D_x(x_I));
    }

    std::cout << "\n Testing integrate_qag with f(x) = dWS_func(x):\n";
    double integralExact_F = WS_func(xmax_F) - WS_func(xmin_F); // ∫V'(x) dx = V(x_max) - V(x_min).
    double integralQag_F = integrate_qag(dWS_func, xmin_F, xmax_F);
    double integralSimpson_F = integrate_simpson(f_F1D_x, x_F1D_x);
    double integralTrapezoidal_F = integrate_trapezoidal(f_F1D_x, x_F1D_x);
    std::cout << std::setw(30) << "integral_exact" << std::setw(15) << integralExact_F << "\n";
    std::cout << std::setw(30) << "integrate_qag" << std::setw(15) << integralQag_F;
    std::cout << std::setw(30) << "error_qag" << std::setw(15) << integralQag_F - integralExact_F << "\n";
    std::cout << std::setw(30) << "integrate_simpson" << std::setw(15) << integralSimpson_F;
    std::cout << std::setw(30) << "error_simpson" << std::setw(15) << integralSimpson_F - integralExact_F << "\n";
    std::cout << std::setw(30) << "integrate_trapezoidal" << std::setw(15) << integralTrapezoidal_F;
    std::cout << std::setw(30) << "error_trapezoidal" << std::setw(15) << integralTrapezoidal_F - integralExact_F << "\n";

    std::cout << std::right << "\n integrand of [Ylm * Yl'm']:\n";
    for (int l1_I = 3; l1_I <= 6; ++l1_I) {
        for (int m1_I = 0; m1_I <= l1_I / 2; ++m1_I) {
            for (int l2_I = 3; l2_I <= 6; ++l2_I) {
                for (int m2_I = 0; m2_I <= l2_I / 2; ++m2_I) {
                    RealRealToRealFunc f_Func = std::bind(Ylm_product, l1_I, m1_I, l2_I, m2_I, std::placeholders::_1, std::placeholders::_2);
                    int Nphi_I = std::max(2 * (m1_I + m2_I) + 1, 2);
                    int Ntheta_I = (l1_I + l2_I) / 2 + 1;
                    double integral_F = integrate_spherical(f_Func, Ntheta_I, Nphi_I);
                    std::cout << std::setw(10) << std::to_string(integral_F);
                }
            }
            std::cout << "\n";
        }
    }
}
