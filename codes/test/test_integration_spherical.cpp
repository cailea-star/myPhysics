/**
 * @file    test_integration_spherical.cpp
 * @author  cailea
 * @date    2026-08-18
 * @brief   Tests for spherical integration utilities.
 */

#include "integration.hpp"

#include <Eigen/Dense>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <gsl/gsl_sf_legendre.h>

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
 * @brief  Test spherical Gauss-Legendre quadrature.
 * @math   M_ij = ∫Y_i(Ω)Y_j(Ω) dΩ
 * @output Labeled overlap matrices and an acceptance assertion.
 */
int main() {
    // M_ij = ∫Y_i(Ω)Y_j(Ω) dΩ.
    int Nlm_I = 12;
    int i_I = 0;
    int j_I = 0;
    Eigen::MatrixXd M_F2D_i_j = Eigen::MatrixXd::Zero(Nlm_I, Nlm_I);
    Eigen::MatrixXd Mref_F2D_i_j = Eigen::MatrixXd::Zero(Nlm_I, Nlm_I);
    for (int l1_I = 3; l1_I <= 6; ++l1_I) {
        for (int m1_I = 0; m1_I <= l1_I / 2; ++m1_I) {
            j_I = 0;
            for (int l2_I = 3; l2_I <= 6; ++l2_I) {
                for (int m2_I = 0; m2_I <= l2_I / 2; ++m2_I) {
                    RealRealToRealFunc f_Func = std::bind(Ylm_product, l1_I, m1_I, l2_I, m2_I, std::placeholders::_1, std::placeholders::_2);
                    int Nphi_I = std::max(2 * (m1_I + m2_I) + 1, 2);
                    int Ntheta_I = (l1_I + l2_I) / 2 + 1;
                    M_F2D_i_j(i_I, j_I) = integrate_spherical(f_Func, Ntheta_I, Nphi_I);
                    Mref_F2D_i_j(i_I, j_I) = l1_I == l2_I && m1_I == m2_I ? (m1_I == 0 ? 1.0 : 0.5) : 0.0;
                    ++j_I;
                }
            }
            ++i_I;
        }
    }
    double tol_F = 1.0e-4;

    // Input → reference → computed.
    std::cout << std::scientific << std::setprecision(4) << std::left;
    std::cout << "[Input] spherical states = " << Nlm_I << "\n";
    std::cout << "[Reference] spherical overlap M\n";
    for (int i_I = 0; i_I < 5; ++i_I) {std::cout << Mref_F2D_i_j.block(i_I, 0, 1, 5) << " ... " << Mref_F2D_i_j.block(i_I, Nlm_I - 5, 1, 5) << "\n";}
    std::cout << "...\n";
    for (int i_I = Nlm_I - 5; i_I < Nlm_I; ++i_I) {std::cout << Mref_F2D_i_j.block(i_I, 0, 1, 5) << " ... " << Mref_F2D_i_j.block(i_I, Nlm_I - 5, 1, 5) << "\n";}
    std::cout << "[Computed] spherical overlap M\n";
    for (int i_I = 0; i_I < 5; ++i_I) {std::cout << M_F2D_i_j.block(i_I, 0, 1, 5) << " ... " << M_F2D_i_j.block(i_I, Nlm_I - 5, 1, 5) << "\n";}
    std::cout << "...\n";
    for (int i_I = Nlm_I - 5; i_I < Nlm_I; ++i_I) {std::cout << M_F2D_i_j.block(i_I, 0, 1, 5) << " ... " << M_F2D_i_j.block(i_I, Nlm_I - 5, 1, 5) << "\n";}

    // max|M - M_ref| < tol.
    assert((M_F2D_i_j - Mref_F2D_i_j).cwiseAbs().maxCoeff() < tol_F);
    return 0;
}
