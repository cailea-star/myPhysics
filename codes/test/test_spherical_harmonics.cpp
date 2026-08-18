/**
 * @file    test_spherical_harmonics.cpp
 * @author  cailea
 * @date    2026-08-17
 * @brief   Tests for real spherical-harmonic expansion utilities.
 */

#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <Eigen/Core>
#include <gsl/gsl_sf_legendre.h>
#include "spherical_harmonics.hpp"

/**
 * @brief  Evaluate Re Y_lm using the GSL normalized associated Legendre function.
 * @math   Re Y_lm(θ, φ) = P_l^|m|(cos θ) cos(mφ)
 * @output Re Y_lm(θ, φ).
 * @note   Uses the real cosine sector.
 */
double ReYlm(int l_I, int m_I, double theta_F, double phi_F) {
    return gsl_sf_legendre_sphPlm(l_I, std::abs(m_I), std::cos(theta_F)) * std::cos(m_I * phi_F);
}

/**
 * @brief  Test Gauss-Legendre/trapezoidal projection and direct reconstruction.
 * @math   P_ij = ⟨Re Y_i, Re Y_j⟩, R_ij = ⟨Re Y_i, reconstruct(P_ij Re Y_j)⟩
 * @output Labeled matrix corners with acceptance assertions.
 * @note   Uses 4 ≤ l ≤ 6 and 0 ≤ m ≤ l/2 in the real cosine sector.
 */
int main() {
    const int Nlm_I = 10;
    const double tol_F = 1.0e-12;
    Eigen::MatrixXd P_F2D_i_j(Nlm_I, Nlm_I);
    Eigen::MatrixXd R_F2D_i_j(Nlm_I, Nlm_I);
    Eigen::MatrixXd Ref_F2D_i_j = Eigen::MatrixXd::Zero(Nlm_I, Nlm_I);

    int i_I = 0;
    for (int l1_I = 4; l1_I <= 6; ++l1_I) {
        for (int m1_I = 0; m1_I <= l1_I / 2; ++m1_I) {
            auto Y1_Func = [&](double theta_F, double phi_F) {
                return ReYlm(l1_I, m1_I, theta_F, phi_F);
            };
            int j_I = 0;
            for (int l2_I = 4; l2_I <= 6; ++l2_I) {
                for (int m2_I = 0; m2_I <= l2_I / 2; ++m2_I) {
                    int Ntheta_I = (l1_I + l2_I) / 2 + 1;
                    int Nphi_I = 2 * (m1_I + m2_I) + 1;
                    P_F2D_i_j(i_I, j_I) = spherical_Omega2lm(Y1_Func, l2_I, m2_I, Ntheta_I, Nphi_I);
                    auto Yrec_Func = [&](double theta_F, double phi_F) {
                        return spherical_lm2Omega(P_F2D_i_j(i_I, j_I), l2_I, m2_I, theta_F, phi_F);
                    };
                    R_F2D_i_j(i_I, j_I) = spherical_Omega2lm(Yrec_Func, l1_I, m1_I, Ntheta_I, Nphi_I);
                    if (l1_I == l2_I && m1_I == m2_I) {Ref_F2D_i_j(i_I, j_I) = (m1_I == 0) ? 1.0 : 0.5;}
                    ++j_I;
                }
            }
            ++i_I;
        }
    }

    std::cout << std::scientific << std::setprecision(3) << std::right;
    std::cout << "\nInput: Re Y_lm, 4 <= l <= 6, 0 <= m <= l/2\n";
    std::cout << "Reference:\n" << Ref_F2D_i_j << "\n";
    std::cout << "Computed P:\n" << P_F2D_i_j << "\n";
    std::cout << "Computed R:\n" << R_F2D_i_j << "\n";

    assert((P_F2D_i_j - Ref_F2D_i_j).cwiseAbs().maxCoeff() < tol_F);
    assert((R_F2D_i_j - Ref_F2D_i_j).cwiseAbs().maxCoeff() < tol_F);
    return 0;
}
