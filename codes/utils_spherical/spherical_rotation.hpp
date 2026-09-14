/**
 * @file    spherical_rotation.hpp
 * @author  cailea
 * @date    2026-09-14
 * @brief   Spherical harmonic-oscillator rotation matrices.
 */

#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <vector>
#include <Eigen/Core>
#include <unsupported/Eigen/CXX11/Tensor>

#include "spherical_setting.hpp"
#include "group_so3_su2.hpp"

using doubleC = std::complex<double>;

class SphericalRotation {
public:
    bool isBuilt_B = false;
    int Nalpha_I = 1;
    int Nbeta_I = 1;
    int Ngamma_I = 1;

    Eigen::VectorXd alpha_F1D_alpha{};  // α_i ∈ [0,2π).
    Eigen::VectorXd weight_F1D_alpha{}; // w_α,i = 2π/N_α.
    Eigen::VectorXd beta_F1D_beta{};    // β_i = acos x_i ∈ (0,π).
    Eigen::VectorXd weight_F1D_beta{};  // ∫₀^π sinβ f(β)dβ.
    Eigen::VectorXd gamma_F1D_gamma{};  // γ_i ∈ [0,4π).
    Eigen::VectorXd weight_F1D_gamma{}; // w_γ,i = 4π/N_γ.

    std::vector<SphericalSPLabel> labels_S1D_sp{}; // m > 0.

    // Rows and columns: (+m,-m); partners are standard Φ_nlj,-m states.
    Eigen::Tensor<doubleC, 3, Eigen::ColMajor> Ry_C3D_2sp_2sp_beta{};
    Eigen::Tensor<doubleC, 3, Eigen::ColMajor> Rz_C3D_2sp_2sp_alpha{};
    Eigen::Tensor<doubleC, 3, Eigen::ColMajor> Rz_C3D_2sp_2sp_gamma{};

    /**
     * @brief  Initialize uniform azimuthal and Gauss-Legendre polar quadrature.
     * @math   x_i = cosβ_i; dΩ = dα sinβ dβ dγ.
     * @output Labels and angle grids; rotation tensors remain empty.
     * @note   N_α,N_β,N_γ > 0; complete m > 0 multiplet representatives.
     * @note   Weights exclude projection normalization factors.
     */
    SphericalRotation(const std::vector<SphericalSPLabel>& labels_S1D_sp_, int Nalpha_I_, int Nbeta_I_, int Ngamma_I_);

    /**
     * @brief  Build y rotations using Jy eigendecomposition.
     * @math   R_y,ab(β) = δ_Na,Nb δ_ja,jb d^ja_ma,mb(β).
     * @output Complex rotation tensor, ordered (+m,-m).
     */
    void build_Ry();

    /**
     * @brief  Build diagonal z rotations using trigonometric phases.
     * @math   R_z,ab(θ) = δ_ab exp(-im_a θ).
     * @output Complex rotation tensors for α and γ.
     */
    void build_Rz();

    /**
     * @brief  Build y and z rotation caches.
     * @math   {α,β,γ} → {R_y,R_z}.
     * @output Complete rotation caches and readiness.
     */
    void build();
};

inline void SphericalRotation::build() {
    // {R_y,R_z} → ready.
    isBuilt_B = false;
    build_Ry();
    build_Rz();
    isBuilt_B = true;
}

inline void SphericalRotation::build_Rz() {
    // α ∈ [0,2π); γ ∈ [0,4π).
    const Eigen::Index Nsp_I = static_cast<Eigen::Index>(labels_S1D_sp.size());
    assert(Nalpha_I > 0 && Ngamma_I > 0);
    assert(alpha_F1D_alpha.size() == Nalpha_I && gamma_F1D_gamma.size() == Ngamma_I);
    assert(alpha_F1D_alpha.allFinite() && gamma_F1D_gamma.allFinite());

    // {R_z(α),R_z(γ)} → allocated storage → 0.
    Rz_C3D_2sp_2sp_alpha.resize(2 * Nsp_I, 2 * Nsp_I, Nalpha_I);
    Rz_C3D_2sp_2sp_gamma.resize(2 * Nsp_I, 2 * Nsp_I, Ngamma_I);
    Rz_C3D_2sp_2sp_alpha.setZero();
    Rz_C3D_2sp_2sp_gamma.setZero();

    // R_z(α) = diag(exp(-imα),exp(+imα)).
    for (int alpha_I = 0; alpha_I < Nalpha_I; ++alpha_I) {
        const double alpha_F = alpha_F1D_alpha(alpha_I);
        for (Eigen::Index sp_I = 0; sp_I < Nsp_I; ++sp_I) {
            const double m_F = 0.5 * labels_S1D_sp[sp_I].twom_I;
            const double cosMAngle_F = std::cos(m_F * alpha_F);
            const double sinMAngle_F = std::sin(m_F * alpha_F);

            // (+m,-m) → (exp(-imα),exp(+imα)).
            Rz_C3D_2sp_2sp_alpha(sp_I, sp_I, alpha_I) = doubleC(cosMAngle_F, -sinMAngle_F);
            Rz_C3D_2sp_2sp_alpha(sp_I + Nsp_I, sp_I + Nsp_I, alpha_I) = doubleC(cosMAngle_F, sinMAngle_F);
        }
    }

    // R_z(γ) = diag(exp(-imγ),exp(+imγ)).
    for (int gamma_I = 0; gamma_I < Ngamma_I; ++gamma_I) {
        const double gamma_F = gamma_F1D_gamma(gamma_I);
        for (Eigen::Index sp_I = 0; sp_I < Nsp_I; ++sp_I) {
            const double m_F = 0.5 * labels_S1D_sp[sp_I].twom_I;
            const double cosMAngle_F = std::cos(m_F * gamma_F);
            const double sinMAngle_F = std::sin(m_F * gamma_F);

            // (+m,-m) → (exp(-imγ),exp(+imγ)).
            Rz_C3D_2sp_2sp_gamma(sp_I, sp_I, gamma_I) = doubleC(cosMAngle_F, -sinMAngle_F);
            Rz_C3D_2sp_2sp_gamma(sp_I + Nsp_I, sp_I + Nsp_I, gamma_I) = doubleC(cosMAngle_F, sinMAngle_F);
        }
    }
}

inline void SphericalRotation::build_Ry() {
    // {m > 0,β} → (N_sp,2j_max).
    const Eigen::Index Nsp_I = static_cast<Eigen::Index>(labels_S1D_sp.size());
    assert(Nbeta_I > 0 && beta_F1D_beta.size() == Nbeta_I);
    assert(beta_F1D_beta.allFinite());
    int twojMax_I = 0;
    for (const SphericalSPLabel& label : labels_S1D_sp) {
        assert(label.twom_I > 0 && label.twom_I <= label.twoj_I);
        twojMax_I = std::max(twojMax_I, label.twoj_I);
    }

    // j → {sp}; 2j = 1,3,...,2j_max.
    std::vector<std::vector<int>> indices_I2D_j_jsp((twojMax_I + 1) / 2);
    for (int sp_I = 0; sp_I < Nsp_I; ++sp_I) {
        indices_I2D_j_jsp[(labels_S1D_sp[sp_I].twoj_I - 1) / 2].push_back(sp_I);
    }

    // R_y: (+m,-m) × (+m,-m).
    Ry_C3D_2sp_2sp_beta.resize(2 * Nsp_I, 2 * Nsp_I, Nbeta_I);
    Ry_C3D_2sp_2sp_beta.setZero();
    for (int twoj_I = 1; twoj_I <= twojMax_I; twoj_I += 2) {
        RepresentationSpin representation(twoj_I);
        const auto& indices_I1D_jsp = indices_I2D_j_jsp[(twoj_I - 1) / 2];
        for (int beta_I = 0; beta_I < Nbeta_I; ++beta_I) {
            const Eigen::MatrixXcd& d_C2D_m_m = representation.calc_Ry(beta_F1D_beta(beta_I));
            assert(d_C2D_m_m.imag().cwiseAbs().maxCoeff() < 1.0e-12);
            for (int sp2_I : indices_I1D_jsp) {
                for (int sp1_I : indices_I1D_jsp) {
                    // δ_N1,N2 δ_j1,j2.
                    const SphericalSPLabel& label1 = labels_S1D_sp[sp1_I];
                    const SphericalSPLabel& label2 = labels_S1D_sp[sp2_I];
                    const double delta_F = static_cast<double>(label1.N_I == label2.N_I && label1.twoj_I == label2.twoj_I);

                    // M = -j,...,j → indices for ±m.
                    const int m1Pos_I = (twoj_I + label1.twom_I) / 2;
                    const int m1Neg_I = (twoj_I - label1.twom_I) / 2;
                    const int m2Pos_I = (twoj_I + label2.twom_I) / 2;
                    const int m2Neg_I = (twoj_I - label2.twom_I) / 2;

                    // R_y^(s₁s₂) = δ_N1,N2 δ_j1,j2 d^j_(s₁m₁,s₂m₂).
                    Ry_C3D_2sp_2sp_beta(sp1_I, sp2_I, beta_I) = delta_F * d_C2D_m_m(m1Pos_I, m2Pos_I);
                    Ry_C3D_2sp_2sp_beta(sp1_I, sp2_I + Nsp_I, beta_I) = delta_F * d_C2D_m_m(m1Pos_I, m2Neg_I);
                    Ry_C3D_2sp_2sp_beta(sp1_I + Nsp_I, sp2_I, beta_I) = delta_F * d_C2D_m_m(m1Neg_I, m2Pos_I);
                    Ry_C3D_2sp_2sp_beta(sp1_I + Nsp_I, sp2_I + Nsp_I, beta_I) = delta_F * d_C2D_m_m(m1Neg_I, m2Neg_I);
                }
            }
        }
    }
}
