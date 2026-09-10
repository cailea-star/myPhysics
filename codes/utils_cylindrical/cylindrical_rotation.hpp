/**
 * @file    cylindrical_rotation.hpp
 * @author  cailea
 * @date    2026-09-06
 * @brief   Axial harmonic-oscillator rotation matrices.
 */

#pragma once

#include <algorithm>
#include <vector>
#include <cassert>
#include <cmath>
#include <complex>
#include <utility>
#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/LU>
#include <unsupported/Eigen/CXX11/Tensor>
#include <gsl/gsl_sf_gamma.h>

#include "cylindrical_config.hpp"
#include "integration_gauss.hpp"

using doubleC = std::complex<double>;

class CylindricalRotation {
public:
    bool isBuilt_B = false; // True after successful build().
    int Nalpha_I = 1; // Number of α quadrature nodes.
    int Nbeta_I = 1; // Number of β quadrature nodes.
    int Ngamma_I = 1; // Number of γ quadrature nodes.
    Eigen::VectorXd beta_F1D_beta{}; // β_i = acos x_i ∈ (0,π) [rad].
    Eigen::VectorXd weight_F1D_beta{}; // Σ_i w_i f(β_i) ≈ ∫₀^π sinβ f(β)dβ.
    double bz_F = 0.0; // Axial oscillator length.
    double br_F = 0.0; // Radial oscillator length.
    std::vector<CylindricalSPLabel> labels_S1D_sp{}; // Positive-Ω labels.

    // Rows and columns: (+Ω,-Ω).
    Eigen::VectorXcd detRySimplexP_C1D_beta{}; // det d₊(β); det d₋ = (det d₊)*.
    Eigen::Tensor<double, 3, Eigen::ColMajor> Ry_F3D_2sp_2sp_beta{}; // ⟨φ_a|R_y(β_i)|φ_b⟩.
    Eigen::Tensor<doubleC, 3, Eigen::ColMajor> Rz_C3D_2sp_2sp_alpha{}; // δ_ab exp(-iΩ_a α_i).
    Eigen::Tensor<doubleC, 3, Eigen::ColMajor> Rz_C3D_2sp_2sp_gamma{}; // δ_ab exp(-iΩ_a γ_i).

private:
    // Ω partners: (Λ,Σ) ↔ (-Λ,-Σ).
    Eigen::MatrixXd COmegaPos_F2D_sp_ny{}; // C(n_r,Λ,n_y).
    Eigen::MatrixXd COmegaNeg_F2D_sp_ny{}; // C(n_r,-Λ,n_y).

    // K(t;β) = K₀(β) exp(½tᵀH(β)t); t = (t_x1,t_z1,t_x2,t_z2)ᵀ.
    Eigen::VectorXd K0_F1D_beta{};
    Eigen::Tensor<double, 3, Eigen::ColMajor> H_F3D_x1z1x2z2_x1z1x2z2_beta{};

public:
    /**
     * @brief  Cache Gauss-Legendre quadrature, binomial coefficients, and LDLT geometry.
     * @math   x_i = cosβ_i; ∫₀^π sinβ f(β)dβ ≈ Σ_i w_i f(acos x_i).
     * @output Quadrature and caches; rotation tensors remain empty.
     * @note   N_α,N_β,N_γ > 0.
     * @note   Weights exclude projection normalization factors.
     */
    CylindricalRotation(const CylindricalSetting& cylindricalsetting_, const std::vector<CylindricalSPLabel>& labels_S1D_sp_, int Nalpha_I_, int Nbeta_I_, int Ngamma_I_) {
        bz_F = cylindricalsetting_.bz_F;
        br_F = cylindricalsetting_.br_F;
        labels_S1D_sp = labels_S1D_sp_;
        Nalpha_I = Nalpha_I_;
        Nbeta_I = Nbeta_I_;
        Ngamma_I = Ngamma_I_;

        // b_r,b_z > 0; N_α,N_β,N_γ > 0.
        assert(std::isfinite(br_F) && br_F > 0.0);
        assert(std::isfinite(bz_F) && bz_F > 0.0);
        assert(Nalpha_I > 0 && Nbeta_I > 0 && Ngamma_I > 0);

        // Gauss-Legendre: x_i = cosβ_i; Σ_i w_i = 2.
        const GaussLegendreMeshes legendre_meshes(Nbeta_I);
        beta_F1D_beta = legendre_meshes.x_F1D_x.array().acos().matrix();
        weight_F1D_beta = legendre_meshes.w_F1D_x;

        // ν_sp = 2n_r+Λ; 0 ≤ n_y ≤ ν_sp.
        const int Nsp_I = static_cast<int>(labels_S1D_sp.size());
        Eigen::VectorXi nu_I1D_sp{};
        nu_I1D_sp.resize(Nsp_I);
        int nuMax_I = 0;
        for (int sp_I = 0; sp_I < Nsp_I; ++sp_I) {
            const CylindricalSPLabel& label_ = labels_S1D_sp[sp_I];
            assert(is_valid(label_.nz_I, label_.nr_I, label_.Lambda_I, label_.twoOmega_I, label_.twoSigma_I, label_.isParityPositive_B));
            nu_I1D_sp(sp_I) = 2 * label_.nr_I + label_.Lambda_I;
            nuMax_I = std::max(nuMax_I, nu_I1D_sp(sp_I));
        }

        // {C⁺,C⁻,K₀,H} → allocated storage.
        COmegaPos_F2D_sp_ny.resize(Nsp_I, nuMax_I + 1);
        COmegaNeg_F2D_sp_ny.resize(Nsp_I, nuMax_I + 1);
        K0_F1D_beta.resize(Nbeta_I);
        H_F3D_x1z1x2z2_x1z1x2z2_beta.resize(4, 4, Nbeta_I);

        // {C⁺,C⁻,K₀,H} → 0.
        COmegaPos_F2D_sp_ny.setZero();
        COmegaNeg_F2D_sp_ny.setZero();
        K0_F1D_beta.setZero();
        H_F3D_x1z1x2z2_x1z1x2z2_beta.setZero();

        // C_sp,ny^± = C(n_r,±Λ,n_y).
        for (int sp_I = 0; sp_I < Nsp_I; ++sp_I) {
            const CylindricalSPLabel& label_ = labels_S1D_sp[sp_I];
            for (int ny_I = 0; ny_I <= nu_I1D_sp(sp_I); ++ny_I) {
                COmegaPos_F2D_sp_ny(sp_I, ny_I) = calc_polar_cartesian_coeff(label_.nr_I, label_.Lambda_I, ny_I);
                COmegaNeg_F2D_sp_ny(sp_I, ny_I) = calc_polar_cartesian_coeff(label_.nr_I, -label_.Lambda_I, ny_I);
            }
        }

        // (b_r,b_z,β) → cached (H,K₀).
        for (int beta_I = 0; beta_I < Nbeta_I; ++beta_I) {
            const auto [H_F2D_t_t, K0_F] = calc_H_K0(bz_F, br_F, beta_F1D_beta(beta_I));
            Eigen::Map<Eigen::Matrix4d>(H_F3D_x1z1x2z2_x1z1x2z2_beta.data() + 16 * static_cast<Eigen::Index>(beta_I)) = H_F2D_t_t;
            K0_F1D_beta(beta_I) = K0_F;
        }
    }

    /**
     * @brief  Build rotations using generating functions.
     * @math   d_ab(β) = (C† R_Cartesian C)_ab d_spin(β)
     * @output R_y(β) and det d₊(β), ordered (+Ω,-Ω).
     * @note   Deformed truncation generally breaks matrix unitarity.
     * @note   Changed configuration or angles require reconstructing CylindricalRotation.
     */
    void build_Ry();

    /**
     * @brief  Build diagonal z rotations using trigonometric phases.
     * @math   R_z(θ)_ab = δ_ab exp(-iΩ_a θ).
     * @output Cached Ω-basis z-rotation matrices.
     * @note   α_i = 2πi/N_α; γ_i = 4πi/N_γ; right endpoints excluded.
     */
    void build_Rz();

    /**
     * @brief  Build y and z rotation caches.
     * @math   {β,N_α,N_γ} → {R_y,R_z}.
     * @output Complete rotation caches and readiness.
     */
    void build();

    /**
     * @brief  Evaluate Gaussian geometry using LDLT factorization.
     * @math   H = 2BᵀM⁻¹B-I; K₀ = 2/(b_r b_z √det M).
     * @output Pair {H,K₀}.
     */
    static std::pair<Eigen::Matrix4d, double> calc_H_K0(double bz_F, double br_F, double beta_F);

    /**
     * @brief  Evaluate Cartesian rotations by cached Gaussian recurrence.
     * @math   d = ⟨n_x1,n_z1|exp(-iβL_y/ℏ)|n_x2,n_z2⟩.
     * @output Real orbital rotation matrix element.
     * @note   The omitted y factor is δ_ny1,ny2.
     */
    double calc_Ry_cartesian(int nx1_I, int nz1_I, int nx2_I, int nz2_I, int beta_I) const;


    /**
     * @brief  Evaluate polar-to-Cartesian coefficients by binomial expansion.
     * @math   ⟨n_x,n_y|n_r,Λ⟩ = (-i)^n_y C.
     * @output Real expansion coefficient.
     */
    static double calc_polar_cartesian_coeff(int nr_I, int Lambda_I, int ny_I);
};

inline void CylindricalRotation::build() {
    // Ready only after both rotation caches complete.
    isBuilt_B = false;
    build_Ry();
    build_Rz();
    isBuilt_B = true;
}

inline void CylindricalRotation::build_Rz() {
    // α ∈ [0,2π); γ ∈ [0,4π).
    const Eigen::Index Nsp_I = static_cast<Eigen::Index>(labels_S1D_sp.size());
    const double pi_F = std::acos(-1.0);
    assert(Nalpha_I > 0 && Ngamma_I > 0);

    // {R_z(α),R_z(γ)} → allocated storage → 0.
    Rz_C3D_2sp_2sp_alpha.resize(2 * Nsp_I, 2 * Nsp_I, Nalpha_I);
    Rz_C3D_2sp_2sp_gamma.resize(2 * Nsp_I, 2 * Nsp_I, Ngamma_I);
    Rz_C3D_2sp_2sp_alpha.setZero();
    Rz_C3D_2sp_2sp_gamma.setZero();

    // R_z(α) = diag(exp(-iΩα),exp(+iΩα)).
    for (int alpha_I = 0; alpha_I < Nalpha_I; ++alpha_I) {
        const double alpha_F = 2.0 * pi_F * alpha_I / Nalpha_I;
        for (Eigen::Index sp_I = 0; sp_I < Nsp_I; ++sp_I) {
            const double Omega_F = 0.5 * labels_S1D_sp[sp_I].twoOmega_I;
            const double cosOmegaAngle_F = std::cos(Omega_F * alpha_F);
            const double sinOmegaAngle_F = std::sin(Omega_F * alpha_F);

            // (+Ω,-Ω) → (exp(-iΩα),exp(+iΩα)).
            Rz_C3D_2sp_2sp_alpha(sp_I, sp_I, alpha_I) = doubleC(cosOmegaAngle_F, -sinOmegaAngle_F);
            Rz_C3D_2sp_2sp_alpha(sp_I + Nsp_I, sp_I + Nsp_I, alpha_I) = doubleC(cosOmegaAngle_F, sinOmegaAngle_F);
        }
    }

    // R_z(γ) = diag(exp(-iΩγ),exp(+iΩγ)).
    for (int gamma_I = 0; gamma_I < Ngamma_I; ++gamma_I) {
        const double gamma_F = 4.0 * pi_F * gamma_I / Ngamma_I;
        for (Eigen::Index sp_I = 0; sp_I < Nsp_I; ++sp_I) {
            const double Omega_F = 0.5 * labels_S1D_sp[sp_I].twoOmega_I;
            const double cosOmegaAngle_F = std::cos(Omega_F * gamma_F);
            const double sinOmegaAngle_F = std::sin(Omega_F * gamma_F);

            // (+Ω,-Ω) → (exp(-iΩγ),exp(+iΩγ)).
            Rz_C3D_2sp_2sp_gamma(sp_I, sp_I, gamma_I) = doubleC(cosOmegaAngle_F, -sinOmegaAngle_F);
            Rz_C3D_2sp_2sp_gamma(sp_I + Nsp_I, sp_I + Nsp_I, gamma_I) = doubleC(cosOmegaAngle_F, sinOmegaAngle_F);
        }
    }
}

inline void CylindricalRotation::build_Ry() {
    // {Ω > 0,β} → (N_sp,N_β).
    const int Nsp_I = static_cast<int>(labels_S1D_sp.size());
    const int Nbeta_I = static_cast<int>(beta_F1D_beta.size());
    assert(std::isfinite(br_F) && br_F > 0.0);
    assert(std::isfinite(bz_F) && bz_F > 0.0);
    assert(Nbeta_I > 0);
    assert(beta_F1D_beta.allFinite());
    assert((beta_F1D_beta.array() > 0.0).all() && (beta_F1D_beta.array() < std::acos(-1.0)).all());
    assert(K0_F1D_beta.size() == Nbeta_I);
    assert(H_F3D_x1z1x2z2_x1z1x2z2_beta.dimension(0) == 4 && H_F3D_x1z1x2z2_x1z1x2z2_beta.dimension(1) == 4 && H_F3D_x1z1x2z2_x1z1x2z2_beta.dimension(2) == Nbeta_I);

    // ν_sp = 2n_r+Λ; 0 ≤ n_y ≤ ν_sp.
    Eigen::VectorXi nu_I1D_sp{};
    nu_I1D_sp.resize(Nsp_I);
    int nuMax_I = 0;
    for (int sp_I = 0; sp_I < Nsp_I; ++sp_I) {
        const CylindricalSPLabel& label_ = labels_S1D_sp[sp_I];
        assert(is_valid(label_.nz_I, label_.nr_I, label_.Lambda_I, label_.twoOmega_I, label_.twoSigma_I, label_.isParityPositive_B));
        nu_I1D_sp(sp_I) = 2 * label_.nr_I + label_.Lambda_I;
        nuMax_I = std::max(nuMax_I, nu_I1D_sp(sp_I));
    }

    assert(COmegaPos_F2D_sp_ny.rows() == Nsp_I && COmegaPos_F2D_sp_ny.cols() == nuMax_I + 1);
    assert(COmegaNeg_F2D_sp_ny.rows() == Nsp_I && COmegaNeg_F2D_sp_ny.cols() == nuMax_I + 1);

    // R_y: (+Ω,-Ω) × (+Ω,-Ω).
    Eigen::MatrixXd Ry_F2D_2sp_2sp{};
    Ry_F2D_2sp_2sp.resize(2 * static_cast<Eigen::Index>(Nsp_I), 2 * static_cast<Eigen::Index>(Nsp_I));
    Ry_F3D_2sp_2sp_beta.resize(2 * static_cast<Eigen::Index>(Nsp_I), 2 * static_cast<Eigen::Index>(Nsp_I), Nbeta_I);
    Ry_F2D_2sp_2sp.setZero();
    Ry_F3D_2sp_2sp_beta.setZero();

    for (int beta_I = 0; beta_I < Nbeta_I; ++beta_I) {
        // d_spin: (+½,-½) × (+½,-½).
        const double beta_F = beta_F1D_beta(beta_I);
        const double cosHalfBeta_F = std::cos(0.5 * beta_F);
        const double sinHalfBeta_F = std::sin(0.5 * beta_F);
        Eigen::Matrix2d RySpin_F2D_spin_spin{};
        RySpin_F2D_spin_spin(0, 0) = cosHalfBeta_F;
        RySpin_F2D_spin_spin(0, 1) = -sinHalfBeta_F;
        RySpin_F2D_spin_spin(1, 0) = sinHalfBeta_F;
        RySpin_F2D_spin_spin(1, 1) = cosHalfBeta_F;

        for (int sp2_I = 0; sp2_I < Nsp_I; ++sp2_I) {
            for (int sp1_I = 0; sp1_I < Nsp_I; ++sp1_I) {
                // (sp1,sp2) → labels, spin indices, shared n_y range.
                const CylindricalSPLabel& label1_ = labels_S1D_sp[sp1_I];
                const CylindricalSPLabel& label2_ = labels_S1D_sp[sp2_I];
                const int spin1_I = (1 - label1_.twoSigma_I) / 2;
                const int spin2_I = (1 - label2_.twoSigma_I) / 2;
                const int nyMax_I = std::min(nu_I1D_sp(sp1_I), nu_I1D_sp(sp2_I));

                // R^{st} = Σ_ny C_1^s R_Cartesian C_2^t.
                double RyOrbitalPosPos_F = 0.0;
                double RyOrbitalPosNeg_F = 0.0;
                double RyOrbitalNegPos_F = 0.0;
                double RyOrbitalNegNeg_F = 0.0;
                for (int ny_I = 0; ny_I <= nyMax_I; ++ny_I) {
                    const int nx1_I = nu_I1D_sp(sp1_I) - ny_I;
                    const int nx2_I = nu_I1D_sp(sp2_I) - ny_I;

                    // One Cartesian element serves all four sign channels.
                    const double RyCartesian_F = calc_Ry_cartesian(nx1_I, label1_.nz_I, nx2_I, label2_.nz_I, beta_I);
                    RyOrbitalPosPos_F += COmegaPos_F2D_sp_ny(sp1_I, ny_I) * RyCartesian_F * COmegaPos_F2D_sp_ny(sp2_I, ny_I);
                    RyOrbitalPosNeg_F += COmegaPos_F2D_sp_ny(sp1_I, ny_I) * RyCartesian_F * COmegaNeg_F2D_sp_ny(sp2_I, ny_I);
                    RyOrbitalNegPos_F += COmegaNeg_F2D_sp_ny(sp1_I, ny_I) * RyCartesian_F * COmegaPos_F2D_sp_ny(sp2_I, ny_I);
                    RyOrbitalNegNeg_F += COmegaNeg_F2D_sp_ny(sp1_I, ny_I) * RyCartesian_F * COmegaNeg_F2D_sp_ny(sp2_I, ny_I);
                }

                // (Λ,Σ) → (-Λ,-Σ); d = R_orbital d_spin.
                Ry_F2D_2sp_2sp(sp1_I, sp2_I) = RyOrbitalPosPos_F * RySpin_F2D_spin_spin(spin1_I, spin2_I);
                Ry_F2D_2sp_2sp(sp1_I, sp2_I + Nsp_I) = RyOrbitalPosNeg_F * RySpin_F2D_spin_spin(spin1_I, 1 - spin2_I);
                Ry_F2D_2sp_2sp(sp1_I + Nsp_I, sp2_I) = RyOrbitalNegPos_F * RySpin_F2D_spin_spin(1 - spin1_I, spin2_I);
                Ry_F2D_2sp_2sp(sp1_I + Nsp_I, sp2_I + Nsp_I) = RyOrbitalNegNeg_F * RySpin_F2D_spin_spin(1 - spin1_I, 1 - spin2_I);

            }
        }
        // ColMajor slices: offset = (2N_sp)² i_β.
        Eigen::Map<Eigen::MatrixXd>(Ry_F3D_2sp_2sp_beta.data() + 4 * static_cast<Eigen::Index>(Nsp_I) * Nsp_I * beta_I, 2 * static_cast<Eigen::Index>(Nsp_I), 2 * static_cast<Eigen::Index>(Nsp_I)) = Ry_F2D_2sp_2sp;
    }

    // T = [I,I;-iE,iE]/√2; E_aa = 2Σ_a.
    Eigen::MatrixXcd T_C2D_Omega2sp_Simplex2sp{};
    Eigen::MatrixXcd RySimplex_C2D_2sp_2sp{};
    T_C2D_Omega2sp_Simplex2sp.resize(2 * Nsp_I, 2 * Nsp_I);
    RySimplex_C2D_2sp_2sp.resize(2 * Nsp_I, 2 * Nsp_I);
    detRySimplexP_C1D_beta.resize(Nbeta_I);
    T_C2D_Omega2sp_Simplex2sp.setZero();

    const double invSqrt2_F = 1.0 / std::sqrt(2.0);
    for (int sp_I = 0; sp_I < Nsp_I; ++sp_I) {
        // T = [I,I;-iE,iE]/√2.
        T_C2D_Omega2sp_Simplex2sp(sp_I, sp_I) = invSqrt2_F;
        T_C2D_Omega2sp_Simplex2sp(sp_I, sp_I + Nsp_I) = invSqrt2_F;
        T_C2D_Omega2sp_Simplex2sp(sp_I + Nsp_I, sp_I) = doubleC(0.0, -labels_S1D_sp[sp_I].twoSigma_I * invSqrt2_F);
        T_C2D_Omega2sp_Simplex2sp(sp_I + Nsp_I, sp_I + Nsp_I) = doubleC(0.0, labels_S1D_sp[sp_I].twoSigma_I * invSqrt2_F);
    }

    for (int beta_I = 0; beta_I < Nbeta_I; ++beta_I) {
        // R_y(β): contiguous ColMajor slice.
        const Eigen::Map<const Eigen::MatrixXd> Ry_F2D_2sp_2sp(Ry_F3D_2sp_2sp_beta.data() + 4 * static_cast<Eigen::Index>(Nsp_I) * Nsp_I * beta_I, 2 * Nsp_I, 2 * Nsp_I);

        // T†R_yT = diag(d₊,d₋).
        RySimplex_C2D_2sp_2sp.noalias() = T_C2D_Omega2sp_Simplex2sp.adjoint() * Ry_F2D_2sp_2sp.cast<doubleC>() * T_C2D_Omega2sp_Simplex2sp;
        assert(RySimplex_C2D_2sp_2sp.allFinite());
        detRySimplexP_C1D_beta(beta_I) = RySimplex_C2D_2sp_2sp.topLeftCorner(Nsp_I, Nsp_I).determinant();
    }
}

inline double CylindricalRotation::calc_Ry_cartesian(int nx1_I, int nz1_I, int nx2_I, int nz2_I, int beta_I) const {
    // {n_x1,n_z1,n_x2,n_z2} ≥ 0; 0 ≤ i_β < N_β.
    assert(nx1_I >= 0 && nz1_I >= 0 && nx2_I >= 0 && nz2_I >= 0);
    assert(beta_I >= 0 && beta_I < K0_F1D_beta.size());
    assert(K0_F1D_beta.size() == beta_F1D_beta.size());
    assert(H_F3D_x1z1x2z2_x1z1x2z2_beta.dimension(0) == 4 && H_F3D_x1z1x2z2_x1z1x2z2_beta.dimension(1) == 4 && H_F3D_x1z1x2z2_x1z1x2z2_beta.dimension(2) == K0_F1D_beta.size());

    // K(t;β) = K₀(β) exp(½tᵀH(β)t).
    const Eigen::Map<const Eigen::Matrix4d> H_F2D_t_t(H_F3D_x1z1x2z2_x1z1x2z2_beta.data() + 16 * static_cast<Eigen::Index>(beta_I));
    const double K0_F = K0_F1D_beta(beta_I);

    // d₀ = K₀; target domain includes every predecessor.
    Eigen::Tensor<double, 4, Eigen::ColMajor> d_F4D_nx1_nz1_nx2_nz2{};
    d_F4D_nx1_nz1_nx2_nz2.resize(nx1_I + 1, nz1_I + 1, nx2_I + 1, nz2_I + 1);
    d_F4D_nx1_nz1_nx2_nz2.setZero();
    d_F4D_nx1_nz1_nx2_nz2(0, 0, 0, 0) = K0_F;

    // √n_i d_n = Σ_j H_ij √(n_j-δ_ij) d_(n-e_i-e_j).
    const auto calc_next_Func = [&](int nx1_I_, int nz1_I_, int nx2_I_, int nz2_I_, int axis_I) {
        Eigen::Vector4i n_I1D_t(nx1_I_, nz1_I_, nx2_I_, nz2_I_);
        const double denominator_F = std::sqrt(static_cast<double>(n_I1D_t(axis_I)));
        --n_I1D_t(axis_I);
        double sum_F = 0.0;
        for (int t_I = 0; t_I < 4; ++t_I) {
            Eigen::Vector4i lower_I1D_t = n_I1D_t;
            lower_I1D_t(t_I) = std::max(0, n_I1D_t(t_I) - 1);
            sum_F += H_F2D_t_t(axis_I, t_I) * std::sqrt(static_cast<double>(n_I1D_t(t_I))) * d_F4D_nx1_nz1_nx2_nz2(lower_I1D_t(0), lower_I1D_t(1), lower_I1D_t(2), lower_I1D_t(3));
        }
        return sum_F / denominator_F;
    };

    // n_x1=n_z1=n_x2=0; increase n_z2.
    for (int nz2_I_ = 1; nz2_I_ <= nz2_I; ++nz2_I_) {d_F4D_nx1_nz1_nx2_nz2(0, 0, 0, nz2_I_) = calc_next_Func(0, 0, 0, nz2_I_, 3);}

    // n_x1=n_z1=0; increase n_x2.
    for (int nx2_I_ = 1; nx2_I_ <= nx2_I; ++nx2_I_) {
        for (int nz2_I_ = 0; nz2_I_ <= nz2_I; ++nz2_I_) {d_F4D_nx1_nz1_nx2_nz2(0, 0, nx2_I_, nz2_I_) = calc_next_Func(0, 0, nx2_I_, nz2_I_, 2);}
    }

    // n_x1=0; increase n_z1.
    for (int nz1_I_ = 1; nz1_I_ <= nz1_I; ++nz1_I_) {
        for (int nx2_I_ = 0; nx2_I_ <= nx2_I; ++nx2_I_) {
            for (int nz2_I_ = 0; nz2_I_ <= nz2_I; ++nz2_I_) {d_F4D_nx1_nz1_nx2_nz2(0, nz1_I_, nx2_I_, nz2_I_) = calc_next_Func(0, nz1_I_, nx2_I_, nz2_I_, 1);}
        }
    }

    // Increase n_x1 over the remaining domain.
    for (int nx1_I_ = 1; nx1_I_ <= nx1_I; ++nx1_I_) {
        for (int nz1_I_ = 0; nz1_I_ <= nz1_I; ++nz1_I_) {
            for (int nx2_I_ = 0; nx2_I_ <= nx2_I; ++nx2_I_) {
                for (int nz2_I_ = 0; nz2_I_ <= nz2_I; ++nz2_I_) {d_F4D_nx1_nz1_nx2_nz2(nx1_I_, nz1_I_, nx2_I_, nz2_I_) = calc_next_Func(nx1_I_, nz1_I_, nx2_I_, nz2_I_, 0);}
            }
        }
    }

    return d_F4D_nx1_nz1_nx2_nz2(nx1_I, nz1_I, nx2_I, nz2_I);
}

inline std::pair<Eigen::Matrix4d, double> CylindricalRotation::calc_H_K0(double bz_F, double br_F, double beta_F) {
    // b_r,b_z > 0; β ∈ ℝ.
    assert(std::isfinite(br_F) && br_F > 0.0);
    assert(std::isfinite(bz_F) && bz_F > 0.0);
    assert(std::isfinite(beta_F));

    // (b_r,b_z) → (b_r⁻¹,b_z⁻¹,b_r⁻²,b_z⁻²).
    const double brInv_F = 1.0 / br_F;
    const double bzInv_F = 1.0 / bz_F;
    const double brInv2_F = brInv_F * brInv_F;
    const double bzInv2_F = bzInv_F * bzInv_F;

    // R_y⁻¹: x′ = cx-sz, z′ = sx+cz.
    const double cosBeta_F = std::cos(beta_F);
    const double sinBeta_F = std::sin(beta_F);

    // Gaussian exponent: -½(x,z)M(x,z)ᵀ + √2(x,z)Bt.
    Eigen::Matrix2d M_F2D_xz_xz{};
    M_F2D_xz_xz(0, 0) = (1.0 + cosBeta_F * cosBeta_F) * brInv2_F + sinBeta_F * sinBeta_F * bzInv2_F;
    M_F2D_xz_xz(0, 1) = cosBeta_F * sinBeta_F * (bzInv2_F - brInv2_F);
    M_F2D_xz_xz(1, 0) = M_F2D_xz_xz(0, 1);
    M_F2D_xz_xz(1, 1) = sinBeta_F * sinBeta_F * brInv2_F + (1.0 + cosBeta_F * cosBeta_F) * bzInv2_F;

    // t = (t_x1,t_z1,t_x2,t_z2)ᵀ.
    Eigen::Matrix<double, 2, 4> B_F2D_xz_t{};
    B_F2D_xz_t(0, 0) = brInv_F;
    B_F2D_xz_t(0, 1) = 0.0;
    B_F2D_xz_t(0, 2) = cosBeta_F * brInv_F;
    B_F2D_xz_t(0, 3) = sinBeta_F * bzInv_F;
    B_F2D_xz_t(1, 0) = 0.0;
    B_F2D_xz_t(1, 1) = bzInv_F;
    B_F2D_xz_t(1, 2) = -sinBeta_F * brInv_F;
    B_F2D_xz_t(1, 3) = cosBeta_F * bzInv_F;

    // H = 2BᵀM⁻¹B-I; K₀ = 2/(b_r b_z √det M).
    const Eigen::Matrix<double, 2, 4> MinvB_F2D_xz_t = M_F2D_xz_xz.ldlt().solve(B_F2D_xz_t);
    Eigen::Matrix4d H_F2D_t_t = 2.0 * B_F2D_xz_t.transpose() * MinvB_F2D_xz_t;
    H_F2D_t_t.diagonal().array() -= 1.0;
    const double K0_F = 2.0 * brInv_F * bzInv_F / std::sqrt(M_F2D_xz_xz.determinant());
    return {H_F2D_t_t, K0_F};
}

inline double CylindricalRotation::calc_polar_cartesian_coeff(int nr_I, int Lambda_I, int ny_I) {
    // (n_r,Λ,n_y) → (n_x,q_min,q_max,C).
    assert(ny_I >= 0 && ny_I <= 2 * nr_I + std::abs(Lambda_I));
    const int nx_I = 2 * nr_I + std::abs(Lambda_I) - ny_I;
    const int qMin_I = std::max(0, nr_I + (std::abs(Lambda_I) - Lambda_I) / 2 - nx_I);
    const int qMax_I = std::min(ny_I, nr_I + (std::abs(Lambda_I) - Lambda_I) / 2);
    const double prefactor_F = std::pow(-1, nr_I) * std::pow(2.0, -nr_I - 0.5 * std::abs(Lambda_I)) * std::exp(0.5 * (gsl_sf_lnfact(nr_I + std::abs(Lambda_I)) + gsl_sf_lnfact(nr_I) - gsl_sf_lnfact(nx_I) - gsl_sf_lnfact(ny_I)));
    double sum_F = 0.0;
    for (int q_I = qMin_I; q_I <= qMax_I; ++q_I) {
        const int kx_I = nr_I - q_I + (std::abs(Lambda_I) - Lambda_I) / 2;
        sum_F += gsl_sf_choose(static_cast<unsigned int>(nx_I), static_cast<unsigned int>(kx_I)) * gsl_sf_choose(static_cast<unsigned int>(ny_I), static_cast<unsigned int>(q_I)) * std::pow(-1, ny_I - q_I);
    }
    return prefactor_F * sum_F;
}
