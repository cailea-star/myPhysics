/**
 * @file    hfb_projection.hpp
 * @author  cailea
 * @date    2026-09-11
 * @brief   Angular-momentum projection of HFB configurations.
 */

#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <numbers>
#include <vector>

#include <gsl/gsl_sf_coupling.h>

#include <Eigen/Core>
#include <unsupported/Eigen/CXX11/Tensor>

#include "group_so3_su2.hpp"
#include "hfb_projection_pnp.hpp"

/**
 * @brief  Project HFB configurations onto particle number and angular momentum.
 * @math   N = TargetN; 2J = TargetTwoJ; K = -J+k, k = 0,…,2J.
 * @note   Rotations and U,V share one ordered single-particle basis.
 * @note   Weights exclude Wigner factors and projector normalization.
 * @note   All calc methods require nonzero vacuum overlaps at quadrature nodes.
 * @note   Returned references alias reusable output tensors.
 */
class HFBProjection : public HFBProjectionPNP {
public:
    int TargetTwoI_I = 0;

    int Nalpha_I = 0;
    int Nbeta_I = 0;
    int Ngamma_I = 0;

    Eigen::Tensor<doubleC, 3, Eigen::ColMajor> Rz_C3D_sp_sp_alpha{};
    Eigen::Tensor<doubleC, 3, Eigen::ColMajor> Ry_C3D_sp_sp_beta{};
    Eigen::Tensor<doubleC, 3, Eigen::ColMajor> Rz_C3D_sp_sp_gamma{};

    Eigen::VectorXd alpha_F1D_alpha{};
    Eigen::VectorXd beta_F1D_beta{};
    Eigen::VectorXd gamma_F1D_gamma{};

    Eigen::VectorXd weight_F1D_alpha{};
    Eigen::VectorXd weight_F1D_beta{};
    Eigen::VectorXd weight_F1D_gamma{};

private:
    RepresentationSpin representation_spin;

    Eigen::Tensor<doubleC, 4, Eigen::ColMajor> result_C4D_cfg1_cfg2_K_K{};

public:
    /**
     * @brief  Allocate projection and Pfaffian workspaces.
     * @math   φ_k = 2πk/Nphi; Δφ/(2π) = 1/Nphi.
     * @output Stored dimensions, target quantum numbers, and allocated workspaces.
     * @note   Initialize projection data; left/right configuration lists may differ.
     */
    HFBProjection(int TargetN_I_, int TargetTwoI_I_, int Nsp_I_, const std::vector<std::vector<int>>& config1_I2D_cfg1_cqp1_, const std::vector<std::vector<int>>& config2_I2D_cfg2_cqp2_, int Nphi_I_, int Nalpha_I_, int Nbeta_I_, int Ngamma_I_)
    : HFBProjectionPNP(TargetN_I_, Nsp_I_, config1_I2D_cfg1_cqp1_, config2_I2D_cfg2_cqp2_, Nphi_I_), representation_spin(TargetTwoI_I_) {
        assert(TargetTwoI_I_ >= 0);
        assert(Nalpha_I_ > 0 && Nbeta_I_ > 0 && Ngamma_I_ > 0);

        // 2I → angular-momentum target.
        TargetTwoI_I = TargetTwoI_I_;

        // Nα,Nβ,Nγ → Euler mesh sizes.
        Nalpha_I = Nalpha_I_;
        Nbeta_I = Nbeta_I_;
        Ngamma_I = Ngamma_I_;

        // Rz(α),Ry(β),Rz(γ) ∈ ℂ^{Nsp×Nsp}.
        Rz_C3D_sp_sp_alpha.resize(Nsp_I, Nsp_I, Nalpha_I);
        Ry_C3D_sp_sp_beta.resize(Nsp_I, Nsp_I, Nbeta_I);
        Rz_C3D_sp_sp_gamma.resize(Nsp_I, Nsp_I, Ngamma_I);

        // {α,wα}, {β,wβ}, {γ,wγ}.
        alpha_F1D_alpha.resize(Nalpha_I);
        beta_F1D_beta.resize(Nbeta_I);
        gamma_F1D_gamma.resize(Ngamma_I);
        weight_F1D_alpha.resize(Nalpha_I);
        weight_F1D_beta.resize(Nbeta_I);
        weight_F1D_gamma.resize(Ngamma_I);

        // result ∈ ℂ^{Ncfg1×Ncfg2×(2J+1)×(2J+1)}.
        result_C4D_cfg1_cfg2_K_K.resize(Ncfg1_I, Ncfg2_I, TargetTwoI_I + 1, TargetTwoI_I + 1);
    }

    /**
     * @brief  Store alpha rotations, nodes, and quadrature weights.
     * @math   Rz ∈ ℂ^{Nsp×Nsp×Nalpha}; weights represent dα.
     * @output Updated Rz(α), alpha nodes, and weights.
     */
    void update_alpha(const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& Rz_C3D_sp_sp_alpha_, const Eigen::VectorXd& alpha_F1D_alpha_, const Eigen::VectorXd& weight_F1D_alpha_);

    /**
     * @brief  Store beta rotations, nodes, and quadrature weights.
     * @math   Ry ∈ ℂ^{Nsp×Nsp×Nbeta}; weights represent sinβ dβ.
     * @output Updated Ry(β), beta nodes, and weights.
     */
    void update_beta(const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& Ry_C3D_sp_sp_beta_, const Eigen::VectorXd& beta_F1D_beta_, const Eigen::VectorXd& weight_F1D_beta_);

    /**
     * @brief  Store gamma rotations, nodes, and quadrature weights.
     * @math   Rz ∈ ℂ^{Nsp×Nsp×Ngamma}; weights represent dγ.
     * @output Updated Rz(γ), gamma nodes, and weights.
     */
    void update_gamma(const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& Rz_C3D_sp_sp_gamma_, const Eigen::VectorXd& gamma_F1D_gamma_, const Eigen::VectorXd& weight_F1D_gamma_);

    /**
     * @brief  Cache rotation using ordered matrix products.
     * @math   D = Rz(α) Ry(β) Rz(γ).
     * @output Updated inherited D_C2D_sp_sp.
     * @note   Requires initialized rotation tables.
     */
    void prepare_rotation(int alpha_I, int beta_I, int gamma_I);

    /**
     * @brief  Integrate overlaps using Euler-angle and uniform gauge quadrature.
     * @math   result_abk₁k₂ = ⟨Φ₁;a|Pᴺ Pᴶ_{K₁K₂}|Φ₂;b⟩.
     * @output Updated result_C4D_cfg1_cfg2_K_K and its const reference.
     * @note   VΩ = Σwα Σwβ Σwγ for complete integration domains.
     * @note   Integer J permits SO(3); half-integer J requires SU(2).
     */
    const Eigen::Tensor<doubleC, 4, Eigen::ColMajor>& calc_overlap();

    /**
     * @brief  Integrate one-body kernels using Euler-angle and uniform gauge quadrature.
     * @math   H¹_ab(g) = Σ_ij OneBody_ij OBTD_ab(i,j;g).
     * @math   result_abk₁k₂ = ⟨Φ₁;a|H¹ Pᴺ Pᴶ_{K₁K₂}|Φ₂;b⟩.
     * @output Updated result_C4D_cfg1_cfg2_K_K and its const reference.
     * @note   Requires rotational invariance and particle-number conservation.
     */
    const Eigen::Tensor<doubleC, 4, Eigen::ColMajor>& calc_one_body(const Eigen::MatrixXd& OneBody_F2D_sp_sp);

    /**
     * @brief  Integrate two-body kernels using Euler-angle and uniform gauge quadrature.
     * @math   H²_ab(g) = ½Σ_ijkl TwoBody_ijkl TBTD_ab(i,j,k,l;g).
     * @math   result_abk₁k₂ = ⟨Φ₁;a|H² Pᴺ Pᴶ_{K₁K₂}|Φ₂;b⟩.
     * @output Updated result_C4D_cfg1_cfg2_K_K and its const reference.
     * @note   TwoBody contains unsymmetrized matrix elements ⟨ij|v|kl⟩.
     * @note   Requires rotational invariance and particle-number conservation.
     */
    const Eigen::Tensor<doubleC, 4, Eigen::ColMajor>& calc_two_body(const Eigen::Tensor<double, 4, Eigen::ColMajor>& TwoBody_F4D_sp_sp_sp_sp);

    /**
     * @brief  Compute reduced multipoles using quadrature and Clebsch–Gordan coupling.
     * @math   multipole_abk₁k₂ = ⟨Φ₁^{J₁;N};K₁,a‖Q_λ‖Φ₂^{J₂;N};K₂,b⟩.
     * @math   dim = Ncfg1 × Ncfg2 × (TwoJ1+1) × (TwoJ2+1).
     * @output Filled caller-owned tensor: cfg1, cfg2, K1, K2.
     * @note   Output dimensions must match Ncfg1,Ncfg2,TwoI1+1,TwoI2+1.
     * @note   twoLambda_I = 2λ; μ = -λ,…,λ.
     * @note   Number-conserving spherical tensor; unnormalized, unmixed configuration states.
     */
    void calc_multipole(int TwoI1_I, int TwoI2_I, int twoLambda_I, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& Multipole_C3D_sp_sp_mu, Eigen::Tensor<doubleC, 4, Eigen::ColMajor>& multipole_C4D_cfg1_cfg2_K1_K2);
};

inline void HFBProjection::update_alpha(const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& Rz_C3D_sp_sp_alpha_, const Eigen::VectorXd& alpha_F1D_alpha_, const Eigen::VectorXd& weight_F1D_alpha_) {
    assert(Rz_C3D_sp_sp_alpha_.dimension(0) == Nsp_I && Rz_C3D_sp_sp_alpha_.dimension(1) == Nsp_I && Rz_C3D_sp_sp_alpha_.dimension(2) == Nalpha_I);
    assert(alpha_F1D_alpha_.size() == Nalpha_I && weight_F1D_alpha_.size() == Nalpha_I);

    // {Rz(α),α,wα} → stored quadrature data.
    Rz_C3D_sp_sp_alpha = Rz_C3D_sp_sp_alpha_;
    alpha_F1D_alpha = alpha_F1D_alpha_;
    weight_F1D_alpha = weight_F1D_alpha_;
}

inline void HFBProjection::update_beta(const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& Ry_C3D_sp_sp_beta_, const Eigen::VectorXd& beta_F1D_beta_, const Eigen::VectorXd& weight_F1D_beta_) {
    assert(Ry_C3D_sp_sp_beta_.dimension(0) == Nsp_I && Ry_C3D_sp_sp_beta_.dimension(1) == Nsp_I && Ry_C3D_sp_sp_beta_.dimension(2) == Nbeta_I);
    assert(beta_F1D_beta_.size() == Nbeta_I && weight_F1D_beta_.size() == Nbeta_I);

    // {Ry(β),β,wβ} → stored quadrature data.
    Ry_C3D_sp_sp_beta = Ry_C3D_sp_sp_beta_;
    beta_F1D_beta = beta_F1D_beta_;
    weight_F1D_beta = weight_F1D_beta_;
}

inline void HFBProjection::update_gamma(const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& Rz_C3D_sp_sp_gamma_, const Eigen::VectorXd& gamma_F1D_gamma_, const Eigen::VectorXd& weight_F1D_gamma_) {
    assert(Rz_C3D_sp_sp_gamma_.dimension(0) == Nsp_I && Rz_C3D_sp_sp_gamma_.dimension(1) == Nsp_I && Rz_C3D_sp_sp_gamma_.dimension(2) == Ngamma_I);
    assert(gamma_F1D_gamma_.size() == Ngamma_I && weight_F1D_gamma_.size() == Ngamma_I);

    // {Rz(γ),γ,wγ} → stored quadrature data.
    Rz_C3D_sp_sp_gamma = Rz_C3D_sp_sp_gamma_;
    gamma_F1D_gamma = gamma_F1D_gamma_;
    weight_F1D_gamma = weight_F1D_gamma_;
}

inline void HFBProjection::prepare_rotation(int alpha_I, int beta_I, int gamma_I) {
    assert(alpha_I >= 0 && alpha_I < Nalpha_I);
    assert(beta_I >= 0 && beta_I < Nbeta_I);
    assert(gamma_I >= 0 && gamma_I < Ngamma_I);

    // D = Rz(α) Ry(β) Rz(γ).
    const Eigen::Map<const Eigen::MatrixXcd> RzAlpha_C2D_sp_sp(Rz_C3D_sp_sp_alpha.data() + static_cast<Eigen::Index>(alpha_I) * Nsp_I * Nsp_I, Nsp_I, Nsp_I);
    const Eigen::Map<const Eigen::MatrixXcd> Ry_C2D_sp_sp(Ry_C3D_sp_sp_beta.data() + static_cast<Eigen::Index>(beta_I) * Nsp_I * Nsp_I, Nsp_I, Nsp_I);
    const Eigen::Map<const Eigen::MatrixXcd> RzGamma_C2D_sp_sp(Rz_C3D_sp_sp_gamma.data() + static_cast<Eigen::Index>(gamma_I) * Nsp_I * Nsp_I, Nsp_I, Nsp_I);
    D_C2D_sp_sp.noalias() = RzAlpha_C2D_sp_sp * Ry_C2D_sp_sp * RzGamma_C2D_sp_sp;
}

inline const Eigen::Tensor<doubleC, 4, Eigen::ColMajor>& HFBProjection::calc_overlap() {
    // (2J+1)/VΩ; weights contain sinβ dβ.
    const double volume_F = weight_F1D_alpha.sum() * weight_F1D_beta.sum() * weight_F1D_gamma.sum();
    assert(std::isfinite(volume_F) && volume_F > 0.0);
    const double normalization_F = (TargetTwoI_I + 1.0) / volume_F;
    result_C4D_cfg1_cfg2_K_K.setZero();

    // dᴶ(β) is reused over α, γ, φ.
    for (int beta_I = 0; beta_I < Nbeta_I; ++beta_I) {
        const auto& RyJ_C2D_K_K = representation_spin.calc_Ry(beta_F1D_beta(beta_I));
        for (int alpha_I = 0; alpha_I < Nalpha_I; ++alpha_I) {
            for (int gamma_I = 0; gamma_I < Ngamma_I; ++gamma_I) {
                prepare_rotation(alpha_I, beta_I, gamma_I);
                // w = (2J+1) wα wβ wγ/VΩ.
                const double weight_F = normalization_F * weight_F1D_alpha(alpha_I) * weight_F1D_beta(beta_I) * weight_F1D_gamma(gamma_I);
                const auto& tmp_C2D_cfg1_cfg2 = calc_overlap_pnp();

                for (int K2_I = 0; K2_I <= TargetTwoI_I; ++K2_I) {
                    const double K2_F = K2_I - 0.5 * TargetTwoI_I;
                    for (int K1_I = 0; K1_I <= TargetTwoI_I; ++K1_I) {
                        // Dᴶ* = exp(iK₁α) dᴶ* exp(iK₂γ).
                        const double K1_F = K1_I - 0.5 * TargetTwoI_I;
                        const doubleC factorOmega_C = std::conj(RyJ_C2D_K_K(K1_I, K2_I)) * std::exp(doubleC(0.0, K1_F * alpha_F1D_alpha(alpha_I) + K2_F * gamma_F1D_gamma(gamma_I)));
                        const doubleC weight_factor_C = weight_F * factorOmega_C;
                        Eigen::Map<Eigen::MatrixXcd> result_C2D_cfg_cfg(result_C4D_cfg1_cfg2_K_K.data() + (static_cast<Eigen::Index>(K2_I) * (TargetTwoI_I + 1) + K1_I) * Ncfg1_I * Ncfg2_I, Ncfg1_I, Ncfg2_I);
                        result_C2D_cfg_cfg += weight_factor_C * tmp_C2D_cfg1_cfg2;
                    }
                }
            }
        }
    }
    return result_C4D_cfg1_cfg2_K_K;
}

inline const Eigen::Tensor<doubleC, 4, Eigen::ColMajor>& HFBProjection::calc_one_body(const Eigen::MatrixXd& OneBody_F2D_sp_sp) {
    assert(OneBody_F2D_sp_sp.rows() == Nsp_I && OneBody_F2D_sp_sp.cols() == Nsp_I && OneBody_F2D_sp_sp.allFinite());

    const Eigen::MatrixXcd OneBody_C2D_sp_sp = OneBody_F2D_sp_sp.cast<doubleC>();

    // (2J+1)/VΩ; weights contain sinβ dβ.
    const double volume_F = weight_F1D_alpha.sum() * weight_F1D_beta.sum() * weight_F1D_gamma.sum();
    assert(std::isfinite(volume_F) && volume_F > 0.0);
    const double normalization_F = (TargetTwoI_I + 1.0) / volume_F;
    result_C4D_cfg1_cfg2_K_K.setZero();

    // dᴶ(β) is reused over α, γ, φ.
    for (int beta_I = 0; beta_I < Nbeta_I; ++beta_I) {
        const auto& RyJ_C2D_K_K = representation_spin.calc_Ry(beta_F1D_beta(beta_I));
        for (int alpha_I = 0; alpha_I < Nalpha_I; ++alpha_I) {
            for (int gamma_I = 0; gamma_I < Ngamma_I; ++gamma_I) {
                prepare_rotation(alpha_I, beta_I, gamma_I);
                // w = (2J+1) wα wβ wγ/VΩ.
                const double weight_F = normalization_F * weight_F1D_alpha(alpha_I) * weight_F1D_beta(beta_I) * weight_F1D_gamma(gamma_I);
                const auto& tmp_C2D_cfg1_cfg2 = calc_one_body_pnp(OneBody_C2D_sp_sp);

                for (int K2_I = 0; K2_I <= TargetTwoI_I; ++K2_I) {
                    const double K2_F = K2_I - 0.5 * TargetTwoI_I;
                    for (int K1_I = 0; K1_I <= TargetTwoI_I; ++K1_I) {
                        // Dᴶ* = exp(iK₁α) dᴶ* exp(iK₂γ).
                        const double K1_F = K1_I - 0.5 * TargetTwoI_I;
                        const doubleC factorOmega_C = std::conj(RyJ_C2D_K_K(K1_I, K2_I)) * std::exp(doubleC(0.0, K1_F * alpha_F1D_alpha(alpha_I) + K2_F * gamma_F1D_gamma(gamma_I)));
                        const doubleC weight_factor_C = weight_F * factorOmega_C;
                        Eigen::Map<Eigen::MatrixXcd> result_C2D_cfg_cfg(result_C4D_cfg1_cfg2_K_K.data() + (static_cast<Eigen::Index>(K2_I) * (TargetTwoI_I + 1) + K1_I) * Ncfg1_I * Ncfg2_I, Ncfg1_I, Ncfg2_I);
                        result_C2D_cfg_cfg += weight_factor_C * tmp_C2D_cfg1_cfg2;
                    }
                }
            }
        }
    }
    return result_C4D_cfg1_cfg2_K_K;
}

inline const Eigen::Tensor<doubleC, 4, Eigen::ColMajor>& HFBProjection::calc_two_body(const Eigen::Tensor<double, 4, Eigen::ColMajor>& TwoBody_F4D_sp_sp_sp_sp) {
    assert(TwoBody_F4D_sp_sp_sp_sp.dimension(0) == Nsp_I && TwoBody_F4D_sp_sp_sp_sp.dimension(1) == Nsp_I && TwoBody_F4D_sp_sp_sp_sp.dimension(2) == Nsp_I && TwoBody_F4D_sp_sp_sp_sp.dimension(3) == Nsp_I);
    assert(Eigen::Map<const Eigen::VectorXd>(TwoBody_F4D_sp_sp_sp_sp.data(), TwoBody_F4D_sp_sp_sp_sp.size()).allFinite());

    // (2J+1)/VΩ; weights contain sinβ dβ.
    const double volume_F = weight_F1D_alpha.sum() * weight_F1D_beta.sum() * weight_F1D_gamma.sum();
    assert(std::isfinite(volume_F) && volume_F > 0.0);
    const double normalization_F = (TargetTwoI_I + 1.0) / volume_F;
    result_C4D_cfg1_cfg2_K_K.setZero();

    // dᴶ(β) is reused over α, γ, φ.
    for (int beta_I = 0; beta_I < Nbeta_I; ++beta_I) {
        const auto& RyJ_C2D_K_K = representation_spin.calc_Ry(beta_F1D_beta(beta_I));
        for (int alpha_I = 0; alpha_I < Nalpha_I; ++alpha_I) {
            for (int gamma_I = 0; gamma_I < Ngamma_I; ++gamma_I) {
                prepare_rotation(alpha_I, beta_I, gamma_I);
                // w = (2J+1) wα wβ wγ/VΩ.
                const double weight_F = normalization_F * weight_F1D_alpha(alpha_I) * weight_F1D_beta(beta_I) * weight_F1D_gamma(gamma_I);
                const auto& tmp_C2D_cfg1_cfg2 = calc_two_body_pnp(TwoBody_F4D_sp_sp_sp_sp);

                for (int K2_I = 0; K2_I <= TargetTwoI_I; ++K2_I) {
                    const double K2_F = K2_I - 0.5 * TargetTwoI_I;
                    for (int K1_I = 0; K1_I <= TargetTwoI_I; ++K1_I) {
                        // Dᴶ* = exp(iK₁α) dᴶ* exp(iK₂γ).
                        const double K1_F = K1_I - 0.5 * TargetTwoI_I;
                        const doubleC factorOmega_C = std::conj(RyJ_C2D_K_K(K1_I, K2_I)) * std::exp(doubleC(0.0, K1_F * alpha_F1D_alpha(alpha_I) + K2_F * gamma_F1D_gamma(gamma_I)));
                        const doubleC weight_factor_C = weight_F * factorOmega_C;
                        Eigen::Map<Eigen::MatrixXcd> result_C2D_cfg_cfg(result_C4D_cfg1_cfg2_K_K.data() + (static_cast<Eigen::Index>(K2_I) * (TargetTwoI_I + 1) + K1_I) * Ncfg1_I * Ncfg2_I, Ncfg1_I, Ncfg2_I);
                        result_C2D_cfg_cfg += weight_factor_C * tmp_C2D_cfg1_cfg2;
                    }
                }
            }
        }
    }
    return result_C4D_cfg1_cfg2_K_K;
}

inline void HFBProjection::calc_multipole(int TwoI1_I, int TwoI2_I, int twoLambda_I, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& Multipole_C3D_sp_sp_mu, Eigen::Tensor<doubleC, 4, Eigen::ColMajor>& multipole_C4D_cfg1_cfg2_K1_K2) {
    assert(TwoI1_I >= 0 && TwoI2_I == TargetTwoI_I && twoLambda_I >= 0 && twoLambda_I % 2 == 0);
    assert(Multipole_C3D_sp_sp_mu.dimension(0) == Nsp_I && Multipole_C3D_sp_sp_mu.dimension(1) == Nsp_I && Multipole_C3D_sp_sp_mu.dimension(2) == twoLambda_I + 1);
    assert(Eigen::Map<const Eigen::VectorXcd>(Multipole_C3D_sp_sp_mu.data(), Multipole_C3D_sp_sp_mu.size()).allFinite());
    assert(multipole_C4D_cfg1_cfg2_K1_K2.dimension(0) == Ncfg1_I && multipole_C4D_cfg1_cfg2_K1_K2.dimension(1) == Ncfg2_I && multipole_C4D_cfg1_cfg2_K1_K2.dimension(2) == TwoI1_I + 1 && multipole_C4D_cfg1_cfg2_K1_K2.dimension(3) == TwoI2_I + 1);
    multipole_C4D_cfg1_cfg2_K1_K2.setZero();
    // |J₁-J₂| ≤ λ ≤ J₁+J₂; J₁-J₂ ∈ ℤ.
    if (std::abs(TwoI1_I - TwoI2_I) > twoLambda_I || TwoI1_I + TwoI2_I < twoLambda_I || (TwoI1_I + TwoI2_I) % 2 != 0) { return; }

    Eigen::MatrixXd CG_F2D_K1_mu{};
    CG_F2D_K1_mu.resize(TwoI1_I + 1, twoLambda_I + 1);
    Eigen::Tensor<doubleC, 3, Eigen::ColMajor> Multipole_C3D_cfg_cfg_mu{};
    Multipole_C3D_cfg_cfg_mu.resize(Ncfg1_I, Ncfg2_I, twoLambda_I + 1);
    // √(2J₁+1) CG = (-1)^(J₂-λ+K₁) (2J₁+1) (J₂ λ J₁; K₁-μ μ -K₁).
    for (int K1_I = 0; K1_I <= TwoI1_I; ++K1_I) {
        const int TwoK1_I = 2 * K1_I - TwoI1_I;
        const int TwoMuMin_I = std::max(-twoLambda_I, TwoK1_I - TwoI2_I);
        const int TwoMuMax_I = std::min(twoLambda_I, TwoK1_I + TwoI2_I);
        const double sign_F = 1.0 - 2.0 * (std::abs((TwoI2_I + TwoK1_I - twoLambda_I) / 2) % 2);
        for (int TwoMu_I = TwoMuMin_I; TwoMu_I <= TwoMuMax_I; TwoMu_I += 2) {
            const int mu_I = (TwoMu_I + twoLambda_I) / 2;
            CG_F2D_K1_mu(K1_I, mu_I) = sign_F * (TwoI1_I + 1.0) * gsl_sf_coupling_3j(TwoI2_I, twoLambda_I, TwoI1_I, TwoK1_I - TwoMu_I, TwoMu_I, -TwoK1_I);
        }
    }

    // (2J+1)/(VΩ Nφ); weights contain sinβ dβ.
    const double volume_F = weight_F1D_alpha.sum() * weight_F1D_beta.sum() * weight_F1D_gamma.sum();
    assert(std::isfinite(volume_F) && volume_F > 0.0);
    const double normalization_F = (TargetTwoI_I + 1.0) / (volume_F * Nphi_I);

    // dᴶ(β) is reused over α, γ, φ.
    for (int beta_I = 0; beta_I < Nbeta_I; ++beta_I) {
        const auto& RyJ_C2D_K_K = representation_spin.calc_Ry(beta_F1D_beta(beta_I));
        for (int alpha_I = 0; alpha_I < Nalpha_I; ++alpha_I) {
            for (int gamma_I = 0; gamma_I < Ngamma_I; ++gamma_I) {
                prepare_rotation(alpha_I, beta_I, gamma_I);
                // w = (2J+1) wα wβ wγ/(VΩ Nφ).
                const double weight_F = normalization_F * weight_F1D_alpha(alpha_I) * weight_F1D_beta(beta_I) * weight_F1D_gamma(gamma_I);
                // Σφ exp(iNφ) O(φ,Ω).
                Multipole_C3D_cfg_cfg_mu.setZero();
                prepare_DUV();
                for (int phi_I = 0; phi_I < Nphi_I; ++phi_I) {
                    // U₂g = exp(-iφ) R U₂; V₂g = exp(iφ) R* V₂.
                    const double phi_F = 2.0 * std::numbers::pi * phi_I / Nphi_I;
                    prepare_contractions(phi_I);
                    const doubleC factorPhi_C = std::exp(doubleC(0.0, TargetN_I * phi_F));
                    // Q_abμ(g) = Σ_ij (Q_λμ)_ij OBTD_ab(i,j;g).
                    for (int sp2_I = 0; sp2_I < Nsp_I; ++sp2_I) {
                        for (int sp1_I = 0; sp1_I < Nsp_I; ++sp1_I) {
                            const auto& OBTD_C2D_cfg_cfg = hfb_pfaffian.calc_obtd(sp1_I, sp2_I);
                            for (int TwoMu_I = -twoLambda_I; TwoMu_I <= twoLambda_I; TwoMu_I += 2) {
                                const int mu_I = (TwoMu_I + twoLambda_I) / 2;
                                Eigen::Map<Eigen::MatrixXcd> Multipole_C2D_cfg_cfg(Multipole_C3D_cfg_cfg_mu.data() + static_cast<Eigen::Index>(mu_I) * Ncfg1_I * Ncfg2_I, Ncfg1_I, Ncfg2_I);
                                Multipole_C2D_cfg_cfg += factorPhi_C * Multipole_C3D_sp_sp_mu(sp1_I, sp2_I, mu_I) * OBTD_C2D_cfg_cfg;
                            }
                        }
                    }
                }

                for (int K2_I = 0; K2_I <= TargetTwoI_I; ++K2_I) {
                    const double K2_F = K2_I - 0.5 * TargetTwoI_I;
                    for (int K1_I = 0; K1_I <= TwoI1_I; ++K1_I) {
                        const int TwoK1_I = 2 * K1_I - TwoI1_I;
                        const int TwoMuMin_I = std::max(-twoLambda_I, TwoK1_I - TwoI2_I);
                        const int TwoMuMax_I = std::min(twoLambda_I, TwoK1_I + TwoI2_I);
                        Eigen::Map<Eigen::MatrixXcd> result_C2D_cfg_cfg(multipole_C4D_cfg1_cfg2_K1_K2.data() + (static_cast<Eigen::Index>(K2_I) * (TwoI1_I + 1) + K1_I) * Ncfg1_I * Ncfg2_I, Ncfg1_I, Ncfg2_I);
                        for (int TwoMu_I = TwoMuMin_I; TwoMu_I <= TwoMuMax_I; TwoMu_I += 2) {
                            const int mu_I = (TwoMu_I + twoLambda_I) / 2;
                            // M = K₁-μ; Dᴶ₂*_(M,K₂) = exp(iMα) dᴶ₂*_(M,K₂) exp(iK₂γ).
                            const int M_I = (TwoK1_I - TwoMu_I + TwoI2_I) / 2;
                            const double M_F = 0.5 * (TwoK1_I - TwoMu_I);
                            const doubleC factorOmega_C = std::conj(RyJ_C2D_K_K(M_I, K2_I)) * std::exp(doubleC(0.0, M_F * alpha_F1D_alpha(alpha_I) + K2_F * gamma_F1D_gamma(gamma_I)));
                            const doubleC weight_factor_C = weight_F * factorOmega_C * CG_F2D_K1_mu(K1_I, mu_I);
                            const Eigen::Map<const Eigen::MatrixXcd> Multipole_C2D_cfg_cfg(Multipole_C3D_cfg_cfg_mu.data() + static_cast<Eigen::Index>(mu_I) * Ncfg1_I * Ncfg2_I, Ncfg1_I, Ncfg2_I);
                            result_C2D_cfg_cfg += weight_factor_C * Multipole_C2D_cfg_cfg;
                        }
                    }
                }
            }
        }
    }
}
