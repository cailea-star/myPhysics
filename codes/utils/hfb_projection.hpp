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

#include <gsl/gsl_sf_coupling.h>

#include <Eigen/Core>
#include <unsupported/Eigen/CXX11/Tensor>

#include "group_so3_su2.hpp"
#include "hfb_pfaffian.hpp"

/**
 * @brief  Project HFB configurations onto particle number and angular momentum.
 * @math   N = TargetN; 2J = TargetTwoJ; K = -J+k, k = 0,…,2J.
 * @note   Rotations and U,V share one ordered single-particle basis.
 * @note   Weights exclude Wigner factors and projector normalization.
 * @note   All calc methods require nonzero vacuum overlaps at quadrature nodes.
 * @note   Returned references alias reusable output tensors.
 */
class HFBProjection {
public:
    int TargetN_I = 0;
    int TargetTwoJ_I = 0;

    int Nsp_I = 0;
    int Ncfg_I = 0;
    int NcqpMax_I = 0;

    int Nphi_I = 0;
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

    Eigen::MatrixXcd U1_C2D_sp_qp1{};
    Eigen::MatrixXcd V1_C2D_sp_qp1{};
    Eigen::MatrixXcd U2_C2D_sp_qp2{};
    Eigen::MatrixXcd V2_C2D_sp_qp2{};

    HFBPfaffian hfb_pfaffian;
    RepresentationSpin representation_spin;

    Eigen::Tensor<doubleC, 4, Eigen::ColMajor> result_C4D_cfg_cfg_K_K{};
    Eigen::Tensor<doubleC, 4, Eigen::ColMajor> multipole_C4D_cfg_cfg_K1_K2{};

private:
    Eigen::MatrixXcd RzRy_C2D_sp_sp{};
    Eigen::MatrixXcd RzRyRz_C2D_sp_sp{};
    Eigen::MatrixXcd U2g_C2D_sp_qp2{};
    Eigen::MatrixXcd V2g_C2D_sp_qp2{};

public:
    /**
     * @brief  Allocate projection and Pfaffian workspaces.
     * @math   φ_k = 2πk/Nphi; Δφ/(2π) = 1/Nphi.
     * @output Stored dimensions, target quantum numbers, and allocated workspaces.
     * @note   Rotations, Euler meshes, weights, and vacua require updates.
     * @note   NcqpMax selects configuration parity and maximum quasiparticle count.
     */
    HFBProjection(int TargetN_I_, int TargetTwoJ_I_, int Nsp_I_, int NcqpMax_I_, int Nphi_I_, int Nalpha_I_, int Nbeta_I_, int Ngamma_I_)
    : hfb_pfaffian(Nsp_I_, NcqpMax_I_, NcqpMax_I_), representation_spin(TargetTwoJ_I_) {
        assert(TargetN_I_ >= 0 && TargetN_I_ <= Nsp_I_ && TargetTwoJ_I_ >= 0);
        assert(Nphi_I_ > 0 && Nalpha_I_ > 0 && Nbeta_I_ > 0 && Ngamma_I_ > 0);

        // (N,J,Nsp,NcqpMax) → fixed projection parameters.
        TargetN_I = TargetN_I_;
        TargetTwoJ_I = TargetTwoJ_I_;
        Nsp_I = Nsp_I_;
        Ncfg_I = static_cast<int>(hfb_pfaffian.config1_I2D_cfg1_cqp1.size());
        NcqpMax_I = NcqpMax_I_;

        // Nφ,Nα,Nβ,Nγ → mesh sizes.
        Nphi_I = Nphi_I_;
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

        // U₁,V₁,U₂,V₂ ∈ ℂ^{Nsp×Nsp}.
        U1_C2D_sp_qp1.resize(Nsp_I, Nsp_I);
        V1_C2D_sp_qp1.resize(Nsp_I, Nsp_I);
        U2_C2D_sp_qp2.resize(Nsp_I, Nsp_I);
        V2_C2D_sp_qp2.resize(Nsp_I, Nsp_I);

        // result ∈ ℂ^{Ncfg×Ncfg×(2J+1)×(2J+1)}.
        result_C4D_cfg_cfg_K_K.resize(Ncfg_I, Ncfg_I, TargetTwoJ_I + 1, TargetTwoJ_I + 1);

        // RzRy,RzRyRz,U₂g,V₂g ∈ ℂ^{Nsp×Nsp}.
        RzRy_C2D_sp_sp.resize(Nsp_I, Nsp_I);
        RzRyRz_C2D_sp_sp.resize(Nsp_I, Nsp_I);
        U2g_C2D_sp_qp2.resize(Nsp_I, Nsp_I);
        V2g_C2D_sp_qp2.resize(Nsp_I, Nsp_I);
    }

    /**
     * @brief  Store left and right Bogoliubov matrices.
     * @math   U₁,V₁,U₂,V₂ ∈ ℝ^{Nsp×Nsp} → complex caches.
     * @output Updated U₁,V₁,U₂,V₂.
     * @note   Requires canonical pairs and invertible U₁,U₂.
     */
    void update_UV(const Eigen::MatrixXd& U1_F2D_sp_qp1_, const Eigen::MatrixXd& V1_F2D_sp_qp1_, const Eigen::MatrixXd& U2_F2D_sp_qp2_, const Eigen::MatrixXd& V2_F2D_sp_qp2_);

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
     * @brief  Integrate overlaps using Euler-angle and uniform gauge quadrature.
     * @math   result_abk₁k₂ = ⟨Φ₁;a|Pᴺ Pᴶ_{K₁K₂}|Φ₂;b⟩.
     * @output Updated result_C4D_cfg_cfg_K_K and its const reference.
     * @note   VΩ = Σwα Σwβ Σwγ for complete integration domains.
     * @note   Integer J permits SO(3); half-integer J requires SU(2).
     */
    const Eigen::Tensor<doubleC, 4, Eigen::ColMajor>& calc_overlap();

    /**
     * @brief  Integrate one-body kernels using Euler-angle and uniform gauge quadrature.
     * @math   H¹_ab(g) = Σ_ij OneBody_ij OBTD_ab(i,j;g).
     * @math   result_abk₁k₂ = ⟨Φ₁;a|H¹ Pᴺ Pᴶ_{K₁K₂}|Φ₂;b⟩.
     * @output Updated result_C4D_cfg_cfg_K_K and its const reference.
     * @note   Requires rotational invariance and particle-number conservation.
     */
    const Eigen::Tensor<doubleC, 4, Eigen::ColMajor>& calc_one_body(const Eigen::MatrixXd& OneBody_F2D_sp_sp);

    /**
     * @brief  Integrate two-body kernels using Euler-angle and uniform gauge quadrature.
     * @math   H²_ab(g) = ½Σ_ijkl TwoBody_ijkl TBTD_ab(i,j,k,l;g).
     * @math   result_abk₁k₂ = ⟨Φ₁;a|H² Pᴺ Pᴶ_{K₁K₂}|Φ₂;b⟩.
     * @output Updated result_C4D_cfg_cfg_K_K and its const reference.
     * @note   TwoBody contains unsymmetrized matrix elements ⟨ij|v|kl⟩.
     * @note   Requires rotational invariance and particle-number conservation.
     */
    const Eigen::Tensor<doubleC, 4, Eigen::ColMajor>& calc_two_body(const Eigen::Tensor<double, 4, Eigen::ColMajor>& TwoBody_F4D_sp_sp_sp_sp);

    /**
     * @brief  Compute reduced multipoles using quadrature and Clebsch–Gordan coupling.
     * @math   multipole_abk₁k₂ = ⟨Φ₁^{J₁;N};K₁,a‖Q_λ‖Φ₂^{J₂;N};K₂,b⟩.
     * @math   dim = Ncfg × Ncfg × (TwoJ1+1) × (TwoJ2+1).
     * @output Updated multipole_C4D_cfg_cfg_K1_K2 and its const reference.
     * @note   twoLambda_I = 2λ; μ = -λ,…,λ.
     * @note   Number-conserving spherical tensor; unnormalized, unmixed configuration states.
     */
    const Eigen::Tensor<doubleC, 4, Eigen::ColMajor>& calc_multipole(int TwoJ1_I, int TwoJ2_I, int twoLambda_I, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& Multipole_C3D_sp_sp_mu);
};

inline void HFBProjection::update_UV(const Eigen::MatrixXd& U1_F2D_sp_qp1_, const Eigen::MatrixXd& V1_F2D_sp_qp1_, const Eigen::MatrixXd& U2_F2D_sp_qp2_, const Eigen::MatrixXd& V2_F2D_sp_qp2_) {
    assert(U1_F2D_sp_qp1_.rows() == Nsp_I && U1_F2D_sp_qp1_.cols() == Nsp_I);
    assert(V1_F2D_sp_qp1_.rows() == Nsp_I && V1_F2D_sp_qp1_.cols() == Nsp_I);
    assert(U2_F2D_sp_qp2_.rows() == Nsp_I && U2_F2D_sp_qp2_.cols() == Nsp_I);
    assert(V2_F2D_sp_qp2_.rows() == Nsp_I && V2_F2D_sp_qp2_.cols() == Nsp_I);
    // (U₁,V₁,U₂,V₂) ∈ ℝ → ℂ.
    U1_C2D_sp_qp1 = U1_F2D_sp_qp1_.cast<doubleC>();
    V1_C2D_sp_qp1 = V1_F2D_sp_qp1_.cast<doubleC>();
    U2_C2D_sp_qp2 = U2_F2D_sp_qp2_.cast<doubleC>();
    V2_C2D_sp_qp2 = V2_F2D_sp_qp2_.cast<doubleC>();
}

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

inline const Eigen::Tensor<doubleC, 4, Eigen::ColMajor>& HFBProjection::calc_overlap() {
    // (2J+1)/(VΩ Nφ); weights contain sinβ dβ.
    const double volume_F = weight_F1D_alpha.sum() * weight_F1D_beta.sum() * weight_F1D_gamma.sum();
    assert(std::isfinite(volume_F) && volume_F > 0.0);
    const double normalization_F = (TargetTwoJ_I + 1.0) / (volume_F * Nphi_I);
    result_C4D_cfg_cfg_K_K.setZero();

    // dᴶ(β) is reused over α, γ, φ.
    for (int beta_I = 0; beta_I < Nbeta_I; ++beta_I) {
        const auto& RyJ_C2D_K_K = representation_spin.calc_Ry(beta_F1D_beta(beta_I));
        const Eigen::Map<const Eigen::MatrixXcd> Ry_C2D_sp_sp(Ry_C3D_sp_sp_beta.data() + static_cast<Eigen::Index>(beta_I) * Nsp_I * Nsp_I, Nsp_I, Nsp_I);
        for (int alpha_I = 0; alpha_I < Nalpha_I; ++alpha_I) {
            const Eigen::Map<const Eigen::MatrixXcd> RzAlpha_C2D_sp_sp(Rz_C3D_sp_sp_alpha.data() + static_cast<Eigen::Index>(alpha_I) * Nsp_I * Nsp_I, Nsp_I, Nsp_I);
            RzRy_C2D_sp_sp.noalias() = RzAlpha_C2D_sp_sp * Ry_C2D_sp_sp;
            for (int gamma_I = 0; gamma_I < Ngamma_I; ++gamma_I) {
                const Eigen::Map<const Eigen::MatrixXcd> RzGamma_C2D_sp_sp(Rz_C3D_sp_sp_gamma.data() + static_cast<Eigen::Index>(gamma_I) * Nsp_I * Nsp_I, Nsp_I, Nsp_I);
                RzRyRz_C2D_sp_sp.noalias() = RzRy_C2D_sp_sp * RzGamma_C2D_sp_sp;
                // w = (2J+1) wα wβ wγ/(VΩ Nφ).
                const double weight_F = normalization_F * weight_F1D_alpha(alpha_I) * weight_F1D_beta(beta_I) * weight_F1D_gamma(gamma_I);
                for (int phi_I = 0; phi_I < Nphi_I; ++phi_I) {
                    // U₂g = exp(-iφ) R U₂; V₂g = exp(iφ) R* V₂.
                    const double phi_F = 2.0 * std::numbers::pi * phi_I / Nphi_I;
                    const doubleC gauge_C = std::exp(doubleC(0.0, -phi_F));
                    U2g_C2D_sp_qp2.noalias() = gauge_C * RzRyRz_C2D_sp_sp * U2_C2D_sp_qp2;
                    V2g_C2D_sp_qp2.noalias() = std::conj(gauge_C) * RzRyRz_C2D_sp_sp.conjugate() * V2_C2D_sp_qp2;
                    hfb_pfaffian.update_contractions(U1_C2D_sp_qp1, V1_C2D_sp_qp1, U2g_C2D_sp_qp2, V2g_C2D_sp_qp2, doubleC(1.0, 0.0));
                    const auto& overlap_C2D_cfg_cfg = hfb_pfaffian.calc_overlap();

                    const doubleC factorPhi_C = std::exp(doubleC(0.0, TargetN_I * phi_F));
                    for (int K2_I = 0; K2_I <= TargetTwoJ_I; ++K2_I) {
                        const double K2_F = K2_I - 0.5 * TargetTwoJ_I;
                        for (int K1_I = 0; K1_I <= TargetTwoJ_I; ++K1_I) {
                            // Dᴶ* = exp(iK₁α) dᴶ* exp(iK₂γ).
                            const double K1_F = K1_I - 0.5 * TargetTwoJ_I;
                            const doubleC factorOmega_C = std::conj(RyJ_C2D_K_K(K1_I, K2_I)) * std::exp(doubleC(0.0, K1_F * alpha_F1D_alpha(alpha_I) + K2_F * gamma_F1D_gamma(gamma_I)));
                            const doubleC weight_factor_C = weight_F * factorPhi_C * factorOmega_C;
                            Eigen::Map<Eigen::MatrixXcd> result_C2D_cfg_cfg(result_C4D_cfg_cfg_K_K.data() + (static_cast<Eigen::Index>(K2_I) * (TargetTwoJ_I + 1) + K1_I) * Ncfg_I * Ncfg_I, Ncfg_I, Ncfg_I);
                            result_C2D_cfg_cfg += weight_factor_C * overlap_C2D_cfg_cfg;
                        }
                    }
                }
            }
        }
    }
    return result_C4D_cfg_cfg_K_K;
}

inline const Eigen::Tensor<doubleC, 4, Eigen::ColMajor>& HFBProjection::calc_one_body(const Eigen::MatrixXd& OneBody_F2D_sp_sp) {
    assert(OneBody_F2D_sp_sp.rows() == Nsp_I && OneBody_F2D_sp_sp.cols() == Nsp_I && OneBody_F2D_sp_sp.allFinite());
    Eigen::MatrixXcd OneBody_C2D_cfg_cfg{};
    OneBody_C2D_cfg_cfg.resize(Ncfg_I, Ncfg_I);

    // (2J+1)/(VΩ Nφ); weights contain sinβ dβ.
    const double volume_F = weight_F1D_alpha.sum() * weight_F1D_beta.sum() * weight_F1D_gamma.sum();
    assert(std::isfinite(volume_F) && volume_F > 0.0);
    const double normalization_F = (TargetTwoJ_I + 1.0) / (volume_F * Nphi_I);
    result_C4D_cfg_cfg_K_K.setZero();

    // dᴶ(β) is reused over α, γ, φ.
    for (int beta_I = 0; beta_I < Nbeta_I; ++beta_I) {
        const auto& RyJ_C2D_K_K = representation_spin.calc_Ry(beta_F1D_beta(beta_I));
        const Eigen::Map<const Eigen::MatrixXcd> Ry_C2D_sp_sp(Ry_C3D_sp_sp_beta.data() + static_cast<Eigen::Index>(beta_I) * Nsp_I * Nsp_I, Nsp_I, Nsp_I);
        for (int alpha_I = 0; alpha_I < Nalpha_I; ++alpha_I) {
            const Eigen::Map<const Eigen::MatrixXcd> RzAlpha_C2D_sp_sp(Rz_C3D_sp_sp_alpha.data() + static_cast<Eigen::Index>(alpha_I) * Nsp_I * Nsp_I, Nsp_I, Nsp_I);
            RzRy_C2D_sp_sp.noalias() = RzAlpha_C2D_sp_sp * Ry_C2D_sp_sp;
            for (int gamma_I = 0; gamma_I < Ngamma_I; ++gamma_I) {
                const Eigen::Map<const Eigen::MatrixXcd> RzGamma_C2D_sp_sp(Rz_C3D_sp_sp_gamma.data() + static_cast<Eigen::Index>(gamma_I) * Nsp_I * Nsp_I, Nsp_I, Nsp_I);
                RzRyRz_C2D_sp_sp.noalias() = RzRy_C2D_sp_sp * RzGamma_C2D_sp_sp;
                // w = (2J+1) wα wβ wγ/(VΩ Nφ).
                const double weight_F = normalization_F * weight_F1D_alpha(alpha_I) * weight_F1D_beta(beta_I) * weight_F1D_gamma(gamma_I);
                for (int phi_I = 0; phi_I < Nphi_I; ++phi_I) {
                    // U₂g = exp(-iφ) R U₂; V₂g = exp(iφ) R* V₂.
                    const double phi_F = 2.0 * std::numbers::pi * phi_I / Nphi_I;
                    const doubleC gauge_C = std::exp(doubleC(0.0, -phi_F));
                    U2g_C2D_sp_qp2.noalias() = gauge_C * RzRyRz_C2D_sp_sp * U2_C2D_sp_qp2;
                    V2g_C2D_sp_qp2.noalias() = std::conj(gauge_C) * RzRyRz_C2D_sp_sp.conjugate() * V2_C2D_sp_qp2;
                    hfb_pfaffian.update_contractions(U1_C2D_sp_qp1, V1_C2D_sp_qp1, U2g_C2D_sp_qp2, V2g_C2D_sp_qp2, doubleC(1.0, 0.0));
                    // H¹_ab(g) = Σ_ij h_ij OBTD_ab(i,j;g).
                    OneBody_C2D_cfg_cfg.setZero();
                    for (int sp2_I = 0; sp2_I < Nsp_I; ++sp2_I) {
                        for (int sp1_I = 0; sp1_I < Nsp_I; ++sp1_I) {
                            OneBody_C2D_cfg_cfg += OneBody_F2D_sp_sp(sp1_I, sp2_I) * hfb_pfaffian.calc_obtd(sp1_I, sp2_I);
                        }
                    }

                    const doubleC factorPhi_C = std::exp(doubleC(0.0, TargetN_I * phi_F));
                    for (int K2_I = 0; K2_I <= TargetTwoJ_I; ++K2_I) {
                        const double K2_F = K2_I - 0.5 * TargetTwoJ_I;
                        for (int K1_I = 0; K1_I <= TargetTwoJ_I; ++K1_I) {
                            // Dᴶ* = exp(iK₁α) dᴶ* exp(iK₂γ).
                            const double K1_F = K1_I - 0.5 * TargetTwoJ_I;
                            const doubleC factorOmega_C = std::conj(RyJ_C2D_K_K(K1_I, K2_I)) * std::exp(doubleC(0.0, K1_F * alpha_F1D_alpha(alpha_I) + K2_F * gamma_F1D_gamma(gamma_I)));
                            const doubleC weight_factor_C = weight_F * factorPhi_C * factorOmega_C;
                            Eigen::Map<Eigen::MatrixXcd> result_C2D_cfg_cfg(result_C4D_cfg_cfg_K_K.data() + (static_cast<Eigen::Index>(K2_I) * (TargetTwoJ_I + 1) + K1_I) * Ncfg_I * Ncfg_I, Ncfg_I, Ncfg_I);
                            result_C2D_cfg_cfg += weight_factor_C * OneBody_C2D_cfg_cfg;
                        }
                    }
                }
            }
        }
    }
    return result_C4D_cfg_cfg_K_K;
}

inline const Eigen::Tensor<doubleC, 4, Eigen::ColMajor>& HFBProjection::calc_two_body(const Eigen::Tensor<double, 4, Eigen::ColMajor>& TwoBody_F4D_sp_sp_sp_sp) {
    assert(TwoBody_F4D_sp_sp_sp_sp.dimension(0) == Nsp_I && TwoBody_F4D_sp_sp_sp_sp.dimension(1) == Nsp_I && TwoBody_F4D_sp_sp_sp_sp.dimension(2) == Nsp_I && TwoBody_F4D_sp_sp_sp_sp.dimension(3) == Nsp_I);
    assert(Eigen::Map<const Eigen::VectorXd>(TwoBody_F4D_sp_sp_sp_sp.data(), TwoBody_F4D_sp_sp_sp_sp.size()).allFinite());
    Eigen::MatrixXcd TwoBody_C2D_cfg_cfg{};
    TwoBody_C2D_cfg_cfg.resize(Ncfg_I, Ncfg_I);

    // (2J+1)/(VΩ Nφ); weights contain sinβ dβ.
    const double volume_F = weight_F1D_alpha.sum() * weight_F1D_beta.sum() * weight_F1D_gamma.sum();
    assert(std::isfinite(volume_F) && volume_F > 0.0);
    const double normalization_F = (TargetTwoJ_I + 1.0) / (volume_F * Nphi_I);
    result_C4D_cfg_cfg_K_K.setZero();

    // dᴶ(β) is reused over α, γ, φ.
    for (int beta_I = 0; beta_I < Nbeta_I; ++beta_I) {
        const auto& RyJ_C2D_K_K = representation_spin.calc_Ry(beta_F1D_beta(beta_I));
        const Eigen::Map<const Eigen::MatrixXcd> Ry_C2D_sp_sp(Ry_C3D_sp_sp_beta.data() + static_cast<Eigen::Index>(beta_I) * Nsp_I * Nsp_I, Nsp_I, Nsp_I);
        for (int alpha_I = 0; alpha_I < Nalpha_I; ++alpha_I) {
            const Eigen::Map<const Eigen::MatrixXcd> RzAlpha_C2D_sp_sp(Rz_C3D_sp_sp_alpha.data() + static_cast<Eigen::Index>(alpha_I) * Nsp_I * Nsp_I, Nsp_I, Nsp_I);
            RzRy_C2D_sp_sp.noalias() = RzAlpha_C2D_sp_sp * Ry_C2D_sp_sp;
            for (int gamma_I = 0; gamma_I < Ngamma_I; ++gamma_I) {
                const Eigen::Map<const Eigen::MatrixXcd> RzGamma_C2D_sp_sp(Rz_C3D_sp_sp_gamma.data() + static_cast<Eigen::Index>(gamma_I) * Nsp_I * Nsp_I, Nsp_I, Nsp_I);
                RzRyRz_C2D_sp_sp.noalias() = RzRy_C2D_sp_sp * RzGamma_C2D_sp_sp;
                // w = (2J+1) wα wβ wγ/(VΩ Nφ).
                const double weight_F = normalization_F * weight_F1D_alpha(alpha_I) * weight_F1D_beta(beta_I) * weight_F1D_gamma(gamma_I);
                for (int phi_I = 0; phi_I < Nphi_I; ++phi_I) {
                    // U₂g = exp(-iφ) R U₂; V₂g = exp(iφ) R* V₂.
                    const double phi_F = 2.0 * std::numbers::pi * phi_I / Nphi_I;
                    const doubleC gauge_C = std::exp(doubleC(0.0, -phi_F));
                    U2g_C2D_sp_qp2.noalias() = gauge_C * RzRyRz_C2D_sp_sp * U2_C2D_sp_qp2;
                    V2g_C2D_sp_qp2.noalias() = std::conj(gauge_C) * RzRyRz_C2D_sp_sp.conjugate() * V2_C2D_sp_qp2;
                    hfb_pfaffian.update_contractions(U1_C2D_sp_qp1, V1_C2D_sp_qp1, U2g_C2D_sp_qp2, V2g_C2D_sp_qp2, doubleC(1.0, 0.0));
                    // H²_ab(g) = ½Σ_ijkl v_ijkl TBTD_ab(i,j,k,l;g).
                    TwoBody_C2D_cfg_cfg.setZero();
                    for (int sp4_I = 0; sp4_I < Nsp_I; ++sp4_I) {
                        for (int sp3_I = 0; sp3_I < Nsp_I; ++sp3_I) {
                            for (int sp2_I = 0; sp2_I < Nsp_I; ++sp2_I) {
                                for (int sp1_I = 0; sp1_I < Nsp_I; ++sp1_I) {
                                    TwoBody_C2D_cfg_cfg += 0.5 * TwoBody_F4D_sp_sp_sp_sp(sp1_I, sp2_I, sp3_I, sp4_I) * hfb_pfaffian.calc_tbtd(sp1_I, sp2_I, sp3_I, sp4_I);
                                }
                            }
                        }
                    }

                    const doubleC factorPhi_C = std::exp(doubleC(0.0, TargetN_I * phi_F));
                    for (int K2_I = 0; K2_I <= TargetTwoJ_I; ++K2_I) {
                        const double K2_F = K2_I - 0.5 * TargetTwoJ_I;
                        for (int K1_I = 0; K1_I <= TargetTwoJ_I; ++K1_I) {
                            // Dᴶ* = exp(iK₁α) dᴶ* exp(iK₂γ).
                            const double K1_F = K1_I - 0.5 * TargetTwoJ_I;
                            const doubleC factorOmega_C = std::conj(RyJ_C2D_K_K(K1_I, K2_I)) * std::exp(doubleC(0.0, K1_F * alpha_F1D_alpha(alpha_I) + K2_F * gamma_F1D_gamma(gamma_I)));
                            const doubleC weight_factor_C = weight_F * factorPhi_C * factorOmega_C;
                            Eigen::Map<Eigen::MatrixXcd> result_C2D_cfg_cfg(result_C4D_cfg_cfg_K_K.data() + (static_cast<Eigen::Index>(K2_I) * (TargetTwoJ_I + 1) + K1_I) * Ncfg_I * Ncfg_I, Ncfg_I, Ncfg_I);
                            result_C2D_cfg_cfg += weight_factor_C * TwoBody_C2D_cfg_cfg;
                        }
                    }
                }
            }
        }
    }
    return result_C4D_cfg_cfg_K_K;
}

inline const Eigen::Tensor<doubleC, 4, Eigen::ColMajor>& HFBProjection::calc_multipole(int TwoJ1_I, int TwoJ2_I, int twoLambda_I, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& Multipole_C3D_sp_sp_mu) {
    assert(TwoJ1_I >= 0 && TwoJ2_I == TargetTwoJ_I && twoLambda_I >= 0 && twoLambda_I % 2 == 0);
    assert(Multipole_C3D_sp_sp_mu.dimension(0) == Nsp_I && Multipole_C3D_sp_sp_mu.dimension(1) == Nsp_I && Multipole_C3D_sp_sp_mu.dimension(2) == twoLambda_I + 1);
    assert(Eigen::Map<const Eigen::VectorXcd>(Multipole_C3D_sp_sp_mu.data(), Multipole_C3D_sp_sp_mu.size()).allFinite());
    multipole_C4D_cfg_cfg_K1_K2.resize(Ncfg_I, Ncfg_I, TwoJ1_I + 1, TwoJ2_I + 1);
    multipole_C4D_cfg_cfg_K1_K2.setZero();
    // |J₁-J₂| ≤ λ ≤ J₁+J₂; J₁-J₂ ∈ ℤ.
    if (std::abs(TwoJ1_I - TwoJ2_I) > twoLambda_I || TwoJ1_I + TwoJ2_I < twoLambda_I || (TwoJ1_I + TwoJ2_I) % 2 != 0) { return multipole_C4D_cfg_cfg_K1_K2; }

    Eigen::MatrixXd CG_F2D_K1_mu{};
    CG_F2D_K1_mu.resize(TwoJ1_I + 1, twoLambda_I + 1);
    Eigen::Tensor<doubleC, 3, Eigen::ColMajor> Multipole_C3D_cfg_cfg_mu{};
    Multipole_C3D_cfg_cfg_mu.resize(Ncfg_I, Ncfg_I, twoLambda_I + 1);
    // √(2J₁+1) CG = (-1)^(J₂-λ+K₁) (2J₁+1) (J₂ λ J₁; K₁-μ μ -K₁).
    for (int K1_I = 0; K1_I <= TwoJ1_I; ++K1_I) {
        const int TwoK1_I = 2 * K1_I - TwoJ1_I;
        const int TwoMuMin_I = std::max(-twoLambda_I, TwoK1_I - TwoJ2_I);
        const int TwoMuMax_I = std::min(twoLambda_I, TwoK1_I + TwoJ2_I);
        const double sign_F = 1.0 - 2.0 * (std::abs((TwoJ2_I + TwoK1_I - twoLambda_I) / 2) % 2);
        for (int TwoMu_I = TwoMuMin_I; TwoMu_I <= TwoMuMax_I; TwoMu_I += 2) {
            const int mu_I = (TwoMu_I + twoLambda_I) / 2;
            CG_F2D_K1_mu(K1_I, mu_I) = sign_F * (TwoJ1_I + 1.0) * gsl_sf_coupling_3j(TwoJ2_I, twoLambda_I, TwoJ1_I, TwoK1_I - TwoMu_I, TwoMu_I, -TwoK1_I);
        }
    }

    // (2J+1)/(VΩ Nφ); weights contain sinβ dβ.
    const double volume_F = weight_F1D_alpha.sum() * weight_F1D_beta.sum() * weight_F1D_gamma.sum();
    assert(std::isfinite(volume_F) && volume_F > 0.0);
    const double normalization_F = (TargetTwoJ_I + 1.0) / (volume_F * Nphi_I);


    // dᴶ(β) is reused over α, γ, φ.
    for (int beta_I = 0; beta_I < Nbeta_I; ++beta_I) {
        const auto& RyJ_C2D_K_K = representation_spin.calc_Ry(beta_F1D_beta(beta_I));
        const Eigen::Map<const Eigen::MatrixXcd> Ry_C2D_sp_sp(Ry_C3D_sp_sp_beta.data() + static_cast<Eigen::Index>(beta_I) * Nsp_I * Nsp_I, Nsp_I, Nsp_I);
        for (int alpha_I = 0; alpha_I < Nalpha_I; ++alpha_I) {
            const Eigen::Map<const Eigen::MatrixXcd> RzAlpha_C2D_sp_sp(Rz_C3D_sp_sp_alpha.data() + static_cast<Eigen::Index>(alpha_I) * Nsp_I * Nsp_I, Nsp_I, Nsp_I);
            RzRy_C2D_sp_sp.noalias() = RzAlpha_C2D_sp_sp * Ry_C2D_sp_sp;
            for (int gamma_I = 0; gamma_I < Ngamma_I; ++gamma_I) {
                const Eigen::Map<const Eigen::MatrixXcd> RzGamma_C2D_sp_sp(Rz_C3D_sp_sp_gamma.data() + static_cast<Eigen::Index>(gamma_I) * Nsp_I * Nsp_I, Nsp_I, Nsp_I);
                RzRyRz_C2D_sp_sp.noalias() = RzRy_C2D_sp_sp * RzGamma_C2D_sp_sp;
                // w = (2J+1) wα wβ wγ/(VΩ Nφ).
                const double weight_F = normalization_F * weight_F1D_alpha(alpha_I) * weight_F1D_beta(beta_I) * weight_F1D_gamma(gamma_I);
                for (int phi_I = 0; phi_I < Nphi_I; ++phi_I) {
                    // U₂g = exp(-iφ) R U₂; V₂g = exp(iφ) R* V₂.
                    const double phi_F = 2.0 * std::numbers::pi * phi_I / Nphi_I;
                    const doubleC gauge_C = std::exp(doubleC(0.0, -phi_F));
                    U2g_C2D_sp_qp2.noalias() = gauge_C * RzRyRz_C2D_sp_sp * U2_C2D_sp_qp2;
                    V2g_C2D_sp_qp2.noalias() = std::conj(gauge_C) * RzRyRz_C2D_sp_sp.conjugate() * V2_C2D_sp_qp2;
                    hfb_pfaffian.update_contractions(U1_C2D_sp_qp1, V1_C2D_sp_qp1, U2g_C2D_sp_qp2, V2g_C2D_sp_qp2, doubleC(1.0, 0.0));
                    // Q_abμ(g) = Σ_ij (Q_λμ)_ij OBTD_ab(i,j;g).
                    Multipole_C3D_cfg_cfg_mu.setZero();
                    for (int sp2_I = 0; sp2_I < Nsp_I; ++sp2_I) {
                        for (int sp1_I = 0; sp1_I < Nsp_I; ++sp1_I) {
                            const auto& OBTD_C2D_cfg_cfg = hfb_pfaffian.calc_obtd(sp1_I, sp2_I);
                            for (int TwoMu_I = -twoLambda_I; TwoMu_I <= twoLambda_I; TwoMu_I += 2) {
                                const int mu_I = (TwoMu_I + twoLambda_I) / 2;
                                Eigen::Map<Eigen::MatrixXcd> Multipole_C2D_cfg_cfg(Multipole_C3D_cfg_cfg_mu.data() + static_cast<Eigen::Index>(mu_I) * Ncfg_I * Ncfg_I, Ncfg_I, Ncfg_I);
                                Multipole_C2D_cfg_cfg += Multipole_C3D_sp_sp_mu(sp1_I, sp2_I, mu_I) * OBTD_C2D_cfg_cfg;
                            }
                        }
                    }

                    const doubleC factorPhi_C = std::exp(doubleC(0.0, TargetN_I * phi_F));
                    for (int K2_I = 0; K2_I <= TargetTwoJ_I; ++K2_I) {
                        const double K2_F = K2_I - 0.5 * TargetTwoJ_I;
                        for (int K1_I = 0; K1_I <= TwoJ1_I; ++K1_I) {
                            const int TwoK1_I = 2 * K1_I - TwoJ1_I;
                            const int TwoMuMin_I = std::max(-twoLambda_I, TwoK1_I - TwoJ2_I);
                            const int TwoMuMax_I = std::min(twoLambda_I, TwoK1_I + TwoJ2_I);
                            Eigen::Map<Eigen::MatrixXcd> result_C2D_cfg_cfg(multipole_C4D_cfg_cfg_K1_K2.data() + (static_cast<Eigen::Index>(K2_I) * (TwoJ1_I + 1) + K1_I) * Ncfg_I * Ncfg_I, Ncfg_I, Ncfg_I);
                            for (int TwoMu_I = TwoMuMin_I; TwoMu_I <= TwoMuMax_I; TwoMu_I += 2) {
                                const int mu_I = (TwoMu_I + twoLambda_I) / 2;
                                // M = K₁-μ; Dᴶ₂*_(M,K₂) = exp(iMα) dᴶ₂*_(M,K₂) exp(iK₂γ).
                                const int M_I = (TwoK1_I - TwoMu_I + TwoJ2_I) / 2;
                                const double M_F = 0.5 * (TwoK1_I - TwoMu_I);
                                const doubleC factorOmega_C = std::conj(RyJ_C2D_K_K(M_I, K2_I)) * std::exp(doubleC(0.0, M_F * alpha_F1D_alpha(alpha_I) + K2_F * gamma_F1D_gamma(gamma_I)));
                                const doubleC weight_factor_C = weight_F * factorPhi_C * factorOmega_C * CG_F2D_K1_mu(K1_I, mu_I);
                                const Eigen::Map<const Eigen::MatrixXcd> Multipole_C2D_cfg_cfg(Multipole_C3D_cfg_cfg_mu.data() + static_cast<Eigen::Index>(mu_I) * Ncfg_I * Ncfg_I, Ncfg_I, Ncfg_I);
                                result_C2D_cfg_cfg += weight_factor_C * Multipole_C2D_cfg_cfg;
                            }
                        }
                    }
                }
            }
        }
    }
    return multipole_C4D_cfg_cfg_K1_K2;
}
