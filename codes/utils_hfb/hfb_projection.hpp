/**
 * @file    hfb_projection.hpp
 * @author  cailea
 * @date    2026-09-11
 * @brief   Angular-momentum projection of HFB configurations.
 */

#pragma once

#include <cassert>
#include <cmath>
#include <vector>

#include <Eigen/Core>
#include <unsupported/Eigen/CXX11/Tensor>

#include "group_so3_su2.hpp"
#include "hfb_projection_pnp.hpp"

/**
 * @brief  Project HFB configurations onto particle number and angular momentum.
 * @math   N = TargetN; 2I = TargetTwoI; K = -I+k, k = 0,…,2I.
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

        // result ∈ ℂ^{Ncfg1×Ncfg2×(2I+1)×(2I+1)}.
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
     * @brief  Integrate overlaps using Euler-angle and uniform gauge quadrature.
     * @math   result_abk₁k₂ = ⟨Φ₁;a|Pᴺ Pᴵ_{K₁K₂}|Φ₂;b⟩.
     * @output Updated result_C4D_cfg1_cfg2_K_K and its const reference.
     * @note   VΩ = Σwα Σwβ Σwγ for complete integration domains.
     * @note   Integer I permits SO(3); half-integer I requires SU(2).
     */
    const Eigen::Tensor<doubleC, 4, Eigen::ColMajor>& calc_overlap();

    /**
     * @brief  Integrate one-body kernels using Euler-angle and uniform gauge quadrature.
     * @math   H¹_ab(g) = Σ_ij OneBody_ij OBTD_ab(i,j;g).
     * @math   result_abk₁k₂ = ⟨Φ₁;a|H¹ Pᴺ Pᴵ_{K₁K₂}|Φ₂;b⟩.
     * @output Updated result_C4D_cfg1_cfg2_K_K and its const reference.
     * @note   Requires rotational invariance and particle-number conservation.
     */
    const Eigen::Tensor<doubleC, 4, Eigen::ColMajor>& calc_one_body(const Eigen::MatrixXd& OneBody_F2D_sp_sp);

    /**
     * @brief  Integrate two-body kernels using Euler-angle and uniform gauge quadrature.
     * @math   H²_ab(g) = ½Σ_ijkl TwoBody_ijkl TBTD_ab(i,j,k,l;g).
     * @math   result_abk₁k₂ = ⟨Φ₁;a|H² Pᴺ Pᴵ_{K₁K₂}|Φ₂;b⟩.
     * @output Updated result_C4D_cfg1_cfg2_K_K and its const reference.
     * @note   TwoBody contains unsymmetrized matrix elements ⟨ij|v|kl⟩.
     * @note   Requires rotational invariance and particle-number conservation.
     */
    const Eigen::Tensor<doubleC, 4, Eigen::ColMajor>& calc_two_body(const Eigen::Tensor<double, 4, Eigen::ColMajor>& TwoBody_F4D_sp_sp_sp_sp);

    /**
     * @brief  Integrate separable kernels using Pfaffians and quadrature.
     * @math   O = Σμ[Q̂λμ†Q̂λμ − Σαδ(Qλμ†Qλμ)αδ cα†cδ].
     * @math   Q ∈ ℝ^{Nsp × Nsp × (2λ+1)}; twoLambda_I = 2λ ≥ 0, λ ∈ ℤ; μ = −λ,…,λ.
     * @output Updated result_C4D_cfg1_cfg2_K_K and its const reference.
     * @note   Complete spherical tensor; excludes coupling strength and extra ½.
     */
    const Eigen::Tensor<doubleC, 4, Eigen::ColMajor>& calc_two_body_Q(int twoLambda_I, const Eigen::Tensor<double, 3, Eigen::ColMajor>& Q_F3D_sp_sp_mu);

    /**
     * @brief  Integrate pairing kernels using Pfaffians and quadrature.
     * @math   P̂λμ† = ½Σab Pab c_a†c_b†; O = Σμ P̂λμ†P̂λμ.
     * @math   P ∈ ℝ^{Nsp × Nsp × (2λ+1)}; Pμᵀ = −Pμ; twoLambda_I = 2λ ≥ 0, λ ∈ ℤ; μ = −λ,…,λ.
     * @output Updated result_C4D_cfg1_cfg2_K_K and its const reference.
     * @note   Complete spherical tensor; excludes coupling strength and attraction sign.
     */
    const Eigen::Tensor<doubleC, 4, Eigen::ColMajor>& calc_two_body_P(int twoLambda_I, const Eigen::Tensor<double, 3, Eigen::ColMajor>& P_F3D_sp_sp_mu);

private:
    /**
     * @brief  Integrate scalar kernels using Euler quadrature.
     * @math   result = (2I+1)/VΩ ΣΩ wΩ Dᴵ*(Ω) K(Ω).
     * @output Updated result_C4D_cfg1_cfg2_K_K and its const reference.
     * @note   Callback returns the PNP kernel for current D.
     * @note   Callback preserves angular data and nuclear target.
     */
    template <typename KernelFunc>
    const Eigen::Tensor<doubleC, 4, Eigen::ColMajor>& integrate_amp(const KernelFunc& kernel_Func);
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

template <typename KernelFunc>
inline const Eigen::Tensor<doubleC, 4, Eigen::ColMajor>& HFBProjection::integrate_amp(const KernelFunc& kernel_Func) {
    // (2I+1)/VΩ; weights contain sinβ dβ.
    const double volume_F = weight_F1D_alpha.sum() * weight_F1D_beta.sum() * weight_F1D_gamma.sum();
    assert(std::isfinite(volume_F) && volume_F > 0.0);
    const double normalization_F = (TargetTwoI_I + 1.0) / volume_F;
    result_C4D_cfg1_cfg2_K_K.setZero();

    // dᴵ(β) is reused over α, γ, φ.
    for (int beta_I = 0; beta_I < Nbeta_I; ++beta_I) {
        const auto& RyI_C2D_K_K = representation_spin.calc_Ry(beta_F1D_beta(beta_I));
        for (int alpha_I = 0; alpha_I < Nalpha_I; ++alpha_I) {
            for (int gamma_I = 0; gamma_I < Ngamma_I; ++gamma_I) {
                // D = Rz(α) Ry(β) Rz(γ).
                const Eigen::Map<const Eigen::MatrixXcd> RzAlpha_C2D_sp_sp(Rz_C3D_sp_sp_alpha.data() + static_cast<Eigen::Index>(alpha_I) * Nsp_I * Nsp_I, Nsp_I, Nsp_I);
                const Eigen::Map<const Eigen::MatrixXcd> Ry_C2D_sp_sp(Ry_C3D_sp_sp_beta.data() + static_cast<Eigen::Index>(beta_I) * Nsp_I * Nsp_I, Nsp_I, Nsp_I);
                const Eigen::Map<const Eigen::MatrixXcd> RzGamma_C2D_sp_sp(Rz_C3D_sp_sp_gamma.data() + static_cast<Eigen::Index>(gamma_I) * Nsp_I * Nsp_I, Nsp_I, Nsp_I);
                D_C2D_sp_sp.noalias() = RzAlpha_C2D_sp_sp * Ry_C2D_sp_sp * RzGamma_C2D_sp_sp;
                // w = (2I+1) wα wβ wγ/VΩ.
                const double weight_F = normalization_F * weight_F1D_alpha(alpha_I) * weight_F1D_beta(beta_I) * weight_F1D_gamma(gamma_I);
                const auto& tmp_C2D_cfg1_cfg2 = kernel_Func();

                for (int K2_I = 0; K2_I <= TargetTwoI_I; ++K2_I) {
                    const double K2_F = K2_I - 0.5 * TargetTwoI_I;
                    for (int K1_I = 0; K1_I <= TargetTwoI_I; ++K1_I) {
                        // Dᴵ* = exp(iK₁α) dᴵ* exp(iK₂γ).
                        const double K1_F = K1_I - 0.5 * TargetTwoI_I;
                        const doubleC factorOmega_C = std::conj(RyI_C2D_K_K(K1_I, K2_I)) * std::exp(doubleC(0.0, K1_F * alpha_F1D_alpha(alpha_I) + K2_F * gamma_F1D_gamma(gamma_I)));
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

inline const Eigen::Tensor<doubleC, 4, Eigen::ColMajor>& HFBProjection::calc_overlap() {
    return integrate_amp([&]() -> const Eigen::MatrixXcd& {
        return calc_overlap_pnp();
    });
}

inline const Eigen::Tensor<doubleC, 4, Eigen::ColMajor>& HFBProjection::calc_one_body(const Eigen::MatrixXd& OneBody_F2D_sp_sp) {
    assert(OneBody_F2D_sp_sp.rows() == Nsp_I && OneBody_F2D_sp_sp.cols() == Nsp_I && OneBody_F2D_sp_sp.allFinite());

    const Eigen::MatrixXcd OneBody_C2D_sp_sp = OneBody_F2D_sp_sp.cast<doubleC>();
    return integrate_amp([&]() -> const Eigen::MatrixXcd& {
        return calc_one_body_pnp(OneBody_C2D_sp_sp);
    });
}

inline const Eigen::Tensor<doubleC, 4, Eigen::ColMajor>& HFBProjection::calc_two_body(const Eigen::Tensor<double, 4, Eigen::ColMajor>& TwoBody_F4D_sp_sp_sp_sp) {
    assert(TwoBody_F4D_sp_sp_sp_sp.dimension(0) == Nsp_I && TwoBody_F4D_sp_sp_sp_sp.dimension(1) == Nsp_I && TwoBody_F4D_sp_sp_sp_sp.dimension(2) == Nsp_I && TwoBody_F4D_sp_sp_sp_sp.dimension(3) == Nsp_I);
    assert(Eigen::Map<const Eigen::VectorXd>(TwoBody_F4D_sp_sp_sp_sp.data(), TwoBody_F4D_sp_sp_sp_sp.size()).allFinite());

    return integrate_amp([&]() -> const Eigen::MatrixXcd& {
        return calc_two_body_pnp(TwoBody_F4D_sp_sp_sp_sp);
    });
}

inline const Eigen::Tensor<doubleC, 4, Eigen::ColMajor>& HFBProjection::calc_two_body_Q(int twoLambda_I, const Eigen::Tensor<double, 3, Eigen::ColMajor>& Q_F3D_sp_sp_mu) {
    return integrate_amp([&]() -> const Eigen::MatrixXcd& {
        return calc_two_body_Q_pnp(twoLambda_I, Q_F3D_sp_sp_mu);
    });
}

inline const Eigen::Tensor<doubleC, 4, Eigen::ColMajor>& HFBProjection::calc_two_body_P(int twoLambda_I, const Eigen::Tensor<double, 3, Eigen::ColMajor>& P_F3D_sp_sp_mu) {
    return integrate_amp([&]() -> const Eigen::MatrixXcd& {
        return calc_two_body_P_pnp(twoLambda_I, P_F3D_sp_sp_mu);
    });
}
