/**
 * @file    hfb_projection.hpp
 * @author  cailea
 * @date    2026-09-11
 * @brief   Angular-momentum projection of HFB configurations.
 */

#pragma once

#include <cassert>

#include <Eigen/Core>
#include <unsupported/Eigen/CXX11/Tensor>

#include "hfb_pfaffian.hpp"

/**
 * @brief  Project HFB configurations onto particle number and angular momentum.
 * @math   N = TargetN; 2J = TargetTwoJ; K = -J+k, k = 0,…,2J.
 * @note   Rotations and U,V share one ordered single-particle basis.
 * @note   Weights exclude Wigner factors and projector normalization.
 * @note   All calc methods require nonzero vacuum overlaps at quadrature nodes.
 * @note   All calc methods overwrite the shared result tensor.
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

    Eigen::MatrixXcd U0_C2D_sp_qp{};
    Eigen::MatrixXcd V0_C2D_sp_qp{};

    HFBPfaffian hfb_pfaffian;

    Eigen::Tensor<doubleC, 4, Eigen::ColMajor> result_C4D_cfg_cfg_K_K{};

private:
    Eigen::MatrixXcd RzRy_C2D_sp_sp{};
    Eigen::MatrixXcd RzRyRz_C2D_sp_sp{};
    Eigen::MatrixXcd Ug_C2D_sp_qp{};
    Eigen::MatrixXcd Vg_C2D_sp_qp{};

public:
    /**
     * @brief  Allocate projection and Pfaffian workspaces.
     * @math   φ_k = 2πk/Nphi; Δφ/(2π) = 1/Nphi.
     * @output Stored dimensions, target quantum numbers, and allocated workspaces.
     * @note   Rotations, Euler meshes, weights, and vacuum require updates.
     * @note   NcqpMax selects configuration parity and maximum quasiparticle count.
     */
    HFBProjection(int TargetN_I_, int TargetTwoJ_I_, int Nsp_I_, int NcqpMax_I_, int Nphi_I_, int Nalpha_I_, int Nbeta_I_, int Ngamma_I_)
    : hfb_pfaffian(Nsp_I_, NcqpMax_I_, NcqpMax_I_) {
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

        // U₀,V₀ ∈ ℂ^{Nsp×Nsp}.
        U0_C2D_sp_qp.resize(Nsp_I, Nsp_I);
        V0_C2D_sp_qp.resize(Nsp_I, Nsp_I);

        // result ∈ ℂ^{Ncfg×Ncfg×(2J+1)×(2J+1)}.
        result_C4D_cfg_cfg_K_K.resize(Ncfg_I, Ncfg_I, TargetTwoJ_I + 1, TargetTwoJ_I + 1);

        // RzRy,RzRyRz,Ug,Vg ∈ ℂ^{Nsp×Nsp}.
        RzRy_C2D_sp_sp.resize(Nsp_I, Nsp_I);
        RzRyRz_C2D_sp_sp.resize(Nsp_I, Nsp_I);
        Ug_C2D_sp_qp.resize(Nsp_I, Nsp_I);
        Vg_C2D_sp_qp.resize(Nsp_I, Nsp_I);
    }

    /**
     * @brief  Store the reference Bogoliubov matrices.
     * @math   U0,V0 ∈ ℂ^{Nsp×Nsp}.
     * @output Updated U0 and V0.
     * @note   Requires canonical U0,V0 and invertible U0.
     */
    void update_UV(const Eigen::MatrixXcd& U0_C2D_sp_qp_, const Eigen::MatrixXcd& V0_C2D_sp_qp_);

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
     * @math   result_abk₁k₂ = ⟨Φ_a|Pᴺ Pᴶ_{K₁K₂}|Φ_b⟩.
     * @output Updated result_C4D_cfg_cfg_K_K and its const reference.
     * @note   VΩ = Σwα Σwβ Σwγ for complete integration domains.
     * @note   Integer J permits SO(3); half-integer J requires SU(2).
     */
    const Eigen::Tensor<doubleC, 4, Eigen::ColMajor>& calc_overlap();

    /**
     * @brief  Integrate one-body kernels using Euler-angle and uniform gauge quadrature.
     * @math   H¹_ab(g) = Σ_ij OneBody_ij OBTD_ab(i,j;g).
     * @math   result_abk₁k₂ = ⟨Φ_a|H¹ Pᴺ Pᴶ_{K₁K₂}|Φ_b⟩.
     * @output Updated result_C4D_cfg_cfg_K_K and its const reference.
     * @note   Requires rotational invariance and particle-number conservation.
     */
    const Eigen::Tensor<doubleC, 4, Eigen::ColMajor>& calc_one_body(const Eigen::MatrixXd& OneBody_F2D_sp_sp);

    /**
     * @brief  Integrate two-body kernels using Euler-angle and uniform gauge quadrature.
     * @math   H²_ab(g) = ½Σ_ijkl TwoBody_ijkl TBTD_ab(i,j,k,l;g).
     * @math   result_abk₁k₂ = ⟨Φ_a|H² Pᴺ Pᴶ_{K₁K₂}|Φ_b⟩.
     * @output Updated result_C4D_cfg_cfg_K_K and its const reference.
     * @note   TwoBody contains unsymmetrized matrix elements ⟨ij|v|kl⟩.
     * @note   Requires rotational invariance and particle-number conservation.
     */
    const Eigen::Tensor<doubleC, 4, Eigen::ColMajor>& calc_two_body(const Eigen::Tensor<double, 4, Eigen::ColMajor>& TwoBody_F4D_sp_sp_sp_sp);
};

inline void HFBProjection::update_UV(const Eigen::MatrixXcd& U0_C2D_sp_qp_, const Eigen::MatrixXcd& V0_C2D_sp_qp_) {
    assert(U0_C2D_sp_qp_.rows() == Nsp_I && U0_C2D_sp_qp_.cols() == Nsp_I);
    assert(V0_C2D_sp_qp_.rows() == Nsp_I && V0_C2D_sp_qp_.cols() == Nsp_I);
    U0_C2D_sp_qp = U0_C2D_sp_qp_;
    V0_C2D_sp_qp = V0_C2D_sp_qp_;
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
