/**
 * @file    hfb_kramers.hpp
 * @author  cailea
 * @date    2026-09-10
 * @brief   Time-reversal-invariant axial HFB calculations.
 */

#pragma once

#include <cassert>
#include <cmath>
#include <functional>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Eigenvalues>

class HFBKramersBlock {
public:
    int Nbsp_I = 0;

    // T|α⟩ = ηα|ᾱ⟩; ηα = (-1)^(j - mj)
    Eigen::VectorXd eta_F1D_bsp{};

    Eigen::MatrixXcd h0PosPos_C2D_bsp_bsp{};
    Eigen::MatrixXcd GammaPosPos_C2D_bsp_bsp{};
    Eigen::MatrixXcd DeltaPosNeg_C2D_bsp_bsp{};

    Eigen::VectorXd Eqp_F1D_bqp{};
    Eigen::VectorXd f_F1D_bqp{};

    // X⁺ = [U⁺; V⁻].
    Eigen::MatrixXcd UPos_C2D_bsp_bqp{};
    Eigen::MatrixXcd VNeg_C2D_bsp_bqp{};

    Eigen::MatrixXcd rhoPosPos_C2D_bsp_bsp{};
    Eigen::MatrixXcd kappaPosNeg_C2D_bsp_bsp{};

    Eigen::MatrixXcd HPos_C2D_2bsp_2bsp{};
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> HPos_eigensolver{};

public:
    /**
     * @brief  Allocate one Kramers-representative block.
     * @math   Nbsp = size(η); ηα² = 1.
     * @output Stored phases and allocated workspace.
     */
    HFBKramersBlock(const Eigen::VectorXd& eta_F1D_bsp_) {
        assert(eta_F1D_bsp_.size() > 0);
        assert((eta_F1D_bsp_.array().abs() == 1.0).all());
        Nbsp_I = static_cast<int>(eta_F1D_bsp_.size());
        eta_F1D_bsp = eta_F1D_bsp_;

        // h₀, Γ, Δ, ρ, κ ∈ ℂ^{Nbsp×Nbsp}.
        h0PosPos_C2D_bsp_bsp.resize(Nbsp_I, Nbsp_I);
        h0PosPos_C2D_bsp_bsp.setZero();
        GammaPosPos_C2D_bsp_bsp.resize(Nbsp_I, Nbsp_I);
        GammaPosPos_C2D_bsp_bsp.setZero();
        DeltaPosNeg_C2D_bsp_bsp.resize(Nbsp_I, Nbsp_I);
        DeltaPosNeg_C2D_bsp_bsp.setZero();
        rhoPosPos_C2D_bsp_bsp.resize(Nbsp_I, Nbsp_I);
        rhoPosPos_C2D_bsp_bsp.setZero();
        kappaPosNeg_C2D_bsp_bsp.resize(Nbsp_I, Nbsp_I);
        kappaPosNeg_C2D_bsp_bsp.setZero();

        // E, f ∈ ℝ^{Nbsp}; U⁺, V⁻ ∈ ℂ^{Nbsp×Nbsp}.
        Eqp_F1D_bqp.resize(Nbsp_I);
        f_F1D_bqp.resize(Nbsp_I);
        UPos_C2D_bsp_bqp.resize(Nbsp_I, Nbsp_I);
        VNeg_C2D_bsp_bqp.resize(Nbsp_I, Nbsp_I);
        Eqp_F1D_bqp.setZero();
        f_F1D_bqp.setZero();
        UPos_C2D_bsp_bqp.setZero();
        VNeg_C2D_bsp_bqp.setZero();

        // H⁺ ∈ ℂ^{2Nbsp×2Nbsp}.
        HPos_C2D_2bsp_2bsp.resize(2 * Nbsp_I, 2 * Nbsp_I);
        HPos_C2D_2bsp_2bsp.setZero();
        HPos_eigensolver = Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd>(2 * Nbsp_I);
    }

    /**
     * @brief  Solve thermal HFB using self-adjoint eigendecomposition.
     * @math   D = diag(η); h = h₀⁺⁺ + Γ⁺⁺.
     * @math   H⁺ = [h-λI, Δ⁺⁻; (Δ⁺⁻)†, -DhD+λI].
     * @math   U⁻ = D(U⁺)*; V⁺ = -D(V⁻)*.
     * @output Updated Eqp, UPos, VNeg, f, and densities.
     * @note   Time-reversal-invariant fields; no zero modes; T ≥ 0.
     */
    void update_UV_E_rho_kappa(double lambda_F, double temperature_F);
};

class HFBKramers {
public:
    using OneBodyFunc = std::function<void(int block_I, Eigen::MatrixXcd& h0PosPos_C2D_bsp_bsp)>;

    using TwoBodyFunc = std::function<void(int block_I, const std::vector<HFBKramersBlock>& blocks, Eigen::MatrixXcd& GammaPosPos_C2D_bsp_bsp, Eigen::MatrixXcd& DeltaPosNeg_C2D_bsp_bsp)>;

    int Nblock_I = 0;

    std::vector<HFBKramersBlock> blocks{};

    OneBodyFunc build_onebody{};
    TwoBodyFunc build_twobody{};

public:
    /**
     * @brief  Allocate blocks and build one-body matrices.
     * @math   {η_b} → {Nbsp_b, h₀,b⁺⁺}.
     * @output Initialized blocks and callbacks.
     */
    HFBKramers(const std::vector<Eigen::VectorXd>& eta_F2D_block_bsp_, const OneBodyFunc& build_onebody_, const TwoBodyFunc& build_twobody_) {
        assert(build_onebody_);
        assert(build_twobody_);
        Nblock_I = static_cast<int>(eta_F2D_block_bsp_.size());
        build_onebody = build_onebody_;
        build_twobody = build_twobody_;

        // {η_b} → {block_b}.
        blocks.reserve(Nblock_I);
        for (int block_I = 0; block_I < Nblock_I; ++block_I) {
            blocks.emplace_back(eta_F2D_block_bsp_[block_I]);
        }

        // block_b → h₀,b⁺⁺.
        for (int block_I = 0; block_I < Nblock_I; ++block_I) {
            build_onebody(block_I, blocks[block_I].h0PosPos_C2D_bsp_bsp);
        }
    }

    /**
     * @brief  Update fields using all representative densities.
     * @math   {ρ_b⁺⁺,κ_b⁺⁻,η_b} → {Γ_b⁺⁺,Δ_b⁺⁻}.
     * @output Overwritten GammaPosPos and DeltaPosNeg.
     */
    void update_Gamma_Delta();

    /**
     * @brief  Solve all blocks using self-adjoint eigendecomposition.
     * @math   {H_b⁺(λ),T} → {E_b,U_b⁺,V_b⁻,f_b,ρ_b⁺⁺,κ_b⁺⁻}.
     * @output Updated representative solutions and densities.
     */
    void update_UV_E_rho_kappa(double lambda_F, double temperature_F);
};

/**
 * @brief  Solve thermal HFB using self-adjoint eigendecomposition.
 * @math   H⁺X⁺ = X⁺E; U⁻ = D(U⁺)*; V⁺ = -D(V⁻)*.
 * @output Updated Eqp, UPos, VNeg, f, and densities.
 */
inline void HFBKramersBlock::update_UV_E_rho_kappa(double lambda_F, double temperature_F) {
    assert(temperature_F >= 0.0);
    assert(h0PosPos_C2D_bsp_bsp.isApprox(h0PosPos_C2D_bsp_bsp.adjoint(), 1.0e-12));
    assert(GammaPosPos_C2D_bsp_bsp.isApprox(GammaPosPos_C2D_bsp_bsp.adjoint(), 1.0e-12));
    // Δ = DΔ†D; D = diag(η).
    assert(DeltaPosNeg_C2D_bsp_bsp.isApprox(eta_F1D_bsp.asDiagonal() * DeltaPosNeg_C2D_bsp_bsp.adjoint() * eta_F1D_bsp.asDiagonal(), 1.0e-12));

    // H⁺₁₁ = h₀⁺⁺ + Γ⁺⁺ - λI.
    HPos_C2D_2bsp_2bsp.topLeftCorner(Nbsp_I, Nbsp_I) = h0PosPos_C2D_bsp_bsp + GammaPosPos_C2D_bsp_bsp;
    HPos_C2D_2bsp_2bsp.topLeftCorner(Nbsp_I, Nbsp_I).diagonal().array() -= lambda_F;

    // H⁺₁₂ = Δ; H⁺₂₁ = Δ†; H⁺₂₂ = -DH⁺₁₁D.
    HPos_C2D_2bsp_2bsp.topRightCorner(Nbsp_I, Nbsp_I) = DeltaPosNeg_C2D_bsp_bsp;
    HPos_C2D_2bsp_2bsp.bottomLeftCorner(Nbsp_I, Nbsp_I) = DeltaPosNeg_C2D_bsp_bsp.adjoint();
    HPos_C2D_2bsp_2bsp.bottomRightCorner(Nbsp_I, Nbsp_I).noalias() = -(eta_F1D_bsp.asDiagonal() * HPos_C2D_2bsp_2bsp.topLeftCorner(Nbsp_I, Nbsp_I) * eta_F1D_bsp.asDiagonal());

    // H⁺ → (E, X⁺).
    HPos_eigensolver.compute(HPos_C2D_2bsp_2bsp);
    assert(HPos_eigensolver.info() == Eigen::Success);
    const Eigen::VectorXd& eigenvalues_F1D_state = HPos_eigensolver.eigenvalues();
    const Eigen::MatrixXcd& eigenvectors_C2D_2bsp_state = HPos_eigensolver.eigenvectors();

    // E₁ ≤ ⋯ ≤ E₂Nbsp; retain Nbsp positive-energy modes.
    assert(eigenvalues_F1D_state(Nbsp_I - 1) < 0.0);
    assert(eigenvalues_F1D_state(Nbsp_I) > 0.0);
    Eqp_F1D_bqp = eigenvalues_F1D_state.tail(Nbsp_I);
    UPos_C2D_bsp_bqp = eigenvectors_C2D_2bsp_state.topRightCorner(Nbsp_I, Nbsp_I);
    VNeg_C2D_bsp_bqp = eigenvectors_C2D_2bsp_state.bottomRightCorner(Nbsp_I, Nbsp_I);

    // T ≤ 10⁻¹² → f = 0; otherwise f = e⁻ᴱᐟᵀ/(1+e⁻ᴱᐟᵀ).
    f_F1D_bqp.setZero();
    if (temperature_F > 1.0e-12) {
        for (int bqp_I = 0; bqp_I < Nbsp_I; ++bqp_I) {
            const double expMinusEOverT_F = std::exp(-Eqp_F1D_bqp(bqp_I) / temperature_F);
            f_F1D_bqp(bqp_I) = expMinusEOverT_F / (1.0 + expMinusEOverT_F);
        }
    }

    // ρ⁺⁺ = DV⁻(1-f)(V⁻)†D + U⁺f(U⁺)†.
    rhoPosPos_C2D_bsp_bsp.noalias() = eta_F1D_bsp.asDiagonal() * VNeg_C2D_bsp_bqp * (1.0 - f_F1D_bqp.array()).matrix().asDiagonal() * VNeg_C2D_bsp_bqp.adjoint() * eta_F1D_bsp.asDiagonal();
    rhoPosPos_C2D_bsp_bsp.noalias() += UPos_C2D_bsp_bqp * f_F1D_bqp.asDiagonal() * UPos_C2D_bsp_bqp.adjoint();

    // κ⁺⁻ = -DV⁻(1-f)(U⁺)†D + U⁺f(V⁻)†.
    kappaPosNeg_C2D_bsp_bsp.noalias() = -(eta_F1D_bsp.asDiagonal() * VNeg_C2D_bsp_bqp * (1.0 - f_F1D_bqp.array()).matrix().asDiagonal() * UPos_C2D_bsp_bqp.adjoint() * eta_F1D_bsp.asDiagonal());
    kappaPosNeg_C2D_bsp_bsp.noalias() += UPos_C2D_bsp_bqp * f_F1D_bqp.asDiagonal() * VNeg_C2D_bsp_bqp.adjoint();
}

/**
 * @brief  Update fields using all representative densities.
 * @math   {ρ_b⁺⁺,κ_b⁺⁻,η_b} → {Γ_b⁺⁺,Δ_b⁺⁻}.
 * @output Overwritten GammaPosPos and DeltaPosNeg.
 */
inline void HFBKramers::update_Gamma_Delta() {
    for (int block_I = 0; block_I < Nblock_I; ++block_I) {
        build_twobody(block_I, blocks, blocks[block_I].GammaPosPos_C2D_bsp_bsp, blocks[block_I].DeltaPosNeg_C2D_bsp_bsp);
    }
}

/**
 * @brief  Solve all blocks using self-adjoint eigendecomposition.
 * @math   {H_b⁺(λ),T} → {E_b,U_b⁺,V_b⁻,f_b,ρ_b⁺⁺,κ_b⁺⁻}.
 * @output Updated representative solutions and densities.
 */
inline void HFBKramers::update_UV_E_rho_kappa(double lambda_F, double temperature_F) {
    for (int block_I = 0; block_I < Nblock_I; ++block_I) {
        blocks[block_I].update_UV_E_rho_kappa(lambda_F, temperature_F);
    }
}
