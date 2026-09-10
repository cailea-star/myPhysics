/**
 * @file    hfb_axial.hpp
 * @author  cailea
 * @date    2026-09-10
 * @brief   Axial HFB without time-reversal symmetry.
 */

#pragma once

#include <cassert>
#include <cmath>
#include <functional>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Eigenvalues>

class HFBAxialBlock {
public:
    int Nbsp_I = 0; // Single-particle dimension per branch.

    Eigen::MatrixXcd h0PosPos_C2D_bsp_bsp{};
    Eigen::MatrixXcd h0NegNeg_C2D_bsp_bsp{};

    Eigen::MatrixXcd GammaPosPos_C2D_bsp_bsp{};
    Eigen::MatrixXcd GammaNegNeg_C2D_bsp_bsp{};
    Eigen::MatrixXcd DeltaPosNeg_C2D_bsp_bsp{};

    Eigen::VectorXd EqpPos_F1D_bqp{};
    Eigen::VectorXd EqpNeg_F1D_bqp{};

    // f = 1/(1+exp(E/T)); T = k_B T_phys.
    Eigen::VectorXd fPos_F1D_bqp{};
    Eigen::VectorXd fNeg_F1D_bqp{};

    // X⁺ = [U⁺; V⁻].
    Eigen::MatrixXcd UPos_C2D_bsp_bqp{};
    Eigen::MatrixXcd VNeg_C2D_bsp_bqp{};

    // X⁻ = [U⁻; V⁺].
    Eigen::MatrixXcd UNeg_C2D_bsp_bqp{};
    Eigen::MatrixXcd VPos_C2D_bsp_bqp{};

    Eigen::MatrixXcd rhoPosPos_C2D_bsp_bsp{};
    Eigen::MatrixXcd rhoNegNeg_C2D_bsp_bsp{};
    Eigen::MatrixXcd kappaPosNeg_C2D_bsp_bsp{};

    Eigen::MatrixXcd HPos_C2D_2bsp_2bsp{};
    Eigen::MatrixXcd HNeg_C2D_2bsp_2bsp{};
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> HPos_eigensolver{};
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> HNeg_eigensolver{};

public:
    /**
     * @brief  Allocate one opposite-Ω block pair.
     * @math   dim(bsp⁺) = dim(bsp⁻) = Nbsp.
     * @output Zeroed matrices and allocated eigensolver workspace.
     */
    HFBAxialBlock(int Nbsp_I_) {
        assert(Nbsp_I_ > 0);
        Nbsp_I = Nbsp_I_;

        // h₀, Γ, Δ, ρ, κ ∈ ℂ^{Nbsp×Nbsp}.
        h0PosPos_C2D_bsp_bsp.resize(Nbsp_I, Nbsp_I);
        h0PosPos_C2D_bsp_bsp.setZero();
        h0NegNeg_C2D_bsp_bsp.resize(Nbsp_I, Nbsp_I);
        h0NegNeg_C2D_bsp_bsp.setZero();
        GammaPosPos_C2D_bsp_bsp.resize(Nbsp_I, Nbsp_I);
        GammaPosPos_C2D_bsp_bsp.setZero();
        GammaNegNeg_C2D_bsp_bsp.resize(Nbsp_I, Nbsp_I);
        GammaNegNeg_C2D_bsp_bsp.setZero();
        DeltaPosNeg_C2D_bsp_bsp.resize(Nbsp_I, Nbsp_I);
        DeltaPosNeg_C2D_bsp_bsp.setZero();
        rhoPosPos_C2D_bsp_bsp.resize(Nbsp_I, Nbsp_I);
        rhoPosPos_C2D_bsp_bsp.setZero();
        rhoNegNeg_C2D_bsp_bsp.resize(Nbsp_I, Nbsp_I);
        rhoNegNeg_C2D_bsp_bsp.setZero();
        kappaPosNeg_C2D_bsp_bsp.resize(Nbsp_I, Nbsp_I);
        kappaPosNeg_C2D_bsp_bsp.setZero();

        // Positive-energy column counts follow diagonalization.
        UPos_C2D_bsp_bqp.resize(Nbsp_I, 0);
        VNeg_C2D_bsp_bqp.resize(Nbsp_I, 0);
        UNeg_C2D_bsp_bqp.resize(Nbsp_I, 0);
        VPos_C2D_bsp_bqp.resize(Nbsp_I, 0);

        // H⁺, H⁻ ∈ ℂ^{2Nbsp×2Nbsp}.
        HPos_C2D_2bsp_2bsp.resize(2 * Nbsp_I, 2 * Nbsp_I);
        HNeg_C2D_2bsp_2bsp.resize(2 * Nbsp_I, 2 * Nbsp_I);
        HPos_C2D_2bsp_2bsp.setZero();
        HNeg_C2D_2bsp_2bsp.setZero();
        HPos_eigensolver = Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd>(2 * Nbsp_I);
        HNeg_eigensolver = Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd>(2 * Nbsp_I);
    }

    /**
     * @brief  Solve thermal branches using self-adjoint eigendecomposition.
     * @math   H⁺ = [h⁺⁺-λI, Δ⁺⁻; (Δ⁺⁻)†, -(h⁻⁻)*+λI].
     * @math   H⁻ = [h⁻⁻-λI, -(Δ⁺⁻)ᵀ; -(Δ⁺⁻)*, -(h⁺⁺)*+λI].
     * @output Updated Eqp, U, V, f, and densities.
     * @note   No zero modes; T = k_B T_phys ≥ 0.
     */
    void update_UV_E_rho_kappa(double lambda_F, double temperature_F);
};

class HFBAxial {
public:
    using OneBodyFunc = std::function<void(int block_I, Eigen::MatrixXcd& h0PosPos_C2D_bsp_bsp, Eigen::MatrixXcd& h0NegNeg_C2D_bsp_bsp)>;

    using TwoBodyFunc = std::function<void(int block_I, const std::vector<HFBAxialBlock>& blocks, Eigen::MatrixXcd& GammaPosPos_C2D_bsp_bsp, Eigen::MatrixXcd& GammaNegNeg_C2D_bsp_bsp, Eigen::MatrixXcd& DeltaPosNeg_C2D_bsp_bsp)>;

    int Nblock_I = 0;

    std::vector<HFBAxialBlock> blocks{};

    OneBodyFunc build_onebody{};
    TwoBodyFunc build_twobody{};

public:
    /**
     * @brief  Allocate blocks and build one-body matrices.
     * @math   {Nbsp_b} → {h₀,b⁺⁺, h₀,b⁻⁻}.
     * @output Initialized blocks and callbacks.
     */
    HFBAxial(const std::vector<int>& Nbsp_I1D_block_, const OneBodyFunc& build_onebody_, const TwoBodyFunc& build_twobody_) {
        assert(build_onebody_);
        assert(build_twobody_);
        Nblock_I = static_cast<int>(Nbsp_I1D_block_.size());
        build_onebody = build_onebody_;
        build_twobody = build_twobody_;

        // {Nbsp_b} → {block_b}.
        blocks.reserve(Nblock_I);
        for (int block_I = 0; block_I < Nblock_I; ++block_I) {
            blocks.emplace_back(Nbsp_I1D_block_[block_I]);
        }

        // block_b → (h₀,b⁺⁺, h₀,b⁻⁻).
        for (int block_I = 0; block_I < Nblock_I; ++block_I) {
            build_onebody(block_I, blocks[block_I].h0PosPos_C2D_bsp_bsp, blocks[block_I].h0NegNeg_C2D_bsp_bsp);
        }
    }

    /**
     * @brief  Update fields using all block densities.
     * @math   {ρ_b,κ_b} → {Γ_b,Δ_b}.
     * @output Overwritten Gamma and Delta in all blocks.
     */
    void update_Gamma_Delta();

    /**
     * @brief  Solve all blocks using self-adjoint eigendecomposition.
     * @math   {H_b(λ),T} → {E_b,U_b,V_b,f_b,ρ_b,κ_b}.
     * @output Updated quasiparticle solutions and densities.
     */
    void update_UV_E_rho_kappa(double lambda_F, double temperature_F);
};


inline void HFBAxialBlock::update_UV_E_rho_kappa(double lambda_F, double temperature_F) {
    assert(temperature_F >= 0.0);
    assert(h0PosPos_C2D_bsp_bsp.isApprox(h0PosPos_C2D_bsp_bsp.adjoint(), 1.0e-12));
    assert(h0NegNeg_C2D_bsp_bsp.isApprox(h0NegNeg_C2D_bsp_bsp.adjoint(), 1.0e-12));
    assert(GammaPosPos_C2D_bsp_bsp.isApprox(GammaPosPos_C2D_bsp_bsp.adjoint(), 1.0e-12));
    assert(GammaNegNeg_C2D_bsp_bsp.isApprox(GammaNegNeg_C2D_bsp_bsp.adjoint(), 1.0e-12));

    // H⁺₁₁ = h⁺⁺ - λI; H⁻₁₁ = h⁻⁻ - λI.
    HPos_C2D_2bsp_2bsp.topLeftCorner(Nbsp_I, Nbsp_I) = h0PosPos_C2D_bsp_bsp + GammaPosPos_C2D_bsp_bsp;
    HNeg_C2D_2bsp_2bsp.topLeftCorner(Nbsp_I, Nbsp_I) = h0NegNeg_C2D_bsp_bsp + GammaNegNeg_C2D_bsp_bsp;
    HPos_C2D_2bsp_2bsp.topLeftCorner(Nbsp_I, Nbsp_I).diagonal().array() -= lambda_F;
    HNeg_C2D_2bsp_2bsp.topLeftCorner(Nbsp_I, Nbsp_I).diagonal().array() -= lambda_F;

    // H⁺ = [h⁺⁺-λI, Δ; Δ†, -(h⁻⁻)*+λI].
    HPos_C2D_2bsp_2bsp.topRightCorner(Nbsp_I, Nbsp_I) = DeltaPosNeg_C2D_bsp_bsp;
    HPos_C2D_2bsp_2bsp.bottomLeftCorner(Nbsp_I, Nbsp_I) = DeltaPosNeg_C2D_bsp_bsp.adjoint();
    HPos_C2D_2bsp_2bsp.bottomRightCorner(Nbsp_I, Nbsp_I) = -HNeg_C2D_2bsp_2bsp.topLeftCorner(Nbsp_I, Nbsp_I).conjugate();

    // H⁻ = [h⁻⁻-λI, -Δᵀ; -Δ*, -(h⁺⁺)*+λI].
    HNeg_C2D_2bsp_2bsp.topRightCorner(Nbsp_I, Nbsp_I) = -DeltaPosNeg_C2D_bsp_bsp.transpose();
    HNeg_C2D_2bsp_2bsp.bottomLeftCorner(Nbsp_I, Nbsp_I) = -DeltaPosNeg_C2D_bsp_bsp.conjugate();
    HNeg_C2D_2bsp_2bsp.bottomRightCorner(Nbsp_I, Nbsp_I) = -HPos_C2D_2bsp_2bsp.topLeftCorner(Nbsp_I, Nbsp_I).conjugate();

    // H± → (E±, X±).
    HPos_eigensolver.compute(HPos_C2D_2bsp_2bsp);
    HNeg_eigensolver.compute(HNeg_C2D_2bsp_2bsp);
    assert(HPos_eigensolver.info() == Eigen::Success);
    assert(HNeg_eigensolver.info() == Eigen::Success);
    const Eigen::VectorXd& eigenvaluesPos_F1D_state = HPos_eigensolver.eigenvalues();
    const Eigen::VectorXd& eigenvaluesNeg_F1D_state = HNeg_eigensolver.eigenvalues();
    const Eigen::MatrixXcd& eigenvectorsPos_C2D_2bsp_state = HPos_eigensolver.eigenvectors();
    const Eigen::MatrixXcd& eigenvectorsNeg_C2D_2bsp_state = HNeg_eigensolver.eigenvectors();
    assert((eigenvaluesPos_F1D_state.array() != 0.0).all());
    assert((eigenvaluesNeg_F1D_state.array() != 0.0).all());

    // E± > 0; Nqp⁺ + Nqp⁻ = 2Nbsp.
    EqpPos_F1D_bqp.resize((eigenvaluesPos_F1D_state.array() > 0.0).count());
    EqpNeg_F1D_bqp.resize((eigenvaluesNeg_F1D_state.array() > 0.0).count());
    assert(EqpPos_F1D_bqp.size() + EqpNeg_F1D_bqp.size() == 2 * Nbsp_I);
    EqpPos_F1D_bqp = eigenvaluesPos_F1D_state.tail(EqpPos_F1D_bqp.size());
    EqpNeg_F1D_bqp = eigenvaluesNeg_F1D_state.tail(EqpNeg_F1D_bqp.size());

    // X⁺ = [U⁺; V⁻]; X⁻ = [U⁻; V⁺].
    UPos_C2D_bsp_bqp = eigenvectorsPos_C2D_2bsp_state.topRightCorner(Nbsp_I, EqpPos_F1D_bqp.size());
    VNeg_C2D_bsp_bqp = eigenvectorsPos_C2D_2bsp_state.bottomRightCorner(Nbsp_I, EqpPos_F1D_bqp.size());
    UNeg_C2D_bsp_bqp = eigenvectorsNeg_C2D_2bsp_state.topRightCorner(Nbsp_I, EqpNeg_F1D_bqp.size());
    VPos_C2D_bsp_bqp = eigenvectorsNeg_C2D_2bsp_state.bottomRightCorner(Nbsp_I, EqpNeg_F1D_bqp.size());

    // T ≤ 10⁻¹² → f± = 0.
    fPos_F1D_bqp.resize(EqpPos_F1D_bqp.size());
    fNeg_F1D_bqp.resize(EqpNeg_F1D_bqp.size());
    fPos_F1D_bqp.setZero();
    fNeg_F1D_bqp.setZero();
    if (temperature_F > 1.0e-12) {
        // f± = e⁻ᴱ±ᐟᵀ/(1+e⁻ᴱ±ᐟᵀ).
        for (int bqp_I = 0; bqp_I < EqpPos_F1D_bqp.size(); ++bqp_I) {
            const double expMinusEOverT_F = std::exp(-EqpPos_F1D_bqp(bqp_I) / temperature_F);
            fPos_F1D_bqp(bqp_I) = expMinusEOverT_F / (1.0 + expMinusEOverT_F);
        }
        for (int bqp_I = 0; bqp_I < EqpNeg_F1D_bqp.size(); ++bqp_I) {
            const double expMinusEOverT_F = std::exp(-EqpNeg_F1D_bqp(bqp_I) / temperature_F);
            fNeg_F1D_bqp(bqp_I) = expMinusEOverT_F / (1.0 + expMinusEOverT_F);
        }
    }

    // ρ⁺⁺ = (V⁺)*(1-f⁻)(V⁺)ᵀ + U⁺f⁺(U⁺)†.
    rhoPosPos_C2D_bsp_bsp.noalias() = VPos_C2D_bsp_bqp.conjugate() * (1.0 - fNeg_F1D_bqp.array()).matrix().asDiagonal() * VPos_C2D_bsp_bqp.transpose();
    rhoPosPos_C2D_bsp_bsp.noalias() += UPos_C2D_bsp_bqp * fPos_F1D_bqp.asDiagonal() * UPos_C2D_bsp_bqp.adjoint();

    // ρ⁻⁻ = (V⁻)*(1-f⁺)(V⁻)ᵀ + U⁻f⁻(U⁻)†.
    rhoNegNeg_C2D_bsp_bsp.noalias() = VNeg_C2D_bsp_bqp.conjugate() * (1.0 - fPos_F1D_bqp.array()).matrix().asDiagonal() * VNeg_C2D_bsp_bqp.transpose();
    rhoNegNeg_C2D_bsp_bsp.noalias() += UNeg_C2D_bsp_bqp * fNeg_F1D_bqp.asDiagonal() * UNeg_C2D_bsp_bqp.adjoint();

    // κ⁺⁻ = (V⁺)*(1-f⁻)(U⁻)ᵀ + U⁺f⁺(V⁻)†.
    kappaPosNeg_C2D_bsp_bsp.noalias() = VPos_C2D_bsp_bqp.conjugate() * (1.0 - fNeg_F1D_bqp.array()).matrix().asDiagonal() * UNeg_C2D_bsp_bqp.transpose();
    kappaPosNeg_C2D_bsp_bsp.noalias() += UPos_C2D_bsp_bqp * fPos_F1D_bqp.asDiagonal() * VNeg_C2D_bsp_bqp.adjoint();
}


inline void HFBAxial::update_Gamma_Delta() {
    for (int block_I = 0; block_I < Nblock_I; ++block_I) {
        build_twobody(block_I, blocks, blocks[block_I].GammaPosPos_C2D_bsp_bsp, blocks[block_I].GammaNegNeg_C2D_bsp_bsp, blocks[block_I].DeltaPosNeg_C2D_bsp_bsp);
    }
}


inline void HFBAxial::update_UV_E_rho_kappa(double lambda_F, double temperature_F) {
    for (int block_I = 0; block_I < Nblock_I; ++block_I) {
        blocks[block_I].update_UV_E_rho_kappa(lambda_F, temperature_F);
    }
}
