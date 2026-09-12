/**
 * @file    hfb_kramers.hpp
 * @author  cailea
 * @date    2026-09-10
 * @brief   Time-reversal-invariant axial HFB calculations.
 */

#pragma once

#include <cassert>
#include <cmath>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Eigenvalues>

struct HFBKramersBlockSolution {
    // T|α⟩ = ηα|ᾱ⟩; ηα = (-1)^(j - mj).
    Eigen::VectorXd eta_F1D_bsp{};

    Eigen::VectorXd Eqp_F1D_bqp{};
    Eigen::VectorXd f_F1D_bqp{};

    // X⁺ = [U⁺; V⁻].
    Eigen::MatrixXd UPos_F2D_bsp_bqp{};
    Eigen::MatrixXd VNeg_F2D_bsp_bqp{};

    Eigen::MatrixXd rhoPosPos_F2D_bsp_bsp{};
    Eigen::MatrixXd kappaPosNeg_F2D_bsp_bsp{};
};

struct HFBKramersBlockField {
    Eigen::MatrixXd h0PosPos_F2D_bsp_bsp{};
    Eigen::MatrixXd GammaPosPos_F2D_bsp_bsp{};
    Eigen::MatrixXd DeltaPosNeg_F2D_bsp_bsp{};
};

class HFBKramers {
public:
    double lambda_F = -7.0; // Fermi energy [MeV].
    double lambda2_F = 0.0; // Lipkin-Nogami λ₂ [MeV].
    double temperature_F = 0.0;
    double EspCut_F = 60.0; // Equivalent single-particle energy cutoff [MeV].

    int Nblock_I = 0;
    std::vector<int> Nbsp_I1D_block{};

    std::vector<HFBKramersBlockSolution> solutions{};
    std::vector<HFBKramersBlockField> fields{};

    std::vector<Eigen::MatrixXd> HPos_F3D_block_2bsp_2bsp{};
    std::vector<Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>> HPos_eigensolvers{};

public:
    /**
     * @brief  Allocate single-species Kramers-representative blocks.
     * @math   Nbsp_b = size(η_b); ηα² = 1.
     * @output Stored phases, zeroed fields, solutions, allocated workspaces.
     */
    HFBKramers(const std::vector<Eigen::VectorXd>& eta_F2D_block_bsp_) {
        Nblock_I = static_cast<int>(eta_F2D_block_bsp_.size());

        // {η_b} → {field_b, solution_b, workspace_b}.
        Nbsp_I1D_block.resize(Nblock_I);
        solutions.resize(Nblock_I);
        fields.resize(Nblock_I);
        HPos_F3D_block_2bsp_2bsp.resize(Nblock_I);
        HPos_eigensolvers.resize(Nblock_I);

        for (int block_I = 0; block_I < Nblock_I; ++block_I) {
            const int Nbsp_I = static_cast<int>(eta_F2D_block_bsp_[block_I].size());
            assert(Nbsp_I > 0);
            assert((eta_F2D_block_bsp_[block_I].array().abs() == 1.0).all());
            Nbsp_I1D_block[block_I] = Nbsp_I;
            HFBKramersBlockField& field = fields[block_I];
            HFBKramersBlockSolution& solution = solutions[block_I];
            solution.eta_F1D_bsp = eta_F2D_block_bsp_[block_I];

            // h₀, Γ, Δ ∈ ℝ^{Nbsp×Nbsp}.
            field.h0PosPos_F2D_bsp_bsp.resize(Nbsp_I, Nbsp_I);
            field.h0PosPos_F2D_bsp_bsp.setZero();
            field.GammaPosPos_F2D_bsp_bsp.resize(Nbsp_I, Nbsp_I);
            field.GammaPosPos_F2D_bsp_bsp.setZero();
            field.DeltaPosNeg_F2D_bsp_bsp.resize(Nbsp_I, Nbsp_I);
            field.DeltaPosNeg_F2D_bsp_bsp.setZero();

            // E, f ∈ ℝ^{Nbsp}.
            solution.Eqp_F1D_bqp.resize(Nbsp_I);
            solution.Eqp_F1D_bqp.setZero();
            solution.f_F1D_bqp.resize(Nbsp_I);
            solution.f_F1D_bqp.setZero();

            // U⁺, V⁻, ρ⁺⁺, κ⁺⁻ ∈ ℝ^{Nbsp×Nbsp}.
            solution.UPos_F2D_bsp_bqp.resize(Nbsp_I, Nbsp_I);
            solution.UPos_F2D_bsp_bqp.setZero();
            solution.VNeg_F2D_bsp_bqp.resize(Nbsp_I, Nbsp_I);
            solution.VNeg_F2D_bsp_bqp.setZero();
            solution.rhoPosPos_F2D_bsp_bsp.resize(Nbsp_I, Nbsp_I);
            solution.rhoPosPos_F2D_bsp_bsp.setZero();
            solution.kappaPosNeg_F2D_bsp_bsp.resize(Nbsp_I, Nbsp_I);
            solution.kappaPosNeg_F2D_bsp_bsp.setZero();

            // H⁺ ∈ ℝ^{2Nbsp×2Nbsp}.
            HPos_F3D_block_2bsp_2bsp[block_I].resize(2 * Nbsp_I, 2 * Nbsp_I);
            HPos_F3D_block_2bsp_2bsp[block_I].setZero();
            HPos_eigensolvers[block_I] = Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>(2 * Nbsp_I);
        }
    }

    /**
     * @brief  Solve thermal HFB using real symmetric eigendecomposition.
     * @math   h = h₀⁺⁺ + Γ⁺⁺ + 4λ₂ρ⁺⁺ - 2λ₂I.
     * @math   H⁺ = [h-λI, Δ⁺⁻; (Δ⁺⁻)ᵀ, -diag(η)h diag(η)+λI].
     * @output Updated representative solutions and densities.
     * @output Mean particle number N = 2Σ_b Tr(ρ_b⁺⁺).
     * @note   Real, time-reversal-invariant fields; no zero modes; T ≥ 0.
     * @note   ε = λ + E(1-2‖V‖²); cutoff affects ρ and κ.
     * @note   EspCut_F = ∞ disables the cutoff.
     * @note   h_LN = h₀ + Γ + 4λ₂ρ - 2λ₂I; pre-solve ρ.
     */
    double update_UV_E_rho_kappa();
};

class HFBKramersNucleus {
public:
    HFBKramers hfb_neutron;
    HFBKramers hfb_proton;

public:
    /**
     * @brief  Initialize species dimensions and workspaces.
     * @math   {η_q,b} → {Nbsp_q,b}; q ∈ {n,p}.
     * @output Allocated solutions, fields, and workspaces.
     */
    HFBKramersNucleus(const std::vector<Eigen::VectorXd>& etaN_F2D_block_bsp_, const std::vector<Eigen::VectorXd>& etaP_F2D_block_bsp_)
    : hfb_neutron(etaN_F2D_block_bsp_), hfb_proton(etaP_F2D_block_bsp_) {}

    /**
     * @brief  Initialize one-body fields in the derived model.
     * @math   h₀,n, h₀,p.
     * @output Initialized neutron and proton one-body fields.
     */
    virtual void initialize_h0() = 0;

    /**
     * @brief  Initialize HFB fields in the derived model.
     * @math   (N,Z) → (Γ_n,Δ_n,Γ_p,Δ_p)_initial.
     * @output Initialized Gamma and Delta for both species.
     */
    virtual void initialize_GammaDelta(int TargetN_I, int TargetZ_I) = 0;

    /**
     * @brief  Update both species using all representative densities.
     * @math   {ρ_q,b⁺⁺,κ_q,b⁺⁻,η_q,b} → {Γ_q,b⁺⁺,Δ_q,b⁺⁻}.
     * @output Overwritten neutron and proton Gamma and Delta.
     */
    virtual void update_Gamma_Delta() = 0;

    /**
     * @brief  Iterate the unblocked HFB equations.
     * @math   (N,Z) → HFB_converged.
     * @output Updated converged solver state.
     */
    void iterate(int TargetN_I, int TargetZ_I);
};

/**
 * @brief  Solve thermal HFB using real symmetric eigendecomposition.
 * @math   H⁺X⁺ = X⁺E; U⁻ = diag(η)U⁺; V⁺ = -diag(η)V⁻.
 * @output Updated representative solutions and densities.
 * @output Mean particle number N = 2Σ_b Tr(ρ_b⁺⁺).
 */
inline double HFBKramers::update_UV_E_rho_kappa() {
    assert(temperature_F >= 0.0);
    assert(static_cast<int>(fields.size()) == Nblock_I);
    assert(static_cast<int>(solutions.size()) == Nblock_I);
    double N_F = 0.0;

    // block_b → (E_b, U_b⁺, V_b⁻, f_b, ρ_b⁺⁺, κ_b⁺⁻).
    for (int block_I = 0; block_I < Nblock_I; ++block_I) {
        const int Nbsp_I = Nbsp_I1D_block[block_I];
        const HFBKramersBlockField& field = fields[block_I];
        HFBKramersBlockSolution& solution = solutions[block_I];
        Eigen::MatrixXd& HPos_F2D_2bsp_2bsp = HPos_F3D_block_2bsp_2bsp[block_I];
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>& HPos_eigensolver = HPos_eigensolvers[block_I];

        assert(field.h0PosPos_F2D_bsp_bsp.isApprox(field.h0PosPos_F2D_bsp_bsp.transpose(), 1.0e-12));
        assert(field.GammaPosPos_F2D_bsp_bsp.isApprox(field.GammaPosPos_F2D_bsp_bsp.transpose(), 1.0e-12));
        // Δ = diag(η)Δᵀdiag(η).
        assert(field.DeltaPosNeg_F2D_bsp_bsp.isApprox(solution.eta_F1D_bsp.asDiagonal() * field.DeltaPosNeg_F2D_bsp_bsp.transpose() * solution.eta_F1D_bsp.asDiagonal(), 1.0e-12));

        // H⁺₁₁ = h₀⁺⁺ + Γ⁺⁺ + 4λ₂ρ⁺⁺ - (λ+2λ₂)I.
        HPos_F2D_2bsp_2bsp.topLeftCorner(Nbsp_I, Nbsp_I) = field.h0PosPos_F2D_bsp_bsp + field.GammaPosPos_F2D_bsp_bsp + 4.0 * lambda2_F * solution.rhoPosPos_F2D_bsp_bsp;
        HPos_F2D_2bsp_2bsp.topLeftCorner(Nbsp_I, Nbsp_I).diagonal().array() -= lambda_F + 2.0 * lambda2_F;

        // H⁺₁₂ = Δ; H⁺₂₁ = Δᵀ; H⁺₂₂ = -diag(η)H⁺₁₁diag(η).
        HPos_F2D_2bsp_2bsp.topRightCorner(Nbsp_I, Nbsp_I) = field.DeltaPosNeg_F2D_bsp_bsp;
        HPos_F2D_2bsp_2bsp.bottomLeftCorner(Nbsp_I, Nbsp_I) = field.DeltaPosNeg_F2D_bsp_bsp.transpose();
        HPos_F2D_2bsp_2bsp.bottomRightCorner(Nbsp_I, Nbsp_I).noalias() = -(solution.eta_F1D_bsp.asDiagonal() * HPos_F2D_2bsp_2bsp.topLeftCorner(Nbsp_I, Nbsp_I) * solution.eta_F1D_bsp.asDiagonal());

        // H⁺ → (E, X⁺).
        HPos_eigensolver.compute(HPos_F2D_2bsp_2bsp);
        assert(HPos_eigensolver.info() == Eigen::Success);
        const Eigen::VectorXd& eigenvalues_F1D_state = HPos_eigensolver.eigenvalues();
        const Eigen::MatrixXd& eigenvectors_F2D_2bsp_state = HPos_eigensolver.eigenvectors();

        // E₁ ≤ ⋯ ≤ E₂Nbsp; retain Nbsp positive-energy modes.
        assert(eigenvalues_F1D_state(Nbsp_I - 1) < 0.0);
        assert(eigenvalues_F1D_state(Nbsp_I) > 0.0);
        solution.Eqp_F1D_bqp = eigenvalues_F1D_state.tail(Nbsp_I);
        solution.UPos_F2D_bsp_bqp = eigenvectors_F2D_2bsp_state.topRightCorner(Nbsp_I, Nbsp_I);
        solution.VNeg_F2D_bsp_bqp = eigenvectors_F2D_2bsp_state.bottomRightCorner(Nbsp_I, Nbsp_I);

        // ε = λ + E(1-2‖V‖²); ε > Ecut + tail → factors = 0.
        const double EspCutTail_F = std::log(1.0 / 1.0e-6 - 1.0) / 100.0;
        Eigen::VectorXd factorU_F1D_bqp{};
        Eigen::VectorXd factorV_F1D_bqp{};
        factorU_F1D_bqp.resize(Nbsp_I);
        factorV_F1D_bqp.resize(Nbsp_I);
        factorU_F1D_bqp.setZero();
        factorV_F1D_bqp.setZero();

        // f = (1-tanh(E/(2T)))/2; active factors = (f,1-f).
        for (int bqp_I = 0; bqp_I < Nbsp_I; ++bqp_I) {
            const double Eqp_F = solution.Eqp_F1D_bqp(bqp_I);
            solution.f_F1D_bqp(bqp_I) = temperature_F > 1.0e-12 ? 0.5 * (1.0 - std::tanh(0.5 * Eqp_F / temperature_F)) : 0.0;
            const double Esp_F = lambda_F + Eqp_F * (1.0 - 2.0 * solution.VNeg_F2D_bsp_bqp.col(bqp_I).squaredNorm());
            if (Esp_F > EspCut_F + EspCutTail_F) {continue;}
            factorU_F1D_bqp(bqp_I) = solution.f_F1D_bqp(bqp_I);
            factorV_F1D_bqp(bqp_I) = 1.0 - solution.f_F1D_bqp(bqp_I);
        }

        // ρ⁺⁺ = diag(η)V⁻(1-f)(V⁻)ᵀdiag(η) + U⁺f(U⁺)ᵀ.
        solution.rhoPosPos_F2D_bsp_bsp.noalias() = solution.eta_F1D_bsp.asDiagonal() * solution.VNeg_F2D_bsp_bqp * factorV_F1D_bqp.asDiagonal() * solution.VNeg_F2D_bsp_bqp.transpose() * solution.eta_F1D_bsp.asDiagonal();
        solution.rhoPosPos_F2D_bsp_bsp.noalias() += solution.UPos_F2D_bsp_bqp * factorU_F1D_bqp.asDiagonal() * solution.UPos_F2D_bsp_bqp.transpose();

        // κ⁺⁻ = -diag(η)V⁻(1-f)(U⁺)ᵀdiag(η) + U⁺f(V⁻)ᵀ.
        solution.kappaPosNeg_F2D_bsp_bsp.noalias() = -(solution.eta_F1D_bsp.asDiagonal() * solution.VNeg_F2D_bsp_bqp * factorV_F1D_bqp.asDiagonal() * solution.UPos_F2D_bsp_bqp.transpose() * solution.eta_F1D_bsp.asDiagonal());
        solution.kappaPosNeg_F2D_bsp_bsp.noalias() += solution.UPos_F2D_bsp_bqp * factorU_F1D_bqp.asDiagonal() * solution.VNeg_F2D_bsp_bqp.transpose();

        // N = 2Σ_b Tr(ρ_b⁺⁺).
        N_F += 2.0 * solution.rhoPosPos_F2D_bsp_bsp.trace();
    }
    return N_F;
}
