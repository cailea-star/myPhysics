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
     * @math   h = h₀⁺⁺ + Γ⁺⁺.
     * @math   H⁺ = [h-λI, Δ⁺⁻; (Δ⁺⁻)ᵀ, -diag(η)h diag(η)+λI].
     * @output Updated representative solutions and densities.
     * @output Mean particle number N = 2Σ_b Tr(ρ_b⁺⁺).
     * @note   Real, time-reversal-invariant fields; no zero modes; T ≥ 0.
     */
    double update_UV_E_rho_kappa(double lambda_F, double temperature_F);
};

class HFBKramersNucleus {
public:
    using OneBodyFunc = std::function<void(int block_I, Eigen::MatrixXd& h0PosPos_F2D_bsp_bsp)>;

    using TwoBodyFunc = std::function<void(const std::vector<HFBKramersBlockSolution>& solutions_n, const std::vector<HFBKramersBlockSolution>& solutions_p, std::vector<HFBKramersBlockField>& fields_n, std::vector<HFBKramersBlockField>& fields_p)>;

    HFBKramers hfb_neutron;
    HFBKramers hfb_proton;

    OneBodyFunc build_onebody_neutron{};
    OneBodyFunc build_onebody_proton{};
    TwoBodyFunc build_twobody{};

public:
    /**
     * @brief  Initialize species and build their one-body matrices.
     * @math   {η_q,b} → {Nbsp_q,b, h₀,q,b⁺⁺}; q ∈ {n,p}.
     * @output Initialized species and stored field callbacks.
     */
    HFBKramersNucleus(const std::vector<Eigen::VectorXd>& etaN_F2D_block_bsp_, const std::vector<Eigen::VectorXd>& etaP_F2D_block_bsp_, const OneBodyFunc& build_onebodyN_, const OneBodyFunc& build_onebodyP_, const TwoBodyFunc& build_twobody_)
    : hfb_neutron(etaN_F2D_block_bsp_), hfb_proton(etaP_F2D_block_bsp_) {
        assert(build_onebodyN_);
        assert(build_onebodyP_);
        assert(build_twobody_);
        build_onebody_neutron = build_onebodyN_;
        build_onebody_proton = build_onebodyP_;
        build_twobody = build_twobody_;

        // block_n → h₀,n⁺⁺.
        for (int block_I = 0; block_I < hfb_neutron.Nblock_I; ++block_I) {
            build_onebody_neutron(block_I, hfb_neutron.fields[block_I].h0PosPos_F2D_bsp_bsp);
        }

        // block_p → h₀,p⁺⁺.
        for (int block_I = 0; block_I < hfb_proton.Nblock_I; ++block_I) {
            build_onebody_proton(block_I, hfb_proton.fields[block_I].h0PosPos_F2D_bsp_bsp);
        }
    }

    /**
     * @brief  Update both species using all representative densities.
     * @math   {ρ_q,b⁺⁺,κ_q,b⁺⁻,η_q,b} → {Γ_q,b⁺⁺,Δ_q,b⁺⁻}.
     * @output Overwritten neutron and proton Gamma and Delta.
     */
    void update_Gamma_Delta();
};

/**
 * @brief  Solve thermal HFB using real symmetric eigendecomposition.
 * @math   H⁺X⁺ = X⁺E; U⁻ = diag(η)U⁺; V⁺ = -diag(η)V⁻.
 * @output Updated representative solutions and densities.
 * @output Mean particle number N = 2Σ_b Tr(ρ_b⁺⁺).
 */
inline double HFBKramers::update_UV_E_rho_kappa(double lambda_F, double temperature_F) {
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

        // H⁺₁₁ = h₀⁺⁺ + Γ⁺⁺ - λI.
        HPos_F2D_2bsp_2bsp.topLeftCorner(Nbsp_I, Nbsp_I) = field.h0PosPos_F2D_bsp_bsp + field.GammaPosPos_F2D_bsp_bsp;
        HPos_F2D_2bsp_2bsp.topLeftCorner(Nbsp_I, Nbsp_I).diagonal().array() -= lambda_F;

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

        // T ≤ 10⁻¹² → f = 0; otherwise f = e⁻ᴱᐟᵀ/(1+e⁻ᴱᐟᵀ).
        solution.f_F1D_bqp.setZero();
        if (temperature_F > 1.0e-12) {
            for (int bqp_I = 0; bqp_I < Nbsp_I; ++bqp_I) {
                const double expMinusEOverT_F = std::exp(-solution.Eqp_F1D_bqp(bqp_I) / temperature_F);
                solution.f_F1D_bqp(bqp_I) = expMinusEOverT_F / (1.0 + expMinusEOverT_F);
            }
        }

        // ρ⁺⁺ = diag(η)V⁻(1-f)(V⁻)ᵀdiag(η) + U⁺f(U⁺)ᵀ.
        solution.rhoPosPos_F2D_bsp_bsp.noalias() = solution.eta_F1D_bsp.asDiagonal() * solution.VNeg_F2D_bsp_bqp * (1.0 - solution.f_F1D_bqp.array()).matrix().asDiagonal() * solution.VNeg_F2D_bsp_bqp.transpose() * solution.eta_F1D_bsp.asDiagonal();
        solution.rhoPosPos_F2D_bsp_bsp.noalias() += solution.UPos_F2D_bsp_bqp * solution.f_F1D_bqp.asDiagonal() * solution.UPos_F2D_bsp_bqp.transpose();

        // κ⁺⁻ = -diag(η)V⁻(1-f)(U⁺)ᵀdiag(η) + U⁺f(V⁻)ᵀ.
        solution.kappaPosNeg_F2D_bsp_bsp.noalias() = -(solution.eta_F1D_bsp.asDiagonal() * solution.VNeg_F2D_bsp_bqp * (1.0 - solution.f_F1D_bqp.array()).matrix().asDiagonal() * solution.UPos_F2D_bsp_bqp.transpose() * solution.eta_F1D_bsp.asDiagonal());
        solution.kappaPosNeg_F2D_bsp_bsp.noalias() += solution.UPos_F2D_bsp_bqp * solution.f_F1D_bqp.asDiagonal() * solution.VNeg_F2D_bsp_bqp.transpose();

        // N = 2Σ_b Tr(ρ_b⁺⁺).
        N_F += 2.0 * solution.rhoPosPos_F2D_bsp_bsp.trace();
    }
    return N_F;
}

/**
 * @brief  Update both species using all representative densities.
 * @math   {ρ_q,b⁺⁺,κ_q,b⁺⁻,η_q,b} → {Γ_q,b⁺⁺,Δ_q,b⁺⁻}.
 * @output Overwritten neutron and proton Gamma and Delta.
 */
inline void HFBKramersNucleus::update_Gamma_Delta() {
    build_twobody(hfb_neutron.solutions, hfb_proton.solutions, hfb_neutron.fields, hfb_proton.fields);
}
