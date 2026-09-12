/**
 * @file    hfb_kramers.hpp
 * @author  cailea
 * @date    2026-09-10
 * @brief   Time-reversal-invariant axial HFB calculations.
 */

#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <map>
#include <numeric>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Eigenvalues>

#include "root.hpp"

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
    struct Element {
        double vPosPosPosPos_F = 0.0;
        double vPosNegPosNeg_F = 0.0;
    };

    using BlockingFunc = std::function<double(std::vector<HFBKramersBlockSolution>& solutions, bool updateTracking_B)>;
    // v̄_{12;34}: Γ uses block13, block24; Δ uses block12, block34.
    // bsp indices follow 1,2,3,4.
    using GammaElementFunc = std::function<Element(int block13_I, int block24_I, int bsp1_I, int bsp2_I, int bsp3_I, int bsp4_I)>;
    using DeltaElementFunc = std::function<Element(int block12_I, int block34_I, int bsp1_I, int bsp2_I, int bsp3_I, int bsp4_I)>;

    int TargetN_I = 0; // Target particle number.
    double lambda_F = -7.0; // Fermi energy [MeV].
    double lambda2_F = 0.0; // Lipkin-Nogami λ₂ [MeV].
    double temperature_F = 0.0;
    double EspCut_F = 60.0; // Equivalent single-particle energy cutoff [MeV].
    double ELipkinNogami_F = 0.0; // Lipkin-Nogami energy [MeV].

    int Nblock_I = 0;
    std::vector<int> Nbsp_I1D_block{};

    std::vector<HFBKramersBlockSolution> solutions{};
    BlockingFunc blocking_Func{}; // updateTracking_B: commit blocking trackers.
    std::vector<HFBKramersBlockField> fields{};

    std::vector<Eigen::MatrixXd> HPos_F3D_block_2bsp_2bsp{};
    std::vector<Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>> HPos_eigensolvers{};

public:
    /**
     * @brief  Allocate single-species Kramers-representative blocks.
     * @math   Nbsp_b = size(η_b); ηα² = 1.
     * @output Stored phases; allocated fields, solutions, and workspaces.
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
            field.GammaPosPos_F2D_bsp_bsp.resize(Nbsp_I, Nbsp_I);
            field.DeltaPosNeg_F2D_bsp_bsp.resize(Nbsp_I, Nbsp_I);

            // E, f ∈ ℝ^{Nbsp}.
            solution.Eqp_F1D_bqp.resize(Nbsp_I);
            solution.f_F1D_bqp.resize(Nbsp_I);

            // U⁺, V⁻, ρ⁺⁺, κ⁺⁻ ∈ ℝ^{Nbsp×Nbsp}.
            solution.UPos_F2D_bsp_bqp.resize(Nbsp_I, Nbsp_I);
            solution.VNeg_F2D_bsp_bqp.resize(Nbsp_I, Nbsp_I);
            solution.rhoPosPos_F2D_bsp_bsp.resize(Nbsp_I, Nbsp_I);
            solution.kappaPosNeg_F2D_bsp_bsp.resize(Nbsp_I, Nbsp_I);

            // H⁺ ∈ ℝ^{2Nbsp×2Nbsp}.
            HPos_F3D_block_2bsp_2bsp[block_I].resize(2 * Nbsp_I, 2 * Nbsp_I);
            HPos_eigensolvers[block_I] = Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>(2 * Nbsp_I);
        }
    }

    /**
     * @brief  Solve thermal HFB using real symmetric eigendecomposition.
     * @math   h = h₀⁺⁺ + Γ⁺⁺.
     * @math   H⁺ = [h-λI, Δ⁺⁻; (Δ⁺⁻)ᵀ, -diag(η)h diag(η)+λI].
     * @output Updated representative solutions and densities.
     * @output Mean particle number N = 2Σ_b Tr(ρ_b⁺⁺).
     * @note   ε = λ+E(1-2‖V‖²) cuts ρ,κ; Ecut = ∞ disables.
     */
    double update_UV_E_rho_kappa();

    /**
     * @brief  Search chemical potential using bracket expansion and Brent.
     * @math   N_blocked(λ) = TargetN.
     * @output Updated chemical potential, blocked solutions, and densities.
     */
    void search_lambda(double lambdaTolerance_F);

    /**
     * @brief Accumulate particle-hole fields by direct matrix-element contraction.
     * @math Γ_{13} += Σ_{24}v̄_{12;34}ρ_42.
     * @output Accumulated Gamma from this species' density.
     * @note Requires time-reversal-invariant matrix elements and density.
     * @note Requires thread-safe matrix-element callbacks.
     */
    void add_Gamma_from_Element(const GammaElementFunc& read_element_Func);

    /**
     * @brief Accumulate pairing fields by direct matrix-element contraction.
     * @math Δ⁺⁻_{12} += Σ_{34}v̄⁺⁻⁺⁻_{12;34}κ⁺⁻_{34}.
     * @output Accumulated Delta from this species' pairing tensor.
     * @note Requires time-reversal-invariant matrix elements.
     * @note Requires thread-safe matrix-element callbacks.
     */
    void add_Delta_from_Element(const DeltaElementFunc& read_element_Func);
};

inline void HFBKramers::add_Gamma_from_Element(const GammaElementFunc& read_element_Func) {
    assert(read_element_Func);

    // Γ⁺⁺_{13} += Σ_{24}(v̄⁺⁺⁺⁺ρ⁺⁺_{42}+v̄⁺⁻⁺⁻ρ⁻⁻_{42}).
    const auto add_Gamma_Func = [&](int block13_I, int bsp1_I, int bsp3_I) {
        double Gamma13PosPos_F = 0.0;
        for (int block24_I = 0; block24_I < Nblock_I; ++block24_I) {
            const auto& solution = solutions[block24_I];
            for (int bsp2_I = 0; bsp2_I < Nbsp_I1D_block[block24_I]; ++bsp2_I) {
                for (int bsp4_I = 0; bsp4_I < Nbsp_I1D_block[block24_I]; ++bsp4_I) {
                    const Element element = read_element_Func(block13_I, block24_I, bsp1_I, bsp2_I, bsp3_I, bsp4_I);
                    const double rho42PosPos_F = solution.rhoPosPos_F2D_bsp_bsp(bsp4_I, bsp2_I);
                    // ρ⁻⁻_{42} = η₄η₂ρ⁺⁺_{42} (real).
                    const double rho42NegNeg_F = solution.eta_F1D_bsp(bsp4_I) * solution.eta_F1D_bsp(bsp2_I) * rho42PosPos_F;
                    Gamma13PosPos_F += element.vPosPosPosPos_F * rho42PosPos_F + element.vPosNegPosNeg_F * rho42NegNeg_F;
                }
            }
        }
        fields[block13_I].GammaPosPos_F2D_bsp_bsp(bsp1_I, bsp3_I) += Gamma13PosPos_F;
    };

    // #pragma omp parallel
    {
        // #pragma omp single
        for (int block13_I = 0; block13_I < Nblock_I; ++block13_I) {
            for (int bsp1_I = 0; bsp1_I < Nbsp_I1D_block[block13_I]; ++bsp1_I) {
                for (int bsp3_I = 0; bsp3_I < Nbsp_I1D_block[block13_I]; ++bsp3_I) {
                    // #pragma omp task firstprivate(block13_I, bsp1_I, bsp3_I)
                    {
                        add_Gamma_Func(block13_I, bsp1_I, bsp3_I);
                    }
                }
            }
        }
    }

}

inline void HFBKramers::add_Delta_from_Element(const DeltaElementFunc& read_element_Func) {
    assert(read_element_Func);

    // κ⁻⁺ = −(κ⁺⁻)ᵀ and v̄_{12;43} = −v̄_{12;34} cancel ½.
    const auto add_Delta_Func = [&](int block12_I, int bsp1_I, int bsp2_I) {
        double Delta12PosNeg_F = 0.0;
        for (int block34_I = 0; block34_I < Nblock_I; ++block34_I) {
            const auto& solution = solutions[block34_I];
            for (int bsp3_I = 0; bsp3_I < Nbsp_I1D_block[block34_I]; ++bsp3_I) {
                for (int bsp4_I = 0; bsp4_I < Nbsp_I1D_block[block34_I]; ++bsp4_I) {
                    const Element element = read_element_Func(block12_I, block34_I, bsp1_I, bsp2_I, bsp3_I, bsp4_I);
                    Delta12PosNeg_F += element.vPosNegPosNeg_F * solution.kappaPosNeg_F2D_bsp_bsp(bsp3_I, bsp4_I);
                }
            }
        }
        fields[block12_I].DeltaPosNeg_F2D_bsp_bsp(bsp1_I, bsp2_I) += Delta12PosNeg_F;
    };

    // #pragma omp parallel
    {
        // #pragma omp single
        for (int block12_I = 0; block12_I < Nblock_I; ++block12_I) {
            for (int bsp1_I = 0; bsp1_I < Nbsp_I1D_block[block12_I]; ++bsp1_I) {
                for (int bsp2_I = 0; bsp2_I < Nbsp_I1D_block[block12_I]; ++bsp2_I) {
                    // #pragma omp task firstprivate(block12_I, bsp1_I, bsp2_I)
                    {
                        add_Delta_Func(block12_I, bsp1_I, bsp2_I);
                    }
                }
            }
        }
    }

}

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

        // H⁺₁₁ = h₀⁺⁺ + Γ⁺⁺ - λI; Γ⁺⁺ includes LN.
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

inline void HFBKramers::search_lambda(double lambdaTolerance_F) {
    assert(TargetN_I >= 0 && TargetN_I <= 2 * std::accumulate(Nbsp_I1D_block.begin(), Nbsp_I1D_block.end(), 0));
    assert(std::isfinite(lambda_F));
    assert(std::isfinite(lambdaTolerance_F) && lambdaTolerance_F > 0.0);
    const double Ntolerance_F = std::min(1.0e-12, lambdaTolerance_F);
    std::map<double, double> NerrorByLambda_Map{};

    // Fixed fields; trial blocking preserves external trackers.
    const auto calc_N_Func = [&](double lambdaTrial_F, bool updateTracking_B) {
        lambda_F = lambdaTrial_F;
        const double Ncalc_F = update_UV_E_rho_kappa();
        if (blocking_Func) {return blocking_Func(solutions, updateTracking_B);}
        return Ncalc_F;
    };

    // λ → N(λ)-TargetN; cache trial residuals.
    const auto calc_Nerror_Func = [&](double lambdaTrial_F) {
        const auto NerrorIt_ = NerrorByLambda_Map.find(lambdaTrial_F);
        if (NerrorIt_ != NerrorByLambda_Map.end()) {return NerrorIt_->second;}
        const double Nerror_F = calc_N_Func(lambdaTrial_F, false) - static_cast<double>(TargetN_I);
        assert(std::isfinite(Nerror_F));
        NerrorByLambda_Map.emplace(lambdaTrial_F, Nerror_F);
        return Nerror_F;
    };

    // λ_initial → N(λ_initial); commit the accepted blocking state.
    const double lambdaInitial_F = lambda_F;
    const double NerrorInitial_F = calc_Nerror_Func(lambdaInitial_F);
    if (std::abs(NerrorInitial_F) <= Ntolerance_F) {
        calc_N_Func(lambdaInitial_F, true);
        return;
    }

    // [λ_min,λ_max] = [λ_initial-2,λ_initial+2].
    double lambdaMax_F = lambdaInitial_F + 2.0;
    double lambdaMin_F = lambdaInitial_F - 2.0;
    double NerrorMin_F = calc_Nerror_Func(lambdaMin_F);
    double NerrorMax_F = calc_Nerror_Func(lambdaMax_F);
    const double NerrorSlope_F = (NerrorMax_F - NerrorMin_F) / (lambdaMax_F - lambdaMin_F);

    // f(λ_min)f(λ_max)>0 → expand one boundary.
    const bool shouldShiftLambdaMin_B = NerrorMin_F * NerrorSlope_F > 0.0;
    double& lambdaBound_F = shouldShiftLambdaMin_B ? lambdaMin_F : lambdaMax_F;
    double& NerrorBound_F = shouldShiftLambdaMin_B ? NerrorMin_F : NerrorMax_F;
    const double lambdaStep_F = shouldShiftLambdaMin_B ? -5.0 : 5.0;
    constexpr int NexpandMax_I = 100;
    int Nexpand_I = 0;
    while (NerrorMin_F * NerrorMax_F > 0.0 && NerrorSlope_F != 0.0 && Nexpand_I < NexpandMax_I) {
        lambdaBound_F += lambdaStep_F;
        NerrorBound_F = calc_Nerror_Func(lambdaBound_F);
        ++Nexpand_I;
    }
    assert(NerrorMin_F * NerrorMax_F <= 0.0);

    // Recompute at λ_root and commit blocking trackers.
    const double lambdaRoot_F = root_brent(calc_Nerror_Func, lambdaMin_F, lambdaMax_F, lambdaTolerance_F);
    assert(std::isfinite(lambdaRoot_F));
    calc_N_Func(lambdaRoot_F, true);
}
