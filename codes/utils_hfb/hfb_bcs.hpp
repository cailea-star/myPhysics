/**
 * @file    hfb_bcs.hpp
 * @author  cailea
 * @date    2026-09-17
 * @brief   Single-particle diagonalization and BCS pairing.
 */

#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <limits>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Eigenvalues>

#include "root.hpp"

struct HFBCSBlockSolution {
    // T|α⟩ = ηα|ᾱ⟩.
    Eigen::VectorXd eta_F1D_bsp{};

    // h f = f diag(Ehf).
    Eigen::VectorXd Ehf_F1D_bhf{};
    Eigen::MatrixXd fhf_F2D_sp_bhf{};

    // bhf ↔ bqp.
    Eigen::VectorXd uPos_F1D_bhf{};
    Eigen::VectorXd vNeg_F1D_bhf{};
    Eigen::VectorXd Eqp_F1D_bhf{};
    Eigen::VectorXd rhoPosPos_F1D_bhf{};
    Eigen::VectorXd kappaPosNeg_F1D_bhf{};

    Eigen::MatrixXd rhoPosPos_F2D_bsp_bsp{};
    Eigen::MatrixXd kappaPosNeg_F2D_bsp_bsp{};
};

struct HFBCSBlockField {
    // h = h₀ + Γ.
    Eigen::MatrixXd h0PosPos_F2D_bsp_bsp{};
    Eigen::MatrixXd GammaPosPos_F2D_bsp_bsp{};
};

class HFBCS {
public:
    struct Element {
        double vPosPosPosPos_F = 0.0;
        double vPosNegPosNeg_F = 0.0;
    };

    // v̄_{12;34}: Γ uses block13, block24.
    using GammaElementFunc = std::function<Element(int block13_I, int block24_I, int bsp1_I, int bsp2_I, int bsp3_I, int bsp4_I)>;

    double lambda_F = -7.0; // Chemical potential [MeV].
    double G_F = 0.0; // Attractive pairing strength [MeV].
    double Delta_F = 1.0; // Pairing gap [MeV].

    int Nblock_I = 0;
    std::vector<int> Nbsp_I1D_block{};

    std::vector<HFBCSBlockSolution> solutions{};
    std::vector<HFBCSBlockField> fields{};

    std::vector<Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>> hPos_eigensolvers{};

public:
    /**
     * @brief Allocate single-species Kramers-representative blocks.
     * @math Nbsp_b = size(η_b); ηα² = 1.
     * @output Stored phases; allocated fields, solutions, workspaces.
     */
    HFBCS(const std::vector<Eigen::VectorXd>& eta_F2D_block_bsp_) {
        Nblock_I = static_cast<int>(eta_F2D_block_bsp_.size());

        // {η_b} → {field_b, solution_b, workspace_b}.
        Nbsp_I1D_block.resize(Nblock_I);
        solutions.resize(Nblock_I);
        fields.resize(Nblock_I);
        hPos_eigensolvers.resize(Nblock_I);

        for (int block_I = 0; block_I < Nblock_I; ++block_I) {
            const int Nbsp_I = static_cast<int>(eta_F2D_block_bsp_[block_I].size());
            assert(Nbsp_I > 0);
            assert((eta_F2D_block_bsp_[block_I].array().abs() == 1.0).all());
            Nbsp_I1D_block[block_I] = Nbsp_I;
            HFBCSBlockField& field = fields[block_I];
            HFBCSBlockSolution& solution = solutions[block_I];
            solution.eta_F1D_bsp = eta_F2D_block_bsp_[block_I];

            // h₀, Γ, f ∈ ℝ^{Nbsp×Nbsp}.
            field.h0PosPos_F2D_bsp_bsp.resize(Nbsp_I, Nbsp_I);
            field.GammaPosPos_F2D_bsp_bsp.resize(Nbsp_I, Nbsp_I);
            solution.fhf_F2D_sp_bhf.resize(Nbsp_I, Nbsp_I);

            // Ehf, u, v, Eqp, ρ⁺⁺, κ⁺⁻ ∈ ℝ^{Nbsp}.
            solution.Ehf_F1D_bhf.resize(Nbsp_I);
            solution.uPos_F1D_bhf.resize(Nbsp_I);
            solution.vNeg_F1D_bhf.resize(Nbsp_I);
            solution.Eqp_F1D_bhf.resize(Nbsp_I);
            solution.rhoPosPos_F1D_bhf.resize(Nbsp_I);
            solution.kappaPosNeg_F1D_bhf.resize(Nbsp_I);

            // ρ⁺⁺, κ⁺⁻ ∈ ℝ^{Nbsp×Nbsp}.
            solution.rhoPosPos_F2D_bsp_bsp.resize(Nbsp_I, Nbsp_I);
            solution.kappaPosNeg_F2D_bsp_bsp.resize(Nbsp_I, Nbsp_I);

            // h⁺⁺ ∈ ℝ^{Nbsp×Nbsp}.
            hPos_eigensolvers[block_I] = Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>(Nbsp_I);
        }
    }

    /**
     * @brief Evaluate zero-temperature BCS amplitudes and densities.
     * @math E_a = √[(Ehf_a−λ)²+Δ_a²]; Δ_a = Δ inside W, otherwise 0.
     * @output Stored λ,Δ; updated u,v,E,ρ,κ.
     * @output Mean particle number N = 2Σ_b Tr(ρ_b⁺⁺).
     * @note Requires current single-particle eigenvalues and eigenvectors.
     * @note EspCut restricts pairing; particle counting includes all states.
     * @note W = {a: Ehf_a ≤ EspCut}; eigenstate partners use η_a=1.
     */
    double update_UV_E_rho_kappa(double lambda_F_, double Delta_F_, double EspCut_F);

    /**
     * @brief Search zero-temperature pairing gap using Brent bracketing.
     * @math 1 = G Σ_{i∈W} 1/(2√[(Ehf_i−λ)²+Δ²]).
     * @output Updated Delta_F; zero when pairing collapses.
     * @note Requires current single-particle eigenvalues; fixed λ.
     * @note Gaps below accuracy_F are returned as zero.
     */
    void search_Delta(double lambda_F_, double EspCut_F, double accuracy_F);

    /**
     * @brief Search chemical potential using spectral bracketing and Brent.
     * @math N(λ,Δ(λ)) = TargetN.
     * @output Updated chemical potential, gap, solutions, densities.
     * @note Diagonalizes supplied h once before searching λ.
     * @note Each trial λ solves Δ before updating densities.
     */
    void search_lambda(int TargetN_I, double EspCut_F, double accuracy_F);

    /**
     * @brief Accumulate particle-hole fields by direct matrix-element contraction.
     * @math Γ_{13} += Σ_{24}v̄_{12;34}ρ_42.
     * @output Accumulated Gamma from this species' density.
     * @note Requires time-reversal invariance and thread-safe callbacks.
     */
    void add_Gamma_from_Element(const GammaElementFunc& read_element_Func);
};

/**
 * @brief Evaluate analytic zero-temperature BCS amplitudes and densities.
 * @math ρ_hf = v²; κ_hf = −uv; N = 2Σρ_hf.
 * @output Stored λ,Δ; updated E,u,v,ρ,κ; returned N.
 * @note Requires strictly positive quasiparticle energies.
 */
inline double HFBCS::update_UV_E_rho_kappa(double lambda_F_, double Delta_F_, double EspCut_F) {
    assert(std::isfinite(lambda_F_));
    assert(std::isfinite(Delta_F_) && Delta_F_ >= 0.0);
    assert(!std::isnan(EspCut_F));
    lambda_F = lambda_F_;
    Delta_F = Delta_F_;
    double N_F = 0.0;

    for (int block_I = 0; block_I < Nblock_I; ++block_I) {
        HFBCSBlockSolution& solution = solutions[block_I];
        const int Nbsp_I = Nbsp_I1D_block[block_I];
        assert(solution.Ehf_F1D_bhf.size() == Nbsp_I && solution.Ehf_F1D_bhf.allFinite());
        assert(solution.fhf_F2D_sp_bhf.rows() == Nbsp_I && solution.fhf_F2D_sp_bhf.cols() == Nbsp_I && solution.fhf_F2D_sp_bhf.allFinite());

        // Δ_a = Δ [ε_a ≤ Ecut]; E_a = hypot(ε_a−λ,Δ_a).
        for (int bhf_I = 0; bhf_I < Nbsp_I; ++bhf_I) {
            const double Ehf_F = solution.Ehf_F1D_bhf(bhf_I);
            const double DeltaCuted_F = Delta_F * static_cast<double>(Ehf_F <= EspCut_F);
            const double Eqp_F = std::hypot(Ehf_F - lambda_F, DeltaCuted_F);
            assert(Eqp_F > 0.0);
            const double ratio_F = (Ehf_F - lambda_F) / Eqp_F;
            solution.Eqp_F1D_bhf(bhf_I) = Eqp_F;
            solution.uPos_F1D_bhf(bhf_I) = std::sqrt(0.5 * (1.0 + ratio_F));
            solution.vNeg_F1D_bhf(bhf_I) = std::sqrt(0.5 * (1.0 - ratio_F));
            // ρ_i⁺⁺ = v_i²; κ_i⁺⁻ = −u_i v_i.
            solution.rhoPosPos_F1D_bhf(bhf_I) = solution.vNeg_F1D_bhf(bhf_I) * solution.vNeg_F1D_bhf(bhf_I);
            solution.kappaPosNeg_F1D_bhf(bhf_I) = -solution.uPos_F1D_bhf(bhf_I) * solution.vNeg_F1D_bhf(bhf_I);
        }

        // ρ⁺⁺ = f diag(ρ_hf⁺⁺) fᵀ.
        solution.rhoPosPos_F2D_bsp_bsp.noalias() = solution.fhf_F2D_sp_bhf * solution.rhoPosPos_F1D_bhf.asDiagonal() * solution.fhf_F2D_sp_bhf.transpose();

        // κ⁺⁻ = f diag(κ_hf⁺⁻) fᵀ diag(η).
        solution.kappaPosNeg_F2D_bsp_bsp.noalias() = solution.fhf_F2D_sp_bhf * solution.kappaPosNeg_F1D_bhf.asDiagonal() * solution.fhf_F2D_sp_bhf.transpose() * solution.eta_F1D_bsp.asDiagonal();

        // N = 2Σ_b Σ_i ρ_bi⁺⁺.
        N_F += 2.0 * solution.rhoPosPos_F1D_bhf.sum();
    }
    return N_F;
}

/**
 * @brief Solve the zero-temperature gap equation using Brent.
 * @math F(Δ) = 1 − (G/2)Σ_{i∈W}1/hypot(Ehf_i−λ,Δ).
 * @output Updated Delta_F; zero below accuracy_F.
 * @note W = {i: Ehf_i ≤ EspCut}; fixed λ.
 */
inline void HFBCS::search_Delta(double lambda_F_, double EspCut_F, double accuracy_F) {
    assert(std::isfinite(lambda_F_));
    assert(std::isfinite(G_F) && G_F >= 0.0);
    assert(!std::isnan(EspCut_F));
    assert(std::isfinite(accuracy_F) && accuracy_F > 0.0);

    // Npair = Σ_b Σ_i [Ehf_bi ≤ Ecut].
    int Npair_I = 0;
    for (int block_I = 0; block_I < Nblock_I; ++block_I) {
        const HFBCSBlockSolution& solution = solutions[block_I];
        assert(solution.Ehf_F1D_bhf.size() == Nbsp_I1D_block[block_I] && solution.Ehf_F1D_bhf.allFinite());
        for (int bhf_I = 0; bhf_I < Nbsp_I1D_block[block_I]; ++bhf_I) {
            Npair_I += static_cast<int>(solution.Ehf_F1D_bhf(bhf_I) <= EspCut_F);
        }
    }

    // Δ > 0 ⇒ every denominator is positive.
    const auto calc_gap_Func = [&](double DeltaTrial_F) {
        double sum_F = 0.0;
        for (int block_I = 0; block_I < Nblock_I; ++block_I) {
            const HFBCSBlockSolution& solution = solutions[block_I];
            for (int bhf_I = 0; bhf_I < Nbsp_I1D_block[block_I]; ++bhf_I) {
                const double Eqp_F = std::hypot(solution.Ehf_F1D_bhf(bhf_I) - lambda_F_, DeltaTrial_F);
                sum_F += static_cast<double>(solution.Ehf_F1D_bhf(bhf_I) <= EspCut_F) * (0.5 * G_F / Eqp_F);
            }
        }
        return 1.0 - sum_F;
    };

    // F increases; F(accuracy) ≥ 0 ⇒ Δ ≤ accuracy.
    Delta_F = 0.0;
    if (calc_gap_Func(accuracy_F) >= 0.0) {return;}

    // Δmax ≥ G Npair ⇒ F(Δmax) ≥ 1/2.
    const double DeltaMax_F = std::max(G_F * Npair_I, 2.0 * accuracy_F);
    assert(std::isfinite(DeltaMax_F));
    Delta_F = root_brent(calc_gap_Func, accuracy_F, DeltaMax_F, accuracy_F);
    assert(std::isfinite(Delta_F) && Delta_F > 0.0);
}

/**
 * @brief Diagonalize h and solve particle number using Brent.
 * @math h f = f diag(Ehf); N(λ,Δ(λ)) = TargetN.
 * @output Updated λ,Δ and all block solutions.
 * @note Requires supplied real symmetric single-particle fields.
 */
inline void HFBCS::search_lambda(int TargetN_I, double EspCut_F, double accuracy_F) {
    assert(Nblock_I > 0);
    assert(std::isfinite(G_F) && G_F >= 0.0);
    assert(!std::isnan(EspCut_F));
    assert(std::isfinite(accuracy_F) && accuracy_F > 0.0);

    // Nsp counts one representative per Kramers pair.
    int Nsp_I = 0;
    double EhfMin_F = std::numeric_limits<double>::infinity();
    double EhfMax_F = -std::numeric_limits<double>::infinity();
    for (int block_I = 0; block_I < Nblock_I; ++block_I) {
        HFBCSBlockSolution& solution = solutions[block_I];
        assert(Nbsp_I1D_block[block_I] > 0);
        assert(fields[block_I].h0PosPos_F2D_bsp_bsp.rows() == Nbsp_I1D_block[block_I] && fields[block_I].h0PosPos_F2D_bsp_bsp.cols() == Nbsp_I1D_block[block_I]);
        assert(fields[block_I].GammaPosPos_F2D_bsp_bsp.rows() == Nbsp_I1D_block[block_I] && fields[block_I].GammaPosPos_F2D_bsp_bsp.cols() == Nbsp_I1D_block[block_I]);
        assert(fields[block_I].h0PosPos_F2D_bsp_bsp.allFinite() && fields[block_I].GammaPosPos_F2D_bsp_bsp.allFinite());

        // (h₀+Γ) f = f diag(Ehf).
        hPos_eigensolvers[block_I].compute(fields[block_I].h0PosPos_F2D_bsp_bsp + fields[block_I].GammaPosPos_F2D_bsp_bsp);
        assert(hPos_eigensolvers[block_I].info() == Eigen::Success);
        solution.Ehf_F1D_bhf = hPos_eigensolvers[block_I].eigenvalues();
        solution.fhf_F2D_sp_bhf = hPos_eigensolvers[block_I].eigenvectors();

        Nsp_I += Nbsp_I1D_block[block_I];
        EhfMin_F = std::min(EhfMin_F, solution.Ehf_F1D_bhf.minCoeff());
        EhfMax_F = std::max(EhfMax_F, solution.Ehf_F1D_bhf.maxCoeff());
    }
    assert(TargetN_I >= 0 && TargetN_I <= 2 * Nsp_I);

    // |ε−λ| > G Nsp ⇒ Δ=0; endpoints have N=0,2Nsp.
    const double lambdaMargin_F = 1.0 + G_F * Nsp_I;
    const double lambdaMin_F = EhfMin_F - lambdaMargin_F;
    const double lambdaMax_F = EhfMax_F + lambdaMargin_F;
    assert(std::isfinite(lambdaMin_F) && std::isfinite(lambdaMax_F) && lambdaMin_F < lambdaMax_F);

    // λ → Δ(λ) → ρ,κ → N−TargetN.
    const auto calc_Nerror_Func = [&](double lambdaTrial_F) {
        search_Delta(lambdaTrial_F, EspCut_F, accuracy_F);
        const double Nerror_F = update_UV_E_rho_kappa(lambdaTrial_F, Delta_F, EspCut_F) - static_cast<double>(TargetN_I);
        assert(std::isfinite(Nerror_F));
        return Nerror_F;
    };

    // Empty/full spaces accept a finite endpoint chemical potential.
    if (std::abs(calc_Nerror_Func(lambdaMin_F)) <= accuracy_F) {return;}
    if (std::abs(calc_Nerror_Func(lambdaMax_F)) <= accuracy_F) {return;}

    // Restore all solution caches at the accepted root.
    const double lambdaRoot_F = root_brent(calc_Nerror_Func, lambdaMin_F, lambdaMax_F, accuracy_F);
    assert(std::isfinite(lambdaRoot_F));
    const double Nerror_F = calc_Nerror_Func(lambdaRoot_F);
    assert(std::abs(Nerror_F) <= accuracy_F);
}

/**
 * @brief Accumulate particle-hole fields using parallel direct contraction.
 * @math Γ⁺⁺_{13} += Σ_{24}(v̄⁺⁺⁺⁺ρ⁺⁺_{42}+v̄⁺⁻⁺⁻ρ⁻⁻_{42}).
 * @output Accumulated Gamma from this species' density.
 * @note Requires time-reversal invariance and thread-safe callbacks.
 */
inline void HFBCS::add_Gamma_from_Element(const GammaElementFunc& read_element_Func) {
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

    #pragma omp parallel
    {
        #pragma omp single
        for (int block13_I = 0; block13_I < Nblock_I; ++block13_I) {
            for (int bsp1_I = 0; bsp1_I < Nbsp_I1D_block[block13_I]; ++bsp1_I) {
                for (int bsp3_I = 0; bsp3_I < Nbsp_I1D_block[block13_I]; ++bsp3_I) {
                    #pragma omp task firstprivate(block13_I, bsp1_I, bsp3_I)
                    {
                        add_Gamma_Func(block13_I, bsp1_I, bsp3_I);
                    }
                }
            }
        }
    }

}
