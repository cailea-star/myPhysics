/**
 * @file    hfb_kramers_blocking.hpp
 * @author  cailea
 * @date    2026-09-12
 * @brief   Track Kramers quasiparticles for equal-filling blocking.
 */

#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>

#include <Eigen/Core>

#include "hfb_kramers.hpp"

/**
 * @brief Track one Kramers pair for equal-filling blocking.
 */
class HFBKramersBlocking {
public:
    int block_I = -1;
    int bqp_I = -1;
    double overlap_F = 0.0;

    Eigen::VectorXd blockedU_F1D_bsp{};
    Eigen::VectorXd blockedV_F1D_bsp{};

public:
    /**
     * @brief Construct an inactive blocking tracker.
     * @math B = ∅.
     * @output Inactive tracker.
     */
    HFBKramersBlocking() = default;

    /**
     * @brief Initialize tracking from a selected quasiparticle.
     * @math (block,bqp) → (U_ref,V_ref).
     * @output Stored reference amplitudes and selected indices.
     */
    HFBKramersBlocking(const std::vector<HFBKramersBlockSolution>& solutions_, int block_I_, int bqp_I_) {
        block_I = block_I_;
        bqp_I = bqp_I_;
        assert(block_I >= 0 && block_I < static_cast<int>(solutions_.size()));
        const auto& solution = solutions_[block_I];
        assert(bqp_I >= 0 && bqp_I < solution.UPos_F2D_bsp_bqp.cols());
        assert(bqp_I < solution.VNeg_F2D_bsp_bqp.cols());
        blockedU_F1D_bsp = solution.UPos_F2D_bsp_bqp.col(bqp_I);
        blockedV_F1D_bsp = solution.VNeg_F2D_bsp_bqp.col(bqp_I);
    }

    /**
     * @brief Apply equal-filling blocking using maximum-overlap tracking.
     * @math (ρ,κ) → (ρ_EFA,κ_EFA); N = 2Σ_b Trρ_b⁺⁺.
     * @output Updated densities and returned species particle number.
     * @note Commit tracker changes only when updateTracking_B is true.
     * @note Requires freshly computed zero-temperature vacuum densities.
     */
    double apply_blocking(std::vector<HFBKramersBlockSolution>& solutions, bool updateTracking_B);

    /**
     * @brief Rank low-energy blocking candidates by hole character.
     * @math E_μ>0, |E_μ−E_min|≤E_cut → {B_μ}.
     * @output At most Nblocking_I initialized candidate trackers.
     */
    static std::vector<HFBKramersBlocking> list_candidates(const std::vector<HFBKramersBlockSolution>& solutions_, int Nblocking_I, double EblockingCut_F);
};

inline double HFBKramersBlocking::apply_blocking(std::vector<HFBKramersBlockSolution>& solutions, bool updateTracking_B) {
    assert((block_I == -1 && bqp_I == -1) || (block_I >= 0 && bqp_I >= 0));
    if (block_I >= 0) {
        assert(block_I < static_cast<int>(solutions.size()));
        auto& solution = solutions[block_I];
        const int Nbsp_I = static_cast<int>(solution.eta_F1D_bsp.size());
        const int Nbqp_I = static_cast<int>(solution.Eqp_F1D_bqp.size());
        assert(blockedU_F1D_bsp.size() == Nbsp_I && blockedV_F1D_bsp.size() == Nbsp_I);
        assert(solution.UPos_F2D_bsp_bqp.rows() == Nbsp_I && solution.VNeg_F2D_bsp_bqp.rows() == Nbsp_I);
        assert(solution.UPos_F2D_bsp_bqp.cols() == Nbqp_I && solution.VNeg_F2D_bsp_bqp.cols() == Nbqp_I);

        // μ = arg max_ν Σ_i(|U_ref,i U_iν|+|V_ref,i V_iν|).
        int bestBqp_I = -1;
        double bestOverlap_F = -1.0;
        for (int bqpTrial_I = 0; bqpTrial_I < Nbqp_I; ++bqpTrial_I) {
            const double overlapCandidate_F = (blockedU_F1D_bsp.array() * solution.UPos_F2D_bsp_bqp.col(bqpTrial_I).array()).abs().sum() + (blockedV_F1D_bsp.array() * solution.VNeg_F2D_bsp_bqp.col(bqpTrial_I).array()).abs().sum();
            if (!(overlapCandidate_F > bestOverlap_F)) {continue;}
            bestOverlap_F = overlapCandidate_F;
            bestBqp_I = bqpTrial_I;
        }
        assert(bestBqp_I >= 0);
        const Eigen::VectorXd U_F1D_bsp = solution.UPos_F2D_bsp_bqp.col(bestBqp_I);
        const Eigen::VectorXd V_F1D_bsp = solution.VNeg_F2D_bsp_bqp.col(bestBqp_I);

        // δρ⁺⁺=½(uuᵀ-DvvᵀD); D=diag(η).
        solution.rhoPosPos_F2D_bsp_bsp.noalias() += 0.5 * U_F1D_bsp * U_F1D_bsp.transpose();
        solution.rhoPosPos_F2D_bsp_bsp.noalias() -= 0.5 * (solution.eta_F1D_bsp.asDiagonal() * V_F1D_bsp) * V_F1D_bsp.transpose() * solution.eta_F1D_bsp.asDiagonal();
        // δκ⁺⁻=½(uvᵀ+DvuᵀD).
        solution.kappaPosNeg_F2D_bsp_bsp.noalias() += 0.5 * U_F1D_bsp * V_F1D_bsp.transpose();
        solution.kappaPosNeg_F2D_bsp_bsp.noalias() += 0.5 * (solution.eta_F1D_bsp.asDiagonal() * V_F1D_bsp) * U_F1D_bsp.transpose() * solution.eta_F1D_bsp.asDiagonal();
        solution.rhoPosPos_F2D_bsp_bsp = 0.5 * (solution.rhoPosPos_F2D_bsp_bsp + solution.rhoPosPos_F2D_bsp_bsp.transpose()).eval();

        // Trial λ preserves the reference; accepted λ commits it.
        if (updateTracking_B) {
            blockedU_F1D_bsp = U_F1D_bsp;
            blockedV_F1D_bsp = V_F1D_bsp;
            overlap_F = bestOverlap_F;
            bqp_I = bestBqp_I;
        }
    }

    // N=2Σ_b Trρ_b⁺⁺, including unblocked blocks.
    double N_F = 0.0;
    for (const auto& solution : solutions) {N_F += 2.0 * solution.rhoPosPos_F2D_bsp_bsp.trace();}
    return N_F;
}

inline std::vector<HFBKramersBlocking> HFBKramersBlocking::list_candidates(const std::vector<HFBKramersBlockSolution>& solutions_, int Nblocking_I, double EblockingCut_F) {
    std::vector<HFBKramersBlocking> blockings_X1D_candidate{};
    if (Nblocking_I <= 0) {return blockings_X1D_candidate;}

    struct BlockingCandidate {
        int block_I = -1;
        int bqp_I = -1;
        double Eqp_F = 0.0;
        double EqpDifference_F = 0.0;
        double Vnorm2_F = 0.0;
    };

    // E_μ>0, ||V_μ||²>10⁻⁴ → E_{min}.
    double EqpMin_F = std::numeric_limits<double>::infinity();
    for (const HFBKramersBlockSolution& block_ : solutions_) {
        const int Nbqp_I = static_cast<int>(block_.Eqp_F1D_bqp.size());
        for (int bqp_I = 0; bqp_I < Nbqp_I; ++bqp_I) {
            const double Eqp_F = block_.Eqp_F1D_bqp(bqp_I);
            const double Vnorm2_F = block_.VNeg_F2D_bsp_bqp.col(bqp_I).squaredNorm();
            if (Eqp_F <= 0.0 || Vnorm2_F <= 1.0e-4) {continue;}
            EqpMin_F = std::min(EqpMin_F, Eqp_F);
        }
    }
    if (!std::isfinite(EqpMin_F)) {return blockings_X1D_candidate;}

    // |E_μ-E_{min}|≤E_{cut} → candidates.
    std::vector<BlockingCandidate> candidates_X1D_candidate{};
    for (int block_I = 0; block_I < static_cast<int>(solutions_.size()); ++block_I) {
        const HFBKramersBlockSolution& block_ = solutions_[block_I];
        const int Nbqp_I = static_cast<int>(block_.Eqp_F1D_bqp.size());
        for (int bqp_I = 0; bqp_I < Nbqp_I; ++bqp_I) {
            const double Eqp_F = block_.Eqp_F1D_bqp(bqp_I);
            const double Vnorm2_F = block_.VNeg_F2D_bsp_bqp.col(bqp_I).squaredNorm();
            const double EqpDifference_F = std::abs(Eqp_F - EqpMin_F);
            const bool isCandidate_B = Eqp_F > 0.0 && EqpDifference_F <= EblockingCut_F;
            if (isCandidate_B) {candidates_X1D_candidate.push_back({block_I, bqp_I, Eqp_F, EqpDifference_F, Vnorm2_F});}
        }
    }
    if (candidates_X1D_candidate.empty()) {return blockings_X1D_candidate;}

    // (|E-E_{min}|,hole,block,bqp) → ascending rank.
    std::sort(candidates_X1D_candidate.begin(), candidates_X1D_candidate.end(), [](const BlockingCandidate& candidateL_, const BlockingCandidate& candidateR_) {
        constexpr double EqpDifferenceTolerance_F = 1.0e-6;
        if (std::abs(candidateL_.EqpDifference_F - candidateR_.EqpDifference_F) > EqpDifferenceTolerance_F) {return candidateL_.EqpDifference_F < candidateR_.EqpDifference_F;}
        const bool isHoleL_B = candidateL_.Vnorm2_F > 0.5;
        const bool isHoleR_B = candidateR_.Vnorm2_F > 0.5;
        if (isHoleL_B != isHoleR_B) {return isHoleL_B;}
        if (candidateL_.block_I != candidateR_.block_I) {return isHoleL_B ? candidateL_.block_I < candidateR_.block_I : candidateL_.block_I > candidateR_.block_I;}
        return candidateL_.bqp_I < candidateR_.bqp_I;
    });

    if (static_cast<int>(candidates_X1D_candidate.size()) > Nblocking_I) {candidates_X1D_candidate.resize(Nblocking_I);}

    // {B_μ} → state trackers.
    blockings_X1D_candidate.reserve(candidates_X1D_candidate.size());
    for (const BlockingCandidate& candidate_ : candidates_X1D_candidate) {
        HFBKramersBlocking blocking_(solutions_, candidate_.block_I, candidate_.bqp_I);
        blocking_.overlap_F = 1.0;
        blockings_X1D_candidate.push_back(std::move(blocking_));
    }
    return blockings_X1D_candidate;
}
