/**
 * @file    axial_hfb_blocking.cpp
 * @author  cailea
 * @date    2026-05-06
 * @brief   Implement axial HFB blocking.
 */

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <utility>
#include <vector>

#include "axial_hfb.hpp"

/**
 * @brief  Apply equal-filling blocking by maximum overlap.
 * @math   (ρ,κ,U_μ,V_μ) → (ρ',κ')
 * @output Updated blocked densities and state tracker.
 * @note   Requires zero-temperature HFB density factors.
 */
void AxialHFBBlocking::apply_blocking(AxialHFBBlockList& blocklist_) {
    assert(block_I >= 0 && block_I < static_cast<int>(blocklist_.blocks_S1D_block.size()));
    AxialHFBBlock& block_ = blocklist_.blocks_S1D_block[block_I];
    const int Nbsp_I = static_cast<int>(block_.labels_S1D_bsp.size());
    const int Nbqp_I = static_cast<int>(block_.Eqp_F1D_bqp.size());
    assert(blockedV_F1D_bsp.size() == Nbsp_I && blockedU_F1D_bsp.size() == Nbsp_I);
    assert(block_.U_F2D_bsp_bqp.cols() == Nbqp_I && block_.V_F2D_bsp_bqp.cols() == Nbqp_I);

    // μ = arg max_ν(|U_μU_ν|₁+|V_μV_ν|₁).
    int bestBqp_I = -1;
    double bestOverlap_F = -1.0;
    for (int bqp_I = 0; bqp_I < Nbqp_I; ++bqp_I) {
        const double overlapCandidate_F = (blockedU_F1D_bsp.array() * block_.U_F2D_bsp_bqp.col(bqp_I).array()).abs().sum() + (blockedV_F1D_bsp.array() * block_.V_F2D_bsp_bqp.col(bqp_I).array()).abs().sum();
        if (overlapCandidate_F > bestOverlap_F) {
            bestOverlap_F = overlapCandidate_F;
            bestBqp_I = bqp_I;
        }
    }
    assert(bestBqp_I >= 0);

    const Eigen::VectorXd U_F1D_bsp = block_.U_F2D_bsp_bqp.col(bestBqp_I);
    const Eigen::VectorXd V_F1D_bsp = block_.V_F2D_bsp_bqp.col(bestBqp_I);

    // ρ' = ρ+(UUᵀ-VVᵀ)/2.
    block_.rho_F2D_bsp_bsp.noalias() += 0.5 * U_F1D_bsp * U_F1D_bsp.transpose();
    block_.rho_F2D_bsp_bsp.noalias() -= 0.5 * V_F1D_bsp * V_F1D_bsp.transpose();

    // κ' = κ+(VUᵀ+UVᵀ)/2.
    block_.kappa_F2D_bsp_bsp.noalias() += 0.5 * V_F1D_bsp * U_F1D_bsp.transpose();
    block_.kappa_F2D_bsp_bsp.noalias() += 0.5 * U_F1D_bsp * V_F1D_bsp.transpose();

    // (ρ,κ) → ([ρ+ρᵀ]/2,[κ+κᵀ]/2).
    block_.rho_F2D_bsp_bsp = 0.5 * (block_.rho_F2D_bsp_bsp + block_.rho_F2D_bsp_bsp.transpose()).eval();
    block_.kappa_F2D_bsp_bsp = 0.5 * (block_.kappa_F2D_bsp_bsp + block_.kappa_F2D_bsp_bsp.transpose()).eval();

    blockedU_F1D_bsp = U_F1D_bsp;
    blockedV_F1D_bsp = V_F1D_bsp;
    overlap_F = bestOverlap_F;
    bqp_I = bestBqp_I;
}

/**
 * @brief  Rank blocking candidates by energy and hole character.
 * @math   E_μ>0, |E_μ-E_{min}|≤E_{cut} → {B_μ}
 * @output Ranked blocking candidates.
 */
std::vector<AxialHFBBlocking> AxialHFBBlocking::list_candidates(const HFBSettings& hfbsettings_, const AxialHFBBlockList& blocklist_, bool isNeutron_B_) {
    const int Nblocking_I = hfbsettings_.Nblocking_I;
    const double EblockingCut_F = hfbsettings_.EblockingCut_F;
    std::vector<AxialHFBBlocking> blockings_S1D_candidate;
    if (Nblocking_I <= 0) {return blockings_S1D_candidate;}

    struct BlockingCandidate {
        int block_I;
        int bqp_I;
        double Eqp_F;
        double EqpDifference_F;
        double Vnorm2_F;
    };

    // E_μ>0, ||V_μ||²>10⁻⁴ → E_{min}.
    const auto defines_EqpMin_Func = [](const AxialHFBBlock& block_, int bqp_I) {
        const double Eqp_F = block_.Eqp_F1D_bqp(bqp_I);
        const double Vnorm2_F = block_.V_F2D_bsp_bqp.col(bqp_I).squaredNorm();
        return Eqp_F > 0.0 && Vnorm2_F > 1.0e-4;
    };

    double EqpMin_F = std::numeric_limits<double>::infinity();
    for (const AxialHFBBlock& block_ : blocklist_.blocks_S1D_block) {
        const int Nbqp_I = static_cast<int>(block_.Eqp_F1D_bqp.size());
        for (int bqp_I = 0; bqp_I < Nbqp_I; ++bqp_I) {
            const double Eqp_F = block_.Eqp_F1D_bqp(bqp_I);
            if (defines_EqpMin_Func(block_, bqp_I) && Eqp_F < EqpMin_F) {EqpMin_F = Eqp_F;}
        }
    }
    if (!std::isfinite(EqpMin_F)) {return blockings_S1D_candidate;}

    // |E_μ-E_{min}|≤E_{cut} → candidates.
    std::vector<BlockingCandidate> candidates_S1D_candidate;
    for (int block_I = 0; block_I < static_cast<int>(blocklist_.blocks_S1D_block.size()); ++block_I) {
        const AxialHFBBlock& block_ = blocklist_.blocks_S1D_block[block_I];
        const int Nbqp_I = static_cast<int>(block_.Eqp_F1D_bqp.size());
        for (int bqp_I = 0; bqp_I < Nbqp_I; ++bqp_I) {
            const double Eqp_F = block_.Eqp_F1D_bqp(bqp_I);
            const double Vnorm2_F = block_.V_F2D_bsp_bqp.col(bqp_I).squaredNorm();
            const double EqpDifference_F = std::abs(Eqp_F - EqpMin_F);
            const bool isCandidate_B = Eqp_F > 0.0 && EqpDifference_F <= EblockingCut_F;
            if (isCandidate_B) {candidates_S1D_candidate.push_back({block_I, bqp_I, Eqp_F, EqpDifference_F, Vnorm2_F});}
        }
    }
    if (candidates_S1D_candidate.empty()) {return blockings_S1D_candidate;}

    // (|E-E_{min}|,hole,block,bqp) → ascending rank.
    std::sort(candidates_S1D_candidate.begin(), candidates_S1D_candidate.end(), [](const BlockingCandidate& candidateL_, const BlockingCandidate& candidateR_) {
        constexpr double EqpDifferenceTolerance_F = 1.0e-6;
        if (std::abs(candidateL_.EqpDifference_F - candidateR_.EqpDifference_F) > EqpDifferenceTolerance_F) {return candidateL_.EqpDifference_F < candidateR_.EqpDifference_F;}
        const bool isHoleL_B = candidateL_.Vnorm2_F > 0.5;
        const bool isHoleR_B = candidateR_.Vnorm2_F > 0.5;
        if (isHoleL_B != isHoleR_B) {return isHoleL_B;}
        if (candidateL_.block_I != candidateR_.block_I) {return isHoleL_B ? candidateL_.block_I < candidateR_.block_I : candidateL_.block_I > candidateR_.block_I;}
        return candidateL_.bqp_I < candidateR_.bqp_I;
    });

    if (static_cast<int>(candidates_S1D_candidate.size()) > Nblocking_I) {candidates_S1D_candidate.resize(Nblocking_I);}

    // {B_μ} → state trackers.
    blockings_S1D_candidate.reserve(candidates_S1D_candidate.size());
    std::cout << "[AxialHFBBlocking] Blocking candidates: eqpmin=" << EqpMin_F << "MeV, window=" << EblockingCut_F << "MeV, count=" << candidates_S1D_candidate.size() << ".\n";
    for (const BlockingCandidate& candidate_ : candidates_S1D_candidate) {
        const AxialHFBBlock& block_ = blocklist_.blocks_S1D_block[candidate_.block_I];
        const Eigen::VectorXd U_F1D_bsp = block_.U_F2D_bsp_bqp.col(candidate_.bqp_I);
        const Eigen::VectorXd V_F1D_bsp = block_.V_F2D_bsp_bqp.col(candidate_.bqp_I);
        Eigen::Index bspMax_I = 0;
        U_F1D_bsp.cwiseAbs().cwiseMax(V_F1D_bsp.cwiseAbs()).maxCoeff(&bspMax_I);
        const AxialSPLabel& label_ = block_.labels_S1D_bsp[bspMax_I];
        std::cout << "  num=" << std::setw(3) << blockings_S1D_candidate.size() + 1 << " block=" << std::setw(3) << candidate_.block_I + 1 << " qp=" << std::setw(3) << candidate_.bqp_I + 1 << " Eqp=" << std::setw(12) << std::setprecision(6) << candidate_.Eqp_F << " diff=" << std::setw(12) << candidate_.EqpDifference_F << " v2=" << std::setw(10) << candidate_.Vnorm2_F << " label=" << label_.twoOmega_I << "/2" << (label_.isParityPositive_B ? "+" : "-") << "[" << label_.N_I << "," << label_.nz_I << "," << label_.Lambda_I << "]\n";
        AxialHFBBlocking blocking_(blocklist_, isNeutron_B_, candidate_.block_I, candidate_.bqp_I);
        blocking_.overlap_F = 1.0;
        blockings_S1D_candidate.push_back(std::move(blocking_));
    }
    return blockings_S1D_candidate;
}
