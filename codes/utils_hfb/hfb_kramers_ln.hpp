/**
 * @file    hfb_kramers_ln.hpp
 * @author  cailea
 * @date    2026-09-12
 * @brief   Kramers Lipkin-Nogami field and energy corrections.
 */

#pragma once

#include "hfb_kramers.hpp"

/**
 * @brief Calculate effective-seniority Lipkin-Nogami corrections.
 * @math λ₂=−G_eff S_N/S_D; E_LN=−4λ₂Σu²v²; Γ+=4λ₂ρ−2λ₂I.
 * @output Updated lambda2, ELipkinNogami, and Gamma matrices.
 * @note Zero-temperature canonical-pair approximation; rebuild Gamma before calling.
 */
inline void add_Gamma_from_lipkin_nogami(HFBKramers& hfb_);

inline void add_Gamma_from_lipkin_nogami(HFBKramers& hfb_) {
    double Su1v3_F = 0.0;
    double Su2v2_F = 0.0;
    double Su3v1_F = 0.0;
    double Su4v4_F = 0.0;
    double Epair_F = 0.0;
    double DeltaRho_F = 0.0;
    double Nparticle_F = 0.0;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver_;

    // Kramers pairs: N=2Trρ⁺⁺; E_pair=Σ_ab Δ⁺⁻_ab κ⁺⁻_ab.
    for (int block_I = 0; block_I < hfb_.Nblock_I; ++block_I) {
        const auto& solution = hfb_.solutions[block_I];
        const auto& field = hfb_.fields[block_I];
        Nparticle_F += 2.0 * solution.rhoPosPos_F2D_bsp_bsp.trace();
        Epair_F += field.DeltaPosNeg_F2D_bsp_bsp.cwiseProduct(solution.kappaPosNeg_F2D_bsp_bsp).sum();
        // Δρ_weighted=2Σ_ab η_b Δ⁺⁻_ab ρ⁺⁺_ab.
        DeltaRho_F += 2.0 * (field.DeltaPosNeg_F2D_bsp_bsp * solution.eta_F1D_bsp.asDiagonal()).cwiseProduct(solution.rhoPosPos_F2D_bsp_bsp).sum();
        eigensolver_.compute(solution.rhoPosPos_F2D_bsp_bsp, Eigen::EigenvaluesOnly);
        assert(eigensolver_.info() == Eigen::Success);

        // Σ_pairs=Σ_{i∈positive branch}; v_i²=eig(ρ⁺⁺)_i.
        for (int bsp_I = 0; bsp_I < hfb_.Nbsp_I1D_block[block_I]; ++bsp_I) {
            const double v2_F = std::clamp(eigensolver_.eigenvalues()(bsp_I), 0.0, 1.0);
            const double u2_F = 1.0 - v2_F;
            const double u_F = std::sqrt(u2_F);
            const double v_F = std::sqrt(v2_F);
            const double u3_F = u2_F * u_F;
            const double v3_F = v2_F * v_F;
            const double u4v4_F = u2_F * u2_F * v2_F * v2_F;
            Su1v3_F += u_F * v3_F;
            Su2v2_F += u2_F * v2_F;
            Su3v1_F += u3_F * v_F;
            Su4v4_F += u4v4_F;
        }
    }

    // Δ̄=−Δρ_weighted/N; G_eff=Δ̄²/E_pair.
    const double lambda2_numer_F = 8.0 * (Su3v1_F * Su1v3_F - Su4v4_F);
    const double lambda2_denom_F = 32.0 * (Su2v2_F * Su2v2_F - Su4v4_F);
    hfb_.lambda2_F = 0.0;
    if (Nparticle_F > 0.0 && Epair_F < 0.0 && lambda2_denom_F > 0.0) {
        const double DeltaAverage_F = -DeltaRho_F / Nparticle_F;
        const double Geff_F = DeltaAverage_F * DeltaAverage_F / Epair_F;
        hfb_.lambda2_F = -Geff_F * lambda2_numer_F / lambda2_denom_F;
    }
    if (!std::isfinite(hfb_.lambda2_F)) {hfb_.lambda2_F = 0.0;}
    if (hfb_.lambda2_F >= 10.0) {hfb_.lambda2_F = 4.0;}
    hfb_.ELipkinNogami_F = -4.0 * hfb_.lambda2_F * Su2v2_F;
    if (!std::isfinite(hfb_.ELipkinNogami_F)) {hfb_.ELipkinNogami_F = 0.0;}

    // Γ⁺⁺ += 4λ₂ρ⁺⁺−2λ₂I.
    #pragma omp parallel for schedule(static)
    for (int block_I = 0; block_I < hfb_.Nblock_I; ++block_I) {
        hfb_.fields[block_I].GammaPosPos_F2D_bsp_bsp += 4.0 * hfb_.lambda2_F * hfb_.solutions[block_I].rhoPosPos_F2D_bsp_bsp;
        hfb_.fields[block_I].GammaPosPos_F2D_bsp_bsp.diagonal().array() -= 2.0 * hfb_.lambda2_F;
    }
}
