/**
 * @file    spherical_setting.hpp
 * @author  cailea
 * @date    2026-09-14
 * @brief   Define spherical harmonic-oscillator labels and settings.
 */

#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <utility>
#include <vector>

class SphericalSPLabel {
public:
    int N_I = 0;                     // N = 2n + l
    int n_I = 0;                     // n ≥ 0
    int l_I = 0;                     // l ≥ 0
    int twoj_I = 0;                  // 2j = 2l ± 1 > 0
    int twom_I = 0;                  // 2m = -2j, -2j+2, ..., 2j
    bool isParityPositive_B = false; // π = (-1)^N = +1

    /**
     * @brief  Construct a spherical single-particle label.
     * @math   l = j-1/2 + [N-(j-1/2)] mod 2, n = (N-l)/2
     * @output Initialized label.
     */
    SphericalSPLabel(int N_I_, int twoj_I_, int twom_I_) {
        assert(N_I_ >= 0);
        assert((twoj_I_ - 1) / 2 <= N_I_);
        assert(twoj_I_ > 0 && twoj_I_ % 2 == 1);
        assert(twom_I_ >= -twoj_I_ && twom_I_ <= twoj_I_);
        assert((twoj_I_ - twom_I_) % 2 == 0);

        // (N,2j,2m) → (N,n,l,2j,2m,π).
        N_I = N_I_;
        twoj_I = twoj_I_;
        twom_I = twom_I_;
        isParityPositive_B = N_I % 2 == 0;
        l_I = (twoj_I - 1) / 2;
        l_I += static_cast<int>((l_I % 2 == 0) != isParityPositive_B);
        n_I = (N_I - l_I) / 2;
    }
};

class SphericalSetting {
public:
    double b_F = 0.0;                        // η = (r/b)²
    int Nr_I = 0;                            // N_r = 2n_max + l_max + 8
    bool useAxialSym_B = false;              // false: N → j → m; true: m → N → j.
    bool useParity_B = false;                // Retain separate N blocks within each group.
    bool useTimeReversal_B = false;          // m > 0 representatives.
    std::vector<int> Nshell_I1D_Nshell{};           // {N₁,N₂,...}
    std::vector<SphericalSPLabel> labels_S1D_sp{}; // α_sp = (n,l,j,m)_sp
    std::vector<std::vector<SphericalSPLabel>> labels_S2D_block_bsp{}; // α_(block,bsp).
    std::vector<std::vector<int>> indices_I2D_block_bsp{}; // sp(block,bsp).

    /**
     * @brief  Construct a spherical harmonic-oscillator setting.
     * @math   N ∈ {N₁,N₂,...}, N_r = 2n_max + l_max + 8
     * @output Labels, shell groups, indices, and radial quadrature order.
     * @note   Shells must be nonempty, nonnegative, and strictly increasing.
     * @note   useAxialSym_B selects m-based instead of N-based enumeration.
     * @note   useParity_B retains N or (m,N) blocks when enabled.
     * @note   Otherwise, use one total block or separate m blocks.
     * @note   useTimeReversal_B retains only m > 0 representatives.
     * @note   Partner states follow from time reversal.
     * @note   Shell groups need not decouple the Hamiltonian.
     */
    SphericalSetting(double b_F_, const std::vector<int>& Nshell_I1D_Nshell_, bool useAxialSym_B_, bool useParity_B_, bool useTimeReversal_B_) {
        // (b,{N},useAxialSym,useParity,useTimeReversal) → setting.
        b_F = b_F_;
        Nshell_I1D_Nshell = Nshell_I1D_Nshell_;
        useAxialSym_B = useAxialSym_B_;
        useParity_B = useParity_B_;
        useTimeReversal_B = useTimeReversal_B_;

        assert(std::isfinite(b_F) && b_F > 0.0);
        assert(!Nshell_I1D_Nshell.empty());
        for (int N_I : Nshell_I1D_Nshell) {assert(N_I >= 0);}
        assert(std::is_sorted(Nshell_I1D_Nshell.begin(), Nshell_I1D_Nshell.end()));
        assert(std::adjacent_find(Nshell_I1D_Nshell.begin(), Nshell_I1D_Nshell.end()) == Nshell_I1D_Nshell.end());

        // useAxialSym: (N → j → m) → (m → N → j).
        if (useAxialSym_B) {fill_labels_by_m();} else {fill_labels_by_N();}
        assert(!labels_S1D_sp.empty());

        // {α_sp} → (n_max,l_max) → N_r.
        int nMax_I = 0;
        int lMax_I = 0;
        for (const SphericalSPLabel& label_ : labels_S1D_sp) {
            nMax_I = std::max(nMax_I, label_.n_I);
            lMax_I = std::max(lMax_I, label_.l_I);
        }
        Nr_I = 2 * nMax_I + lMax_I + 8;
    }

private:
    /**
     * @brief  Enumerate spherical harmonic-oscillator labels.
     * @math   N → j(n,l) → m
     * @output Updated labels, shell groups, and global indices.
     * @note   Time-reversal reduction applies to all label/index containers.
     * @note   Disabled useParity_B merges all N blocks.
     */
    void fill_labels_by_N();

    /**
     * @brief  Enumerate labels by angular-momentum projection.
     * @math   m → N → j(n,l)
     * @output Updated labels, projection groups, and global indices.
     * @note   useParity_B retains separate (m,N) blocks when enabled.
     * @note   Otherwise, each m forms one block.
     * @note   useTimeReversal_B retains only m > 0 representatives.
     */
    void fill_labels_by_m();
};

inline void SphericalSetting::fill_labels_by_N() {
    // ({α_sp},{α_block},{sp_block}) ← ∅.
    labels_S1D_sp.clear();
    labels_S2D_block_bsp.clear();
    indices_I2D_block_bsp.clear();

    // N → j(n,l) → m.
    for (int N_I : Nshell_I1D_Nshell) {
        // block = N.
        labels_S2D_block_bsp.emplace_back();
        indices_I2D_block_bsp.emplace_back();

        for (int twoj_I = 1; twoj_I <= 2 * N_I + 1; twoj_I += 2) {
            // 2m_min = -2j (full), 1 (time reversal).
            const int twomMin_I = -twoj_I + static_cast<int>(useTimeReversal_B) * (twoj_I + 1);
            for (int twom_I = twomMin_I; twom_I <= twoj_I; twom_I += 2) {
                // (N,2j,2m) → α_sp → (α_block,bsp,sp_block,bsp).
                const int sp_I = static_cast<int>(labels_S1D_sp.size());
                labels_S1D_sp.emplace_back(N_I, twoj_I, twom_I);
                labels_S2D_block_bsp.back().push_back(labels_S1D_sp.back());
                indices_I2D_block_bsp.back().push_back(sp_I);
            }
        }
    }

    // useParity = false: {B_N} → B.
    if (!useParity_B) {
        labels_S2D_block_bsp.clear();
        indices_I2D_block_bsp.clear();
        labels_S2D_block_bsp.push_back(labels_S1D_sp);
        indices_I2D_block_bsp.emplace_back();

        // bsp = sp.
        for (int sp_I = 0; sp_I < static_cast<int>(labels_S1D_sp.size()); ++sp_I) {
            indices_I2D_block_bsp.back().push_back(sp_I);
        }
    }
}

inline void SphericalSetting::fill_labels_by_m() {
    // ({α_sp},{α_block},{sp_block}) ← ∅.
    labels_S1D_sp.clear();
    labels_S2D_block_bsp.clear();
    indices_I2D_block_bsp.clear();

    // 2m_max = 2N_max+1; 2m_min = -2m_max or 1.
    const int twomMax_I = 2 * Nshell_I1D_Nshell.back() + 1;
    const int twomMin_I = -twomMax_I + static_cast<int>(useTimeReversal_B) * (twomMax_I + 1);

    // m → N → j(n,l).
    for (int twom_I = twomMin_I; twom_I <= twomMax_I; twom_I += 2) {
        std::vector<SphericalSPLabel> labelsM_S1D_bsp{};
        std::vector<int> indicesM_I1D_bsp{};

        for (int N_I : Nshell_I1D_Nshell) {
            if (std::abs(twom_I) > 2 * N_I + 1) {continue;}

            // useParity: B_m → B_(m,N).
            std::vector<SphericalSPLabel> labelsN_S1D_bsp{};
            std::vector<int> indicesN_I1D_bsp{};
            auto& labelsTarget_S1D_bsp = useParity_B ? labelsN_S1D_bsp : labelsM_S1D_bsp;
            auto& indicesTarget_I1D_bsp = useParity_B ? indicesN_I1D_bsp : indicesM_I1D_bsp;

            // |m| ≤ j ≤ N+1/2.
            for (int twoj_I = std::abs(twom_I); twoj_I <= 2 * N_I + 1; twoj_I += 2) {
                // (N,2j,2m) → α_sp → (α_block,bsp,sp_block,bsp).
                const int sp_I = static_cast<int>(labels_S1D_sp.size());
                labels_S1D_sp.emplace_back(N_I, twoj_I, twom_I);
                labelsTarget_S1D_bsp.push_back(labels_S1D_sp.back());
                indicesTarget_I1D_bsp.push_back(sp_I);
            }

            // B_(m,N) → blocks.
            if (useParity_B && !labelsN_S1D_bsp.empty()) {labels_S2D_block_bsp.push_back(std::move(labelsN_S1D_bsp));}
            if (useParity_B && !indicesN_I1D_bsp.empty()) {indices_I2D_block_bsp.push_back(std::move(indicesN_I1D_bsp));}
        }

        // B_m → blocks.
        if (!useParity_B && !labelsM_S1D_bsp.empty()) {labels_S2D_block_bsp.push_back(std::move(labelsM_S1D_bsp));}
        if (!useParity_B && !indicesM_I1D_bsp.empty()) {indices_I2D_block_bsp.push_back(std::move(indicesM_I1D_bsp));}
    }
}
