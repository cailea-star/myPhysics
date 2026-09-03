/**
 * @file    axial_config.hpp
 * @author  cailea
 * @date    2026-05-02
 * @brief   Define axial harmonic-oscillator labels and configurations.
 */

#pragma once
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <vector>
#include <utility>

class AxialSPLabel {
public:
    int N_I = 0;                     // N = n_z + 2n_r + Λ
    int nz_I = 0;                    // n_z ≥ 0
    int nr_I = 0;                    // n_r ≥ 0
    int Lambda_I = 0;                // Λ ≥ 0
    int twoOmega_I = 0;              // 2Ω = 2Λ ± 1
    int twoSigma_I = 0;              // 2Σ = ±1
    bool isParityPositive_B = false; // π = (-1)^(n_z + Λ) = +1

    /**
     * @brief  Construct an axial single-particle label.
     * @math   N = n_z + 2n_r + Λ, 2Ω = 2Λ ± 1, π = (-1)^(n_z + Λ)
     * @output Initialized label.
     */
    AxialSPLabel(int nz_I_, int nr_I_, int Lambda_I_, int twoOmega_I_) {
        N_I = nz_I_ + 2 * nr_I_ + Lambda_I_;
        nz_I = nz_I_;
        nr_I = nr_I_;
        Lambda_I = Lambda_I_;
        twoOmega_I = twoOmega_I_;
        twoSigma_I = twoOmega_I_ - 2 * Lambda_I_;
        isParityPositive_B = (nz_I_ + Lambda_I_) % 2 == 0;
    }
};

/**
 * @brief  Check the axial single-particle quantum-number constraints.
 * @math   n_z,n_r,Λ ≥ 0, 2Ω = 2Λ + 2Σ, π = (-1)^(n_z+Λ)
 * @output True for a valid retained label.
 */
inline bool is_valid(int nz_I, int nr_I, int Lambda_I, int twoOmega_I, int twoSigma_I, bool isParityPositive_B) {
    // (n_z,n_r,Λ,2Ω,2Σ,π) → valid(α).
    const bool nzisNonnegative_B = nz_I >= 0;
    const bool nrisNonnegative_B = nr_I >= 0;
    const bool LambdaisNonnegative_B = Lambda_I >= 0;
    const bool twoOmegaisPositive_B = twoOmega_I > 0;
    const bool twoSigmaisValid_B = std::abs(twoSigma_I) == 1;
    const bool twoOmegaisConsistent_B = twoOmega_I == 2 * Lambda_I + twoSigma_I;
    const bool parityisConsistent_B = ((nz_I + Lambda_I) % 2 == 0) == isParityPositive_B;
    return nzisNonnegative_B & nrisNonnegative_B & LambdaisNonnegative_B & twoOmegaisPositive_B & twoSigmaisValid_B & twoOmegaisConsistent_B & parityisConsistent_B;
}

class AxialConfig {
public:
    int Nshell_I = 0;                                           // N_shell ≥ 0
    int Nz_I = 0;                                               // N_z = 2n_z^max + 8
    int Nr_I = 0;                                               // N_r = 2n_r^max + Λ^max + 8
    double bz_F = 0.0;                                          // ζ = z / b_z
    double br_F = 0.0;                                          // η = (r / b_r)^2
    bool useReflection_B = false;                               // z ≥ 0 when true
    std::vector<AxialSPLabel> labels_S1D_sp{};                  // α_sp = (n_z,n_r,Λ,Ω,Σ,π)_sp
    std::vector<std::vector<AxialSPLabel>> labels_S2D_block_bsp{}; // α_(block,bsp): labels grouped by Ω or (Ω,π)
    std::vector<std::vector<int>> indices_I2D_block_bsp{};      // sp(block,bsp): global indices of block labels

    /**
     * @brief  Construct an axial harmonic-oscillator configuration.
     * @math   N_z = 2n_z^max + 8, N_r = 2n_r^max + Λ^max + 8
     * @output Single-particle labels, symmetry blocks, and quadrature orders.
     */
    AxialConfig(double bz_F_, double br_F_, int Nshell_I_, bool useReflection_B_) {
        Nshell_I = Nshell_I_;
        bz_F = bz_F_;
        br_F = br_F_;
        useReflection_B = useReflection_B_;

        assert(std::isfinite(bz_F) && bz_F > 0.0);
        assert(std::isfinite(br_F) && br_F > 0.0);
        assert(Nshell_I >= 0);

        // (N_shell,b_z,b_r,reflection) → ({α_sp},{α_block},{sp_block}).
        fill_labels();
        assert(!labels_S1D_sp.empty());

        // {α_sp} → (n_z^max,n_r^max,Λ^max)
        int nzMax_I = 0;
        int nrMax_I = 0;
        int LambdaMax_I = 0;
        for (const AxialSPLabel& label_ : labels_S1D_sp) {
            nzMax_I = std::max(nzMax_I, label_.nz_I);
            nrMax_I = std::max(nrMax_I, label_.nr_I);
            LambdaMax_I = std::max(LambdaMax_I, label_.Lambda_I);
        }

        Nz_I = 2 * nzMax_I + 8;
        Nr_I = 2 * nrMax_I + LambdaMax_I + 8;
    }

    /**
     * @brief  Map (b₀,β₂₀) to axial oscillator lengths.
     * @math   b_z = b_0 exp[√(5/(16π))β_20], b_r = b_0 exp[-√(5/(16π))β_20/2]
     * @output {b_z,b_r}.
     */
    static std::pair<double, double> b0beta20_to_bzbr(double b0_F, double beta20_F);

    /**
     * @brief  Map axial oscillator lengths to (b₀,β₂₀).
     * @math   b_0 = (b_r²b_z)^(1/3), β_20 = (8/3)√(π/5)log(b_z/b_r)
     * @output {b_0,β_20}.
     */
    static std::pair<double, double> bzbr_to_b0beta20(double bz_F, double br_F);

private:
    /**
     * @brief  Evaluate the reduced axial oscillator energy.
     * @math   ε = (n_z+1/2)/b_z² + (2n_r+Λ+1)/b_r²
     * @output Reduced single-particle energy.
     */
    double calc_epsilon(int nz_I, int nr_I, int Lambda_I) const;

    /**
     * @brief  Determine the energy cutoff by ordered multiplicity accumulation.
     * @math   Σ_(ε≤E_cut) d_ν > N_HO/2
     * @output Reduced cutoff energy E_cut.
     */
    double calc_Ecut() const;

    /**
     * @brief  Enumerate retained positive-Ω labels and symmetry blocks.
     * @math   N = n_z+2n_r+Λ, 2Ω = 2Λ±1, π = (-1)^(n_z+Λ)
     * @output Updated label and index containers.
     */
    void fill_labels();
};

inline std::pair<double, double> AxialConfig::b0beta20_to_bzbr(double b0_F, double beta20_F) {
    // (b_0,β_20) → (b_z,b_r)
    assert(std::isfinite(b0_F) && b0_F > 0.0);
    assert(std::isfinite(beta20_F));
    const double pi_F = std::acos(-1.0);
    const double deformation_F = std::sqrt(5.0 / (16.0 * pi_F)) * beta20_F;
    return {b0_F * std::exp(deformation_F), b0_F * std::exp(-0.5 * deformation_F)};
}

inline std::pair<double, double> AxialConfig::bzbr_to_b0beta20(double bz_F, double br_F) {
    // (b_z,b_r) → (b_0,β_20)
    assert(std::isfinite(bz_F) && bz_F > 0.0);
    assert(std::isfinite(br_F) && br_F > 0.0);
    const double pi_F = std::acos(-1.0);
    const double b0_F = std::cbrt(br_F * br_F * bz_F);
    const double beta20_F = (8.0 / 3.0) * std::sqrt(pi_F / 5.0) * std::log(bz_F / br_F);
    return {b0_F, beta20_F};
}

inline double AxialConfig::calc_epsilon(int nz_I, int nr_I, int Lambda_I) const {
    // (b_z,b_r) → (ε_z,ε_r)
    const double epsilonz_F = 1.0 / (bz_F * bz_F);
    const double epsilonr_F = 1.0 / (br_F * br_F);
    return epsilonz_F * (nz_I + 0.5) + epsilonr_F * (2 * nr_I + Lambda_I + 1);
}

inline double AxialConfig::calc_Ecut() const {
    // N_HO=Σ_{N=0}^{N_shell}(N+1)(N+2)/2; n_z^cut=ν^cut=N_HO/2.
    const int NHO_I = (Nshell_I + 1) * (Nshell_I + 2) * (Nshell_I + 3) / 6;
    const int nzCut_I = NHO_I / 2;
    const int nuCut_I = NHO_I / 2;

    // (n_z,ν,d_ν) candidates, ν = 2n_r + Λ
    std::vector<std::array<int, 3>> candidates_I2D_candidate_field;
    candidates_I2D_candidate_field.reserve((nzCut_I + 1) * (nuCut_I + 1));
    for (int nz_I = 0; nz_I <= nzCut_I; ++nz_I) {
        for (int nrml_I = 0; nrml_I <= nuCut_I; ++nrml_I) {
            const int nrMax_I = static_cast<int>(std::floor(nrml_I / 2.0)); // n_r^max = floor(ν/2).
            candidates_I2D_candidate_field.push_back({nz_I, nrml_I, nrMax_I + 1});
        }
    }

    // {(n_z,ν,d_ν)} → ascending ε
    const auto compareEpsilon_Func = [&](const std::array<int, 3>& candidate1_I1D_field, const std::array<int, 3>& candidate2_I1D_field) {
        return calc_epsilon(candidate1_I1D_field[0], 0, candidate1_I1D_field[1]) < calc_epsilon(candidate2_I1D_field[0], 0, candidate2_I1D_field[1]);
    };
    std::sort(candidates_I2D_candidate_field.begin(), candidates_I2D_candidate_field.end(), compareEpsilon_Func);

    // Σd_ν > N_HO/2
    const int Ntarget_I = NHO_I / 2;
    int Naccum_I = 0;
    int candidate_I = 0;
    for (; candidate_I < static_cast<int>(candidates_I2D_candidate_field.size()); ++candidate_I) {
        Naccum_I += candidates_I2D_candidate_field[candidate_I][2];
        if (Naccum_I > Ntarget_I) {break;}
    }
    assert(candidate_I < static_cast<int>(candidates_I2D_candidate_field.size()));
    return calc_epsilon(candidates_I2D_candidate_field[candidate_I][0], 0, candidates_I2D_candidate_field[candidate_I][1]);
}

inline void AxialConfig::fill_labels() {
    // (N_shell,b_z,b_r) → (bounds,E_cut).
    const int NHO_I = (Nshell_I + 1) * (Nshell_I + 2) * (Nshell_I + 3) / 6;
    const int nzCut_I = NHO_I / 2;
    const int nuCut_I = NHO_I / 2;
    const double Ecut_F = calc_Ecut();

    // ({α_sp},{α_block},{sp_block}) ← ∅.
    labels_S1D_sp.clear();
    labels_S2D_block_bsp.clear();
    indices_I2D_block_bsp.clear();
    labels_S1D_sp.reserve(NHO_I);

    // Ω → π → n_r → n_z → Λ_± → Σ_±
    for (int twoOmega_I = 1; twoOmega_I <= 2 * nuCut_I + 1; ++twoOmega_I) {
        std::vector<AxialSPLabel> labelsOmega_S1D_bsp;
        std::vector<int> indicesOmega_I1D_bsp;
        const int LambdaUp_I = (twoOmega_I - 1) / 2;
        const int LambdaDown_I = (twoOmega_I + 1) / 2;
        const int nrMax_I = (nuCut_I - LambdaUp_I + 1) / 2;
        for (bool isParityPositive_B : {true, false}) {
            std::vector<AxialSPLabel> labelsParity_S1D_bsp;
            std::vector<int> indicesParity_I1D_bsp;
            auto& labelsTarget_S1D_bsp = useReflection_B ? labelsParity_S1D_bsp : labelsOmega_S1D_bsp;
            auto& indicesTarget_I1D_bsp = useReflection_B ? indicesParity_I1D_bsp : indicesOmega_I1D_bsp;
            for (int nr_I = 0; nr_I <= nrMax_I; ++nr_I) {
                for (int nz_I = 0; nz_I <= nzCut_I; ++nz_I) {
                    for (int Lambda_I : {LambdaUp_I, LambdaDown_I}) {
                        for (int twoSigma_I : {1, -1}) {
                            if (!is_valid(nz_I, nr_I, Lambda_I, twoOmega_I, twoSigma_I, isParityPositive_B)) {continue;}
                            if (calc_epsilon(nz_I, nr_I, Lambda_I) > Ecut_F) {continue;}
                            const AxialSPLabel label_(nz_I, nr_I, Lambda_I, twoOmega_I);
                            const int sp_I = static_cast<int>(labels_S1D_sp.size());
                            labels_S1D_sp.push_back(label_);
                            labelsTarget_S1D_bsp.push_back(label_);
                            indicesTarget_I1D_bsp.push_back(sp_I);
                        }
                    }
                }
            }
            if (useReflection_B && !labelsParity_S1D_bsp.empty()) {labels_S2D_block_bsp.push_back(std::move(labelsParity_S1D_bsp));}
            if (useReflection_B && !indicesParity_I1D_bsp.empty()) {indices_I2D_block_bsp.push_back(std::move(indicesParity_I1D_bsp));}
        }
        if (!useReflection_B && !labelsOmega_S1D_bsp.empty()) {labels_S2D_block_bsp.push_back(std::move(labelsOmega_S1D_bsp));}
        if (!useReflection_B && !indicesOmega_I1D_bsp.empty()) {indices_I2D_block_bsp.push_back(std::move(indicesOmega_I1D_bsp));}
    }
}
