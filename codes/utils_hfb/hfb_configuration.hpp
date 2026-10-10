/**
 * @file    hfb_configuration.hpp
 * @author  cailea
 * @date    2026-09-27
 * @brief   Quasiparticle configuration selection.
 */

#pragma once

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <utility>
#include <vector>
#include <Eigen/Core>

class HFBConfiguration {
public:
    int Nqp_I = 0; // Complete quasiparticle space, including time-reversed partners.

    // q = Ncqp(i); cutoffs share index i.
    Eigen::VectorXi Ncqp_I1D_Ncqp{};
    Eigen::VectorXd ECut_F1D_Ncqp{};     // ΣE_qp ≤ ECut [MeV].
    Eigen::VectorXi NCut_I1D_Ncqp{};     // Retained configuration count ≤ NCut.
    Eigen::VectorXi TwoKCut_I1D_Ncqp{};  // |Σ2K_qp| ≤ TwoKCut.

    // Ascending quasiparticle indices; {} denotes vacuum.
    std::vector<std::vector<int>> config_I2D_cfg_cqp{};

public:
    /**
     * @brief Construct empty configuration settings.
     * @math C = ∅.
     * @output Empty cutoff vectors and configuration list.
     */
    HFBConfiguration() = default;

    /**
     * @brief Construct configuration settings from supplied cutoffs.
     * @math q = Ncqp(i); E ≤ ECut(i), 2|K| ≤ TwoKCut(i).
     * @output Stored Nqp and cutoffs; empty configuration list.
     * @note Equal-length nonnegative vectors; Ncqp strictly increasing; Nqp ≥ 0.
     */
    HFBConfiguration(int Nqp_I_, const Eigen::VectorXi& Ncqp_I1D_Ncqp_, const Eigen::VectorXd& ECut_F1D_Ncqp_, const Eigen::VectorXi& NCut_I1D_Ncqp_, const Eigen::VectorXi& TwoKCut_I1D_Ncqp_) {
        assert(Nqp_I_ >= 0);
        assert((Ncqp_I1D_Ncqp_.array() >= 0).all() && (Ncqp_I1D_Ncqp_.array() <= Nqp_I_).all());
        assert(Ncqp_I1D_Ncqp_.size() == ECut_F1D_Ncqp_.size());
        assert(Ncqp_I1D_Ncqp_.size() == NCut_I1D_Ncqp_.size());
        assert(Ncqp_I1D_Ncqp_.size() == TwoKCut_I1D_Ncqp_.size());

        // (Nqp,Ncqp,ECut,NCut,TwoKCut) → stored settings.
        Nqp_I = Nqp_I_;
        Ncqp_I1D_Ncqp = Ncqp_I1D_Ncqp_;
        ECut_F1D_Ncqp = ECut_F1D_Ncqp_;
        NCut_I1D_Ncqp = NCut_I1D_Ncqp_;
        TwoKCut_I1D_Ncqp = TwoKCut_I1D_Ncqp_;
    }

    /**
     * @brief Enumerate zero-to-four quasiparticle configurations using nested loops.
     * @math |c| = q; E_c = ΣE_i; 2K_c = Σ2K_i.
     * @output Rebuilt configurations, ordered by count, energy, then indices.
     * @note Inputs follow U,V columns; empty allowed list permits all orbitals.
     */
    void build_config(const Eigen::VectorXd& Eqp_F1D_qp, const Eigen::VectorXi& TwoK_I1D_qp, const std::vector<int>& qp_I1D_allowed = {});
};

inline void HFBConfiguration::build_config(const Eigen::VectorXd& Eqp_F1D_qp, const Eigen::VectorXi& TwoK_I1D_qp, const std::vector<int>& qp_I1D_allowed) {
    assert(Nqp_I >= 0);
    assert(Eqp_F1D_qp.size() == Nqp_I && TwoK_I1D_qp.size() == Nqp_I);
    for (const int qp_I : qp_I1D_allowed) {assert(qp_I >= 0 && qp_I < Nqp_I);}

    assert(Ncqp_I1D_Ncqp.size() == ECut_F1D_Ncqp.size());
    assert(Ncqp_I1D_Ncqp.size() == NCut_I1D_Ncqp.size());
    assert(Ncqp_I1D_Ncqp.size() == TwoKCut_I1D_Ncqp.size());
    assert((Ncqp_I1D_Ncqp.array() >= 0).all() && (Ncqp_I1D_Ncqp.array() <= 4).all());

    // Empty list → full space; ascending unique qp indices.
    std::vector<int> qp_I1D_candidate = qp_I1D_allowed;
    if (qp_I1D_candidate.empty()) {for (int qp_I = 0; qp_I < Nqp_I; ++qp_I) {qp_I1D_candidate.push_back(qp_I);}}
    std::sort(qp_I1D_candidate.begin(), qp_I1D_candidate.end());
    qp_I1D_candidate.erase(std::unique(qp_I1D_candidate.begin(), qp_I1D_candidate.end()), qp_I1D_candidate.end());
    const int Ncandidate_I = static_cast<int>(qp_I1D_candidate.size());

    // C → ∅; reusable candidate buffer.
    config_I2D_cfg_cqp.clear();
    std::vector<std::pair<double, std::vector<int>>> candidates{};

    // q = 0: E = K = 0.
    const auto build_0qp = [&](int NCut_I, double ECut_F, int twoKCut_I) {
        if (NCut_I > 0 && 0.0 <= ECut_F && 0 <= twoKCut_I) {candidates.emplace_back(0.0, std::vector<int>{});}
    };

    // q = 1: increasing orbital indices.
    const auto build_1qp = [&](int NCut_I, double ECut_F, int twoKCut_I) {
        for (const int qp1_I : qp_I1D_candidate) {
            // E_c = ΣE_i; 2K_c = Σ2K_i.
            const double Econfig_F = Eqp_F1D_qp(qp1_I);
            const long long TwoKconfig_I = static_cast<long long>(TwoK_I1D_qp(qp1_I));
            // E_c ≤ Emax; |2K_c| ≤ 2Kmax.
            if (Econfig_F <= ECut_F && std::abs(TwoKconfig_I) <= twoKCut_I) {candidates.emplace_back(Econfig_F, std::vector<int>{qp1_I});}
        }

        std::sort(candidates.begin(), candidates.end());
        candidates.resize(std::min(candidates.size(), static_cast<std::size_t>(NCut_I)));
    };

    // q = 2: increasing orbital indices.
    const auto build_2qp = [&](int NCut_I, double ECut_F, int twoKCut_I) {
        for (int i1_I = 0; i1_I < Ncandidate_I - 1; ++i1_I) {
            const int qp1_I = qp_I1D_candidate[i1_I];
            for (int i2_I = i1_I + 1; i2_I < Ncandidate_I; ++i2_I) {
                const int qp2_I = qp_I1D_candidate[i2_I];
                // E_c = ΣE_i; 2K_c = Σ2K_i.
                const double Econfig_F = Eqp_F1D_qp(qp1_I) + Eqp_F1D_qp(qp2_I);
                const long long TwoKconfig_I = static_cast<long long>(TwoK_I1D_qp(qp1_I)) + TwoK_I1D_qp(qp2_I);
                // E_c ≤ Emax; |2K_c| ≤ 2Kmax.
                if (Econfig_F <= ECut_F && std::abs(TwoKconfig_I) <= twoKCut_I) {candidates.emplace_back(Econfig_F, std::vector<int>{qp1_I, qp2_I});}
            }
        }

        std::sort(candidates.begin(), candidates.end());
        candidates.resize(std::min(candidates.size(), static_cast<std::size_t>(NCut_I)));
    };

    // q = 3: increasing orbital indices.
    const auto build_3qp = [&](int NCut_I, double ECut_F, int twoKCut_I) {
        for (int i1_I = 0; i1_I < Ncandidate_I - 2; ++i1_I) {
            const int qp1_I = qp_I1D_candidate[i1_I];
            for (int i2_I = i1_I + 1; i2_I < Ncandidate_I - 1; ++i2_I) {
                const int qp2_I = qp_I1D_candidate[i2_I];
                for (int i3_I = i2_I + 1; i3_I < Ncandidate_I; ++i3_I) {
                    const int qp3_I = qp_I1D_candidate[i3_I];
                    // E_c = ΣE_i; 2K_c = Σ2K_i.
                    const double Econfig_F = Eqp_F1D_qp(qp1_I) + Eqp_F1D_qp(qp2_I) + Eqp_F1D_qp(qp3_I);
                    const long long TwoKconfig_I = static_cast<long long>(TwoK_I1D_qp(qp1_I)) + TwoK_I1D_qp(qp2_I) + TwoK_I1D_qp(qp3_I);
                    // E_c ≤ Emax; |2K_c| ≤ 2Kmax.
                    if (Econfig_F <= ECut_F && std::abs(TwoKconfig_I) <= twoKCut_I) {candidates.emplace_back(Econfig_F, std::vector<int>{qp1_I, qp2_I, qp3_I});}
                }
            }
        }

        std::sort(candidates.begin(), candidates.end());
        candidates.resize(std::min(candidates.size(), static_cast<std::size_t>(NCut_I)));
    };

    // q = 4: increasing orbital indices.
    const auto build_4qp = [&](int NCut_I, double ECut_F, int twoKCut_I) {
        for (int i1_I = 0; i1_I < Ncandidate_I - 3; ++i1_I) {
            const int qp1_I = qp_I1D_candidate[i1_I];
            for (int i2_I = i1_I + 1; i2_I < Ncandidate_I - 2; ++i2_I) {
                const int qp2_I = qp_I1D_candidate[i2_I];
                for (int i3_I = i2_I + 1; i3_I < Ncandidate_I - 1; ++i3_I) {
                    const int qp3_I = qp_I1D_candidate[i3_I];
                    for (int i4_I = i3_I + 1; i4_I < Ncandidate_I; ++i4_I) {
                        const int qp4_I = qp_I1D_candidate[i4_I];
                        // E_c = ΣE_i; 2K_c = Σ2K_i.
                        const double Econfig_F = Eqp_F1D_qp(qp1_I) + Eqp_F1D_qp(qp2_I) + Eqp_F1D_qp(qp3_I) + Eqp_F1D_qp(qp4_I);
                        const long long TwoKconfig_I = static_cast<long long>(TwoK_I1D_qp(qp1_I)) + TwoK_I1D_qp(qp2_I) + TwoK_I1D_qp(qp3_I) + TwoK_I1D_qp(qp4_I);
                        // E_c ≤ Emax; |2K_c| ≤ 2Kmax.
                        if (Econfig_F <= ECut_F && std::abs(TwoKconfig_I) <= twoKCut_I) {candidates.emplace_back(Econfig_F, std::vector<int>{qp1_I, qp2_I, qp3_I, qp4_I});}
                    }
                }
            }
        }

        std::sort(candidates.begin(), candidates.end());
        candidates.resize(std::min(candidates.size(), static_cast<std::size_t>(NCut_I)));
    };

    for (Eigen::Index iNcqp_I = 0; iNcqp_I < Ncqp_I1D_Ncqp.size(); ++iNcqp_I) {
        // (q, Emax, Nmax, 2Kmax).
        const int Ncqp_I = Ncqp_I1D_Ncqp(iNcqp_I);
        const double ECut_F = ECut_F1D_Ncqp(iNcqp_I);
        const int NCut_I = NCut_I1D_Ncqp(iNcqp_I);
        const int twoKCut_I = TwoKCut_I1D_Ncqp(iNcqp_I);
        if (NCut_I == 0 || Ncqp_I > Ncandidate_I) {continue;}

        // Reuse buffers for fixed q.
        candidates.clear();
        switch (Ncqp_I) {
            case 0: build_0qp(NCut_I, ECut_F, twoKCut_I); break;
            case 1: build_1qp(NCut_I, ECut_F, twoKCut_I); break;
            case 2: build_2qp(NCut_I, ECut_F, twoKCut_I); break;
            case 3: build_3qp(NCut_I, ECut_F, twoKCut_I); break;
            case 4: build_4qp(NCut_I, ECut_F, twoKCut_I); break;
        }

        const std::size_t Ncfg_I = candidates.size();

        // Append the lowest Ncfg configurations.
        for (std::size_t cfg_I = 0; cfg_I < Ncfg_I; ++cfg_I) {
            config_I2D_cfg_cqp.push_back(std::move(candidates[cfg_I].second));
        }
    }
}
