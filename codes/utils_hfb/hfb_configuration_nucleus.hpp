/**
 * @file    hfb_configuration_nucleus.hpp
 * @author  cailea
 * @date    2026-10-06
 * @brief   Neutron-proton quasiparticle configuration selection.
 */

#pragma once

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <utility>
#include <vector>
#include <Eigen/Core>

class HFBConfigurationNucleus {
public:
    // Ncqp(i,n/p); each row has independent E,N,2K cuts.
    Eigen::MatrixXi Ncqp_I2D_Ncqp_np{};
    Eigen::VectorXd ECut_F1D_Ncqp{};     // ΣE_qp ≤ ECut [MeV].
    Eigen::VectorXi NCut_I1D_Ncqp{};     // Retained configuration count ≤ NCut.
    Eigen::VectorXi TwoKCut_I1D_Ncqp{};  // |Σ2K_qp| ≤ TwoKCut.

    // Φ_cfg = Φn[confign(cfg)] ⊗ Φp[configp(cfg)].
    std::vector<std::vector<int>> confign_I2D_cfgnp_cqpn{};
    std::vector<std::vector<int>> configp_I2D_cfgnp_cqpp{};

    // Φn ⊗ Φp: unique species factors.
    std::vector<std::vector<int>> confign_I2D_tensorcfgn_cqpn{};
    std::vector<std::vector<int>> configp_I2D_tensorcfgp_cqpp{};

public:
    /**
     * @brief Construct empty neutron-proton configuration settings.
     * @math Cn = Cp = ∅.
     * @output Empty nuclear cutoffs and paired neutron/proton lists.
     */
    HFBConfigurationNucleus() = default;

    /**
     * @brief Store supplied nuclear configuration cutoffs.
     * @math Ncqp(i,n/p), Emax(i), Nmax(i), 2Kmax(i) ≥ 0.
     * @output Stored cutoffs; empty paired neutron/proton configuration lists.
     * @note Ncqp has two columns; all cutoff lengths match rows.
     */
    HFBConfigurationNucleus(const Eigen::MatrixXi& Ncqp_I2D_Ncqp_np_, const Eigen::VectorXd& ECut_F1D_Ncqp_, const Eigen::VectorXi& NCut_I1D_Ncqp_, const Eigen::VectorXi& TwoKCut_I1D_Ncqp_) {
        assert(Ncqp_I2D_Ncqp_np_.cols() == 2);
        assert(Ncqp_I2D_Ncqp_np_.rows() == ECut_F1D_Ncqp_.size() && Ncqp_I2D_Ncqp_np_.rows() == NCut_I1D_Ncqp_.size() && Ncqp_I2D_Ncqp_np_.rows() == TwoKCut_I1D_Ncqp_.size());
        assert((Ncqp_I2D_Ncqp_np_.array() >= 0).all());
        assert((ECut_F1D_Ncqp_.array() >= 0.0).all() && (NCut_I1D_Ncqp_.array() >= 0).all() && (TwoKCut_I1D_Ncqp_.array() >= 0).all());

        // (Ncqp,Emax,Nmax,2Kmax) → owned settings.
        Ncqp_I2D_Ncqp_np = Ncqp_I2D_Ncqp_np_;
        ECut_F1D_Ncqp = ECut_F1D_Ncqp_;
        NCut_I1D_Ncqp = NCut_I1D_Ncqp_;
        TwoKCut_I1D_Ncqp = TwoKCut_I1D_Ncqp_;
    }

    /**
     * @brief Select EE nuclear configurations using standard PSM windows.
     * @math Φcfg = Φn,cn ⊗ Φp,cp; Ecfg = ΣE_qp.
     * @output Mixed configurations and unique tensor factors.
     * @note Nilsson-ordered representatives; n2p2 requires n2p0 and n0p2 cuts.
     */
    void build_config_kramers_ee(const Eigen::VectorXd& Eqpn_F1D_qpn, const Eigen::VectorXi& TwoKn_I1D_qpn, const std::vector<int>& qpn_I1D_representative, const Eigen::VectorXd& Eqpp_F1D_qpp, const Eigen::VectorXi& TwoKp_I1D_qpp, const std::vector<int>& qpp_I1D_representative);

    /**
     * @brief Select EO configurations using standard PSM windows.
     * @math Φcfg = Φn,0⊗Φp,1 or Φn,2⊗Φp,1.
     * @output Mixed configurations and unique tensor factors.
     * @note n2p0 supplies auxiliary pair cuts; n0p1 selects odd orbits.
     */
    void build_config_kramers_eo(const Eigen::VectorXd& Eqpn_F1D_qpn, const Eigen::VectorXi& TwoKn_I1D_qpn, const std::vector<int>& qpn_I1D_representative, const Eigen::VectorXd& Eqpp_F1D_qpp, const Eigen::VectorXi& TwoKp_I1D_qpp, const std::vector<int>& qpp_I1D_representative);

    /**
     * @brief Select OE configurations using standard PSM windows.
     * @math Φcfg = Φn,1⊗Φp,0 or Φn,1⊗Φp,2.
     * @output Mixed configurations and unique tensor factors.
     * @note n0p2 supplies auxiliary pair cuts; n1p0 selects odd orbits.
     */
    void build_config_kramers_oe(const Eigen::VectorXd& Eqpn_F1D_qpn, const Eigen::VectorXi& TwoKn_I1D_qpn, const std::vector<int>& qpn_I1D_representative, const Eigen::VectorXd& Eqpp_F1D_qpp, const Eigen::VectorXi& TwoKp_I1D_qpp, const std::vector<int>& qpp_I1D_representative);

    /**
     * @brief Select OO configurations using standard PSM windows.
     * @math Φcfg = βn†βp†Φ or βn†βp̄†Φ.
     * @output Mixed configurations and unique tensor factors.
     * @note n1p0 and n0p1 supply auxiliary single-orbit windows.
     */
    void build_config_kramers_oo(const Eigen::VectorXd& Eqpn_F1D_qpn, const Eigen::VectorXi& TwoKn_I1D_qpn, const std::vector<int>& qpn_I1D_representative, const Eigen::VectorXd& Eqpp_F1D_qpp, const Eigen::VectorXi& TwoKp_I1D_qpp, const std::vector<int>& qpp_I1D_representative);

    /**
     * @brief Find tensor-factor indices by configuration equality.
     * @math cfg → (tensorcfgn,tensorcfgp).
     * @output Neutron and proton tensor-factor indices.
     */
    std::pair<int, int> find_tensorcfg(int cfg_I) const;

private:
    /**
     * @brief Stable-sort mixed configurations and deduplicate species factors.
     * @math (Ncqp_n,Ncqp_p) ↑; cfg → tensorcfgn,tensorcfgp.
     * @output Ordered mixed configurations and unique tensor-factor tables.
     */
    void build_tensor_config();

    /**
     * @brief Select two-quasiparticle branches using standard PSM windows.
     * @math NCutHalf = NCut/2; K = Ki+Kj or Ki−Kj.
     * @output Separately energy-sorted same/partner pair branches.
     * @note cutSingle_B selects individual cuts; partners differ by Nqp/2.
     */
    static std::vector<std::vector<int>> calc_config_2qp(int NCut_I, double ECut_F, int TwoKCut_I, const Eigen::VectorXd& Eqp_F1D_qp, const Eigen::VectorXi& TwoK_I1D_qp, const std::vector<int>& qp_I1D_representative, bool cutSingle_B = false);

};

inline std::pair<int, int> HFBConfigurationNucleus::find_tensorcfg(int cfg_I) const {
    assert(cfg_I >= 0 && cfg_I < static_cast<int>(confign_I2D_cfgnp_cqpn.size()) && confign_I2D_cfgnp_cqpn.size() == configp_I2D_cfgnp_cqpp.size());
    const auto confign_iter = std::find(confign_I2D_tensorcfgn_cqpn.begin(), confign_I2D_tensorcfgn_cqpn.end(), confign_I2D_cfgnp_cqpn[cfg_I]);
    const auto configp_iter = std::find(configp_I2D_tensorcfgp_cqpp.begin(), configp_I2D_tensorcfgp_cqpp.end(), configp_I2D_cfgnp_cqpp[cfg_I]);
    assert(confign_iter != confign_I2D_tensorcfgn_cqpn.end() && configp_iter != configp_I2D_tensorcfgp_cqpp.end());
    return {static_cast<int>(confign_iter - confign_I2D_tensorcfgn_cqpn.begin()), static_cast<int>(configp_iter - configp_I2D_tensorcfgp_cqpp.begin())};
}

inline void HFBConfigurationNucleus::build_tensor_config() {
    assert(confign_I2D_cfgnp_cqpn.size() == configp_I2D_cfgnp_cqpp.size());
    // (Ncqp_n,Ncqp_p) ↑; within-class order unchanged.
    std::vector<int> cfg_I1D_order{};
    for (int cfg_I = 0; cfg_I < static_cast<int>(confign_I2D_cfgnp_cqpn.size()); ++cfg_I) {cfg_I1D_order.push_back(cfg_I);}
    std::stable_sort(cfg_I1D_order.begin(), cfg_I1D_order.end(), [&](int cfg1_I, int cfg2_I) {return std::make_pair(confign_I2D_cfgnp_cqpn[cfg1_I].size(), configp_I2D_cfgnp_cqpp[cfg1_I].size()) < std::make_pair(confign_I2D_cfgnp_cqpn[cfg2_I].size(), configp_I2D_cfgnp_cqpp[cfg2_I].size());});
    const auto confign_I2D_unsortedcfg_cqpn = confign_I2D_cfgnp_cqpn;
    const auto configp_I2D_unsortedcfg_cqpp = configp_I2D_cfgnp_cqpp;
    for (std::size_t cfg_I = 0; cfg_I < cfg_I1D_order.size(); ++cfg_I) {
        confign_I2D_cfgnp_cqpn[cfg_I] = confign_I2D_unsortedcfg_cqpn[cfg_I1D_order[cfg_I]];
        configp_I2D_cfgnp_cqpp[cfg_I] = configp_I2D_unsortedcfg_cqpp[cfg_I1D_order[cfg_I]];
    }

    // cfg → unique tensorcfgn,tensorcfgp.
    const auto select_config = [](const auto& configs_I2D_cfg_cqp, auto& unique_I2D_cfg_cqp) {
        unique_I2D_cfg_cqp.clear();
        for (const auto& config_I1D_cqp : configs_I2D_cfg_cqp) {
            if (std::find(unique_I2D_cfg_cqp.begin(), unique_I2D_cfg_cqp.end(), config_I1D_cqp) == unique_I2D_cfg_cqp.end()) {unique_I2D_cfg_cqp.push_back(config_I1D_cqp);}
        }
    };
    select_config(confign_I2D_cfgnp_cqpn, confign_I2D_tensorcfgn_cqpn);
    select_config(configp_I2D_cfgnp_cqpp, configp_I2D_tensorcfgp_cqpp);
}

inline void HFBConfigurationNucleus::build_config_kramers_ee(const Eigen::VectorXd& Eqpn_F1D_qpn, const Eigen::VectorXi& TwoKn_I1D_qpn, const std::vector<int>& qpn_I1D_representative, const Eigen::VectorXd& Eqpp_F1D_qpp, const Eigen::VectorXi& TwoKp_I1D_qpp, const std::vector<int>& qpp_I1D_representative) {
    assert(Eqpn_F1D_qpn.size() > 0 && Eqpn_F1D_qpn.size() % 2 == 0 && TwoKn_I1D_qpn.size() == Eqpn_F1D_qpn.size());
    assert(Eqpp_F1D_qpp.size() > 0 && Eqpp_F1D_qpp.size() % 2 == 0 && TwoKp_I1D_qpp.size() == Eqpp_F1D_qpp.size());
    assert(Ncqp_I2D_Ncqp_np.cols() == 2 && Ncqp_I2D_Ncqp_np.rows() == ECut_F1D_Ncqp.size() && Ncqp_I2D_Ncqp_np.rows() == NCut_I1D_Ncqp.size() && Ncqp_I2D_Ncqp_np.rows() == TwoKCut_I1D_Ncqp.size());
    assert(((Ncqp_I2D_Ncqp_np.array() == 0) || (Ncqp_I2D_Ncqp_np.array() == 2)).all());
    assert((ECut_F1D_Ncqp.array() >= 0.0).all() && (NCut_I1D_Ncqp.array() >= 0).all() && (TwoKCut_I1D_Ncqp.array() >= 0).all());
    for (const int qpn_I : qpn_I1D_representative) {assert(qpn_I >= 0 && qpn_I < Eqpn_F1D_qpn.size());}
    for (const int qpp_I : qpp_I1D_representative) {assert(qpp_I >= 0 && qpp_I < Eqpp_F1D_qpp.size());}

    // Φ0, Φn,2⊗Φp,0, Φn,0⊗Φp,2.
    const auto calc_config_kramers_n0p0 = [&](int NCut_I, double ECut_F, int TwoKCut_I) {
        assert(NCut_I >= 0 && ECut_F >= 0.0 && TwoKCut_I >= 0);
        std::vector<std::vector<int>> confign_I2D_cfgnp_cqpn{};
        std::vector<std::vector<int>> configp_I2D_cfgnp_cqpp{};
        if (NCut_I > 0 && ECut_F >= 0.0 && TwoKCut_I >= 0) {
            confign_I2D_cfgnp_cqpn.emplace_back();
            configp_I2D_cfgnp_cqpp.emplace_back();
        }
        return std::make_pair(confign_I2D_cfgnp_cqpn, configp_I2D_cfgnp_cqpp);
    };

    const auto calc_config_kramers_n2p0 = [&](int NCut_I, double ECut_F, int TwoKCut_I) {
        const auto confign_I2D_cfgnp_cqpn = calc_config_2qp(NCut_I, ECut_F, TwoKCut_I, Eqpn_F1D_qpn, TwoKn_I1D_qpn, qpn_I1D_representative);
        const std::vector<std::vector<int>> configp_I2D_cfgnp_cqpp(confign_I2D_cfgnp_cqpn.size());
        return std::make_pair(confign_I2D_cfgnp_cqpn, configp_I2D_cfgnp_cqpp);
    };

    const auto calc_config_kramers_n0p2 = [&](int NCut_I, double ECut_F, int TwoKCut_I) {
        const auto configp_I2D_cfgnp_cqpp = calc_config_2qp(NCut_I, ECut_F, TwoKCut_I, Eqpp_F1D_qpp, TwoKp_I1D_qpp, qpp_I1D_representative);
        const std::vector<std::vector<int>> confign_I2D_cfgnp_cqpn(configp_I2D_cfgnp_cqpp.size());
        return std::make_pair(confign_I2D_cfgnp_cqpn, configp_I2D_cfgnp_cqpp);
    };

    // Pure-species cuts determine the n2p2 Cartesian factors.
    Eigen::Index iN2P0_I = -1;
    Eigen::Index iN0P2_I = -1;
    for (Eigen::Index iNcqp_I = 0; iNcqp_I < Ncqp_I2D_Ncqp_np.rows(); ++iNcqp_I) {
        if (Ncqp_I2D_Ncqp_np(iNcqp_I, 0) == 2 && Ncqp_I2D_Ncqp_np(iNcqp_I, 1) == 0) {
            assert(iN2P0_I == -1);
            iN2P0_I = iNcqp_I;
        }
        if (Ncqp_I2D_Ncqp_np(iNcqp_I, 0) == 0 && Ncqp_I2D_Ncqp_np(iNcqp_I, 1) == 2) {
            assert(iN0P2_I == -1);
            iN0P2_I = iNcqp_I;
        }
    }
    std::pair<std::vector<std::vector<int>>, std::vector<std::vector<int>>> config_n2p0{};
    std::pair<std::vector<std::vector<int>>, std::vector<std::vector<int>>> config_n0p2{};
    if (iN2P0_I >= 0) {config_n2p0 = calc_config_kramers_n2p0(NCut_I1D_Ncqp(iN2P0_I), ECut_F1D_Ncqp(iN2P0_I), TwoKCut_I1D_Ncqp(iN2P0_I));}
    if (iN0P2_I >= 0) {config_n0p2 = calc_config_kramers_n0p2(NCut_I1D_Ncqp(iN0P2_I), ECut_F1D_Ncqp(iN0P2_I), TwoKCut_I1D_Ncqp(iN0P2_I));}

    // Φn,2⊗Φp,2: n outer, p inner; independent cuts.
    const auto calc_config_kramers_n2p2 = [&](int NCut_I, double ECut_F, int TwoKCut_I) {
        assert(NCut_I >= 0 && ECut_F >= 0.0 && TwoKCut_I >= 0);
        assert(Eqpn_F1D_qpn.size() == TwoKn_I1D_qpn.size() && Eqpp_F1D_qpp.size() == TwoKp_I1D_qpp.size());
        for (const auto& config_I1D_cqp : config_n2p0.first) {assert(config_I1D_cqp.size() == 2);}
        for (const auto& config_I1D_cqp : config_n0p2.second) {assert(config_I1D_cqp.size() == 2);}
        std::vector<std::vector<int>> confign_I2D_cfgnp_cqpn{};
        std::vector<std::vector<int>> configp_I2D_cfgnp_cqpp{};
        for (const auto& confign_I1D_cqpn : config_n2p0.first) {
            for (const auto& configp_I1D_cqpp : config_n0p2.second) {
                if (confign_I2D_cfgnp_cqpn.size() >= static_cast<std::size_t>(NCut_I)) {return std::make_pair(confign_I2D_cfgnp_cqpn, configp_I2D_cfgnp_cqpp);}
                const double Econfig_F = Eqpn_F1D_qpn(confign_I1D_cqpn[0]) + Eqpn_F1D_qpn(confign_I1D_cqpn[1]) + Eqpp_F1D_qpp(configp_I1D_cqpp[0]) + Eqpp_F1D_qpp(configp_I1D_cqpp[1]);
                const long long TwoKconfig_I = static_cast<long long>(TwoKn_I1D_qpn(confign_I1D_cqpn[0])) + TwoKn_I1D_qpn(confign_I1D_cqpn[1]) + TwoKp_I1D_qpp(configp_I1D_cqpp[0]) + TwoKp_I1D_qpp(configp_I1D_cqpp[1]);
                if (Econfig_F <= ECut_F && std::abs(TwoKconfig_I) <= TwoKCut_I) {
                    confign_I2D_cfgnp_cqpn.push_back(confign_I1D_cqpn);
                    configp_I2D_cfgnp_cqpp.push_back(configp_I1D_cqpp);
                }
            }
        }
        return std::make_pair(confign_I2D_cfgnp_cqpn, configp_I2D_cfgnp_cqpp);
    };

    // Row order → paired neutron/proton output order.
    confign_I2D_cfgnp_cqpn.clear();
    configp_I2D_cfgnp_cqpp.clear();
    for (Eigen::Index iNcqp_I = 0; iNcqp_I < Ncqp_I2D_Ncqp_np.rows(); ++iNcqp_I) {
        const int Ncqpn_I = Ncqp_I2D_Ncqp_np(iNcqp_I, 0);
        const int Ncqpp_I = Ncqp_I2D_Ncqp_np(iNcqp_I, 1);
        const int NCut_I = NCut_I1D_Ncqp(iNcqp_I);
        const double ECut_F = ECut_F1D_Ncqp(iNcqp_I);
        const int TwoKCut_I = TwoKCut_I1D_Ncqp(iNcqp_I);
        if (NCut_I == 0) {continue;}
        std::pair<std::vector<std::vector<int>>, std::vector<std::vector<int>>> config{};
        if (Ncqpn_I == 0 && Ncqpp_I == 0) {config = calc_config_kramers_n0p0(NCut_I, ECut_F, TwoKCut_I);}
        if (Ncqpn_I == 2 && Ncqpp_I == 0) {config = config_n2p0;}
        if (Ncqpn_I == 0 && Ncqpp_I == 2) {config = config_n0p2;}
        if (Ncqpn_I == 2 && Ncqpp_I == 2) {assert(iN2P0_I >= 0 && iN0P2_I >= 0); config = calc_config_kramers_n2p2(NCut_I, ECut_F, TwoKCut_I);}
        confign_I2D_cfgnp_cqpn.insert(confign_I2D_cfgnp_cqpn.end(), config.first.begin(), config.first.end());
        configp_I2D_cfgnp_cqpp.insert(configp_I2D_cfgnp_cqpp.end(), config.second.begin(), config.second.end());
    }
    build_tensor_config();
}

inline void HFBConfigurationNucleus::build_config_kramers_eo(const Eigen::VectorXd& Eqpn_F1D_qpn, const Eigen::VectorXi& TwoKn_I1D_qpn, const std::vector<int>& qpn_I1D_representative, const Eigen::VectorXd& Eqpp_F1D_qpp, const Eigen::VectorXi& TwoKp_I1D_qpp, const std::vector<int>& qpp_I1D_representative) {
    assert(Ncqp_I2D_Ncqp_np.cols() == 2 && Ncqp_I2D_Ncqp_np.rows() == ECut_F1D_Ncqp.size() && Ncqp_I2D_Ncqp_np.rows() == NCut_I1D_Ncqp.size() && Ncqp_I2D_Ncqp_np.rows() == TwoKCut_I1D_Ncqp.size());
    Eigen::Index iN2P0_I = -1;
    Eigen::Index iN0P1_I = -1;
    Eigen::Index iN2P1_I = -1;
    // (2,0): auxiliary even pairs; (0,1),(2,1): basis.
    for (Eigen::Index iNcqp_I = 0; iNcqp_I < Ncqp_I2D_Ncqp_np.rows(); ++iNcqp_I) {
        const int Ncqpn_I = Ncqp_I2D_Ncqp_np(iNcqp_I, 0);
        const int Ncqpp_I = Ncqp_I2D_Ncqp_np(iNcqp_I, 1);
        assert((Ncqpn_I == 2 && Ncqpp_I == 0) || (Ncqpn_I == 0 && Ncqpp_I == 1) || (Ncqpn_I == 2 && Ncqpp_I == 1));
        if (Ncqpn_I == 2 && Ncqpp_I == 0) {assert(iN2P0_I == -1); iN2P0_I = iNcqp_I;}
        if (Ncqpn_I == 0 && Ncqpp_I == 1) {assert(iN0P1_I == -1); iN0P1_I = iNcqp_I;}
        if (Ncqpn_I == 2 && Ncqpp_I == 1) {assert(iN2P1_I == -1); iN2P1_I = iNcqp_I;}
    }
    assert(iN0P1_I >= 0 && (iN2P1_I < 0 || iN2P0_I >= 0));
    assert((ECut_F1D_Ncqp.array() >= 0.0).all() && (NCut_I1D_Ncqp.array() >= 0).all() && (TwoKCut_I1D_Ncqp.array() >= 0).all());

    // Odd auxiliary pairs → relevant 1qp representatives, Nilsson order.
    const auto configp_I2D_pair_cqpp = calc_config_2qp(NCut_I1D_Ncqp(iN0P1_I), ECut_F1D_Ncqp(iN0P1_I), TwoKCut_I1D_Ncqp(iN0P1_I), Eqpp_F1D_qpp, TwoKp_I1D_qpp, qpp_I1D_representative, true);
    const int NqppHalf_I = static_cast<int>(Eqpp_F1D_qpp.size()) / 2;
    std::vector<int> qpp_I1D_used{};
    for (const int qpp_I : qpp_I1D_representative) {
        const bool used_B = std::any_of(configp_I2D_pair_cqpp.begin(), configp_I2D_pair_cqpp.end(), [&](const auto& config_I1D_cqp) {return config_I1D_cqp[0] % NqppHalf_I == qpp_I % NqppHalf_I || config_I1D_cqp[1] % NqppHalf_I == qpp_I % NqppHalf_I;});
        if (used_B) {qpp_I1D_used.push_back(qpp_I);}
    }
    std::vector<std::vector<int>> confign_I2D_pair_cqpn{};
    if (iN2P0_I >= 0) {confign_I2D_pair_cqpn = calc_config_2qp(NCut_I1D_Ncqp(iN2P0_I), ECut_F1D_Ncqp(iN2P0_I), TwoKCut_I1D_Ncqp(iN2P0_I), Eqpn_F1D_qpn, TwoKn_I1D_qpn, qpn_I1D_representative);}

    confign_I2D_cfgnp_cqpn.clear();
    configp_I2D_cfgnp_cqpp.clear();
    // Row order retained; (2,0) is not an odd-nucleus basis state.
    for (Eigen::Index iNcqp_I = 0; iNcqp_I < Ncqp_I2D_Ncqp_np.rows(); ++iNcqp_I) {
        if (iNcqp_I == iN0P1_I) {
            for (const int qpp_I : qpp_I1D_used) {
                confign_I2D_cfgnp_cqpn.emplace_back();
                configp_I2D_cfgnp_cqpp.push_back({qpp_I});
            }
        }
        if (iNcqp_I != iN2P1_I) {continue;}
        const int NCut_I = NCut_I1D_Ncqp(iNcqp_I);
        const double ECut_F = ECut_F1D_Ncqp(iNcqp_I);
        const int TwoKCut_I = TwoKCut_I1D_Ncqp(iNcqp_I);
        const int NpairHalf_I = static_cast<int>(confign_I2D_pair_cqpn.size()) / 2;
        const int NqpnHalf_I = static_cast<int>(Eqpn_F1D_qpn.size()) / 2;
        int Nconfig_I = 0;
        // BASE01: first NTQ(n)/2 pairs; representative then partner.
        for (int branch_I = 0; branch_I < 2; ++branch_I) {
            const int Nlimit_I = (branch_I + 1) * NCut_I / 2;
            for (int pair_I = 0; pair_I < NpairHalf_I && Nconfig_I < Nlimit_I; ++pair_I) {
                const auto& confign_I1D_cqpn = confign_I2D_pair_cqpn[pair_I];
                if (branch_I == 1 && confign_I1D_cqpn[0] % NqpnHalf_I == confign_I1D_cqpn[1] % NqpnHalf_I) {continue;}
                for (const int qpp_I : qpp_I1D_used) {
                    if (Nconfig_I >= Nlimit_I) {break;}
                    const int qppSelected_I = (qpp_I + branch_I * NqppHalf_I) % Eqpp_F1D_qpp.size();
                    const double Econfig_F = Eqpn_F1D_qpn(confign_I1D_cqpn[0]) + Eqpn_F1D_qpn(confign_I1D_cqpn[1]) + Eqpp_F1D_qpp(qpp_I);
                    const long long TwoKconfig_I = static_cast<long long>(TwoKn_I1D_qpn(confign_I1D_cqpn[0])) + TwoKn_I1D_qpn(confign_I1D_cqpn[1]) + (1 - 2 * branch_I) * TwoKp_I1D_qpp(qpp_I);
                    if (Econfig_F > ECut_F || std::abs(TwoKconfig_I) > TwoKCut_I) {continue;}
                    confign_I2D_cfgnp_cqpn.push_back(confign_I1D_cqpn);
                    configp_I2D_cfgnp_cqpp.push_back({qppSelected_I});
                    ++Nconfig_I;
                }
            }
        }
    }
    build_tensor_config();
}

inline void HFBConfigurationNucleus::build_config_kramers_oe(const Eigen::VectorXd& Eqpn_F1D_qpn, const Eigen::VectorXi& TwoKn_I1D_qpn, const std::vector<int>& qpn_I1D_representative, const Eigen::VectorXd& Eqpp_F1D_qpp, const Eigen::VectorXi& TwoKp_I1D_qpp, const std::vector<int>& qpp_I1D_representative) {
    // BASE10 = BASE01 with n ↔ p.
    Eigen::MatrixXi Ncqp_I2D_Ncqp_pn = Ncqp_I2D_Ncqp_np.rowwise().reverse();
    HFBConfigurationNucleus configuration_eo(Ncqp_I2D_Ncqp_pn, ECut_F1D_Ncqp, NCut_I1D_Ncqp, TwoKCut_I1D_Ncqp);
    configuration_eo.build_config_kramers_eo(Eqpp_F1D_qpp, TwoKp_I1D_qpp, qpp_I1D_representative, Eqpn_F1D_qpn, TwoKn_I1D_qpn, qpn_I1D_representative);
    confign_I2D_cfgnp_cqpn = configuration_eo.configp_I2D_cfgnp_cqpp;
    configp_I2D_cfgnp_cqpp = configuration_eo.confign_I2D_cfgnp_cqpn;
    build_tensor_config();
}

inline void HFBConfigurationNucleus::build_config_kramers_oo(const Eigen::VectorXd& Eqpn_F1D_qpn, const Eigen::VectorXi& TwoKn_I1D_qpn, const std::vector<int>& qpn_I1D_representative, const Eigen::VectorXd& Eqpp_F1D_qpp, const Eigen::VectorXi& TwoKp_I1D_qpp, const std::vector<int>& qpp_I1D_representative) {
    assert(Ncqp_I2D_Ncqp_np.cols() == 2 && Ncqp_I2D_Ncqp_np.rows() == ECut_F1D_Ncqp.size() && Ncqp_I2D_Ncqp_np.rows() == NCut_I1D_Ncqp.size() && Ncqp_I2D_Ncqp_np.rows() == TwoKCut_I1D_Ncqp.size());
    assert(Eqpn_F1D_qpn.size() > 0 && Eqpn_F1D_qpn.size() % 2 == 0 && TwoKn_I1D_qpn.size() == Eqpn_F1D_qpn.size());
    assert(Eqpp_F1D_qpp.size() > 0 && Eqpp_F1D_qpp.size() % 2 == 0 && TwoKp_I1D_qpp.size() == Eqpp_F1D_qpp.size());
    assert((ECut_F1D_Ncqp.array() >= 0.0).all() && (NCut_I1D_Ncqp.array() >= 0).all() && (TwoKCut_I1D_Ncqp.array() >= 0).all());
    Eigen::Index iN1P0_I = -1;
    Eigen::Index iN0P1_I = -1;
    Eigen::Index iN1P1_I = -1;
    // (1,0),(0,1): auxiliary windows; (1,1): basis.
    for (Eigen::Index iNcqp_I = 0; iNcqp_I < Ncqp_I2D_Ncqp_np.rows(); ++iNcqp_I) {
        const int Ncqpn_I = Ncqp_I2D_Ncqp_np(iNcqp_I, 0);
        const int Ncqpp_I = Ncqp_I2D_Ncqp_np(iNcqp_I, 1);
        assert((Ncqpn_I == 1 && Ncqpp_I == 0) || (Ncqpn_I == 0 && Ncqpp_I == 1) || (Ncqpn_I == 1 && Ncqpp_I == 1));
        if (Ncqpn_I == 1 && Ncqpp_I == 0) {assert(iN1P0_I == -1); iN1P0_I = iNcqp_I;}
        if (Ncqpn_I == 0 && Ncqpp_I == 1) {assert(iN0P1_I == -1); iN0P1_I = iNcqp_I;}
        if (Ncqpn_I == 1 && Ncqpp_I == 1) {assert(iN1P1_I == -1); iN1P1_I = iNcqp_I;}
    }
    assert(iN1P0_I >= 0 && iN0P1_I >= 0 && iN1P1_I >= 0);
    confign_I2D_cfgnp_cqpn.clear();
    configp_I2D_cfgnp_cqpp.clear();
    const int NCut_I = NCut_I1D_Ncqp(iN1P1_I);
    if (NCut_I == 0 || NCut_I1D_Ncqp(iN1P0_I) == 0 || NCut_I1D_Ncqp(iN0P1_I) == 0) {build_tensor_config(); return;}

    // BASE11: NH = single-species NCut, without halving.
    const auto calc_window = [](int Nwindow_I, const Eigen::VectorXd& Eqp_F1D_qp, const std::vector<int>& qp_I1D_representative) {
        const int Nallowed_I = static_cast<int>(qp_I1D_representative.size());
        assert(Nwindow_I > 0 && Nwindow_I <= Nallowed_I);
        for (const int qp_I : qp_I1D_representative) {assert(qp_I >= 0 && qp_I < Eqp_F1D_qp.size());}
        const auto min_iter = std::min_element(qp_I1D_representative.begin(), qp_I1D_representative.end(), [&](int qp1_I, int qp2_I) {return Eqp_F1D_qp(qp1_I) < Eqp_F1D_qp(qp2_I);});
        const int iMin_I = static_cast<int>(min_iter - qp_I1D_representative.begin());
        const int iBegin_I = std::min(std::max(0, iMin_I - Nwindow_I / 2), Nallowed_I - Nwindow_I);
        return std::make_pair(iBegin_I, iBegin_I + Nwindow_I);
    };
    const auto windown = calc_window(NCut_I1D_Ncqp(iN1P0_I), Eqpn_F1D_qpn, qpn_I1D_representative);
    const auto windowp = calc_window(NCut_I1D_Ncqp(iN0P1_I), Eqpp_F1D_qpp, qpp_I1D_representative);
    int Nconfig_I = 0;
    // Representative/partner branches; n outer, p inner.
    for (int branch_I = 0; branch_I < 2; ++branch_I) {
        const int Nlimit_I = (branch_I + 1) * NCut_I / 2;
        std::vector<double> Econfig_F1D_cfg{};
        std::pair<std::vector<std::vector<int>>, std::vector<std::vector<int>>> candidate{};
        for (int in_I = windown.first; in_I < windown.second && Nconfig_I < Nlimit_I; ++in_I) {
            const int qpn_I = qpn_I1D_representative[in_I];
            if (Eqpn_F1D_qpn(qpn_I) > ECut_F1D_Ncqp(iN1P0_I) || std::abs(TwoKn_I1D_qpn(qpn_I)) > TwoKCut_I1D_Ncqp(iN1P0_I)) {continue;}
            for (int ip_I = windowp.first; ip_I < windowp.second && Nconfig_I < Nlimit_I; ++ip_I) {
                const int qpp_I = qpp_I1D_representative[ip_I];
                if (Eqpp_F1D_qpp(qpp_I) > ECut_F1D_Ncqp(iN0P1_I) || std::abs(TwoKp_I1D_qpp(qpp_I)) > TwoKCut_I1D_Ncqp(iN0P1_I)) {continue;}
                const double Econfig_F = Eqpn_F1D_qpn(qpn_I) + Eqpp_F1D_qpp(qpp_I);
                const long long TwoKconfig_I = static_cast<long long>(TwoKn_I1D_qpn(qpn_I)) + (1 - 2 * branch_I) * TwoKp_I1D_qpp(qpp_I);
                if (Econfig_F > ECut_F1D_Ncqp(iN1P1_I) || std::abs(TwoKconfig_I) > TwoKCut_I1D_Ncqp(iN1P1_I)) {continue;}
                const int qppSelected_I = (qpp_I + branch_I * (Eqpp_F1D_qpp.size() / 2)) % Eqpp_F1D_qpp.size();
                Econfig_F1D_cfg.push_back(Econfig_F);
                candidate.first.push_back({qpn_I});
                candidate.second.push_back({qppSelected_I});
                ++Nconfig_I;
            }
        }
        // ORD2QP: E,n,p swap together; branches sort separately.
        const std::size_t Ncfg_I = Econfig_F1D_cfg.size();
        for (std::size_t cfg_I = 0; cfg_I < Ncfg_I; ++cfg_I) {
            std::size_t cfgMin_I = cfg_I;
            for (std::size_t cfgNext_I = cfg_I + 1; cfgNext_I < Ncfg_I; ++cfgNext_I) {
                if (Econfig_F1D_cfg[cfgNext_I] < Econfig_F1D_cfg[cfgMin_I]) {cfgMin_I = cfgNext_I;}
            }
            std::swap(Econfig_F1D_cfg[cfg_I], Econfig_F1D_cfg[cfgMin_I]);
            std::swap(candidate.first[cfg_I], candidate.first[cfgMin_I]);
            std::swap(candidate.second[cfg_I], candidate.second[cfgMin_I]);
            confign_I2D_cfgnp_cqpn.push_back(candidate.first[cfg_I]);
            configp_I2D_cfgnp_cqpp.push_back(candidate.second[cfg_I]);
        }
    }
    build_tensor_config();
}

inline std::vector<std::vector<int>> HFBConfigurationNucleus::calc_config_2qp(int NCut_I, double ECut_F, int TwoKCut_I, const Eigen::VectorXd& Eqp_F1D_qp, const Eigen::VectorXi& TwoK_I1D_qp, const std::vector<int>& qp_I1D_representative, bool cutSingle_B) {
    assert(NCut_I >= 0 && ECut_F >= 0.0 && TwoKCut_I >= 0);
    assert(Eqp_F1D_qp.size() > 0 && Eqp_F1D_qp.size() % 2 == 0 && Eqp_F1D_qp.size() == TwoK_I1D_qp.size());
    for (const int qp_I : qp_I1D_representative) {assert(qp_I >= 0 && qp_I < Eqp_F1D_qp.size());}

    // ORD2QP: E and qp indices swap together.
    const auto sort_append_config = [](std::vector<std::vector<int>>& result_I2D_cfg_cqp, std::vector<double>& Econfig_F1D_cfg, std::vector<std::vector<int>>& candidate_I2D_cfg_cqp) {
        const std::size_t Ncfg_I = candidate_I2D_cfg_cqp.size();
        for (std::size_t cfg_I = 0; cfg_I < Ncfg_I; ++cfg_I) {
            std::size_t cfgMin_I = cfg_I;
            for (std::size_t cfgNext_I = cfg_I + 1; cfgNext_I < Ncfg_I; ++cfgNext_I) {
                if (Econfig_F1D_cfg[cfgNext_I] < Econfig_F1D_cfg[cfgMin_I]) {cfgMin_I = cfgNext_I;}
            }
            std::swap(Econfig_F1D_cfg[cfg_I], Econfig_F1D_cfg[cfgMin_I]);
            std::swap(candidate_I2D_cfg_cqp[cfg_I], candidate_I2D_cfg_cqp[cfgMin_I]);
            result_I2D_cfg_cqp.push_back(candidate_I2D_cfg_cqp[cfg_I]);
        }
    };

    std::vector<std::vector<int>> result_I2D_cfg_cqp{};
    if (NCut_I == 0) {return result_I2D_cfg_cqp;}
    const int Nqp_I = static_cast<int>(Eqp_F1D_qp.size());
    // NH = NCut/2; MI = max(0,LOW−NH/2); MF = MI+NH.
    const int NCutHalf_I = NCut_I / 2;
    const int Nallowed_I = static_cast<int>(qp_I1D_representative.size());
    assert(NCutHalf_I > 0 && NCutHalf_I <= Nallowed_I);
    const auto min_iter = std::min_element(qp_I1D_representative.begin(), qp_I1D_representative.end(), [&](int qp1_I, int qp2_I) {return Eqp_F1D_qp(qp1_I) < Eqp_F1D_qp(qp2_I);});
    const int iMin_I = static_cast<int>(min_iter - qp_I1D_representative.begin());
    const int iBegin_I = std::min(std::max(0, iMin_I - NCutHalf_I / 2), Nallowed_I - NCutHalf_I);
    const int iEnd_I = iBegin_I + NCutHalf_I;
    std::vector<double> EconfigSame_F1D_cfg{};
    std::vector<double> EconfigPartner_F1D_cfg{};
    std::vector<std::vector<int>> configSame_I2D_cfg_cqp{};
    std::vector<std::vector<int>> configPartner_I2D_cfg_cqp{};

    // βi†βj†: i<j; retain first NH accepted pairs.
    for (int i1_I = iBegin_I; i1_I < iEnd_I && configSame_I2D_cfg_cqp.size() < static_cast<std::size_t>(NCutHalf_I); ++i1_I) {
        for (int i2_I = i1_I + 1; i2_I < iEnd_I && configSame_I2D_cfg_cqp.size() < static_cast<std::size_t>(NCutHalf_I); ++i2_I) {
            const int qp1_I = qp_I1D_representative[i1_I];
            const int qp2_I = qp_I1D_representative[i2_I];
            const double Econfig_F = Eqp_F1D_qp(qp1_I) + Eqp_F1D_qp(qp2_I);
            const long long TwoKconfig_I = static_cast<long long>(TwoK_I1D_qp(qp1_I)) + TwoK_I1D_qp(qp2_I);
            // Odd species: individual cuts; even species: pair cuts.
            bool accept_B = Econfig_F <= ECut_F && std::abs(TwoKconfig_I) <= TwoKCut_I;
            if (cutSingle_B) {accept_B = Eqp_F1D_qp(qp1_I) <= ECut_F && Eqp_F1D_qp(qp2_I) <= ECut_F && std::abs(TwoK_I1D_qp(qp1_I)) <= TwoKCut_I && std::abs(TwoK_I1D_qp(qp2_I)) <= TwoKCut_I;}
            if (accept_B) {
                EconfigSame_F1D_cfg.push_back(Econfig_F);
                configSame_I2D_cfg_cqp.push_back({std::min(qp1_I, qp2_I), std::max(qp1_I, qp2_I)});
            }
        }
    }

    // βi†βj̄†: i≤j; retain up to NCut total pairs.
    for (int i1_I = iBegin_I; i1_I < iEnd_I && configSame_I2D_cfg_cqp.size() + configPartner_I2D_cfg_cqp.size() < static_cast<std::size_t>(NCut_I); ++i1_I) {
        for (int i2_I = i1_I; i2_I < iEnd_I && configSame_I2D_cfg_cqp.size() + configPartner_I2D_cfg_cqp.size() < static_cast<std::size_t>(NCut_I); ++i2_I) {
            const int qp1_I = qp_I1D_representative[i1_I];
            const int qp2_I = qp_I1D_representative[i2_I];
            const int qp2Partner_I = (qp2_I + Nqp_I / 2) % Nqp_I;
            const double Econfig_F = Eqp_F1D_qp(qp1_I) + Eqp_F1D_qp(qp2_I);
            const long long TwoKconfig_I = static_cast<long long>(TwoK_I1D_qp(qp1_I)) - TwoK_I1D_qp(qp2_I);
            // Odd species: individual cuts; even species: pair cuts.
            bool accept_B = Econfig_F <= ECut_F && std::abs(TwoKconfig_I) <= TwoKCut_I;
            if (cutSingle_B) {accept_B = Eqp_F1D_qp(qp1_I) <= ECut_F && Eqp_F1D_qp(qp2_I) <= ECut_F && std::abs(TwoK_I1D_qp(qp1_I)) <= TwoKCut_I && std::abs(TwoK_I1D_qp(qp2_I)) <= TwoKCut_I;}
            if (accept_B) {
                EconfigPartner_F1D_cfg.push_back(Econfig_F);
                configPartner_I2D_cfg_cqp.push_back({std::min(qp1_I, qp2Partner_I), std::max(qp1_I, qp2Partner_I)});
            }
        }
    }

    sort_append_config(result_I2D_cfg_cqp, EconfigSame_F1D_cfg, configSame_I2D_cfg_cqp);
    sort_append_config(result_I2D_cfg_cqp, EconfigPartner_F1D_cfg, configPartner_I2D_cfg_cqp);
    return result_I2D_cfg_cqp;

}
