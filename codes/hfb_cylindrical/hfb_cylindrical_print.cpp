/**
 * @file    hfb_cylindrical_print.cpp
 * @author  cailea
 * @date    2026-09-12
 * @brief   Print HFB iteration summaries and blocking states.
 */

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "hfb_cylindrical_observable.hpp"

namespace {

/**
 * @brief Print one HFB iteration summary.
 * @math (i,ε,α,O) → stdout
 * @output Formatted iteration row.
 */
void print_summary(const HFBCylindricalObservable& observable_, const HFBEDFSetting& hfbsettings_, int iteration_I, double residual_F, double mixing_F) {
    if (iteration_I == 0) {
        std::cout << "  ----------------------------------------------------------------------------------------------------------------------------------\n";
        std::cout << "  i          si     mix   beta      Etot       A       rn      rp         En      Dn      Ep      Dp        Ln      Lp\n";
        std::cout << "  ----------------------------------------------------------------------------------------------------------------------------------" << std::endl;
    }

    // Δ_q ← |Δ_q|+λ_{2,q}.
    double Delta_n_F = std::abs(observable_.Delta_n_F);
    double Delta_p_F = std::abs(observable_.Delta_p_F);
    if (hfbsettings_.useLipkinNogami_B) {
        Delta_n_F += observable_.lambda2_n_F;
        Delta_p_F += observable_.lambda2_p_F;
    }

    // (i,ε,α,O) → stdout.
    std::cout << std::fixed;
    std::cout << std::setw(4) << iteration_I << (iteration_I > 0 ? "B" : "L") << " ";
    std::cout << std::setw(12) << std::setprecision(8) << residual_F;
    std::cout << std::setw(5) << std::setprecision(2) << mixing_F;
    std::cout << std::setw(7) << std::setprecision(3) << observable_.beta2_F;
    std::cout << std::setw(13) << std::setprecision(6) << observable_.Etot_F << " ";
    std::cout << std::setw(6) << std::setprecision(1) << observable_.A_F;
    std::cout << std::setw(8) << std::setprecision(3) << observable_.rms_n_F;
    std::cout << std::setw(8) << std::setprecision(3) << observable_.rms_p_F << " | ";
    std::cout << std::setw(8) << std::setprecision(3) << observable_.Epair_n_F;
    std::cout << std::setw(8) << std::setprecision(3) << Delta_n_F << " | ";
    std::cout << std::setw(8) << std::setprecision(3) << observable_.Epair_p_F;
    std::cout << std::setw(8) << std::setprecision(3) << Delta_p_F << " | ";
    std::cout << std::setw(8) << std::setprecision(3) << observable_.lambda_n_F;
    std::cout << std::setw(8) << std::setprecision(3) << observable_.lambda_p_F;
    std::cout << std::endl;
}

/**
 * @brief Print active quasiparticle blocking states.
 * @math {B_μ} → stdout
 * @output Formatted blocking rows.
 */
void print_blocking(const CylindricalSetting& cylindricalsetting_, const std::vector<HFBKramersBlocking>& activeBlockings_, bool isNeutron_B) {
    for (const HFBKramersBlocking& blocking_ : activeBlockings_) {
        if (blocking_.block_I < 0 || blocking_.block_I >= static_cast<int>(cylindricalsetting_.labels_S2D_block_bsp.size())) {continue;}
        if (blocking_.bqp_I < 0) {continue;}
        const auto& labels_S1D_bsp = cylindricalsetting_.labels_S2D_block_bsp[blocking_.block_I];
        if (blocking_.blockedU_F1D_bsp.size() != labels_S1D_bsp.size()) {continue;}
        if (blocking_.blockedV_F1D_bsp.size() != labels_S1D_bsp.size()) {continue;}

        // max(|U|,|V|) → bsp_max.
        Eigen::Index bspMax_I = 0;
        (blocking_.blockedU_F1D_bsp.cwiseAbs().cwiseMax(blocking_.blockedV_F1D_bsp.cwiseAbs())).maxCoeff(&bspMax_I);
        const CylindricalSPLabel& label_ = labels_S1D_bsp[bspMax_I];

        // (q,block,bqp,overlap,α_bsp) → stdout.
        std::cout << "      blocking " << (isNeutron_B ? "n" : "p");
        std::cout << " block=" << std::setw(3) << (blocking_.block_I + 1);
        std::cout << " qp=" << std::setw(3) << (blocking_.bqp_I + 1);
        std::cout << " ovlp=" << std::setw(8) << std::setprecision(4) << blocking_.overlap_F;
        std::cout << " label=" << label_.twoOmega_I << "/2";
        std::cout << (label_.isParityPositive_B ? "+" : "-") << "[";
        std::cout << label_.N_I << ",";
        std::cout << label_.nz_I << ",";
        std::cout << label_.Lambda_I << "]" << std::endl;
    }
}

} // namespace

void HFBKramersNucleusCylindrical::print_abstract(int iteration_I, double error_F, double mixing_F) const {
    print_abstract(iteration_I, error_F, mixing_F, {}, {});
}

void HFBKramersNucleusCylindrical::print_abstract(int iteration_I, double error_F, double mixing_F, const std::vector<HFBKramersBlocking>& neutronBlockings_, const std::vector<HFBKramersBlocking>& protonBlockings_) const {
    HFBCylindricalObservable observable_;
    observable_.update_observable(*this, neutronBlockings_, protonBlockings_);
    // Unknown blocked indices: report the actual chemical potential.
    if (hfb_neutron.blocking_Func && neutronBlockings_.empty()) {observable_.lambda_n_F = hfb_neutron.lambda_F;}
    if (hfb_proton.blocking_Func && protonBlockings_.empty()) {observable_.lambda_p_F = hfb_proton.lambda_F;}
    print_summary(observable_, hfbedfsetting, iteration_I, error_F, mixing_F);
    print_blocking(cylindricalsetting, neutronBlockings_, true);
    print_blocking(cylindricalsetting, protonBlockings_, false);
}
