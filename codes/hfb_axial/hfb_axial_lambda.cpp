/**
 * @file    hfb_axial_lambda.cpp
 * @author  cailea
 * @date    2026-09-03
 * @brief   Update axial HFB chemical potentials.
 */

#include "hfb_axial.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <map>

#include "root.hpp"

namespace {

/**
 * @brief  Calculate trial particle number from both density sectors.
 * @math   N(λ)=Σ_b[Tr ρ_b⁺⁺+Tr ρ_b⁻⁻]
 * @output Trial block amplitudes, densities, and particle number.
 * @note   Uses block workspaces and the current blocking implementation.
 */
double calc_N_at_lambda(double lambda_F, AxialHFBBlockList& blocklist_, const HFBSettings& hfbsettings_, const std::vector<AxialHFBBlocking>& activeBlockings_) {
    // λ → ℋ⁺ → {U⁺,V⁻,ρ⁺⁺,ρ⁻⁻,κ⁺⁻,κ⁻⁺}; reuse block workspaces.
    blocklist_.lambda_F = lambda_F;
    blocklist_.update_UV_E_rho_kappa(hfbsettings_);

    // Trial blocking leaves the persistent state trackers unchanged.
    for (const AxialHFBBlocking& blocking_ : activeBlockings_) {
        if (blocking_.isNeutron_B != blocklist_.isNeutron_B) {continue;}
        AxialHFBBlocking trialBlocking_ = blocking_;
        trialBlocking_.apply_blocking(blocklist_);
    }

    // N=Σ_b[Tr ρ_b⁺⁺+Tr ρ_b⁻⁻].
    double Ncalc_F = 0.0;
    for (const AxialHFBBlock& block_ : blocklist_.blocks_X1D_block) {
        Ncalc_F += block_.rhoPosPos_F2D_bsp_bsp.trace() + block_.rhoNegNeg_F2D_bsp_bsp.trace();
    }
    return Ncalc_F;
}
} // namespace

void AxialHFB::update_blocklist_lambda(AxialHFBBlockList& blocklist_, int Ntarget_I, std::vector<AxialHFBBlocking>& activeBlockings_, double lambdaTolerance_F) {
    const double Ntolerance_F = std::max(1.0e-8, 1.0e-10 * std::max(1, Ntarget_I));
    std::map<double, double> NerrorByLambda_Map{};

    // {B_μ^q} → (ρ_q,κ_q).
    const auto apply_active_blockings_Func = [&]() {
        for (AxialHFBBlocking& blocking_ : activeBlockings_) {
            if (blocking_.isNeutron_B != blocklist_.isNeutron_B) {continue;}
            blocking_.apply_blocking(blocklist_);
        }
    };

    // λ → N_q(λ)-N_q^{target}.
    const auto calc_Nerror_Func = [&](double lambda_F) {
        const auto NerrorIt_ = NerrorByLambda_Map.find(lambda_F);
        if (NerrorIt_ != NerrorByLambda_Map.end()) {return NerrorIt_->second;}
        const double Ncalc_F = calc_N_at_lambda(lambda_F, blocklist_, hfbsettings, activeBlockings_);
        const double Nerror_F = Ncalc_F - static_cast<double>(Ntarget_I);
        NerrorByLambda_Map.emplace(lambda_F, Nerror_F);
        return Nerror_F;
    };

    // λ_0 ← λ_previous or -7 MeV.
    assert(std::isfinite(blocklist_.lambda_F));
    const double lambdaInitial_F = blocklist_.lambda_F;
    const double NerrorInitial_F = calc_Nerror_Func(lambdaInitial_F);
    if (std::abs(NerrorInitial_F) <= Ntolerance_F) {
        blocklist_.lambda_F = lambdaInitial_F;
        blocklist_.update_UV_E_rho_kappa(hfbsettings);
        apply_active_blockings_Func();
        return;
    }

    // [λ_min,λ_max] ← [λ_0-2,λ_0+2].
    double lambdaMax_F = lambdaInitial_F + 2.0;
    double lambdaMin_F = lambdaInitial_F - 2.0;
    double NerrorMin_F = calc_Nerror_Func(lambdaMin_F);
    double NerrorMax_F = calc_Nerror_Func(lambdaMax_F);
    const double NerrorSlope_F = (NerrorMax_F - NerrorMin_F) / (lambdaMax_F - lambdaMin_F);

    // f(λ_min)f(λ_max)>0 → expand bracket.
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

    // f(λ_q)=0 → (U,V,E,ρ,κ)_q.
    const double lambdaRoot_F = root_brent(calc_Nerror_Func, lambdaMin_F, lambdaMax_F, lambdaTolerance_F);
    blocklist_.lambda_F = lambdaRoot_F;
    blocklist_.update_UV_E_rho_kappa(hfbsettings);
    apply_active_blockings_Func();
}
