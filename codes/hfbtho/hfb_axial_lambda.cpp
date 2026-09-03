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

#include "root.hpp"

void AxialHFB::update_blocklist_lambda(AxialHFBBlockList& blocklist_, int Ntarget_I, std::vector<AxialHFBBlocking>& activeBlockings_, bool isNeutron_B, double lambdaTolerance_F) {
    const double Ntolerance_F = std::max(1.0e-8, 1.0e-10 * std::max(1, Ntarget_I));

    // {B_μ^q} → (ρ_q,κ_q).
    const auto apply_active_blockings_Func = [&]() {
        for (AxialHFBBlocking& blocking_ : activeBlockings_) {
            if (blocking_.isNeutron_B != isNeutron_B) {continue;}
            blocking_.apply_blocking(blocklist_);
        }
    };

    // λ → N_q(λ)-N_q^{target}.
    const auto calc_Nerror_Func = [&](double lambda_F) {
        blocklist_.lambda_F = lambda_F;
        blocklist_.update_UV_E_rho_kappa(hfbsettings);
        for (const AxialHFBBlocking& blocking_ : activeBlockings_) {
            if (blocking_.isNeutron_B != isNeutron_B) {continue;}
            AxialHFBBlocking blockingTrial_ = blocking_;
            blockingTrial_.apply_blocking(blocklist_);
        }

        double Ncalc_F = 0.0;
        for (const AxialHFBBlock& block_ : blocklist_.blocks_X1D_block) {
            Ncalc_F += 2.0 * block_.rho_F2D_bsp_bsp.trace();
        }
        return Ncalc_F - static_cast<double>(Ntarget_I);
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
    while (NerrorMin_F * NerrorMax_F > 0.0 && NerrorSlope_F != 0.0) {
        lambdaBound_F += lambdaStep_F;
        NerrorBound_F = calc_Nerror_Func(lambdaBound_F);
    }

    // f(λ_q)=0 → (U,V,E,ρ,κ)_q.
    const double lambdaRoot_F = root_brent(calc_Nerror_Func, lambdaMin_F, lambdaMax_F, lambdaTolerance_F);
    blocklist_.lambda_F = lambdaRoot_F;
    blocklist_.update_UV_E_rho_kappa(hfbsettings);
    apply_active_blockings_Func();
}
