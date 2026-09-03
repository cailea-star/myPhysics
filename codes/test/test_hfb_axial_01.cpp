/**
 * @file    test_hfb_axial_01.cpp
 * @author  cailea
 * @date    2026-09-02
 * @brief   Validate Ca47 neutron-hole blocking.
 */

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "hfb_axial.hpp"
#include "hfb_axial_observable.hpp"
#include "hfb_edf_skyrme.hpp"
#include "test_hfb_axial_reference.hpp"

/**
 * @brief Validate Ca47 neutron-hole observables against HFBTHO.
 * @math ⁴⁸Ca → ⁴⁷Ca(n⁻¹)
 * @output Passed reference assertions.
 */
int main() {
    std::cout << std::scientific << std::setprecision(12) << std::right;

    // (b₀,N_shell,N_r,N_z) ← run_000001.
    const int Nshell_I = 10;
    const double b0_F = 1.75000371403579;
    AxialConfig axialconfig_(b0_F, b0_F, Nshell_I, true);
    axialconfig_.Nr_I = 40;
    axialconfig_.Nz_I = 40;

    // ⁴⁸Ca → ⁴⁷Ca(n⁻¹).
    const int Ncore_I = 28;
    const int Zcore_I = 20;
    const int Nfinal_I = 27;
    const int Zfinal_I = 20;

    // P_HFB ← run_000001.
    HFBSettings hfbsettings_;
    hfbsettings_.useLipkinNogami_B = false;
    hfbsettings_.useCmCorrection_B = true;
    hfbsettings_.Nblocking_I = 6;
    hfbsettings_.EspCut_F = 60.0;
    hfbsettings_.accuracy_F = 1.0e-5;

    std::cout << "[Input] N, Z = " << Nfinal_I << ", " << Zfinal_I << "\n";
    std::cout << "[Input] b_0, N_shell = " << b0_F << ", " << Nshell_I << "\n";

    // ⁴⁸Ca → {B_μ | μ=1,…,6}.
    EDFParamsSkyrme edf_skyrme_ = HFBfunctionals::SKMstar();
    AxialHFB hfb_(axialconfig_, hfbsettings_, edf_skyrme_);
    hfb_.initialize_WS_field(Ncore_I, Zcore_I);
    hfb_.iterate(Ncore_I, Zcore_I);
    const AxialHFB hfbCore_ = hfb_;
    std::vector<AxialHFBBlocking> blockings_S1D_candidate = AxialHFBBlocking::list_candidates(hfbsettings_, hfb_.blocklist_n, true);

    // O_ref(μ=1,…,6) ← HFBTHO run_000001.
    const std::vector<AxialHFBObservableReference> references_S1D_candidate{{-409.353490, -10.181604, -13.510115, -0.000000, -0.000000, 0.000070, -0.000014, 3.573004, 3.442663}, {-409.309844, -10.180494, -13.628915, -0.000000, -0.000000, 0.000054, -0.000010, 3.572602, 3.442358}, {-409.288385, -10.222255, -13.743649, -0.000000, -0.000000, 0.000059, -0.000009, 3.572603, 3.442320}, {-409.437538, -10.222038, -13.407578, -0.000000, -0.000000, 0.000071, -0.000023, 3.573810, 3.443265}, {-407.779143, -9.536332, -13.272025, -5.809006, -0.000000, 1.310967, -0.000024, 3.589751, 3.448737}, {-407.776599, -9.536052, -13.300916, -5.806758, -0.000000, 1.310750, -0.000023, 3.589762, 3.448738}};
    assert(blockings_S1D_candidate.size() == references_S1D_candidate.size());

    for (int candidate_I = 0; candidate_I < static_cast<int>(blockings_S1D_candidate.size()); ++candidate_I) {
        // (HFB_core,B_μ) → O_μ.
        AxialHFB hfbBlocked_ = hfbCore_;
        AxialHFBBlocking& activeBlocking_ = blockings_S1D_candidate[candidate_I];
        hfbBlocked_.iterate(Nfinal_I, Zfinal_I, activeBlocking_);
        AxialHFBObservable observable_;
        observable_.update_observable(hfbBlocked_, activeBlocking_);
        const AxialHFBObservableReference& reference_ = references_S1D_candidate[candidate_I];

        // ε = max|O-O_ref| by units.
        const double energyError_F = std::max({std::abs(observable_.Etot_F - reference_.Etot_F), std::abs(observable_.Epair_n_F - reference_.Epair_n_F), std::abs(observable_.Epair_p_F - reference_.Epair_p_F)});
        const double lambdaError_F = std::max(std::abs(observable_.lambda_n_F - reference_.lambda_n_F), std::abs(observable_.lambda_p_F - reference_.lambda_p_F));
        const double DeltaError_F = std::max(std::abs(observable_.Delta_n_F - reference_.Delta_n_F), std::abs(observable_.Delta_p_F - reference_.Delta_p_F));
        const double rmsError_F = std::max(std::abs(observable_.rms_n_F - reference_.rms_n_F), std::abs(observable_.rms_p_F - reference_.rms_p_F));

        std::cout << "[Input] blocking candidate = " << candidate_I + 1 << "\n";
        std::cout << "[Reference] Etot, lambda_n, lambda_p, Epair_n, Epair_p, Delta_n, Delta_p, rms_n, rms_p = " << reference_.Etot_F << ", " << reference_.lambda_n_F << ", " << reference_.lambda_p_F << ", " << reference_.Epair_n_F << ", " << reference_.Epair_p_F << ", " << reference_.Delta_n_F << ", " << reference_.Delta_p_F << ", " << reference_.rms_n_F << ", " << reference_.rms_p_F << "\n";
        std::cout << "[Computed] Etot, lambda_n, lambda_p, Epair_n, Epair_p, Delta_n, Delta_p, rms_n, rms_p = " << observable_.Etot_F << ", " << observable_.lambda_n_F << ", " << observable_.lambda_p_F << ", " << observable_.Epair_n_F << ", " << observable_.Epair_p_F << ", " << observable_.Delta_n_F << ", " << observable_.Delta_p_F << ", " << observable_.rms_n_F << ", " << observable_.rms_p_F << "\n";
        std::cout << "[Error] energy, lambda, Delta, rms = " << energyError_F << ", " << lambdaError_F << ", " << DeltaError_F << ", " << rmsError_F << "\n";

        // (ε_E,ε_λ,ε_Δ,ε_rms) ≤ (10⁻³,10⁻³,10⁻³,10⁻⁴).
        assert(energyError_F <= 1.0e-3);
        assert(lambdaError_F <= 1.0e-3);
        assert(DeltaError_F <= 1.0e-3);
        assert(rmsError_F <= 1.0e-4);
    }

    return 0;
}
