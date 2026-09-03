/**
 * @file    test_hfb_axial_02.cpp
 * @author  cailea
 * @date    2026-09-02
 * @brief   Validate Ca49 neutron-particle blocking.
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
 * @brief Validate Ca49 neutron-particle observables against HFBTHO.
 * @math ^{48}Ca → ^{49}Ca(n)
 * @output Passed reference assertions.
 */
int main() {
    std::cout << std::scientific << std::setprecision(12) << std::right;

    // (b₀,N_shell,N_r,N_z) ← run_000002.
    const int Nshell_I = 10;
    const double b0_F = 1.75000371403579;
    AxialConfig axialconfig_(b0_F, b0_F, Nshell_I, true);
    axialconfig_.Nr_I = 40;
    axialconfig_.Nz_I = 40;

    // ⁴⁸Ca → ⁴⁹Ca(n).
    const int Ncore_I = 28;
    const int Zcore_I = 20;
    const int Nfinal_I = 29;
    const int Zfinal_I = 20;

    // P_HFB ← run_000002.
    HFBSettings hfbsettings_;
    hfbsettings_.useLipkinNogami_B = false;
    hfbsettings_.useCmCorrection_B = true;
    hfbsettings_.Nblocking_I = 6;
    hfbsettings_.EspCut_F = 60.0;
    hfbsettings_.accuracy_F = 1.0e-5;

    std::cout << "[Input] N, Z = " << Nfinal_I << ", " << Zfinal_I << "\n";
    std::cout << "[Input] b_0, N_shell = " << b0_F << ", " << Nshell_I << "\n";
    std::cout << "[Input] blocking candidate = 5\n";

    // ⁴⁸Ca → B₅ → ⁴⁹Ca.
    EDFParamsSkyrme edf_skyrme_ = HFBfunctionals::SKMstar();
    AxialHFB hfb_(axialconfig_, hfbsettings_, edf_skyrme_);
    hfb_.initialize_WS_field(Ncore_I, Zcore_I);
    hfb_.iterate(Ncore_I, Zcore_I);
    std::vector<AxialHFBBlocking> blockings_S1D_candidate = AxialHFBBlocking::list_candidates(hfbsettings_, hfb_.blocklist_n, true);
    assert(blockings_S1D_candidate.size() > 4);
    AxialHFBBlocking activeBlocking_ = blockings_S1D_candidate[4];
    hfb_.iterate(Nfinal_I, Zfinal_I, activeBlocking_);
    AxialHFBObservable observable_;
    observable_.update_observable(hfb_, activeBlocking_);

    // O_ref ← HFBTHO run_000002.
    const AxialHFBObservableReference reference_{-425.463058, -10.343656, -15.028470, -0.000000, -0.000000, 0.000047, -0.000003, 3.632518, 3.455211};

    // ε = max|O-O_ref| by units.
    const double energyError_F = std::max({std::abs(observable_.Etot_F - reference_.Etot_F), std::abs(observable_.Epair_n_F - reference_.Epair_n_F), std::abs(observable_.Epair_p_F - reference_.Epair_p_F)});
    const double lambdaError_F = std::max(std::abs(observable_.lambda_n_F - reference_.lambda_n_F), std::abs(observable_.lambda_p_F - reference_.lambda_p_F));
    const double DeltaError_F = std::max(std::abs(observable_.Delta_n_F - reference_.Delta_n_F), std::abs(observable_.Delta_p_F - reference_.Delta_p_F));
    const double rmsError_F = std::max(std::abs(observable_.rms_n_F - reference_.rms_n_F), std::abs(observable_.rms_p_F - reference_.rms_p_F));

    std::cout << "[Reference] Etot, lambda_n, lambda_p, Epair_n, Epair_p, Delta_n, Delta_p, rms_n, rms_p = " << reference_.Etot_F << ", " << reference_.lambda_n_F << ", " << reference_.lambda_p_F << ", " << reference_.Epair_n_F << ", " << reference_.Epair_p_F << ", " << reference_.Delta_n_F << ", " << reference_.Delta_p_F << ", " << reference_.rms_n_F << ", " << reference_.rms_p_F << "\n";
    std::cout << "[Computed] Etot, lambda_n, lambda_p, Epair_n, Epair_p, Delta_n, Delta_p, rms_n, rms_p = " << observable_.Etot_F << ", " << observable_.lambda_n_F << ", " << observable_.lambda_p_F << ", " << observable_.Epair_n_F << ", " << observable_.Epair_p_F << ", " << observable_.Delta_n_F << ", " << observable_.Delta_p_F << ", " << observable_.rms_n_F << ", " << observable_.rms_p_F << "\n";
    std::cout << "[Error] energy, lambda, Delta, rms = " << energyError_F << ", " << lambdaError_F << ", " << DeltaError_F << ", " << rmsError_F << "\n";

    // (ε_E,ε_λ,ε_Δ,ε_rms) ≤ (10⁻³,10⁻³,10⁻³,10⁻⁴).
    assert(energyError_F <= 1.0e-3);
    assert(lambdaError_F <= 1.0e-3);
    assert(DeltaError_F <= 1.0e-3);
    assert(rmsError_F <= 1.0e-4);

    return 0;
}
