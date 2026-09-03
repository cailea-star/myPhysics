/**
 * @file    test_hfb_axial_10.cpp
 * @author  cailea
 * @date    2026-09-02
 * @brief   Validate blocked SkM* HFB against run 000010.
 */

#include "hfb_axial.hpp"
#include "hfb_axial_observable.hpp"
#include "hfb_edf_skyrme.hpp"
#include "test_hfb_axial_reference.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

/**
 * @brief  Validate neutron-proton blocked SkM* observables.
 * @math   O_{C++}^{48K} ↔ O_{000010}^{48K}
 * @output Zero on accepted observable errors.
 */
int main() {
    // (N_{shell},b_0,N_r,N_z) ← run 000010.
    const int Nshell_I = 12;
    const double b0_F = 1.75000371403579;
    AxialConfig axialconfig_(b0_F, b0_F, Nshell_I, true);
    axialconfig_.Nr_I = 40;
    axialconfig_.Nz_I = 40;

    // ^{48}Ca(28,20) → ^{48}K(29,19).
    const int Ncore_I = 28;
    const int Zcore_I = 20;
    const int Nfinal_I = 29;
    const int Zfinal_I = 19;

    // P_{HFB} ← run 000010.
    HFBSettings hfbsettings_ = HFBSettings::setting_skyrme();
    hfbsettings_.Nblocking_I = 7;
    hfbsettings_.EspCut_F = 60.0;
    hfbsettings_.accuracy_F = 1.0e-5;

    // ^{48}Ca → HFB_{core}.
    const EDFParamsSkyrme edf_skyrme_ = HFBfunctionals::SKMstar();
    AxialHFB hfb_(axialconfig_, hfbsettings_, edf_skyrme_);
    hfb_.initialize_WS_field(Ncore_I, Zcore_I);
    hfb_.iterate(Ncore_I, Zcore_I);

    // HFB_{core} → B_n ⊕ B_p.
    const std::vector<AxialHFBBlocking> blockingCandidates_S1D_n = AxialHFBBlocking::list_candidates(hfbsettings_, hfb_.blocklist_n, true);
    const std::vector<AxialHFBBlocking> blockingCandidates_S1D_p = AxialHFBBlocking::list_candidates(hfbsettings_, hfb_.blocklist_p, false);
    assert(blockingCandidates_S1D_n.size() > 4);
    assert(blockingCandidates_S1D_p.size() > 3);
    std::vector<AxialHFBBlocking> activeBlockings_S1D_q;
    activeBlockings_S1D_q.reserve(2);
    activeBlockings_S1D_q.push_back(blockingCandidates_S1D_n[4]);
    activeBlockings_S1D_q.push_back(blockingCandidates_S1D_p[3]);

    // (HFB_{core},B_n,B_p) → O_{C++}^{48K}.
    hfb_.iterate(Nfinal_I, Zfinal_I, activeBlockings_S1D_q);
    AxialHFBObservable observable_;
    observable_.update_observable(hfb_, activeBlockings_S1D_q);

    // (O_{C++}^{48K},O_{000010}^{48K}) → max errors.
    const AxialHFBObservableReference reference_{-407.810403, -9.067451, -16.342254, -0.000000, -0.000000, 0.000001, -0.000006, 3.649260, 3.438629};
    const double energyErrorMax_F = std::max({std::abs(observable_.Etot_F - reference_.Etot_F), std::abs(observable_.Epair_n_F - reference_.Epair_n_F), std::abs(observable_.Epair_p_F - reference_.Epair_p_F)});
    const double lambdaDeltaError_F = std::max({std::abs(observable_.lambda_n_F - reference_.lambda_n_F), std::abs(observable_.lambda_p_F - reference_.lambda_p_F), std::abs(observable_.Delta_n_F - reference_.Delta_n_F), std::abs(observable_.Delta_p_F - reference_.Delta_p_F)});
    const double rmsError_F = std::max(std::abs(observable_.rms_n_F - reference_.rms_n_F), std::abs(observable_.rms_p_F - reference_.rms_p_F));

    std::cout << std::scientific << std::setprecision(12) << std::right;

    // P_{in} → stdout.
    std::cout << "[Input] Ncore = " << Ncore_I << std::endl;
    std::cout << "[Input] Zcore = " << Zcore_I << std::endl;
    std::cout << "[Input] Nfinal = " << Nfinal_I << std::endl;
    std::cout << "[Input] Zfinal = " << Zfinal_I << std::endl;

    // O_{000010}^{48K} → stdout.
    std::cout << "[Reference] Etot = " << reference_.Etot_F << " MeV" << std::endl;
    std::cout << "[Reference] lambda_n = " << reference_.lambda_n_F << " MeV" << std::endl;
    std::cout << "[Reference] lambda_p = " << reference_.lambda_p_F << " MeV" << std::endl;
    std::cout << "[Reference] Epair_n = " << reference_.Epair_n_F << " MeV" << std::endl;
    std::cout << "[Reference] Epair_p = " << reference_.Epair_p_F << " MeV" << std::endl;
    std::cout << "[Reference] Delta_n = " << reference_.Delta_n_F << " MeV" << std::endl;
    std::cout << "[Reference] Delta_p = " << reference_.Delta_p_F << " MeV" << std::endl;
    std::cout << "[Reference] rms_n = " << reference_.rms_n_F << " fm" << std::endl;
    std::cout << "[Reference] rms_p = " << reference_.rms_p_F << " fm" << std::endl;

    // O_{C++}^{48K} → stdout.
    std::cout << "[Computed] Etot = " << observable_.Etot_F << " MeV" << std::endl;
    std::cout << "[Computed] lambda_n = " << observable_.lambda_n_F << " MeV" << std::endl;
    std::cout << "[Computed] lambda_p = " << observable_.lambda_p_F << " MeV" << std::endl;
    std::cout << "[Computed] Epair_n = " << observable_.Epair_n_F << " MeV" << std::endl;
    std::cout << "[Computed] Epair_p = " << observable_.Epair_p_F << " MeV" << std::endl;
    std::cout << "[Computed] Delta_n = " << observable_.Delta_n_F << " MeV" << std::endl;
    std::cout << "[Computed] Delta_p = " << observable_.Delta_p_F << " MeV" << std::endl;
    std::cout << "[Computed] rms_n = " << observable_.rms_n_F << " fm" << std::endl;
    std::cout << "[Computed] rms_p = " << observable_.rms_p_F << " fm" << std::endl;

    // |O_{C++}^{48K}-O_{000010}^{48K}| → stdout.
    std::cout << "[Error] |Etot - Etot_ref| = " << std::abs(observable_.Etot_F - reference_.Etot_F) << " MeV" << std::endl;
    std::cout << "[Error] |Epair_n - Epair_n_ref| = " << std::abs(observable_.Epair_n_F - reference_.Epair_n_F) << " MeV" << std::endl;
    std::cout << "[Error] |Epair_p - Epair_p_ref| = " << std::abs(observable_.Epair_p_F - reference_.Epair_p_F) << " MeV" << std::endl;
    std::cout << "[Error] |lambda_n - lambda_n_ref| = " << std::abs(observable_.lambda_n_F - reference_.lambda_n_F) << " MeV" << std::endl;
    std::cout << "[Error] |lambda_p - lambda_p_ref| = " << std::abs(observable_.lambda_p_F - reference_.lambda_p_F) << " MeV" << std::endl;
    std::cout << "[Error] |Delta_n - Delta_n_ref| = " << std::abs(observable_.Delta_n_F - reference_.Delta_n_F) << " MeV" << std::endl;
    std::cout << "[Error] |Delta_p - Delta_p_ref| = " << std::abs(observable_.Delta_p_F - reference_.Delta_p_F) << " MeV" << std::endl;
    std::cout << "[Error] |rms_n - rms_n_ref| = " << std::abs(observable_.rms_n_F - reference_.rms_n_F) << " fm" << std::endl;
    std::cout << "[Error] |rms_p - rms_p_ref| = " << std::abs(observable_.rms_p_F - reference_.rms_p_F) << " fm" << std::endl;

    // (ε_E,ε_{λΔ},ε_r) ≤ (10^{-3},10^{-3},10^{-4}).
    assert(energyErrorMax_F <= 1.0e-3);
    assert(lambdaDeltaError_F <= 1.0e-3);
    assert(rmsError_F <= 1.0e-4);
    return 0;
}
