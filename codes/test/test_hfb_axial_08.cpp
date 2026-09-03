/**
 * @file    test_hfb_axial_08.cpp
 * @author  cailea
 * @date    2026-09-02
 * @brief   Validate finite-temperature SLY4 against run 000008.
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

/**
 * @brief  Validate finite-temperature SLY4 HFB observables.
 * @math   O_{C++} ↔ O_{000008}
 * @output Zero on accepted observable errors.
 */
int main() {
    // (N,Z,N_{shell},b_0,N_r,N_z) ← run 000008.
    const int Ntarget_I = 26;
    const int Ztarget_I = 24;
    const int Nshell_I = 12;
    const double b0_F = 1.76221461551038;
    AxialConfig axialconfig_(b0_F, b0_F, Nshell_I, true);
    axialconfig_.Nr_I = 40;
    axialconfig_.Nz_I = 40;

    // P_{HFB} ← run 000008.
    HFBSettings hfbsettings_ = HFBSettings::setting_skyrme();
    hfbsettings_.useEspCut_B = true;
    hfbsettings_.useCmCorrection_B = true;
    hfbsettings_.temperature_F = 1.5;
    hfbsettings_.EspCut_F = 60.0;
    hfbsettings_.accuracy_F = 1.0e-9;

    // (C_{axial},P_{HFB},SLY4) → O_{C++}.
    const EDFParamsSkyrme edf_skyrme_ = HFBfunctionals::SLY4();
    AxialHFB hfb_(axialconfig_, hfbsettings_, edf_skyrme_);
    hfb_.initialize_WS_field(Ntarget_I, Ztarget_I);
    hfb_.iterate(Ntarget_I, Ztarget_I);
    AxialHFBObservable observable_;
    observable_.update_observable(hfb_);

    // (O_{C++},O_{000008}) → max errors.
    const AxialHFBObservableReference reference_{-423.459451, -11.806442, -6.886468, -0.000000, -0.000000, 0.000000, -0.000000, 3.591964, 3.594830};
    const double energyErrorMax_F = std::max({std::abs(observable_.Etot_F - reference_.Etot_F), std::abs(observable_.Epair_n_F - reference_.Epair_n_F), std::abs(observable_.Epair_p_F - reference_.Epair_p_F)});
    const double lambdaDeltaError_F = std::max({std::abs(observable_.lambda_n_F - reference_.lambda_n_F), std::abs(observable_.lambda_p_F - reference_.lambda_p_F), std::abs(observable_.Delta_n_F - reference_.Delta_n_F), std::abs(observable_.Delta_p_F - reference_.Delta_p_F)});
    const double rmsError_F = std::max(std::abs(observable_.rms_n_F - reference_.rms_n_F), std::abs(observable_.rms_p_F - reference_.rms_p_F));

    std::cout << std::scientific << std::setprecision(12) << std::right;

    // P_{in} → stdout.
    std::cout << "[Input] N = " << Ntarget_I << std::endl;
    std::cout << "[Input] Z = " << Ztarget_I << std::endl;
    std::cout << "[Input] Nshell = " << Nshell_I << std::endl;
    std::cout << "[Input] T = " << hfbsettings_.temperature_F << " MeV" << std::endl;

    // O_{000008} → stdout.
    std::cout << "[Reference] Etot = " << reference_.Etot_F << " MeV" << std::endl;
    std::cout << "[Reference] lambda_n = " << reference_.lambda_n_F << " MeV" << std::endl;
    std::cout << "[Reference] lambda_p = " << reference_.lambda_p_F << " MeV" << std::endl;
    std::cout << "[Reference] Epair_n = " << reference_.Epair_n_F << " MeV" << std::endl;
    std::cout << "[Reference] Epair_p = " << reference_.Epair_p_F << " MeV" << std::endl;
    std::cout << "[Reference] Delta_n = " << reference_.Delta_n_F << " MeV" << std::endl;
    std::cout << "[Reference] Delta_p = " << reference_.Delta_p_F << " MeV" << std::endl;
    std::cout << "[Reference] rms_n = " << reference_.rms_n_F << " fm" << std::endl;
    std::cout << "[Reference] rms_p = " << reference_.rms_p_F << " fm" << std::endl;

    // O_{C++} → stdout.
    std::cout << "[Computed] Etot = " << observable_.Etot_F << " MeV" << std::endl;
    std::cout << "[Computed] lambda_n = " << observable_.lambda_n_F << " MeV" << std::endl;
    std::cout << "[Computed] lambda_p = " << observable_.lambda_p_F << " MeV" << std::endl;
    std::cout << "[Computed] Epair_n = " << observable_.Epair_n_F << " MeV" << std::endl;
    std::cout << "[Computed] Epair_p = " << observable_.Epair_p_F << " MeV" << std::endl;
    std::cout << "[Computed] Delta_n = " << observable_.Delta_n_F << " MeV" << std::endl;
    std::cout << "[Computed] Delta_p = " << observable_.Delta_p_F << " MeV" << std::endl;
    std::cout << "[Computed] rms_n = " << observable_.rms_n_F << " fm" << std::endl;
    std::cout << "[Computed] rms_p = " << observable_.rms_p_F << " fm" << std::endl;

    // |O_{C++}-O_{000008}| → stdout.
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
