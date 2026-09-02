/**
 * @file    test_hfb_axial_14.cpp
 * @author  cailea
 * @date    2026-09-02
 * @brief   Validate deformed Sn120 D1S observables.
 */

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>

#include "hfb_axial.hpp"
#include "hfb_axial_observable.hpp"
#include "hfb_edf_gogny.hpp"
#include "test_hfb_axial_reference.hpp"

/**
 * @brief  Validate observables using absolute tolerances.
 * @math   O_{C++} ≈ O_{HFBTHO}
 * @output 0 when all observable errors satisfy tolerances.
 */
int main() {
    std::cout << std::scientific << std::setprecision(12) << std::right;

    // (N,Z,N_shell,b_0,β_20,N_z,N_r) → C_axial.
    const int Ntarget_I = 70;
    const int Ztarget_I = 50;
    const int Nshell_I = 20;
    const double b0_F = 2.13396610;
    const double beta20_F = 0.40;
    const auto [bz_F, br_F] = AxialConfig::b0beta20_to_bzbr(b0_F, beta20_F);
    AxialConfig axialconfig_(bz_F, br_F, Nshell_I, true);
    axialconfig_.Nr_I = 50;
    axialconfig_.Nz_I = 50;

    // S_HFB → S_D1S.
    HFBSettings hfbsettings_;
    hfbsettings_.termSwitches = HFBSettings::TermSwitches::gogny();
    hfbsettings_.termSwitches.addLocalCoulomb_B = false;
    hfbsettings_.termSwitches.addFiniteRangeCoulomb_B = false;
    hfbsettings_.useEspCut_B = false;
    hfbsettings_.useCmCorrection_B = true;
    hfbsettings_.useLipkinNogami_B = false;
    hfbsettings_.accuracy_F = 1.0e-8;
    const EDFParamsGogny edf_gogny_ = EDFParamsGogny::D1S();

    // Inputs → stdout.
    std::cout << "[Input] (N,Z,Nshell,b0,beta20,Nz,Nr) = (" << Ntarget_I << "," << Ztarget_I << "," << Nshell_I << "," << b0_F << "," << beta20_F << "," << axialconfig_.Nz_I << "," << axialconfig_.Nr_I << ")\n";
    std::cout << "[Input] force = " << edf_gogny_.forceName_Str << ", Coulomb = false\n";

    // (C_axial,S_HFB,D1S) → HFB → O.
    AxialHFB hfb_(axialconfig_, hfbsettings_, edf_gogny_.make_local_edf(), edf_gogny_);
    hfb_.initialize_WS_field(Ntarget_I, Ztarget_I, beta20_F);
    hfb_.iterate(Ntarget_I, Ztarget_I);
    AxialHFBObservable observable_;
    observable_.update_observable(hfb_);

    // O_ref ← HFBTHO/outputs/run_000014/hfbtho.out.
    const AxialHFBObservableReference reference_{-1385.166912, -7.365802, -24.759824, -20.714699, -0.000000, 2.246873, -0.000000, 4.626410, 4.438686};

    // O_ref → stdout.
    std::cout << "[Reference] Etot [MeV] = " << reference_.Etot_F << "\n";
    std::cout << "[Reference] lambda_n [MeV] = " << reference_.lambda_n_F << "\n";
    std::cout << "[Reference] lambda_p [MeV] = " << reference_.lambda_p_F << "\n";
    std::cout << "[Reference] Epair_n [MeV] = " << reference_.Epair_n_F << "\n";
    std::cout << "[Reference] Epair_p [MeV] = " << reference_.Epair_p_F << "\n";
    std::cout << "[Reference] Delta_n [MeV] = " << reference_.Delta_n_F << "\n";
    std::cout << "[Reference] Delta_p [MeV] = " << reference_.Delta_p_F << "\n";
    std::cout << "[Reference] rms_n [fm] = " << reference_.rms_n_F << "\n";
    std::cout << "[Reference] rms_p [fm] = " << reference_.rms_p_F << "\n";

    // O_C++ → stdout.
    std::cout << "[Computed] Etot [MeV] = " << observable_.Etot_F << "\n";
    std::cout << "[Computed] lambda_n [MeV] = " << observable_.lambda_n_F << "\n";
    std::cout << "[Computed] lambda_p [MeV] = " << observable_.lambda_p_F << "\n";
    std::cout << "[Computed] Epair_n [MeV] = " << observable_.Epair_n_F << "\n";
    std::cout << "[Computed] Epair_p [MeV] = " << observable_.Epair_p_F << "\n";
    std::cout << "[Computed] Delta_n [MeV] = " << observable_.Delta_n_F << "\n";
    std::cout << "[Computed] Delta_p [MeV] = " << observable_.Delta_p_F << "\n";
    std::cout << "[Computed] rms_n [fm] = " << observable_.rms_n_F << "\n";
    std::cout << "[Computed] rms_p [fm] = " << observable_.rms_p_F << "\n";

    // |O_C++-O_ref| → stdout.
    std::cout << "[Error] |Etot-Etot_ref| [MeV] = " << std::abs(observable_.Etot_F - reference_.Etot_F) << "\n";
    std::cout << "[Error] |lambda_n-lambda_n_ref| [MeV] = " << std::abs(observable_.lambda_n_F - reference_.lambda_n_F) << "\n";
    std::cout << "[Error] |lambda_p-lambda_p_ref| [MeV] = " << std::abs(observable_.lambda_p_F - reference_.lambda_p_F) << "\n";
    std::cout << "[Error] |Epair_n-Epair_n_ref| [MeV] = " << std::abs(observable_.Epair_n_F - reference_.Epair_n_F) << "\n";
    std::cout << "[Error] |Epair_p-Epair_p_ref| [MeV] = " << std::abs(observable_.Epair_p_F - reference_.Epair_p_F) << "\n";
    std::cout << "[Error] |Delta_n-Delta_n_ref| [MeV] = " << std::abs(observable_.Delta_n_F - reference_.Delta_n_F) << "\n";
    std::cout << "[Error] |Delta_p-Delta_p_ref| [MeV] = " << std::abs(observable_.Delta_p_F - reference_.Delta_p_F) << "\n";
    std::cout << "[Error] |rms_n-rms_n_ref| [fm] = " << std::abs(observable_.rms_n_F - reference_.rms_n_F) << "\n";
    std::cout << "[Error] |rms_p-rms_p_ref| [fm] = " << std::abs(observable_.rms_p_F - reference_.rms_p_F) << "\n";

    // max|δE|≤10^-3 MeV; max|δr|≤10^-4 fm.
    assert(std::abs(observable_.Etot_F - reference_.Etot_F) <= 1.0e-3);
    assert(std::max(std::abs(observable_.Epair_n_F - reference_.Epair_n_F), std::abs(observable_.Epair_p_F - reference_.Epair_p_F)) <= 1.0e-3);
    assert(std::max(std::abs(observable_.lambda_n_F - reference_.lambda_n_F), std::abs(observable_.lambda_p_F - reference_.lambda_p_F)) <= 1.0e-3);
    assert(std::max(std::abs(observable_.Delta_n_F - reference_.Delta_n_F), std::abs(observable_.Delta_p_F - reference_.Delta_p_F)) <= 1.0e-3);
    assert(std::max(std::abs(observable_.rms_n_F - reference_.rms_n_F), std::abs(observable_.rms_p_F - reference_.rms_p_F)) <= 1.0e-4);
    return 0;
}
