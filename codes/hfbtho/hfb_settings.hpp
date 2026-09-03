/**
 * @file    hfb_settings.hpp
 * @author  cailea
 * @date    2026-09-03
 * @brief   Define HFB iteration controls and term switches.
 */

#pragma once

#include <cassert>
#include <cmath>

#include "hfb_edf_skyrme.hpp"

/**
 * @brief Store enabled HFB energy terms.
 */
class HFBTermSwitches {
public:
    bool addKinetic_B = false; // E_kin = (ℏ²/2m)τ.
    bool addLocalRhoRho_B = false; // E_ρρ = Σ_t C_t^ρρ_t².
    bool addLocalRhoAlpha_B = false; // E_ρα = Σ_t C_t^{ρ,α}ρ_0^αρ_t².
    bool addLocalRhoTau_B = false; // E_ρτ = Σ_t C_t^τρ_tτ_t.
    bool addLocalSurface_B = false; // E_surf = Σ_t(C_t^{Δρ}ρ_tΔρ_t+C_t^{∇ρ}|∇ρ_t|²).
    bool addLocalSpinOrbit_B = false; // E_SO = Σ_t(C_t^{∇J}ρ_t∇·J_t+C_t^{J∇ρ}J_t·∇ρ_t).
    bool addLocalTensor_B = false; // E_tensor = Σ_t(C_t^J J_t²+C_t^{J̄}J_{t,ab}J_{t,ba}).
    bool addLocalCoulomb_B = false; // E_C = E_C^dir+E_C^Slater.
    bool addLocalPair_B = false; // E_pair = Σ_q∫κ_q*Δ_q d³r.
    bool useLocalPairRegularization_B = false; // C_pair → C_pair^reg.
    bool addFiniteRangeGogny_B = false; // V_G = Σ_g(W_g+B_gP_σ-H_gP_τ-M_gP_σP_τ)e^{-r₁₂²/μ_g²}.
    bool addFiniteRangeCoulomb_B = false; // V_C = e²/r₁₂.

public:
    /**
     * @brief Construct disabled HFB term switches.
     * @math ∅ → 0
     * @output Disabled HFB terms.
     */
    HFBTermSwitches() = default;

    /**
     * @brief Build local Skyrme term switches.
     * @math ∅ → S_{Skyrme}
     * @output Standard local-EDF switches.
     */
    static HFBTermSwitches skyrme() {
        HFBTermSwitches termSwitches_;
        termSwitches_.addKinetic_B = true;
        termSwitches_.addLocalRhoRho_B = true;
        termSwitches_.addLocalRhoAlpha_B = true;
        termSwitches_.addLocalRhoTau_B = true;
        termSwitches_.addLocalSurface_B = true;
        termSwitches_.addLocalSpinOrbit_B = true;
        termSwitches_.addLocalCoulomb_B = true;
        termSwitches_.addLocalPair_B = true;
        return termSwitches_;
    }

    /**
     * @brief Build finite-range Gogny term switches.
     * @math ∅ → S_{Gogny}
     * @output Gogny interaction switches.
     */
    static HFBTermSwitches gogny() {
        HFBTermSwitches termSwitches_;
        termSwitches_.addKinetic_B = true;
        termSwitches_.addLocalRhoAlpha_B = true;
        termSwitches_.addLocalSpinOrbit_B = true;
        termSwitches_.addFiniteRangeGogny_B = true;
        termSwitches_.addFiniteRangeCoulomb_B = true;
        return termSwitches_;
    }
};

/**
 * @brief Store HFB iteration controls and term switches.
 */
class HFBSettings {
public:
    int Nblocking_I = 6; // Number of blocking candidates.
    double EblockingCut_F = 1.0; // Blocking window [MeV].
    double temperature_F = 0.0; // Temperature [MeV].
    double accuracy_F = 1.0e-5; // SCF convergence tolerance.
    double EspCut_F = 60.0; // Active quasiparticle cutoff [MeV].
    bool useEspCut_B = false; // E_qp ≤ E_sp^cut.
    bool useLipkinNogami_B = false; // Enable Lipkin-Nogami correction.
    bool useCmCorrection_B = false; // ℏ²/(2m) → [1-1/A]ℏ²/(2m).
    HFBTermSwitches termSwitches{}; // Enabled energy terms.

public:
    /**
     * @brief Construct default HFB controls.
     * @math ∅ → P_{HFB}
     * @output Initialized solver settings.
     */
    HFBSettings() = default;

    /**
     * @brief Build standard Skyrme HFB settings.
     * @math ∅ → P_{HFB}^{Skyrme}
     * @output Skyrme solver settings.
     */
    static HFBSettings setting_skyrme() {
        HFBSettings hfbsettings_;
        hfbsettings_.useEspCut_B = true;
        hfbsettings_.useCmCorrection_B = true;
        hfbsettings_.termSwitches = HFBTermSwitches::skyrme();
        return hfbsettings_;
    }

    /**
     * @brief Build standard Gogny HFB settings.
     * @math ∅ → P_{HFB}^{Gogny}
     * @output Gogny solver settings.
     */
    static HFBSettings setting_gogny() {
        HFBSettings hfbsettings_;
        hfbsettings_.useCmCorrection_B = true;
        hfbsettings_.termSwitches = HFBTermSwitches::gogny();
        return hfbsettings_;
    }

    /**
     * @brief Calculate the spherical oscillator length.
     * @math A → b_0
     * @output Oscillator length b_0 [fm].
     */
    static double calc_b0(int Atarget_I) {
        assert(Atarget_I > 0);
        const double hbzero_F = EDFParamsSkyrme{}.hbzero_F;
        const double r0_F = 1.20;
        const double hbar_omega_F = 41.0 * std::pow(static_cast<double>(Atarget_I), -1.0 / 3.0) * r0_F;
        return std::sqrt(2.0 * hbzero_F / hbar_omega_F);
    }

    /**
     * @brief Apply nucleus-dependent EDF corrections.
     * @math (C,A,P_{HFB}) → C_{active}
     * @output Active Skyrme EDF parameters.
     */
    EDFParamsSkyrme make_active_edf(const EDFParamsSkyrme& base_edf_, int Atarget_I) const {
        assert(Atarget_I > 0);
        EDFParamsSkyrme active_edf_ = base_edf_;
        const double cm_factor_F = 1.0 - static_cast<double>(useCmCorrection_B) / static_cast<double>(Atarget_I);

        // C_{kin} → s_{kin}(1-s_{cm}/A)C_{kin}.
        active_edf_.hbzero_F *= cm_factor_F * static_cast<double>(termSwitches.addKinetic_B);
        active_edf_.hbzeron_F *= cm_factor_F * static_cast<double>(termSwitches.addKinetic_B);
        active_edf_.hbzerop_F *= cm_factor_F * static_cast<double>(termSwitches.addKinetic_B);

        // C_i → s_i C_i.
        active_edf_.Crho_0_F *= static_cast<double>(termSwitches.addLocalRhoRho_B);
        active_edf_.Crho_1_F *= static_cast<double>(termSwitches.addLocalRhoRho_B);
        active_edf_.Cdrho_0_F *= static_cast<double>(termSwitches.addLocalRhoAlpha_B);
        active_edf_.Cdrho_1_F *= static_cast<double>(termSwitches.addLocalRhoAlpha_B);
        active_edf_.Ctau_0_F *= static_cast<double>(termSwitches.addLocalRhoTau_B);
        active_edf_.Ctau_1_F *= static_cast<double>(termSwitches.addLocalRhoTau_B);
        active_edf_.CrDr_0_F *= static_cast<double>(termSwitches.addLocalSurface_B);
        active_edf_.CrDr_1_F *= static_cast<double>(termSwitches.addLocalSurface_B);
        active_edf_.Cnrho_0_F *= static_cast<double>(termSwitches.addLocalSurface_B);
        active_edf_.Cnrho_1_F *= static_cast<double>(termSwitches.addLocalSurface_B);
        active_edf_.CrdJ_0_F *= static_cast<double>(termSwitches.addLocalSpinOrbit_B);
        active_edf_.CrdJ_1_F *= static_cast<double>(termSwitches.addLocalSpinOrbit_B);
        active_edf_.CJdr_0_F *= static_cast<double>(termSwitches.addLocalSpinOrbit_B);
        active_edf_.CJdr_1_F *= static_cast<double>(termSwitches.addLocalSpinOrbit_B);
        active_edf_.CJ_0_F *= static_cast<double>(termSwitches.addLocalTensor_B);
        active_edf_.CJ_1_F *= static_cast<double>(termSwitches.addLocalTensor_B);
        active_edf_.CJbar_0_F *= static_cast<double>(termSwitches.addLocalTensor_B);
        active_edf_.CJbar_1_F *= static_cast<double>(termSwitches.addLocalTensor_B);
        active_edf_.CpV0_0_F *= static_cast<double>(termSwitches.addLocalPair_B);
        active_edf_.CpV0_1_F *= static_cast<double>(termSwitches.addLocalPair_B);
        return active_edf_;
    }
};
