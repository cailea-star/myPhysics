/**
 * @file    hfb_edf_skyrme.hpp
 * @author  cailea
 * @date    2026-05-06
 * @brief   Define Skyrme EDF couplings and conversions.
 */

#pragma once

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

namespace {

constexpr double kfconst_F = 2.4554457015685776; // k_F/ρ^(1/3) = kfconst_F = (3π²/2)^(1/3).
constexpr double CK_F = 3.6175281560109624; // τ/ρ^(5/3) = CK_F = 3kfconst_F²/5.

} // namespace

/**
 * @brief Store Skyrme EDF and nuclear-matter parameters.
 */
struct EDFParamsSkyrme {
    // t ∈ {0,1}; q ∈ {n,p}; ρ_0 = ρ_p + ρ_n; ρ_1 = ρ_n - ρ_p.
    // 𝓔_even(r) = Σ_{t=0,1}[(C_t^ρ + C_t^{ρ,α}ρ_0^α)ρ_t² + C_t^τρ_tτ_t
    //             + C_t^{Δρ}ρ_tΔρ_t + C_t^{∇J}ρ_t∇·J_t
    //             + C_t^JΣ_{ab}J_{t,ab}² + C_t^{J̄}Σ_{ab}J_{t,ab}J_{t,ba}
    //             + C_t^{∇ρ}(∇ρ_t)² + C_t^{J∇ρ}J_t·∇ρ_t].
    double Crho_0_F = 0.0; // C_0^ρ.
    double Crho_1_F = 0.0; // C_1^ρ.
    double Cdrho_0_F = 0.0; // C_0^{ρ,α}.
    double Cdrho_1_F = 0.0; // C_1^{ρ,α}.
    double Ctau_0_F = 0.0; // C_0^τ.
    double Ctau_1_F = 0.0; // C_1^τ.
    double CrDr_0_F = 0.0; // C_0^{Δρ}.
    double CrDr_1_F = 0.0; // C_1^{Δρ}.
    double CrdJ_0_F = 0.0; // C_0^{∇J}.
    double CrdJ_1_F = 0.0; // C_1^{∇J}.
    double CJ_0_F = 0.0; // C_0^J.
    double CJ_1_F = 0.0; // C_1^J.
    double CJbar_0_F = 0.0; // C_0^{J̄}.
    double CJbar_1_F = 0.0; // C_1^{J̄}.
    double Cnrho_0_F = 0.0; // C_0^{∇ρ}.
    double Cnrho_1_F = 0.0; // C_1^{∇ρ}.
    double CJdr_0_F = 0.0; // C_0^{J∇ρ}.
    double CJdr_1_F = 0.0; // C_1^{J∇ρ}.

    // 𝓔_odd(r) = Σ_{t=0,1}[(C_t^s + C_t^{s,α}ρ_0^α)s_t² + C_t^j j_t²
    //            + C_t^{Δs}s_t·Δs_t + C_t^{∇j}s_t·(∇×j_t)
    //            + C_t^T s_t·T_t + C_t^{∇s}(∇·s_t)² + C_t^F s_t·F_t].
    double Cs_0_F = 0.0; // C_0^s.
    double Cs_1_F = 0.0; // C_1^s.
    double Cds_0_F = 0.0; // C_0^{s,α}.
    double Cds_1_F = 0.0; // C_1^{s,α}.
    double Csj_0_F = 0.0; // C_0^j.
    double Csj_1_F = 0.0; // C_1^j.
    double CsDs_0_F = 0.0; // C_0^{Δs}.
    double CsDs_1_F = 0.0; // C_1^{Δs}.
    double CsdJ_0_F = 0.0; // C_0^{∇j}.
    double CsdJ_1_F = 0.0; // C_1^{∇j}.
    double CT_0_F = 0.0; // C_0^T.
    double CT_1_F = 0.0; // C_1^T.
    double CsNabs_0_F = 0.0; // C_0^{∇s}.
    double CsNabs_1_F = 0.0; // C_1^{∇s}.
    double CF_0_F = 0.0; // C_0^F.
    double CF_1_F = 0.0; // C_1^F.

    // V_pair,q(r) ≡ Δ_q(r) = C_q^{V0}[1 - C_q^{V1}ρ_0(r)/ρ_c]κ_q(r).
    // E_pair^{HFBTHO} = Σ_{q=n,p}∫d³r κ_q^*(r)Δ_q(r).
    double CpV0_0_F = 0.0; // Neutron strength [MeV fm³].
    double CpV0_1_F = 0.0; // Proton strength [MeV fm³].
    double CpV1_0_F = 0.5; // Neutron density factor.
    double CpV1_1_F = 0.5; // Proton density factor.

    double sigma_F = 1.0 / 6.0; // Density exponent α.

    // Physical constants.
    double hbzero_F = 20.735530; // ℏ²/(2m) [MeV fm²].
    double hbzeron_F = 20.721247; // ℏ²/(2m_n) [MeV fm²].
    double hbzerop_F = 20.749809; // ℏ²/(2m_p) [MeV fm²].
    double e2charg_F = 1.439978408596513; // e² [MeV fm].
    double CExPar_F = 1.0; // Coulomb-exchange factor.

    // Nuclear-matter observables.
    double rhoNM_F = 0.16; // ρ_sat [fm⁻³].
    double ENM_F = -16.0; // E/A [MeV].
    double KNM_F = 230.0; // K [MeV].
    double SMASSNM_F = 0.9; // m_s*/m.
    double VMASSNM_F = 1.0; // m_v*/m.
    double ASSNM_F = 30.0; // a_sym [MeV].
    double LASSNM_F = 50.0; // L_sym [MeV].
    double PNM_F = 0.0; // P_sat [MeV fm⁻³].

    std::string functionalName_Str{"Unknown"}; // EDF identifier.
};

namespace SkyrmeINM {

/**
 * @brief Store nuclear-matter-derived volume couplings.
 */
struct VolumeCouplings {
    double Crho_0_F = 0.0;
    double Crho_1_F = 0.0;
    double Cdrho_0_F = 0.0;
    double Cdrho_1_F = 0.0;
    double Ctau_0_F = 0.0;
    double Ctau_1_F = 0.0;
    double sigma_F = 0.0;
};

/**
 * @brief  Calculate standard Skyrme volume couplings.
 * @math   (ρ,E,K,m_s^*,a_sym,L_sym,m_v^*) → (C_t^ρ,C_t^{ρ,α},C_t^τ,α)
 * @output Nuclear-matter-derived volume couplings.
 * @note   Excludes DME, three-body, and finite-range terms.
 */
inline VolumeCouplings calc_volume_couplings_from_nuclear_matter(double rho_F, double E_F, double K_F, double SMASS_F, double ASS_F, double LASS_F, double VMASS_F, double hbzero_F) {
    VolumeCouplings couplings_;

    // (ρ,C_K) → τ_c.
    [[maybe_unused]] constexpr double kf_F = kfconst_F;
    constexpr double CK_F = ::CK_F;
    constexpr double c13_F = 1.0 / 3.0;
    constexpr double c23_F = 2.0 / 3.0;
    const double tauc_F = CK_F * std::pow(rho_F, c23_F);

    // α = N_α/D_α.
    const double sigma_numer_F = c13_F * (-K_F + tauc_F * hbzero_F * (-3.0 + 4.0 * SMASS_F) - 9.0 * E_F);
    const double sigma_denom_F = tauc_F * hbzero_F * (-3.0 + 2.0 * SMASS_F) + 3.0 * E_F;
    couplings_.sigma_F = sigma_numer_F / sigma_denom_F;
    const double sigma_F = couplings_.sigma_F;

    // (ρ,E,m_s^*,α) → C_0^ρ,C_0^{ρ,α},C_0^τ.
    couplings_.Crho_0_F = c13_F * (tauc_F * hbzero_F * (-3.0 + (2.0 - 3.0 * sigma_F) * SMASS_F) + 3.0 * (1.0 + sigma_F) * E_F) / (sigma_F * rho_F);
    couplings_.Cdrho_0_F = c13_F * std::pow(rho_F, -1.0 - sigma_F) * (tauc_F * hbzero_F * (3.0 - 2.0 * SMASS_F) - 3.0 * E_F) / sigma_F;
    couplings_.Ctau_0_F = hbzero_F * (SMASS_F - 1.0) / rho_F;

    // (ρ,a_sym,L_sym,m_v^*,α,C_0^τ) → C_1^ρ,C_1^{ρ,α},C_1^τ.
    const double Ctau_0_F = couplings_.Ctau_0_F;
    couplings_.Crho_1_F = (27.0 * ASS_F * (1.0 + sigma_F) - 9.0 * LASS_F + 5.0 * tauc_F * hbzero_F * (5.0 - 6.0 * VMASS_F + 3.0 * sigma_F * (-4.0 + 3.0 * VMASS_F)) + rho_F * 40.0 * tauc_F * Ctau_0_F * (1.0 - 1.5 * sigma_F)) / (27.0 * sigma_F * rho_F);
    couplings_.Cdrho_1_F = -std::pow(rho_F, -1.0 - sigma_F) * (27.0 * ASS_F - 9.0 * LASS_F + 5.0 * tauc_F * hbzero_F * (5.0 - 6.0 * VMASS_F) + 40.0 * tauc_F * rho_F * Ctau_0_F) / (27.0 * sigma_F);
    couplings_.Ctau_1_F = (hbzero_F * (1.0 - VMASS_F) + rho_F * Ctau_0_F) / rho_F;

    return couplings_;
}

} // namespace SkyrmeINM


namespace HFBfunctionals {

/**
 * @brief  Build the UNEDF1 Skyrme functional.
 * @math   NM → (C_t^ρ,C_t^{ρ,α},C_t^τ); P_sat = 0
 * @output UNEDF1 EDF parameters.
 * @note   Matches HFBTHO UNE1 parameters.
 */
inline EDFParamsSkyrme UNEDF1() {
    EDFParamsSkyrme edf_parameters_;
    edf_parameters_.functionalName_Str = "UNEDF1";

    // (ρ,E,K,m_s*,a_sym,L_sym,m_v*) ← UNE1.
    constexpr double rho_F = 0.158706769332587;
    constexpr double E_F = -15.8;
    constexpr double K_F = 220.0;
    constexpr double SMASS_F = 0.992423332283364;
    constexpr double ASS_F = 28.9867890577721;
    constexpr double LASS_F = 40.0047904804136;
    constexpr double VMASS_F = 1.24983857423227;

    // (ℏ²/2m,e²,C_ex) ← UNE1.
    edf_parameters_.hbzero_F = 20.73553;
    edf_parameters_.hbzeron_F = 20.73553;
    edf_parameters_.hbzerop_F = 20.73553;
    edf_parameters_.e2charg_F = 1.439978408596513;
    edf_parameters_.CExPar_F = 1.0;

    // NM → volume couplings.
    const SkyrmeINM::VolumeCouplings volume_couplings_ = SkyrmeINM::calc_volume_couplings_from_nuclear_matter(rho_F, E_F, K_F, SMASS_F, ASS_F, LASS_F, VMASS_F, edf_parameters_.hbzero_F);
    edf_parameters_.sigma_F = volume_couplings_.sigma_F;
    edf_parameters_.Crho_0_F = volume_couplings_.Crho_0_F;
    edf_parameters_.Crho_1_F = volume_couplings_.Crho_1_F;
    edf_parameters_.Cdrho_0_F = volume_couplings_.Cdrho_0_F;
    edf_parameters_.Cdrho_1_F = volume_couplings_.Cdrho_1_F;
    edf_parameters_.Ctau_0_F = volume_couplings_.Ctau_0_F;
    edf_parameters_.Ctau_1_F = volume_couplings_.Ctau_1_F;

    // (C_t^{Δρ},C_t^{∇J},C_t^J) ← UNE1.
    edf_parameters_.CrDr_0_F = -45.1351310222373;
    edf_parameters_.CrDr_1_F = -145.382167908057;
    edf_parameters_.CrdJ_0_F = -74.0263331764599;
    edf_parameters_.CrdJ_1_F = -35.6582611147917;
    edf_parameters_.CJ_0_F = 0.0;
    edf_parameters_.CJ_1_F = 0.0;
    edf_parameters_.CJbar_0_F = 0.0;
    edf_parameters_.CJbar_1_F = 0.0;
    edf_parameters_.Cnrho_0_F = 0.0;
    edf_parameters_.Cnrho_1_F = 0.0;
    edf_parameters_.CJdr_0_F = 0.0;
    edf_parameters_.CJdr_1_F = 0.0;

    // (C_q^{V0},C_q^{V1}) ← UNE1.
    edf_parameters_.CpV0_0_F = -186.065399575124;
    edf_parameters_.CpV0_1_F = -206.579593890106;
    edf_parameters_.CpV1_0_F = 0.5;
    edf_parameters_.CpV1_1_F = 0.5;

    // C_t^j=-C_t^τ; C_t^{∇j}=C_t^{∇J}; C_t^T=C_t^F=0.
    edf_parameters_.Cs_0_F = 0.0;
    edf_parameters_.Cs_1_F = 0.0;
    edf_parameters_.Cds_0_F = 0.0;
    edf_parameters_.Cds_1_F = 0.0;
    edf_parameters_.Csj_0_F = -edf_parameters_.Ctau_0_F;
    edf_parameters_.Csj_1_F = -edf_parameters_.Ctau_1_F;
    edf_parameters_.CsDs_0_F = 0.0;
    edf_parameters_.CsDs_1_F = 0.0;
    edf_parameters_.CsdJ_0_F = edf_parameters_.CrdJ_0_F;
    edf_parameters_.CsdJ_1_F = edf_parameters_.CrdJ_1_F;
    edf_parameters_.CT_0_F = 0.0;
    edf_parameters_.CT_1_F = 0.0;
    edf_parameters_.CsNabs_0_F = 0.0;
    edf_parameters_.CsNabs_1_F = 0.0;
    edf_parameters_.CF_0_F = 0.0;
    edf_parameters_.CF_1_F = 0.0;

    // (ρ,E,K,m_s*,m_v*,a_sym,L_sym,P) → EDF metadata.
    edf_parameters_.rhoNM_F = rho_F;
    edf_parameters_.ENM_F = E_F;
    edf_parameters_.KNM_F = K_F;
    edf_parameters_.SMASSNM_F = SMASS_F;
    edf_parameters_.VMASSNM_F = VMASS_F;
    edf_parameters_.ASSNM_F = ASS_F;
    edf_parameters_.LASSNM_F = LASS_F;
    edf_parameters_.PNM_F = 0.0;

    return edf_parameters_;
}

/**
 * @brief  Build the SkM* Skyrme functional.
 * @math   SkM^* → \{C_t\}
 * @output SkM* EDF parameters.
 */
inline EDFParamsSkyrme SKMstar() {
    EDFParamsSkyrme edf_parameters_;

    // {C_t,α,ℏ²/2m,e²,C_ex} ← SkM*.
    edf_parameters_.functionalName_Str = "SKM*";
    edf_parameters_.Crho_0_F = -991.875;
    edf_parameters_.Crho_1_F = 390.1375;
    edf_parameters_.Cdrho_0_F = 974.6875;
    edf_parameters_.Cdrho_1_F = -324.895833333333314;
    edf_parameters_.Ctau_0_F = 34.6875;
    edf_parameters_.Ctau_1_F = -34.0625;
    edf_parameters_.CrDr_0_F = -68.203125;
    edf_parameters_.CrDr_1_F = 17.109375;
    edf_parameters_.CrdJ_0_F = -97.5;
    edf_parameters_.CrdJ_1_F = -32.5;
    edf_parameters_.CJ_0_F = 34.0625;
    edf_parameters_.CJ_1_F = 34.0625;
    edf_parameters_.CJbar_0_F = 0.0;
    edf_parameters_.CJbar_1_F = 0.0;
    edf_parameters_.Cnrho_0_F = 0.0;
    edf_parameters_.Cnrho_1_F = 0.0;
    edf_parameters_.CJdr_0_F = 0.0;
    edf_parameters_.CJdr_1_F = 0.0;
    edf_parameters_.CpV0_0_F = -280.0;
    edf_parameters_.CpV0_1_F = -280.0;
    edf_parameters_.CpV1_0_F = 0.5;
    edf_parameters_.CpV1_1_F = 0.5;
    edf_parameters_.sigma_F = 1.0 / 6.0;
    edf_parameters_.hbzero_F = 20.73;
    edf_parameters_.hbzeron_F = 20.73;
    edf_parameters_.hbzerop_F = 20.73;
    edf_parameters_.e2charg_F = 1.439978408596513;
    edf_parameters_.CExPar_F = 1.0;

    // (ρ,E,K,m_s*,m_v*,a_sym,L_sym,P) ← SkM*.
    edf_parameters_.rhoNM_F = 0.160318515964671449;
    edf_parameters_.ENM_F = -15.7762359616147627;
    edf_parameters_.KNM_F = 216.657542303868979;
    edf_parameters_.SMASSNM_F = 1.26826090316085582;
    edf_parameters_.VMASSNM_F = 1.53168827653502970;
    edf_parameters_.ASSNM_F = 30.0323865018228560;
    edf_parameters_.LASSNM_F = 45.7703887670942962;
    edf_parameters_.PNM_F = 0.0;

    return edf_parameters_;
}

/**
 * @brief  Build the SLY4 Skyrme functional.
 * @math   SLY4 → \{C_t\}
 * @output SLY4 EDF parameters.
 */
inline EDFParamsSkyrme SLY4() {
    EDFParamsSkyrme edf_parameters_;

    // {C_t,α,ℏ²/2m,e²,C_ex} ← SLY4.
    edf_parameters_.functionalName_Str = "SLY4";
    edf_parameters_.Crho_0_F = -933.342374999999947;
    edf_parameters_.Crho_1_F = 830.052485500000103;
    edf_parameters_.Cdrho_0_F = 861.062500000000000;
    edf_parameters_.Cdrho_1_F = -1064.27324999999996;
    edf_parameters_.Ctau_0_F = 57.1286874999999981;
    edf_parameters_.Ctau_1_F = 24.6567364999999974;
    edf_parameters_.CrDr_0_F = -76.9962031249999939;
    edf_parameters_.CrDr_1_F = 15.6571351249999999;
    edf_parameters_.CrdJ_0_F = -92.2500000000000000;
    edf_parameters_.CrdJ_1_F = -30.7500000000000000;
    edf_parameters_.CJ_0_F = 17.2096115000000012;
    edf_parameters_.CJ_1_F = 64.5758124999999978;
    edf_parameters_.CJbar_0_F = 0.0;
    edf_parameters_.CJbar_1_F = 0.0;
    edf_parameters_.Cnrho_0_F = 0.0;
    edf_parameters_.Cnrho_1_F = 0.0;
    edf_parameters_.CJdr_0_F = 0.0;
    edf_parameters_.CJdr_1_F = 0.0;
    edf_parameters_.CpV0_0_F = -300.0;
    edf_parameters_.CpV0_1_F = -300.0;
    edf_parameters_.CpV1_0_F = 0.5;
    edf_parameters_.CpV1_1_F = 0.5;
    edf_parameters_.sigma_F = 1.0 / 6.0;
    edf_parameters_.hbzero_F = 20.735530;
    edf_parameters_.hbzeron_F = 20.735530;
    edf_parameters_.hbzerop_F = 20.735530;
    edf_parameters_.e2charg_F = 1.43997840859651349;
    edf_parameters_.CExPar_F = 1.0;

    // (ρ,E,K,m_s*,m_v*,a_sym,L_sym,P) ← SLY4.
    edf_parameters_.rhoNM_F = 0.159538756711733343;
    edf_parameters_.ENM_F = -15.9721491414446017;
    edf_parameters_.KNM_F = 229.900964482603513;
    edf_parameters_.SMASSNM_F = 1.43954698897607836;
    edf_parameters_.VMASSNM_F = 1.24983854719625342;
    edf_parameters_.ASSNM_F = 32.0043028150520144;
    edf_parameters_.LASSNM_F = 45.9617514804616079;
    edf_parameters_.PNM_F = 0.0;

    return edf_parameters_;
}

} // namespace HFBfunctionals


namespace SkyrmeConversion {

/**
 * @brief Store Skyrme force parameters.
 */
struct TParameters {
    double t0_F = 0.0;
    double t1_F = 0.0;
    double t2_F = 0.0;
    double t3_F = 0.0;
    double x0_F = 0.0;
    double x1_F = 0.0;
    double x2_F = 0.0;
    double x3_F = 0.0;
    double b4_F = 0.0;
    double b4p_F = 0.0;
    double te_F = 0.0;
    double to_F = 0.0;
};

/**
 * @brief  Convert Skyrme force parameters to EDF couplings.
 * @math   (t_i,x_i,b_4,b_4',t_e,t_o,α) → \{C_t\}
 * @output Skyrme EDF parameters.
 */
inline EDFParamsSkyrme t_to_C(double t0_F, double t1_F, double t2_F, double t3_F, double x0_F, double x1_F, double x2_F, double x3_F, double b4_F, double b4p_F, double te_F, double to_F, double sigma_F, double CpV0_n_F, double CpV0_p_F, double CpV1_n_F = 0.5, double CpV1_p_F = 0.5) {
    EDFParamsSkyrme edf_parameters_;
    edf_parameters_.sigma_F = sigma_F;

    // (t_i,x_i,b_4,b_4',t_e,t_o) → time-even couplings.
    edf_parameters_.Crho_0_F = 3.0 / 8.0 * t0_F;
    edf_parameters_.Crho_1_F = -0.25 * t0_F * (0.5 + x0_F);
    edf_parameters_.Cdrho_0_F = 1.0 / 16.0 * t3_F;
    edf_parameters_.Cdrho_1_F = -1.0 / 24.0 * t3_F * (0.5 + x3_F);
    edf_parameters_.Ctau_0_F = 3.0 / 16.0 * t1_F + 0.25 * t2_F * (1.25 + x2_F);
    edf_parameters_.Ctau_1_F = -0.125 * t1_F * (0.5 + x1_F) + 0.125 * t2_F * (0.5 + x2_F);
    edf_parameters_.CrDr_0_F = 1.0 / 16.0 * t2_F * (1.25 + x2_F) - 9.0 / 64.0 * t1_F;
    edf_parameters_.CrDr_1_F = 3.0 / 32.0 * t1_F * (0.5 + x1_F) + 1.0 / 32.0 * t2_F * (0.5 + x2_F);
    edf_parameters_.CrdJ_0_F = -b4_F - 0.5 * b4p_F;
    edf_parameters_.CrdJ_1_F = -0.5 * b4p_F;
    edf_parameters_.CJ_0_F = -1.0 / 16.0 * (t1_F * (2.0 * x1_F - 1.0) + t2_F * (2.0 * x2_F + 1.0) - 5.0 * te_F - 15.0 * to_F);
    edf_parameters_.CJ_1_F = -1.0 / 16.0 * (t2_F - t1_F + 5.0 * te_F - 5.0 * to_F);
    edf_parameters_.CJbar_0_F = -3.0 / 8.0 * (te_F + 3.0 * to_F);
    edf_parameters_.CJbar_1_F = -3.0 / 8.0 * (te_F - to_F);

    // (t_i,x_i,t_e,t_o) → time-odd couplings.
    edf_parameters_.Cs_0_F = -0.25 * t0_F * (0.5 - x0_F);
    edf_parameters_.Cs_1_F = -0.125 * t0_F;
    edf_parameters_.Cds_0_F = -1.0 / 24.0 * t3_F * (0.5 - x3_F);
    edf_parameters_.Cds_1_F = -1.0 / 48.0 * t3_F;
    edf_parameters_.Csj_0_F = -edf_parameters_.Ctau_0_F;
    edf_parameters_.Csj_1_F = -edf_parameters_.Ctau_1_F;
    edf_parameters_.CsDs_0_F = 3.0 / 32.0 * t1_F * (0.5 - x1_F) + 1.0 / 32.0 * t2_F * (0.5 + x2_F);
    edf_parameters_.CsDs_1_F = 3.0 / 64.0 * t1_F + 1.0 / 64.0 * t2_F;
    edf_parameters_.CsdJ_0_F = edf_parameters_.CrdJ_0_F;
    edf_parameters_.CsdJ_1_F = edf_parameters_.CrdJ_1_F;
    edf_parameters_.CT_0_F = -edf_parameters_.CJ_0_F;
    edf_parameters_.CT_1_F = -edf_parameters_.CJ_1_F;
    edf_parameters_.CsNabs_0_F = -3.0 / 32.0 * (te_F - to_F);
    edf_parameters_.CsNabs_1_F = -1.0 / 32.0 * (3.0 * te_F + to_F);
    edf_parameters_.CF_0_F = -2.0 * edf_parameters_.CJbar_0_F;
    edf_parameters_.CF_1_F = -2.0 * edf_parameters_.CJbar_1_F;

    // (C_n^{V0},C_p^{V0},C_n^{V1},C_p^{V1}) → pairing couplings.
    edf_parameters_.CpV0_0_F = CpV0_n_F;
    edf_parameters_.CpV0_1_F = CpV0_p_F;
    edf_parameters_.CpV1_0_F = CpV1_n_F;
    edf_parameters_.CpV1_1_F = CpV1_p_F;

    return edf_parameters_;
}

/**
 * @brief  Convert EDF couplings to Skyrme force parameters.
 * @math   \{C_t\} → (t_i,x_i,b_4,b_4',t_e,t_o)
 * @output Skyrme force parameters.
 */
inline TParameters C_to_t(const EDFParamsSkyrme& edf_parameters_) {
    TParameters t_parameters_;

    // C → (t_0,t_1,t_2,t_3).
    t_parameters_.t0_F = 8.0 / 3.0 * edf_parameters_.Crho_0_F;
    t_parameters_.t1_F = 4.0 / 3.0 * (edf_parameters_.Ctau_0_F - 4.0 * edf_parameters_.CrDr_0_F);
    t_parameters_.t2_F = 4.0 / 3.0 * (3.0 * edf_parameters_.Ctau_0_F - 6.0 * edf_parameters_.Ctau_1_F + 4.0 * edf_parameters_.CrDr_0_F - 8.0 * edf_parameters_.CrDr_1_F);
    t_parameters_.t3_F = 16.0 * edf_parameters_.Cdrho_0_F;

    // C → (x_0,x_1,x_2,x_3).
    t_parameters_.x0_F = -0.5 * (3.0 * edf_parameters_.Crho_1_F / edf_parameters_.Crho_0_F + 1.0);
    t_parameters_.x1_F = 2.0 * (-edf_parameters_.Ctau_0_F - 3.0 * edf_parameters_.Ctau_1_F + 4.0 * edf_parameters_.CrDr_0_F + 12.0 * edf_parameters_.CrDr_1_F) / t_parameters_.t1_F / 3.0;
    t_parameters_.x2_F = -2.0 * (3.0 * edf_parameters_.Ctau_0_F - 15.0 * edf_parameters_.Ctau_1_F + 4.0 * edf_parameters_.CrDr_0_F - 20.0 * edf_parameters_.CrDr_1_F) / t_parameters_.t2_F / 3.0;
    t_parameters_.x3_F = -0.5 * (3.0 * edf_parameters_.Cdrho_1_F / edf_parameters_.Cdrho_0_F + 1.0);

    // C → (b_4,b_4',t_e,t_o).
    t_parameters_.b4_F = edf_parameters_.CrdJ_1_F - edf_parameters_.CrdJ_0_F;
    t_parameters_.b4p_F = -2.0 * edf_parameters_.CrdJ_1_F;
    t_parameters_.te_F = 4.0 / 15.0 * (3.0 * edf_parameters_.CJ_0_F - 9.0 * edf_parameters_.CJ_1_F - 4.0 * edf_parameters_.CrDr_0_F + 12.0 * edf_parameters_.CrDr_1_F - 2.0 * edf_parameters_.Ctau_0_F + 6.0 * edf_parameters_.Ctau_1_F);
    t_parameters_.to_F = 4.0 / 15.0 * (3.0 * edf_parameters_.CJ_0_F + 3.0 * edf_parameters_.CJ_1_F + 4.0 * edf_parameters_.CrDr_0_F + 4.0 * edf_parameters_.CrDr_1_F);

    return t_parameters_;
}

} // namespace SkyrmeConversion


/**
 * @brief  Print Skyrme EDF and force parameters.
 * @math   EDF → stdout
 * @output Fixed-precision parameter tables.
 */
inline void print_edf(const EDFParamsSkyrme& edf_skyrme_) {
    // EDF identifier → table header.
    std::cout << "\n";
    std::cout << "  " << edf_skyrme_.functionalName_Str << " functional\n";
    std::cout << "  ----------------------------------------\n";

    /**
     * @brief  Print two fixed-precision values.
     * @math   (l_0,v_0,l_1,v_1) → row
     * @output One formatted stream row.
     */
    const auto print_pair_Func = [](std::ostream& output_, const char* label0_Str, double value0_F, const char* label1_Str, double value1_F) {
        const int label_width_I = 10;
        const int value_width_I = 24;
        const int precision_I = 15;
        output_ << "  " << std::left << std::setw(label_width_I) << label0_Str << "= " << std::right << std::setw(value_width_I) << std::fixed << std::setprecision(precision_I) << value0_F;
        output_ << "    ;";
        output_ << "  " << std::left << std::setw(label_width_I) << label1_Str << "= " << std::right << std::setw(value_width_I) << std::fixed << std::setprecision(precision_I) << value1_F;
        output_ << "\n";
    };

    // C_t → coupling table.
    print_pair_Func(std::cout, "Crho(0)", edf_skyrme_.Crho_0_F, "Crho(1)", edf_skyrme_.Crho_1_F);
    print_pair_Func(std::cout, "CDrho(0)", edf_skyrme_.Cdrho_0_F,  "CDrho(1)", edf_skyrme_.Cdrho_1_F);
    print_pair_Func(std::cout, "Ctau(0)", edf_skyrme_.Ctau_0_F, "Ctau(1)", edf_skyrme_.Ctau_1_F);
    print_pair_Func(std::cout, "CrDr(0)", edf_skyrme_.CrDr_0_F, "CrDr(1)", edf_skyrme_.CrDr_1_F);
    print_pair_Func(std::cout, "CrdJ(0)", edf_skyrme_.CrdJ_0_F, "CrdJ(1)", edf_skyrme_.CrdJ_1_F);
    print_pair_Func(std::cout, "CJ(0)", edf_skyrme_.CJ_0_F, "CJ(1)", edf_skyrme_.CJ_1_F);
    print_pair_Func(std::cout, "CpV0(0)", edf_skyrme_.CpV0_0_F, "CpV0(1)", edf_skyrme_.CpV0_1_F);
    print_pair_Func(std::cout, "CpV1(0)", edf_skyrme_.CpV1_0_F, "CpV1(1)", edf_skyrme_.CpV1_1_F);
    std::cout << "  " << std::left << std::setw(10) << "sigma" << "= " << std::right << std::setw(24) << std::fixed << std::setprecision(15) << edf_skyrme_.sigma_F << "    ;" << "  " << std::left << std::setw(10) << "hbzero" << "= " << std::right << std::setw(24) << std::fixed << std::setprecision(15) << edf_skyrme_.hbzero_F << "\n";

    // EDF options → option table.
    std::cout << "  functional has DME couplings: F\n";
    std::cout << "  use 3-Nucleon DME couplings: F\n";
    std::cout << "  e^2 chrg" << "= " << std::right << std::setw(24) << std::fixed << std::setprecision(15) << edf_skyrme_.e2charg_F << "    ;" << "  " << std::left << std::setw(10) << "CExPar" << "= " << std::right << std::setw(24) << std::fixed << std::setprecision(15) << edf_skyrme_.CExPar_F << "\n";
    std::cout << "  c.m. correction: F, chr. density in direct Coul: F\n";
    std::cout << "  use tensor terms: F\n";
    std::cout << "\n";

    // NM observables → nuclear-matter table.
    std::cout << "  Nuclear matter properties\n";
    std::cout << "  ----------------------------------------\n";
    print_pair_Func(std::cout, "E_NM", edf_skyrme_.ENM_F, "K_NM", edf_skyrme_.KNM_F);
    print_pair_Func(std::cout, "P_NM", 0.0, "RHO_NM", edf_skyrme_.rhoNM_F);
    print_pair_Func(std::cout, "ASS_NM", edf_skyrme_.ASSNM_F, "LASS_NM", edf_skyrme_.LASSNM_F);
    print_pair_Func(std::cout, "SMASS_NM", edf_skyrme_.SMASSNM_F, "VMASS_NM", edf_skyrme_.VMASSNM_F);
    std::cout << "\n";

    // C → (t,x) table.
    std::cout << "  Associated (t,x)-coupling constants\n";
    std::cout << "  ----------------------------------------\n";
    const SkyrmeConversion::TParameters t_parameters_ = SkyrmeConversion::C_to_t(edf_skyrme_);
    print_pair_Func(std::cout, "t0", t_parameters_.t0_F, "x0", t_parameters_.x0_F);
    print_pair_Func(std::cout, "t1", t_parameters_.t1_F, "x1", t_parameters_.x1_F);
    print_pair_Func(std::cout, "t2", t_parameters_.t2_F, "x2", t_parameters_.x2_F);
    print_pair_Func(std::cout, "t3", t_parameters_.t3_F, "x3", t_parameters_.x3_F);
    print_pair_Func(std::cout, "b4", t_parameters_.b4_F, "b4p", t_parameters_.b4p_F);
    print_pair_Func(std::cout, "te", t_parameters_.te_F, "to", t_parameters_.to_F);
    std::cout << "  " << std::left << std::setw(10) << "sigma" << "= " << std::right << std::setw(24) << std::fixed << std::setprecision(15) << edf_skyrme_.sigma_F << "    ;" << "  " << std::left << std::setw(10) << "hbzero" << "= " << std::right << std::setw(24) << std::fixed << std::setprecision(15) << edf_skyrme_.hbzero_F << "\n";
    std::cout << std::endl;
}
