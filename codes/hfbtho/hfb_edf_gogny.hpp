/**
 * @file    hfb_edf_gogny.hpp
 * @author  cailea
 * @date    2026-05-22
 * @brief   Define Gogny-force parameters and local EDF couplings.
 */

#pragma once

#include <array>
#include <string>

#include "hfb_edf_skyrme.hpp"

/**
 * @brief Store two-range Gogny-force parameters.
 */
class EDFParamsGogny {
public:
    static constexpr int Ng_I = 2;
    using GaussianValues = std::array<double, Ng_I>;

    std::string forceName_Str{"Unknown"};
    GaussianValues mu_F1D_g{}; // μ_g [fm].
    GaussianValues W_F1D_g{}; // W_g [MeV].
    GaussianValues B_F1D_g{}; // B_g [MeV].
    GaussianValues H_F1D_g{}; // H_g [MeV].
    GaussianValues M_F1D_g{}; // M_g [MeV].
    double t3G_F = 0.0; // t_3^G [MeV fm^{3+3α}].
    double x3_F = 1.0; // x_3.
    double alpha_F = 1.0 / 3.0; // α.
    double WLS_F = 0.0; // W_{LS} [MeV fm⁵].
    double hbzero_F = 20.73667552957479; // ℏ²/(2m) [MeV fm²].

public:
    /**
     * @brief  Construct an empty Gogny parameter set.
     * @math   ∅ → P_G
     * @output Zero-initialized Gogny parameters.
     */
    EDFParamsGogny() = default;

    /**
     * @brief  Construct a two-range Gogny parameter set.
     * @math   {μ_g,W_g,B_g,H_g,M_g,t_3^G,x_3,α,W_{LS}} → P_G
     * @output Initialized Gogny parameters.
     */
    EDFParamsGogny(const std::string& forceName_Str_, const GaussianValues& mu_F1D_g_, const GaussianValues& W_F1D_g_, const GaussianValues& B_F1D_g_, const GaussianValues& H_F1D_g_, const GaussianValues& M_F1D_g_, double t3G_F_, double x3_F_, double alpha_F_, double WLS_F_, double hbzero_F_) {
        // {μ_g,W_g,B_g,H_g,M_g,t_3^G,x_3,α,W_{LS}} → P_G.
        forceName_Str = forceName_Str_;
        mu_F1D_g = mu_F1D_g_;
        W_F1D_g = W_F1D_g_;
        B_F1D_g = B_F1D_g_;
        H_F1D_g = H_F1D_g_;
        M_F1D_g = M_F1D_g_;
        t3G_F = t3G_F_;
        x3_F = x3_F_;
        alpha_F = alpha_F_;
        WLS_F = WLS_F_;
        hbzero_F = hbzero_F_;
    }

    /**
     * @brief  Convert local Gogny terms into EDF couplings.
     * @math   (t_3^G,x_3,α,W_{LS}) → {C_t}
     * @output Local companion EDF.
     */
    EDFParamsSkyrme make_local_edf() const;

    /**
     * @brief  Build the D1 Gogny parameterization.
     * @math   D1 → P_G
     * @output D1 parameters.
     */
    static EDFParamsGogny D1() {
        return EDFParamsGogny("D1", GaussianValues{0.7, 1.2}, GaussianValues{-402.40, -21.30}, GaussianValues{-100.00, -11.77}, GaussianValues{-496.20, 37.27}, GaussianValues{-23.56, -68.81}, 1350.00, 1.0, 1.0 / 3.0, 115.000, 20.73667552957479);
    }

    /**
     * @brief  Build the D1S Gogny parameterization.
     * @math   D1S → P_G
     * @output D1S parameters.
     */
    static EDFParamsGogny D1S() {
        return EDFParamsGogny("D1S", GaussianValues{0.7, 1.2}, GaussianValues{-1720.30, 103.64}, GaussianValues{1300.00, -163.48}, GaussianValues{-1813.53, 162.81}, GaussianValues{1397.60, -223.93}, 1390.600, 1.0, 1.0 / 3.0, 130.000, 20.73667622931579050281);
    }

    /**
     * @brief  Build the D1N Gogny parameterization.
     * @math   D1N → P_G
     * @output D1N parameters.
     */
    static EDFParamsGogny D1N() {
        return EDFParamsGogny("D1N", GaussianValues{0.8, 1.2}, GaussianValues{-2047.61, 293.02}, GaussianValues{1700.00, -300.78}, GaussianValues{-2414.93, 414.59}, GaussianValues{1519.35, -316.84}, 1609.50, 1.0, 1.0 / 3.0, 115.000, 20.73667552957479);
    }
};

/**
 * @brief  Convert local Gogny terms into EDF couplings.
 * @math   (t_3^G,x_3,α,W_{LS}) → {C_t}
 * @output Local companion EDF.
 */
inline EDFParamsSkyrme EDFParamsGogny::make_local_edf() const {
    EDFParamsSkyrme local_edf_;
    // (ℏ²/2m,e²,α,C_ex) ← P_G.
    local_edf_.functionalName_Str = forceName_Str + "-local";
    local_edf_.hbzero_F = hbzero_F;
    local_edf_.hbzeron_F = hbzero_F;
    local_edf_.hbzerop_F = hbzero_F;
    local_edf_.e2charg_F = 1.439978408596513;
    local_edf_.CExPar_F = 1.0;
    local_edf_.sigma_F = alpha_F;

    // t_3=6t_3^G → {C_t^{ρ,α},C_t^{s,α}}.
    const double t3_F = 6.0 * t3G_F;
    local_edf_.Cdrho_0_F = t3_F / 16.0;
    local_edf_.Cdrho_1_F = -t3_F * (0.5 + x3_F) / 24.0;
    local_edf_.Cds_0_F = -t3_F * (0.5 - x3_F) / 24.0;
    local_edf_.Cds_1_F = -t3_F / 48.0;

    // b_4=b'_4=W_{LS}/2 → C_t^{∇J}=C_t^{∇j}.
    const double b4_F = 0.5 * WLS_F;
    const double b4p_F = 0.5 * WLS_F;
    local_edf_.CrdJ_0_F = -b4_F - 0.5 * b4p_F;
    local_edf_.CrdJ_1_F = -0.5 * b4p_F;
    local_edf_.CsdJ_0_F = local_edf_.CrdJ_0_F;
    local_edf_.CsdJ_1_F = local_edf_.CrdJ_1_F;

    // C_q^{V0}=C_q^{V1}=0.
    local_edf_.CpV0_0_F = 0.0;
    local_edf_.CpV0_1_F = 0.0;
    local_edf_.CpV1_0_F = 0.0;
    local_edf_.CpV1_1_F = 0.0;
    return local_edf_;
}
