/**
 * @file    hfb_projection_nucleus.hpp
 * @author  cailea
 * @date    2026-09-17
 * @brief   Joint neutron-proton projection interfaces.
 */

#pragma once

#include <cassert>
#include <functional>
#include <utility>

#include "hfb_projection.hpp"

/**
 * @brief  Project neutron-proton configurations using gauge and Euler quadrature.
 * @math   |Φκ⟩ = |Φν,κν⟩ ⊗ |Φπ,κπ⟩.
 * @note   τ = 0: neutron; τ = 1: proton.
 * @note   Species-local orbital indices start at zero.
 * @note   Both species share the Euler-angle nodes.
 * @note   Returned references alias nuclear output workspaces.
 * @note   Calculation implementations remain pending.
 */
class HFBProjectionNucleus {
public:
    // ⟨τ₁α₁|O|τ₂α₂⟩.
    using OneBodyElementFunc = std::function<doubleC(int tau1_I, int sp1_I, int tau2_I, int sp2_I)>;

    // ⟨τ₁α₁,τ₂α₂|O|τ₃α₃,τ₄α₄⟩: unsymmetrized.
    using TwoBodyElementFunc = std::function<doubleC(int tau1_I, int sp1_I, int tau2_I, int sp2_I, int tau3_I, int sp3_I, int tau4_I, int sp4_I)>;

    // ⟨τ₁α₁|Tλμ|τ₂α₂⟩; μ = -λ,…,λ.
    using MultipoleElementFunc = std::function<doubleC(int tau1_I, int sp1_I, int tau2_I, int sp2_I, int mu_I)>;

    HFBProjection projection_neutron;
    HFBProjection projection_proton;

    Eigen::Tensor<doubleC, 6, Eigen::ColMajor> result_C6D_Ncfg1_Zcfg1_Ncfg2_Zcfg2_K1_K2{};

public:
    /**
     * @brief  Move species projections into owned workspaces.
     * @math   |Φκ⟩ = |Φν,κν⟩ ⊗ |Φπ,κπ⟩.
     * @output Object owning neutron and proton projections.
     */
    HFBProjectionNucleus(HFBProjection projection_neutron_, HFBProjection projection_proton_)
    : projection_neutron(std::move(projection_neutron_)), projection_proton(std::move(projection_proton_)) {
        assert(projection_neutron.TargetTwoI_I == projection_proton.TargetTwoI_I);
        assert(projection_neutron.Nalpha_I == projection_proton.Nalpha_I);
        assert(projection_neutron.Nbeta_I == projection_proton.Nbeta_I);
        assert(projection_neutron.Ngamma_I == projection_proton.Ngamma_I);
    }

    /**
     * @brief  Integrate nuclear overlaps using gauge and Euler quadrature.
     * @math   result_abK₁K₂ = ⟨Φ₁;a|Pνᴺ Pπᶻ Pᴵ_{K₁K₂}|Φ₂;b⟩.
     * @output Tensor indexed by Ncfg1, Zcfg1, Ncfg2, Zcfg2, K1, K2.
     */
    const Eigen::Tensor<doubleC, 6, Eigen::ColMajor>& calc_overlap();

    /**
     * @brief  Integrate scalar one-body kernels using Pfaffians and quadrature.
     * @math   O¹ = Σ₁₂ ⟨1|O|2⟩ c₁†c₂; 1 = (τ₁,α₁).
     * @output Tensor indexed by Ncfg1, Zcfg1, Ncfg2, Zcfg2, K1, K2.
     * @note   Includes species-conserving and charge-exchange channels.
     */
    const Eigen::Tensor<doubleC, 6, Eigen::ColMajor>& calc_one_body(const OneBodyElementFunc& read_element_Func);

    /**
     * @brief  Integrate scalar two-body kernels using Pfaffians and quadrature.
     * @math   O² = ½Σ₁₂₃₄ ⟨12|O|34⟩ c₁†c₂†c₄c₃.
     * @output Tensor indexed by Ncfg1, Zcfg1, Ncfg2, Zcfg2, K1, K2.
     * @note   Callback supplies unsymmetrized elements; integration includes ½.
     */
    const Eigen::Tensor<doubleC, 6, Eigen::ColMajor>& calc_two_body(const TwoBodyElementFunc& read_element_Func);

    /**
     * @brief  Integrate reduced multipoles using quadrature and Clebsch–Gordan coupling.
     * @math   result_abK₁K₂ = ⟨Φ₁^{I₁};K₁,a‖Tλ‖Φ₂^{I₂};K₂,b⟩.
     * @output Tensor indexed by Ncfg1, Zcfg1, Ncfg2, Zcfg2, K1, K2.
     * @note   TwoI1 = 2I₁; TwoI2 = 2I₂; twoLambda = 2λ.
     * @note   Callback mu_I is the magnetic quantum number.
     * @note   States are unnormalized and unmixed.
     */
    const Eigen::Tensor<doubleC, 6, Eigen::ColMajor>& calc_multipole(int TwoI1_I, int TwoI2_I, int twoLambda_I, const MultipoleElementFunc& read_element_Func);
};
