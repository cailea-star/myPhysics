/**
 * @file    hfb_projection_nucleus.hpp
 * @author  cailea
 * @date    2026-09-17
 * @brief   Joint neutron-proton projection interfaces.
 */

#pragma once

#include <cassert>
#include <cmath>
#include <utility>

#include "hfb_projection.hpp"

/**
 * @brief  Project neutron-proton configurations using gauge and Euler quadrature.
 * @math   |Φκ⟩ = |Φν,κν⟩ ⊗ |Φπ,κπ⟩.
 * @note   τ = 0: neutron; τ = 1: proton.
 * @note   Species-local indices start at zero; Euler meshes are shared.
 * @note   Returned references alias nuclear output workspaces.
 * @note   Particle-number targets apply to the right state.
 * @note   Requires even reference vacua and initialized projection data.
 * @note   Vacuum overlaps must remain nonzero at quadrature nodes.
 */
class HFBProjectionNucleus {
public:
    int TargetTwoI_I = 0;

    HFBProjection projection_neutron;
    HFBProjection projection_proton;

private:
    Eigen::Index Ncfg1pcfg2p_I = 0;
    Eigen::Index Ncfg1ncfg2n_I = 0;
    Eigen::Index Ncfg1pcfg2pcfg1ncfg2n_I = 0;

    RepresentationSpin representation_spin;

    Eigen::Tensor<doubleC, 4, Eigen::ColMajor> tmp_C4D_cfg1p_cfg2p_cfg1n_cfg2n{};
    Eigen::Tensor<doubleC, 6, Eigen::ColMajor> result_C6D_cfg1p_cfg2p_cfg1n_cfg2n_K1_K2{};

    public:
    /**
     * @brief  Move species projections into owned workspaces.
     * @math   |Φκ⟩ = |Φν,κν⟩ ⊗ |Φπ,κπ⟩.
     * @output Object owning neutron and proton projections.
     */
    HFBProjectionNucleus(HFBProjection projection_neutron_, HFBProjection projection_proton_)
    : projection_neutron(std::move(projection_neutron_)), projection_proton(std::move(projection_proton_)), representation_spin(projection_neutron.TargetTwoI_I) {
        assert(projection_neutron.TargetTwoI_I == projection_proton.TargetTwoI_I);
        assert(projection_neutron.Nalpha_I == projection_proton.Nalpha_I);
        assert(projection_neutron.Nbeta_I == projection_proton.Nbeta_I);
        assert(projection_neutron.Ngamma_I == projection_proton.Ngamma_I);

        TargetTwoI_I = projection_neutron.TargetTwoI_I;
        Ncfg1pcfg2p_I = static_cast<Eigen::Index>(projection_proton.Ncfg1_I) * projection_proton.Ncfg2_I;
        Ncfg1ncfg2n_I = static_cast<Eigen::Index>(projection_neutron.Ncfg1_I) * projection_neutron.Ncfg2_I;
        Ncfg1pcfg2pcfg1ncfg2n_I = Ncfg1pcfg2p_I * Ncfg1ncfg2n_I;

        // result ∈ ℂ^{cfg1p×cfg2p×cfg1n×cfg2n×(2I+1)×(2I+1)}.
        tmp_C4D_cfg1p_cfg2p_cfg1n_cfg2n.resize(projection_proton.Ncfg1_I, projection_proton.Ncfg2_I, projection_neutron.Ncfg1_I, projection_neutron.Ncfg2_I);
        result_C6D_cfg1p_cfg2p_cfg1n_cfg2n_K1_K2.resize(projection_proton.Ncfg1_I, projection_proton.Ncfg2_I, projection_neutron.Ncfg1_I, projection_neutron.Ncfg2_I, TargetTwoI_I + 1, TargetTwoI_I + 1);
    }

    /**
     * @brief  Integrate nuclear overlaps using gauge and Euler quadrature.
     * @math   result_abK₁K₂ = ⟨Φ₁;a|Pνᴺ Pπᶻ Pᴵ_{K₁K₂}|Φ₂;b⟩.
     * @output Tensor indexed by cfg1p, cfg2p, cfg1n, cfg2n, K1, K2.
     */
    const Eigen::Tensor<doubleC, 6, Eigen::ColMajor>& calc_overlap();

    /**
     * @brief  Integrate NN kernels using Pfaffians and quadrature.
     * @math   O_NN = Σαβ Oαβ cνα† cνβ.
     * @output Tensor indexed by cfg1p, cfg2p, cfg1n, cfg2n, K1, K2.
     * @note   Requires rotationally scalar matrix elements.
     */
    const Eigen::Tensor<doubleC, 6, Eigen::ColMajor>& calc_one_body_NN(const Eigen::MatrixXcd& OneBody_C2D_spn_spn);

    /**
     * @brief  Integrate PP kernels using Pfaffians and quadrature.
     * @math   O_PP = Σαβ Oαβ cπα† cπβ.
     * @output Tensor indexed by cfg1p, cfg2p, cfg1n, cfg2n, K1, K2.
     * @note   Requires rotationally scalar matrix elements.
     */
    const Eigen::Tensor<doubleC, 6, Eigen::ColMajor>& calc_one_body_PP(const Eigen::MatrixXcd& OneBody_C2D_spp_spp);

    /**
     * @brief  Integrate PN kernels using Pfaffians and quadrature.
     * @math   O_PN = Σαβ Oαβ cπα† cνβ.
     * @output Tensor indexed by cfg1p, cfg2p, cfg1n, cfg2n, K1, K2.
     * @note   Requires rotationally scalar matrix elements.
     */
    const Eigen::Tensor<doubleC, 6, Eigen::ColMajor>& calc_one_body_PN(const Eigen::MatrixXcd& OneBody_C2D_spp_spn);

    /**
     * @brief  Integrate NNNN kernels using Pfaffians and quadrature.
     * @math   O_NNNN = ½Σ₁₂₃₄ V₁₂₃₄ cν₁† cν₂† cν₄ cν₃.
     * @output Tensor indexed by cfg1p, cfg2p, cfg1n, cfg2n, K1, K2.
     * @note   V contains unsymmetrized, rotationally scalar matrix elements.
     */
    const Eigen::Tensor<doubleC, 6, Eigen::ColMajor>& calc_two_body_NNNN(const Eigen::Tensor<double, 4, Eigen::ColMajor>& TwoBody_F4D_spn_spn_spn_spn);

    /**
     * @brief  Integrate PPPP kernels using Pfaffians and quadrature.
     * @math   O_PPPP = ½Σ₁₂₃₄ V₁₂₃₄ cπ₁† cπ₂† cπ₄ cπ₃.
     * @output Tensor indexed by cfg1p, cfg2p, cfg1n, cfg2n, K1, K2.
     * @note   V contains unsymmetrized, rotationally scalar matrix elements.
     */
    const Eigen::Tensor<doubleC, 6, Eigen::ColMajor>& calc_two_body_PPPP(const Eigen::Tensor<double, 4, Eigen::ColMajor>& TwoBody_F4D_spp_spp_spp_spp);

    /**
     * @brief  Integrate NPPN kernels using Pfaffians and quadrature.
     * @math   O_NPPN = Σ₁₂₃₄ W₁₂₃₄ cν₁† cπ₂† cπ₃ cν₄.
     * @math   W₁₂₃₄ = ½(Vνπνπ₁₂₄₃ + Vπνπν₂₁₃₄).
     * @output Tensor indexed by cfg1p, cfg2p, cfg1n, cfg2n, K1, K2.
     * @note   Input axes: π₂,π₃,ν₁,ν₄; operator order: ν₁,π₂,π₃,ν₄.
     * @note   Input W combines both orderings; no additional ½.
     * @note   Requires rotationally scalar matrix elements.
     */
    const Eigen::Tensor<doubleC, 6, Eigen::ColMajor>& calc_two_body_NPPN(const Eigen::Tensor<double, 4, Eigen::ColMajor>& TwoBody_F4D_spp2_spp3_spn1_spn4);
};

inline const Eigen::Tensor<doubleC, 6, Eigen::ColMajor>& HFBProjectionNucleus::calc_overlap() {
    auto& neutron = projection_neutron;
    auto& proton = projection_proton;

    // Ων = Ωπ; wν = wπ.
    assert(neutron.alpha_F1D_alpha.isApprox(proton.alpha_F1D_alpha));
    assert(neutron.beta_F1D_beta.isApprox(proton.beta_F1D_beta));
    assert(neutron.gamma_F1D_gamma.isApprox(proton.gamma_F1D_gamma));
    assert(neutron.weight_F1D_alpha.isApprox(proton.weight_F1D_alpha));
    assert(neutron.weight_F1D_beta.isApprox(proton.weight_F1D_beta));
    assert(neutron.weight_F1D_gamma.isApprox(proton.weight_F1D_gamma));

    // VΩ = Σwα Σwβ Σwγ; wβ includes sinβ.
    const double volume_F = neutron.weight_F1D_alpha.sum() * neutron.weight_F1D_beta.sum() * neutron.weight_F1D_gamma.sum();
    assert(std::isfinite(volume_F) && volume_F > 0.0);
    const double normalization_F = (TargetTwoI_I + 1.0) / volume_F;

    result_C6D_cfg1p_cfg2p_cfg1n_cfg2n_K1_K2.setZero();

    // dᴵ(β) is shared by both species.
    for (int beta_I = 0; beta_I < neutron.Nbeta_I; ++beta_I) {
        const auto& RyI_C2D_K_K = representation_spin.calc_Ry(neutron.beta_F1D_beta(beta_I));
        for (int alpha_I = 0; alpha_I < neutron.Nalpha_I; ++alpha_I) {
            for (int gamma_I = 0; gamma_I < neutron.Ngamma_I; ++gamma_I) {
                neutron.prepare_rotation(alpha_I, beta_I, gamma_I);
                proton.prepare_rotation(alpha_I, beta_I, gamma_I);

                const auto& overlap_C2D_cfg1n_cfg2n = neutron.calc_overlap_pnp();
                const auto& overlap_C2D_cfg1p_cfg2p = proton.calc_overlap_pnp();

                // ColMajor: (p₁,p₂) → rows; (n₁,n₂) → columns.
                const Eigen::Map<const Eigen::VectorXcd> overlap_C1D_cfg1ncfg2n(overlap_C2D_cfg1n_cfg2n.data(), Ncfg1ncfg2n_I);
                const Eigen::Map<const Eigen::VectorXcd> overlap_C1D_cfg1pcfg2p(overlap_C2D_cfg1p_cfg2p.data(), Ncfg1pcfg2p_I);
                Eigen::Map<Eigen::MatrixXcd> overlap_C2D_cfg1pcfg2p_cfg1ncfg2n(tmp_C4D_cfg1p_cfg2p_cfg1n_cfg2n.data(), Ncfg1pcfg2p_I, Ncfg1ncfg2n_I);
                overlap_C2D_cfg1pcfg2p_cfg1ncfg2n.noalias() = overlap_C1D_cfg1pcfg2p * overlap_C1D_cfg1ncfg2n.transpose();

                // w = (2I+1) wα wβ wγ / VΩ.
                const double weight_F = normalization_F * neutron.weight_F1D_alpha(alpha_I) * neutron.weight_F1D_beta(beta_I) * neutron.weight_F1D_gamma(gamma_I);
                for (int K2_I = 0; K2_I < TargetTwoI_I + 1; ++K2_I) {
                    const double K2_F = K2_I - 0.5 * TargetTwoI_I;
                    for (int K1_I = 0; K1_I < TargetTwoI_I + 1; ++K1_I) {
                        const double K1_F = K1_I - 0.5 * TargetTwoI_I;

                        // Dᴵ* = exp(iK₁α) dᴵ* exp(iK₂γ).
                        const doubleC factorOmega_C = std::conj(RyI_C2D_K_K(K1_I, K2_I)) * std::exp(doubleC(0.0, K1_F * neutron.alpha_F1D_alpha(alpha_I) + K2_F * neutron.gamma_F1D_gamma(gamma_I)));
                        const doubleC weight_factor_C = weight_F * factorOmega_C;

                        // Each K₁,K₂ block contains all four configuration dimensions.
                        const Eigen::Index offset_I = (static_cast<Eigen::Index>(K2_I) * (TargetTwoI_I + 1) + K1_I) * Ncfg1pcfg2pcfg1ncfg2n_I;
                        Eigen::Map<Eigen::MatrixXcd> result_C2D_cfg1pcfg2p_cfg1ncfg2n(result_C6D_cfg1p_cfg2p_cfg1n_cfg2n_K1_K2.data() + offset_I, Ncfg1pcfg2p_I, Ncfg1ncfg2n_I);
                        result_C2D_cfg1pcfg2p_cfg1ncfg2n += weight_factor_C * overlap_C2D_cfg1pcfg2p_cfg1ncfg2n;
                    }
                }
            }
        }
    }
    return result_C6D_cfg1p_cfg2p_cfg1n_cfg2n_K1_K2;
}

inline const Eigen::Tensor<doubleC, 6, Eigen::ColMajor>& HFBProjectionNucleus::calc_one_body_NN(const Eigen::MatrixXcd& OneBody_C2D_spn_spn) {
    auto& neutron = projection_neutron;
    auto& proton = projection_proton;
    assert(OneBody_C2D_spn_spn.rows() == neutron.Nsp_I && OneBody_C2D_spn_spn.cols() == neutron.Nsp_I && OneBody_C2D_spn_spn.allFinite());

    // Ων = Ωπ; wν = wπ.
    assert(neutron.alpha_F1D_alpha.isApprox(proton.alpha_F1D_alpha));
    assert(neutron.beta_F1D_beta.isApprox(proton.beta_F1D_beta));
    assert(neutron.gamma_F1D_gamma.isApprox(proton.gamma_F1D_gamma));
    assert(neutron.weight_F1D_alpha.isApprox(proton.weight_F1D_alpha));
    assert(neutron.weight_F1D_beta.isApprox(proton.weight_F1D_beta));
    assert(neutron.weight_F1D_gamma.isApprox(proton.weight_F1D_gamma));

    // VΩ = Σwα Σwβ Σwγ; wβ includes sinβ.
    const double volume_F = neutron.weight_F1D_alpha.sum() * neutron.weight_F1D_beta.sum() * neutron.weight_F1D_gamma.sum();
    assert(std::isfinite(volume_F) && volume_F > 0.0);
    const double normalization_F = (TargetTwoI_I + 1.0) / volume_F;

    result_C6D_cfg1p_cfg2p_cfg1n_cfg2n_K1_K2.setZero();

    // dᴵ(β) is shared by both species.
    for (int beta_I = 0; beta_I < neutron.Nbeta_I; ++beta_I) {
        const auto& RyI_C2D_K_K = representation_spin.calc_Ry(neutron.beta_F1D_beta(beta_I));
        for (int alpha_I = 0; alpha_I < neutron.Nalpha_I; ++alpha_I) {
            for (int gamma_I = 0; gamma_I < neutron.Ngamma_I; ++gamma_I) {
                neutron.prepare_rotation(alpha_I, beta_I, gamma_I);
                proton.prepare_rotation(alpha_I, beta_I, gamma_I);

                const auto& overlap_C2D_cfg1p_cfg2p = proton.calc_overlap_pnp();
                const auto& onebody_C2D_cfg1n_cfg2n = neutron.calc_one_body_pnp(OneBody_C2D_spn_spn);

                // ColMajor: (p₁,p₂) → rows; (n₁,n₂) → columns.
                const Eigen::Map<const Eigen::VectorXcd> onebody_C1D_cfg1ncfg2n(onebody_C2D_cfg1n_cfg2n.data(), Ncfg1ncfg2n_I);
                const Eigen::Map<const Eigen::VectorXcd> overlap_C1D_cfg1pcfg2p(overlap_C2D_cfg1p_cfg2p.data(), Ncfg1pcfg2p_I);
                Eigen::Map<Eigen::MatrixXcd> onebody_C2D_cfg1pcfg2p_cfg1ncfg2n(tmp_C4D_cfg1p_cfg2p_cfg1n_cfg2n.data(), Ncfg1pcfg2p_I, Ncfg1ncfg2n_I);
                onebody_C2D_cfg1pcfg2p_cfg1ncfg2n.noalias() = overlap_C1D_cfg1pcfg2p * onebody_C1D_cfg1ncfg2n.transpose();

                // w = (2I+1) wα wβ wγ / VΩ.
                const double weight_F = normalization_F * neutron.weight_F1D_alpha(alpha_I) * neutron.weight_F1D_beta(beta_I) * neutron.weight_F1D_gamma(gamma_I);
                for (int K2_I = 0; K2_I < TargetTwoI_I + 1; ++K2_I) {
                    const double K2_F = K2_I - 0.5 * TargetTwoI_I;
                    for (int K1_I = 0; K1_I < TargetTwoI_I + 1; ++K1_I) {
                        const double K1_F = K1_I - 0.5 * TargetTwoI_I;

                        // Dᴵ* = exp(iK₁α) dᴵ* exp(iK₂γ).
                        const doubleC factorOmega_C = std::conj(RyI_C2D_K_K(K1_I, K2_I)) * std::exp(doubleC(0.0, K1_F * neutron.alpha_F1D_alpha(alpha_I) + K2_F * neutron.gamma_F1D_gamma(gamma_I)));
                        const doubleC weight_factor_C = weight_F * factorOmega_C;

                        // Each K₁,K₂ block contains all four configuration dimensions.
                        const Eigen::Index offset_I = (static_cast<Eigen::Index>(K2_I) * (TargetTwoI_I + 1) + K1_I) * Ncfg1pcfg2pcfg1ncfg2n_I;
                        Eigen::Map<Eigen::MatrixXcd> result_C2D_cfg1pcfg2p_cfg1ncfg2n(result_C6D_cfg1p_cfg2p_cfg1n_cfg2n_K1_K2.data() + offset_I, Ncfg1pcfg2p_I, Ncfg1ncfg2n_I);
                        result_C2D_cfg1pcfg2p_cfg1ncfg2n += weight_factor_C * onebody_C2D_cfg1pcfg2p_cfg1ncfg2n;
                    }
                }
            }
        }
    }
    return result_C6D_cfg1p_cfg2p_cfg1n_cfg2n_K1_K2;
}

inline const Eigen::Tensor<doubleC, 6, Eigen::ColMajor>& HFBProjectionNucleus::calc_one_body_PP(const Eigen::MatrixXcd& OneBody_C2D_spp_spp) {
    auto& neutron = projection_neutron;
    auto& proton = projection_proton;
    assert(OneBody_C2D_spp_spp.rows() == proton.Nsp_I && OneBody_C2D_spp_spp.cols() == proton.Nsp_I && OneBody_C2D_spp_spp.allFinite());

    // Ων = Ωπ; wν = wπ.
    assert(neutron.alpha_F1D_alpha.isApprox(proton.alpha_F1D_alpha));
    assert(neutron.beta_F1D_beta.isApprox(proton.beta_F1D_beta));
    assert(neutron.gamma_F1D_gamma.isApprox(proton.gamma_F1D_gamma));
    assert(neutron.weight_F1D_alpha.isApprox(proton.weight_F1D_alpha));
    assert(neutron.weight_F1D_beta.isApprox(proton.weight_F1D_beta));
    assert(neutron.weight_F1D_gamma.isApprox(proton.weight_F1D_gamma));

    // VΩ = Σwα Σwβ Σwγ; wβ includes sinβ.
    const double volume_F = neutron.weight_F1D_alpha.sum() * neutron.weight_F1D_beta.sum() * neutron.weight_F1D_gamma.sum();
    assert(std::isfinite(volume_F) && volume_F > 0.0);
    const double normalization_F = (TargetTwoI_I + 1.0) / volume_F;

    result_C6D_cfg1p_cfg2p_cfg1n_cfg2n_K1_K2.setZero();

    // dᴵ(β) is shared by both species.
    for (int beta_I = 0; beta_I < neutron.Nbeta_I; ++beta_I) {
        const auto& RyI_C2D_K_K = representation_spin.calc_Ry(neutron.beta_F1D_beta(beta_I));
        for (int alpha_I = 0; alpha_I < neutron.Nalpha_I; ++alpha_I) {
            for (int gamma_I = 0; gamma_I < neutron.Ngamma_I; ++gamma_I) {
                neutron.prepare_rotation(alpha_I, beta_I, gamma_I);
                proton.prepare_rotation(alpha_I, beta_I, gamma_I);

                const auto& overlap_C2D_cfg1n_cfg2n = neutron.calc_overlap_pnp();
                const auto& onebody_C2D_cfg1p_cfg2p = proton.calc_one_body_pnp(OneBody_C2D_spp_spp);

                // ColMajor: (p₁,p₂) → rows; (n₁,n₂) → columns.
                const Eigen::Map<const Eigen::VectorXcd> overlap_C1D_cfg1ncfg2n(overlap_C2D_cfg1n_cfg2n.data(), Ncfg1ncfg2n_I);
                const Eigen::Map<const Eigen::VectorXcd> onebody_C1D_cfg1pcfg2p(onebody_C2D_cfg1p_cfg2p.data(), Ncfg1pcfg2p_I);
                Eigen::Map<Eigen::MatrixXcd> onebody_C2D_cfg1pcfg2p_cfg1ncfg2n(tmp_C4D_cfg1p_cfg2p_cfg1n_cfg2n.data(), Ncfg1pcfg2p_I, Ncfg1ncfg2n_I);
                onebody_C2D_cfg1pcfg2p_cfg1ncfg2n.noalias() = onebody_C1D_cfg1pcfg2p * overlap_C1D_cfg1ncfg2n.transpose();

                // w = (2I+1) wα wβ wγ / VΩ.
                const double weight_F = normalization_F * neutron.weight_F1D_alpha(alpha_I) * neutron.weight_F1D_beta(beta_I) * neutron.weight_F1D_gamma(gamma_I);
                for (int K2_I = 0; K2_I < TargetTwoI_I + 1; ++K2_I) {
                    const double K2_F = K2_I - 0.5 * TargetTwoI_I;
                    for (int K1_I = 0; K1_I < TargetTwoI_I + 1; ++K1_I) {
                        const double K1_F = K1_I - 0.5 * TargetTwoI_I;

                        // Dᴵ* = exp(iK₁α) dᴵ* exp(iK₂γ).
                        const doubleC factorOmega_C = std::conj(RyI_C2D_K_K(K1_I, K2_I)) * std::exp(doubleC(0.0, K1_F * neutron.alpha_F1D_alpha(alpha_I) + K2_F * neutron.gamma_F1D_gamma(gamma_I)));
                        const doubleC weight_factor_C = weight_F * factorOmega_C;

                        // Each K₁,K₂ block contains all four configuration dimensions.
                        const Eigen::Index offset_I = (static_cast<Eigen::Index>(K2_I) * (TargetTwoI_I + 1) + K1_I) * Ncfg1pcfg2pcfg1ncfg2n_I;
                        Eigen::Map<Eigen::MatrixXcd> result_C2D_cfg1pcfg2p_cfg1ncfg2n(result_C6D_cfg1p_cfg2p_cfg1n_cfg2n_K1_K2.data() + offset_I, Ncfg1pcfg2p_I, Ncfg1ncfg2n_I);
                        result_C2D_cfg1pcfg2p_cfg1ncfg2n += weight_factor_C * onebody_C2D_cfg1pcfg2p_cfg1ncfg2n;
                    }
                }
            }
        }
    }
    return result_C6D_cfg1p_cfg2p_cfg1n_cfg2n_K1_K2;
}

inline const Eigen::Tensor<doubleC, 6, Eigen::ColMajor>& HFBProjectionNucleus::calc_one_body_PN(const Eigen::MatrixXcd& OneBody_C2D_spp_spn) {
    auto& neutron = projection_neutron;
    auto& proton = projection_proton;
    assert(OneBody_C2D_spp_spn.rows() == proton.Nsp_I && OneBody_C2D_spp_spn.cols() == neutron.Nsp_I && OneBody_C2D_spp_spn.allFinite());

    // Ων = Ωπ; wν = wπ.
    assert(neutron.alpha_F1D_alpha.isApprox(proton.alpha_F1D_alpha));
    assert(neutron.beta_F1D_beta.isApprox(proton.beta_F1D_beta));
    assert(neutron.gamma_F1D_gamma.isApprox(proton.gamma_F1D_gamma));
    assert(neutron.weight_F1D_alpha.isApprox(proton.weight_F1D_alpha));
    assert(neutron.weight_F1D_beta.isApprox(proton.weight_F1D_beta));
    assert(neutron.weight_F1D_gamma.isApprox(proton.weight_F1D_gamma));

    // VΩ = Σwα Σwβ Σwγ; wβ includes sinβ.
    const double volume_F = neutron.weight_F1D_alpha.sum() * neutron.weight_F1D_beta.sum() * neutron.weight_F1D_gamma.sum();
    assert(std::isfinite(volume_F) && volume_F > 0.0);
    const double normalization_F = (TargetTwoI_I + 1.0) / volume_F;

    result_C6D_cfg1p_cfg2p_cfg1n_cfg2n_K1_K2.setZero();

    // ColMajor: (p₁,p₂) → rows; (n₁,n₂) → columns.
    Eigen::Map<Eigen::MatrixXcd> tmp_C2D_cfg1pcfg2p_cfg1ncfg2n(tmp_C4D_cfg1p_cfg2p_cfg1n_cfg2n.data(), Ncfg1pcfg2p_I, Ncfg1ncfg2n_I);

    // Aν Oᵀ ∈ ℂ^{Ncfg1ncfg2n×Nspp}.
    Eigen::MatrixXcd tmp_C2D_cfg1ncfg2n_spp{};
    tmp_C2D_cfg1ncfg2n_spp.resize(Ncfg1ncfg2n_I, proton.Nsp_I);

    // dᴵ(β) is shared by both species.
    for (int beta_I = 0; beta_I < neutron.Nbeta_I; ++beta_I) {
        const auto& RyI_C2D_K_K = representation_spin.calc_Ry(neutron.beta_F1D_beta(beta_I));
        for (int alpha_I = 0; alpha_I < neutron.Nalpha_I; ++alpha_I) {
            for (int gamma_I = 0; gamma_I < neutron.Ngamma_I; ++gamma_I) {
                neutron.prepare_rotation(alpha_I, beta_I, gamma_I);
                proton.prepare_rotation(alpha_I, beta_I, gamma_I);

                // Aν and Cπ include their gauge normalizations.
                const auto& annihilator_C3D_cfg1n_cfg2n_spn = neutron.calc_annihilator_pnp();
                const auto& creator_C3D_cfg1p_cfg2p_spp = proton.calc_creator_pnp();
                const Eigen::Map<const Eigen::MatrixXcd> annihilator_C2D_cfg1ncfg2n_spn(annihilator_C3D_cfg1n_cfg2n_spn.data(), Ncfg1ncfg2n_I, neutron.Nsp_I);
                const Eigen::Map<const Eigen::MatrixXcd> creator_C2D_cfg1pcfg2p_spp(creator_C3D_cfg1p_cfg2p_spp.data(), Ncfg1pcfg2p_I, proton.Nsp_I);

                // tmp = Aν Oᵀ; no complex conjugation.
                tmp_C2D_cfg1ncfg2n_spp.noalias() = annihilator_C2D_cfg1ncfg2n_spn * OneBody_C2D_spp_spn.transpose();

                // cπ†cν: (-1)^r₁; r₁ counts left neutron quasiparticles.
                for (int cfg2n_I = 0; cfg2n_I < neutron.Ncfg2_I; ++cfg2n_I) {
                    for (int cfg1n_I = 0; cfg1n_I < neutron.Ncfg1_I; ++cfg1n_I) {
                        const double sign_F = 1.0 - 2.0 * (neutron.hfb_pfaffian.config1_I2D_cfg1_cqp1[cfg1n_I].size() % 2);
                        const Eigen::Index row_I = cfg1n_I + static_cast<Eigen::Index>(neutron.Ncfg1_I) * cfg2n_I;
                        tmp_C2D_cfg1ncfg2n_spp.row(row_I) *= sign_F;
                    }
                }

                // tmp = Cπ [(-1)^r₁ Aν Oᵀ]ᵀ.
                tmp_C2D_cfg1pcfg2p_cfg1ncfg2n.noalias() = creator_C2D_cfg1pcfg2p_spp * tmp_C2D_cfg1ncfg2n_spp.transpose();

                // w = (2I+1) wα wβ wγ / VΩ.
                const double weight_F = normalization_F * neutron.weight_F1D_alpha(alpha_I) * neutron.weight_F1D_beta(beta_I) * neutron.weight_F1D_gamma(gamma_I);
                for (int K2_I = 0; K2_I < TargetTwoI_I + 1; ++K2_I) {
                    const double K2_F = K2_I - 0.5 * TargetTwoI_I;
                    for (int K1_I = 0; K1_I < TargetTwoI_I + 1; ++K1_I) {
                        const double K1_F = K1_I - 0.5 * TargetTwoI_I;

                        // Dᴵ* = exp(iK₁α) dᴵ* exp(iK₂γ).
                        const doubleC factorOmega_C = std::conj(RyI_C2D_K_K(K1_I, K2_I)) * std::exp(doubleC(0.0, K1_F * neutron.alpha_F1D_alpha(alpha_I) + K2_F * neutron.gamma_F1D_gamma(gamma_I)));
                        const doubleC weight_factor_C = weight_F * factorOmega_C;

                        // Each K₁,K₂ block contains all four configuration dimensions.
                        const Eigen::Index offset_I = (static_cast<Eigen::Index>(K2_I) * (TargetTwoI_I + 1) + K1_I) * Ncfg1pcfg2pcfg1ncfg2n_I;
                        Eigen::Map<Eigen::MatrixXcd> result_C2D_cfg1pcfg2p_cfg1ncfg2n(result_C6D_cfg1p_cfg2p_cfg1n_cfg2n_K1_K2.data() + offset_I, Ncfg1pcfg2p_I, Ncfg1ncfg2n_I);
                        result_C2D_cfg1pcfg2p_cfg1ncfg2n += weight_factor_C * tmp_C2D_cfg1pcfg2p_cfg1ncfg2n;
                    }
                }
            }
        }
    }
    return result_C6D_cfg1p_cfg2p_cfg1n_cfg2n_K1_K2;
}

inline const Eigen::Tensor<doubleC, 6, Eigen::ColMajor>& HFBProjectionNucleus::calc_two_body_NNNN(const Eigen::Tensor<double, 4, Eigen::ColMajor>& TwoBody_F4D_spn_spn_spn_spn) {
    auto& neutron = projection_neutron;
    auto& proton = projection_proton;

    // V ∈ ℝ^{Nspn×Nspn×Nspn×Nspn}.
    assert(TwoBody_F4D_spn_spn_spn_spn.dimension(0) == neutron.Nsp_I && TwoBody_F4D_spn_spn_spn_spn.dimension(1) == neutron.Nsp_I && TwoBody_F4D_spn_spn_spn_spn.dimension(2) == neutron.Nsp_I && TwoBody_F4D_spn_spn_spn_spn.dimension(3) == neutron.Nsp_I);
    assert(Eigen::Map<const Eigen::VectorXd>(TwoBody_F4D_spn_spn_spn_spn.data(), TwoBody_F4D_spn_spn_spn_spn.size()).allFinite());

    // Ων = Ωπ; wν = wπ.
    assert(neutron.alpha_F1D_alpha.isApprox(proton.alpha_F1D_alpha));
    assert(neutron.beta_F1D_beta.isApprox(proton.beta_F1D_beta));
    assert(neutron.gamma_F1D_gamma.isApprox(proton.gamma_F1D_gamma));
    assert(neutron.weight_F1D_alpha.isApprox(proton.weight_F1D_alpha));
    assert(neutron.weight_F1D_beta.isApprox(proton.weight_F1D_beta));
    assert(neutron.weight_F1D_gamma.isApprox(proton.weight_F1D_gamma));

    // VΩ = Σwα Σwβ Σwγ; wβ includes sinβ.
    const double volume_F = neutron.weight_F1D_alpha.sum() * neutron.weight_F1D_beta.sum() * neutron.weight_F1D_gamma.sum();
    assert(std::isfinite(volume_F) && volume_F > 0.0);
    const double normalization_F = (TargetTwoI_I + 1.0) / volume_F;

    result_C6D_cfg1p_cfg2p_cfg1n_cfg2n_K1_K2.setZero();

    // dᴵ(β) is shared by both species.
    for (int beta_I = 0; beta_I < neutron.Nbeta_I; ++beta_I) {
        const auto& RyI_C2D_K_K = representation_spin.calc_Ry(neutron.beta_F1D_beta(beta_I));
        for (int alpha_I = 0; alpha_I < neutron.Nalpha_I; ++alpha_I) {
            for (int gamma_I = 0; gamma_I < neutron.Ngamma_I; ++gamma_I) {
                neutron.prepare_rotation(alpha_I, beta_I, gamma_I);
                proton.prepare_rotation(alpha_I, beta_I, gamma_I);

                // Hν²(Ω) Nπ(Ω); gauge normalization and ½ included.
                const auto& overlap_C2D_cfg1p_cfg2p = proton.calc_overlap_pnp();
                const auto& twobody_C2D_cfg1n_cfg2n = neutron.calc_two_body_pnp(TwoBody_F4D_spn_spn_spn_spn);

                // ColMajor: (p₁,p₂) → rows; (n₁,n₂) → columns.
                const Eigen::Map<const Eigen::VectorXcd> twobody_C1D_cfg1ncfg2n(twobody_C2D_cfg1n_cfg2n.data(), Ncfg1ncfg2n_I);
                const Eigen::Map<const Eigen::VectorXcd> overlap_C1D_cfg1pcfg2p(overlap_C2D_cfg1p_cfg2p.data(), Ncfg1pcfg2p_I);
                Eigen::Map<Eigen::MatrixXcd> twobody_C2D_cfg1pcfg2p_cfg1ncfg2n(tmp_C4D_cfg1p_cfg2p_cfg1n_cfg2n.data(), Ncfg1pcfg2p_I, Ncfg1ncfg2n_I);
                twobody_C2D_cfg1pcfg2p_cfg1ncfg2n.noalias() = overlap_C1D_cfg1pcfg2p * twobody_C1D_cfg1ncfg2n.transpose();

                // w = (2I+1) wα wβ wγ / VΩ.
                const double weight_F = normalization_F * neutron.weight_F1D_alpha(alpha_I) * neutron.weight_F1D_beta(beta_I) * neutron.weight_F1D_gamma(gamma_I);
                for (int K2_I = 0; K2_I < TargetTwoI_I + 1; ++K2_I) {
                    const double K2_F = K2_I - 0.5 * TargetTwoI_I;
                    for (int K1_I = 0; K1_I < TargetTwoI_I + 1; ++K1_I) {
                        const double K1_F = K1_I - 0.5 * TargetTwoI_I;

                        // Dᴵ* = exp(iK₁α) dᴵ* exp(iK₂γ).
                        const doubleC factorOmega_C = std::conj(RyI_C2D_K_K(K1_I, K2_I)) * std::exp(doubleC(0.0, K1_F * neutron.alpha_F1D_alpha(alpha_I) + K2_F * neutron.gamma_F1D_gamma(gamma_I)));
                        const doubleC weight_factor_C = weight_F * factorOmega_C;

                        // Each K₁,K₂ block contains all four configuration dimensions.
                        const Eigen::Index offset_I = (static_cast<Eigen::Index>(K2_I) * (TargetTwoI_I + 1) + K1_I) * Ncfg1pcfg2pcfg1ncfg2n_I;
                        Eigen::Map<Eigen::MatrixXcd> result_C2D_cfg1pcfg2p_cfg1ncfg2n(result_C6D_cfg1p_cfg2p_cfg1n_cfg2n_K1_K2.data() + offset_I, Ncfg1pcfg2p_I, Ncfg1ncfg2n_I);
                        result_C2D_cfg1pcfg2p_cfg1ncfg2n += weight_factor_C * twobody_C2D_cfg1pcfg2p_cfg1ncfg2n;
                    }
                }
            }
        }
    }
    return result_C6D_cfg1p_cfg2p_cfg1n_cfg2n_K1_K2;
}

inline const Eigen::Tensor<doubleC, 6, Eigen::ColMajor>& HFBProjectionNucleus::calc_two_body_PPPP(const Eigen::Tensor<double, 4, Eigen::ColMajor>& TwoBody_F4D_spp_spp_spp_spp) {
    auto& neutron = projection_neutron;
    auto& proton = projection_proton;

    // V ∈ ℝ^{Nspp×Nspp×Nspp×Nspp}.
    assert(TwoBody_F4D_spp_spp_spp_spp.dimension(0) == proton.Nsp_I && TwoBody_F4D_spp_spp_spp_spp.dimension(1) == proton.Nsp_I && TwoBody_F4D_spp_spp_spp_spp.dimension(2) == proton.Nsp_I && TwoBody_F4D_spp_spp_spp_spp.dimension(3) == proton.Nsp_I);
    assert(Eigen::Map<const Eigen::VectorXd>(TwoBody_F4D_spp_spp_spp_spp.data(), TwoBody_F4D_spp_spp_spp_spp.size()).allFinite());

    // Ων = Ωπ; wν = wπ.
    assert(neutron.alpha_F1D_alpha.isApprox(proton.alpha_F1D_alpha));
    assert(neutron.beta_F1D_beta.isApprox(proton.beta_F1D_beta));
    assert(neutron.gamma_F1D_gamma.isApprox(proton.gamma_F1D_gamma));
    assert(neutron.weight_F1D_alpha.isApprox(proton.weight_F1D_alpha));
    assert(neutron.weight_F1D_beta.isApprox(proton.weight_F1D_beta));
    assert(neutron.weight_F1D_gamma.isApprox(proton.weight_F1D_gamma));

    // VΩ = Σwα Σwβ Σwγ; wβ includes sinβ.
    const double volume_F = neutron.weight_F1D_alpha.sum() * neutron.weight_F1D_beta.sum() * neutron.weight_F1D_gamma.sum();
    assert(std::isfinite(volume_F) && volume_F > 0.0);
    const double normalization_F = (TargetTwoI_I + 1.0) / volume_F;

    result_C6D_cfg1p_cfg2p_cfg1n_cfg2n_K1_K2.setZero();

    // dᴵ(β) is shared by both species.
    for (int beta_I = 0; beta_I < neutron.Nbeta_I; ++beta_I) {
        const auto& RyI_C2D_K_K = representation_spin.calc_Ry(neutron.beta_F1D_beta(beta_I));
        for (int alpha_I = 0; alpha_I < neutron.Nalpha_I; ++alpha_I) {
            for (int gamma_I = 0; gamma_I < neutron.Ngamma_I; ++gamma_I) {
                neutron.prepare_rotation(alpha_I, beta_I, gamma_I);
                proton.prepare_rotation(alpha_I, beta_I, gamma_I);

                // Nν(Ω) Hπ²(Ω); gauge normalization and ½ included.
                const auto& overlap_C2D_cfg1n_cfg2n = neutron.calc_overlap_pnp();
                const auto& twobody_C2D_cfg1p_cfg2p = proton.calc_two_body_pnp(TwoBody_F4D_spp_spp_spp_spp);

                // ColMajor: (p₁,p₂) → rows; (n₁,n₂) → columns.
                const Eigen::Map<const Eigen::VectorXcd> overlap_C1D_cfg1ncfg2n(overlap_C2D_cfg1n_cfg2n.data(), Ncfg1ncfg2n_I);
                const Eigen::Map<const Eigen::VectorXcd> twobody_C1D_cfg1pcfg2p(twobody_C2D_cfg1p_cfg2p.data(), Ncfg1pcfg2p_I);
                Eigen::Map<Eigen::MatrixXcd> twobody_C2D_cfg1pcfg2p_cfg1ncfg2n(tmp_C4D_cfg1p_cfg2p_cfg1n_cfg2n.data(), Ncfg1pcfg2p_I, Ncfg1ncfg2n_I);
                twobody_C2D_cfg1pcfg2p_cfg1ncfg2n.noalias() = twobody_C1D_cfg1pcfg2p * overlap_C1D_cfg1ncfg2n.transpose();

                // w = (2I+1) wα wβ wγ / VΩ.
                const double weight_F = normalization_F * neutron.weight_F1D_alpha(alpha_I) * neutron.weight_F1D_beta(beta_I) * neutron.weight_F1D_gamma(gamma_I);
                for (int K2_I = 0; K2_I < TargetTwoI_I + 1; ++K2_I) {
                    const double K2_F = K2_I - 0.5 * TargetTwoI_I;
                    for (int K1_I = 0; K1_I < TargetTwoI_I + 1; ++K1_I) {
                        const double K1_F = K1_I - 0.5 * TargetTwoI_I;

                        // Dᴵ* = exp(iK₁α) dᴵ* exp(iK₂γ).
                        const doubleC factorOmega_C = std::conj(RyI_C2D_K_K(K1_I, K2_I)) * std::exp(doubleC(0.0, K1_F * neutron.alpha_F1D_alpha(alpha_I) + K2_F * neutron.gamma_F1D_gamma(gamma_I)));
                        const doubleC weight_factor_C = weight_F * factorOmega_C;

                        // Each K₁,K₂ block contains all four configuration dimensions.
                        const Eigen::Index offset_I = (static_cast<Eigen::Index>(K2_I) * (TargetTwoI_I + 1) + K1_I) * Ncfg1pcfg2pcfg1ncfg2n_I;
                        Eigen::Map<Eigen::MatrixXcd> result_C2D_cfg1pcfg2p_cfg1ncfg2n(result_C6D_cfg1p_cfg2p_cfg1n_cfg2n_K1_K2.data() + offset_I, Ncfg1pcfg2p_I, Ncfg1ncfg2n_I);
                        result_C2D_cfg1pcfg2p_cfg1ncfg2n += weight_factor_C * twobody_C2D_cfg1pcfg2p_cfg1ncfg2n;
                    }
                }
            }
        }
    }
    return result_C6D_cfg1p_cfg2p_cfg1n_cfg2n_K1_K2;
}

inline const Eigen::Tensor<doubleC, 6, Eigen::ColMajor>& HFBProjectionNucleus::calc_two_body_NPPN(const Eigen::Tensor<double, 4, Eigen::ColMajor>& TwoBody_F4D_spp2_spp3_spn1_spn4) {
    auto& neutron = projection_neutron;
    auto& proton = projection_proton;

    // W(2,3,1,4) ∈ ℝ^{Nspp×Nspp×Nspn×Nspn}.
    assert(TwoBody_F4D_spp2_spp3_spn1_spn4.dimension(0) == proton.Nsp_I && TwoBody_F4D_spp2_spp3_spn1_spn4.dimension(1) == proton.Nsp_I && TwoBody_F4D_spp2_spp3_spn1_spn4.dimension(2) == neutron.Nsp_I && TwoBody_F4D_spp2_spp3_spn1_spn4.dimension(3) == neutron.Nsp_I);
    assert(Eigen::Map<const Eigen::VectorXd>(TwoBody_F4D_spp2_spp3_spn1_spn4.data(), TwoBody_F4D_spp2_spp3_spn1_spn4.size()).allFinite());

    // Ων = Ωπ; wν = wπ.
    assert(neutron.alpha_F1D_alpha.isApprox(proton.alpha_F1D_alpha));
    assert(neutron.beta_F1D_beta.isApprox(proton.beta_F1D_beta));
    assert(neutron.gamma_F1D_gamma.isApprox(proton.gamma_F1D_gamma));
    assert(neutron.weight_F1D_alpha.isApprox(proton.weight_F1D_alpha));
    assert(neutron.weight_F1D_beta.isApprox(proton.weight_F1D_beta));
    assert(neutron.weight_F1D_gamma.isApprox(proton.weight_F1D_gamma));

    // VΩ = Σwα Σwβ Σwγ; wβ includes sinβ.
    const double volume_F = neutron.weight_F1D_alpha.sum() * neutron.weight_F1D_beta.sum() * neutron.weight_F1D_gamma.sum();
    assert(std::isfinite(volume_F) && volume_F > 0.0);
    const double normalization_F = (TargetTwoI_I + 1.0) / volume_F;

    result_C6D_cfg1p_cfg2p_cfg1n_cfg2n_K1_K2.setZero();

    // ColMajor: (p₁,p₂) → rows; (n₁,n₂) → columns.
    Eigen::Map<Eigen::MatrixXcd> tmp_C2D_cfg1pcfg2p_cfg1ncfg2n(tmp_C4D_cfg1p_cfg2p_cfg1n_cfg2n.data(), Ncfg1pcfg2p_I, Ncfg1ncfg2n_I);

    // Qπ(cfgπ,23), Qν(cfgν,14), W(23,14).
    const Eigen::Index Nspnspn_I = static_cast<Eigen::Index>(neutron.Nsp_I) * neutron.Nsp_I;
    const Eigen::Index Nsppspp_I = static_cast<Eigen::Index>(proton.Nsp_I) * proton.Nsp_I;
    // ColMajor input (2,3,1,4) maps directly to W(23,14).
    const Eigen::Map<const Eigen::MatrixXd> TwoBody_F2D_sppspp_spnspn(TwoBody_F4D_spp2_spp3_spn1_spn4.data(), Nsppspp_I, Nspnspn_I);
    Eigen::MatrixXcd tmp_C2D_cfg1pcfg2p_spn1spn4{};
    tmp_C2D_cfg1pcfg2p_spn1spn4.resize(Ncfg1pcfg2p_I, Nspnspn_I);

    // dᴵ(β) is shared by both species.
    for (int beta_I = 0; beta_I < neutron.Nbeta_I; ++beta_I) {
        const auto& RyI_C2D_K_K = representation_spin.calc_Ry(neutron.beta_F1D_beta(beta_I));
        for (int alpha_I = 0; alpha_I < neutron.Nalpha_I; ++alpha_I) {
            for (int gamma_I = 0; gamma_I < neutron.Ngamma_I; ++gamma_I) {
                neutron.prepare_rotation(alpha_I, beta_I, gamma_I);
                proton.prepare_rotation(alpha_I, beta_I, gamma_I);

                // Qν(14), Qπ(23) include their gauge normalizations.
                const auto& OBTD_C4D_cfg1n_cfg2n_spn1_spn4 = neutron.calc_obtd_pnp();
                const auto& OBTD_C4D_cfg1p_cfg2p_spp2_spp3 = proton.calc_obtd_pnp();
                const Eigen::Map<const Eigen::MatrixXcd> OBTD_C2D_cfg1ncfg2n_spn1spn4(OBTD_C4D_cfg1n_cfg2n_spn1_spn4.data(), Ncfg1ncfg2n_I, Nspnspn_I);
                const Eigen::Map<const Eigen::MatrixXcd> OBTD_C2D_cfg1pcfg2p_spp2spp3(OBTD_C4D_cfg1p_cfg2p_spp2_spp3.data(), Ncfg1pcfg2p_I, Nsppspp_I);

                // cν₁†cπ₂†cπ₃cν₄ = (cν₁†cν₄)(cπ₂†cπ₃).
                // tmp = (Qπ W) Qνᵀ; merged W needs no ½.
                tmp_C2D_cfg1pcfg2p_spn1spn4.noalias() = OBTD_C2D_cfg1pcfg2p_spp2spp3 * TwoBody_F2D_sppspp_spnspn.cast<doubleC>();
                tmp_C2D_cfg1pcfg2p_cfg1ncfg2n.noalias() = tmp_C2D_cfg1pcfg2p_spn1spn4 * OBTD_C2D_cfg1ncfg2n_spn1spn4.transpose();

                // w = (2I+1) wα wβ wγ / VΩ.
                const double weight_F = normalization_F * neutron.weight_F1D_alpha(alpha_I) * neutron.weight_F1D_beta(beta_I) * neutron.weight_F1D_gamma(gamma_I);
                for (int K2_I = 0; K2_I < TargetTwoI_I + 1; ++K2_I) {
                    const double K2_F = K2_I - 0.5 * TargetTwoI_I;
                    for (int K1_I = 0; K1_I < TargetTwoI_I + 1; ++K1_I) {
                        const double K1_F = K1_I - 0.5 * TargetTwoI_I;

                        // Dᴵ* = exp(iK₁α) dᴵ* exp(iK₂γ).
                        const doubleC factorOmega_C = std::conj(RyI_C2D_K_K(K1_I, K2_I)) * std::exp(doubleC(0.0, K1_F * neutron.alpha_F1D_alpha(alpha_I) + K2_F * neutron.gamma_F1D_gamma(gamma_I)));
                        const doubleC weight_factor_C = weight_F * factorOmega_C;

                        // Each K₁,K₂ block contains all four configuration dimensions.
                        const Eigen::Index offset_I = (static_cast<Eigen::Index>(K2_I) * (TargetTwoI_I + 1) + K1_I) * Ncfg1pcfg2pcfg1ncfg2n_I;
                        Eigen::Map<Eigen::MatrixXcd> result_C2D_cfg1pcfg2p_cfg1ncfg2n(result_C6D_cfg1p_cfg2p_cfg1n_cfg2n_K1_K2.data() + offset_I, Ncfg1pcfg2p_I, Ncfg1ncfg2n_I);
                        result_C2D_cfg1pcfg2p_cfg1ncfg2n += weight_factor_C * tmp_C2D_cfg1pcfg2p_cfg1ncfg2n;
                    }
                }
            }
        }
    }
    return result_C6D_cfg1p_cfg2p_cfg1n_cfg2n_K1_K2;
}
