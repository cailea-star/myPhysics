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
     * @brief  Integrate nn kernels using Pfaffians and quadrature.
     * @math   O_nn = Σαβ Oαβ cνα† cνβ.
     * @output Tensor indexed by cfg1p, cfg2p, cfg1n, cfg2n, K1, K2.
     * @note   Requires rotationally scalar matrix elements.
     */
    const Eigen::Tensor<doubleC, 6, Eigen::ColMajor>& calc_one_body_nn(const Eigen::MatrixXcd& OneBody_C2D_spn_spn);

    /**
     * @brief  Integrate pp kernels using Pfaffians and quadrature.
     * @math   O_pp = Σαβ Oαβ cπα† cπβ.
     * @output Tensor indexed by cfg1p, cfg2p, cfg1n, cfg2n, K1, K2.
     * @note   Requires rotationally scalar matrix elements.
     */
    const Eigen::Tensor<doubleC, 6, Eigen::ColMajor>& calc_one_body_pp(const Eigen::MatrixXcd& OneBody_C2D_spp_spp);

    /**
     * @brief  Integrate pn kernels using Pfaffians and quadrature.
     * @math   O_pn = Σαβ Oαβ cπα† cνβ.
     * @output Tensor indexed by cfg1p, cfg2p, cfg1n, cfg2n, K1, K2.
     * @note   Requires rotationally scalar matrix elements.
     */
    const Eigen::Tensor<doubleC, 6, Eigen::ColMajor>& calc_one_body_pn(const Eigen::MatrixXcd& OneBody_C2D_spp_spn);

    /**
     * @brief  Integrate separable nnnn kernels using Pfaffians and quadrature.
     * @math   O = Σμ[Q̂λμ†Q̂λμ − Σαδ(Qλμ†Qλμ)αδ cνα†cνδ].
     * @math   Q ∈ ℝ^{Nspν × Nspν × (2λ+1)}; twoLambda_I = 2λ ≥ 0, λ ∈ ℤ; μ = −λ,…,λ.
     * @output Tensor indexed by cfg1p, cfg2p, cfg1n, cfg2n, K1, K2.
     * @note   Complete spherical tensor; excludes coupling strength and extra ½.
     */
    const Eigen::Tensor<doubleC, 6, Eigen::ColMajor>& calc_two_body_Q_nnnn(int twoLambda_I, const Eigen::Tensor<double, 3, Eigen::ColMajor>& Q_F3D_spn_spn_mu);

    /**
     * @brief  Integrate separable pppp kernels using Pfaffians and quadrature.
     * @math   O = Σμ[Q̂λμ†Q̂λμ − Σαδ(Qλμ†Qλμ)αδ cπα†cπδ].
     * @math   Q ∈ ℝ^{Nspπ × Nspπ × (2λ+1)}; twoLambda_I = 2λ ≥ 0, λ ∈ ℤ; μ = −λ,…,λ.
     * @output Tensor indexed by cfg1p, cfg2p, cfg1n, cfg2n, K1, K2.
     * @note   Complete spherical tensor; excludes coupling strength and extra ½.
     */
    const Eigen::Tensor<doubleC, 6, Eigen::ColMajor>& calc_two_body_Q_pppp(int twoLambda_I, const Eigen::Tensor<double, 3, Eigen::ColMajor>& Q_F3D_spp_spp_mu);

    /**
     * @brief  Integrate mixed QQ kernels using Pfaffians and quadrature.
     * @math   O = Σμ Q̂ν,λμ† Q̂π,λμ; cν†cν cπ†cπ = cν†cπ†cπcν.
     * @math   Qτ ∈ ℝ^{Nspτ × Nspτ × (2λ+1)}; twoLambda_I = 2λ ≥ 0, λ ∈ ℤ; μ = −λ,…,λ.
     * @output Tensor indexed by cfg1p, cfg2p, cfg1n, cfg2n, K1, K2.
     * @note   Consistent complete spherical tensors; coupling strength is external.
     */
    const Eigen::Tensor<doubleC, 6, Eigen::ColMajor>& calc_two_body_Q_nppn(int twoLambda_I, const Eigen::Tensor<double, 3, Eigen::ColMajor>& Q_F3D_spn_spn_mu, const Eigen::Tensor<double, 3, Eigen::ColMajor>& Q_F3D_spp_spp_mu);

    /**
     * @brief  Integrate neutron pairing kernels using Pfaffians and quadrature.
     * @math   P̂ν,λμ† = ½Σab Pab cνa†cνb†; O = Σμ P̂ν,λμ†P̂ν,λμ.
     * @math   P ∈ ℝ^{Nspν × Nspν × (2λ+1)}; Pμᵀ = −Pμ; twoLambda_I = 2λ ≥ 0, λ ∈ ℤ; μ = −λ,…,λ.
     * @output Tensor indexed by cfg1p, cfg2p, cfg1n, cfg2n, K1, K2.
     * @note   Complete spherical tensor; excludes coupling strength and attraction sign.
     */
    const Eigen::Tensor<doubleC, 6, Eigen::ColMajor>& calc_two_body_P_nnnn(int twoLambda_I, const Eigen::Tensor<double, 3, Eigen::ColMajor>& P_F3D_spn_spn_mu);

    /**
     * @brief  Integrate proton pairing kernels using Pfaffians and quadrature.
     * @math   P̂π,λμ† = ½Σab Pab cπa†cπb†; O = Σμ P̂π,λμ†P̂π,λμ.
     * @math   P ∈ ℝ^{Nspπ × Nspπ × (2λ+1)}; Pμᵀ = −Pμ; twoLambda_I = 2λ ≥ 0, λ ∈ ℤ; μ = −λ,…,λ.
     * @output Tensor indexed by cfg1p, cfg2p, cfg1n, cfg2n, K1, K2.
     * @note   Complete spherical tensor; excludes coupling strength and attraction sign.
     */
    const Eigen::Tensor<doubleC, 6, Eigen::ColMajor>& calc_two_body_P_pppp(int twoLambda_I, const Eigen::Tensor<double, 3, Eigen::ColMajor>& P_F3D_spp_spp_mu);

private:
    /**
     * @brief  Integrate nuclear kernels using Euler quadrature.
     * @math   result = (2I+1)/VΩ ΣΩ wΩ Dᴵ*(Ω) K(Ω).
     * @output Updated result_C6D_cfg1p_cfg2p_cfg1n_cfg2n_K1_K2.
     * @note   Callback fills tmp_C4D at current neutron/proton rotations.
     * @note   Callback includes both species' gauge normalizations and overlaps.
     * @note   Callback preserves angular data and result storage.
     */
    template <typename KernelFunc>
    const Eigen::Tensor<doubleC, 6, Eigen::ColMajor>& integrate_amp(const KernelFunc& kernel_Func);
};

template <typename KernelFunc>
inline const Eigen::Tensor<doubleC, 6, Eigen::ColMajor>& HFBProjectionNucleus::integrate_amp(const KernelFunc& kernel_Func) {
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

    Eigen::Map<const Eigen::MatrixXcd> tmp_C2D_cfg1pcfg2p_cfg1ncfg2n(tmp_C4D_cfg1p_cfg2p_cfg1n_cfg2n.data(), Ncfg1pcfg2p_I, Ncfg1ncfg2n_I);

    // dᴵ(β) is shared by both species.
    for (int beta_I = 0; beta_I < neutron.Nbeta_I; ++beta_I) {
        const auto& RyI_C2D_K_K = representation_spin.calc_Ry(neutron.beta_F1D_beta(beta_I));
        for (int alpha_I = 0; alpha_I < neutron.Nalpha_I; ++alpha_I) {
            for (int gamma_I = 0; gamma_I < neutron.Ngamma_I; ++gamma_I) {
                {
                    // D = Rz(α) Ry(β) Rz(γ).
                    const Eigen::Map<const Eigen::MatrixXcd> RzAlpha_C2D_sp_sp(neutron.Rz_C3D_sp_sp_alpha.data() + static_cast<Eigen::Index>(alpha_I) * neutron.Nsp_I * neutron.Nsp_I, neutron.Nsp_I, neutron.Nsp_I);
                    const Eigen::Map<const Eigen::MatrixXcd> Ry_C2D_sp_sp(neutron.Ry_C3D_sp_sp_beta.data() + static_cast<Eigen::Index>(beta_I) * neutron.Nsp_I * neutron.Nsp_I, neutron.Nsp_I, neutron.Nsp_I);
                    const Eigen::Map<const Eigen::MatrixXcd> RzGamma_C2D_sp_sp(neutron.Rz_C3D_sp_sp_gamma.data() + static_cast<Eigen::Index>(gamma_I) * neutron.Nsp_I * neutron.Nsp_I, neutron.Nsp_I, neutron.Nsp_I);
                    neutron.D_C2D_sp_sp.noalias() = RzAlpha_C2D_sp_sp * Ry_C2D_sp_sp * RzGamma_C2D_sp_sp;
                }
                {
                    // D = Rz(α) Ry(β) Rz(γ).
                    const Eigen::Map<const Eigen::MatrixXcd> RzAlpha_C2D_sp_sp(proton.Rz_C3D_sp_sp_alpha.data() + static_cast<Eigen::Index>(alpha_I) * proton.Nsp_I * proton.Nsp_I, proton.Nsp_I, proton.Nsp_I);
                    const Eigen::Map<const Eigen::MatrixXcd> Ry_C2D_sp_sp(proton.Ry_C3D_sp_sp_beta.data() + static_cast<Eigen::Index>(beta_I) * proton.Nsp_I * proton.Nsp_I, proton.Nsp_I, proton.Nsp_I);
                    const Eigen::Map<const Eigen::MatrixXcd> RzGamma_C2D_sp_sp(proton.Rz_C3D_sp_sp_gamma.data() + static_cast<Eigen::Index>(gamma_I) * proton.Nsp_I * proton.Nsp_I, proton.Nsp_I, proton.Nsp_I);
                    proton.D_C2D_sp_sp.noalias() = RzAlpha_C2D_sp_sp * Ry_C2D_sp_sp * RzGamma_C2D_sp_sp;
                }

                kernel_Func();

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

inline const Eigen::Tensor<doubleC, 6, Eigen::ColMajor>& HFBProjectionNucleus::calc_overlap() {
    auto& neutron = projection_neutron;
    auto& proton = projection_proton;

    return integrate_amp([&]() {
        const auto& overlap_C2D_cfg1n_cfg2n = neutron.calc_overlap_pnp();
        const auto& overlap_C2D_cfg1p_cfg2p = proton.calc_overlap_pnp();

        // ColMajor: (p₁,p₂) → rows; (n₁,n₂) → columns.
        const Eigen::Map<const Eigen::VectorXcd> overlap_C1D_cfg1ncfg2n(overlap_C2D_cfg1n_cfg2n.data(), Ncfg1ncfg2n_I);
        const Eigen::Map<const Eigen::VectorXcd> overlap_C1D_cfg1pcfg2p(overlap_C2D_cfg1p_cfg2p.data(), Ncfg1pcfg2p_I);
        Eigen::Map<Eigen::MatrixXcd> overlap_C2D_cfg1pcfg2p_cfg1ncfg2n(tmp_C4D_cfg1p_cfg2p_cfg1n_cfg2n.data(), Ncfg1pcfg2p_I, Ncfg1ncfg2n_I);
        overlap_C2D_cfg1pcfg2p_cfg1ncfg2n.noalias() = overlap_C1D_cfg1pcfg2p * overlap_C1D_cfg1ncfg2n.transpose();
    });
}

inline const Eigen::Tensor<doubleC, 6, Eigen::ColMajor>& HFBProjectionNucleus::calc_one_body_nn(const Eigen::MatrixXcd& OneBody_C2D_spn_spn) {
    auto& neutron = projection_neutron;
    auto& proton = projection_proton;
    assert(OneBody_C2D_spn_spn.rows() == neutron.Nsp_I && OneBody_C2D_spn_spn.cols() == neutron.Nsp_I && OneBody_C2D_spn_spn.allFinite());

    return integrate_amp([&]() {
        const auto& overlap_C2D_cfg1p_cfg2p = proton.calc_overlap_pnp();
        const auto& onebody_C2D_cfg1n_cfg2n = neutron.calc_one_body_pnp(OneBody_C2D_spn_spn);

        // ColMajor: (p₁,p₂) → rows; (n₁,n₂) → columns.
        const Eigen::Map<const Eigen::VectorXcd> onebody_C1D_cfg1ncfg2n(onebody_C2D_cfg1n_cfg2n.data(), Ncfg1ncfg2n_I);
        const Eigen::Map<const Eigen::VectorXcd> overlap_C1D_cfg1pcfg2p(overlap_C2D_cfg1p_cfg2p.data(), Ncfg1pcfg2p_I);
        Eigen::Map<Eigen::MatrixXcd> onebody_C2D_cfg1pcfg2p_cfg1ncfg2n(tmp_C4D_cfg1p_cfg2p_cfg1n_cfg2n.data(), Ncfg1pcfg2p_I, Ncfg1ncfg2n_I);
        onebody_C2D_cfg1pcfg2p_cfg1ncfg2n.noalias() = overlap_C1D_cfg1pcfg2p * onebody_C1D_cfg1ncfg2n.transpose();
    });
}

inline const Eigen::Tensor<doubleC, 6, Eigen::ColMajor>& HFBProjectionNucleus::calc_one_body_pp(const Eigen::MatrixXcd& OneBody_C2D_spp_spp) {
    auto& neutron = projection_neutron;
    auto& proton = projection_proton;
    assert(OneBody_C2D_spp_spp.rows() == proton.Nsp_I && OneBody_C2D_spp_spp.cols() == proton.Nsp_I && OneBody_C2D_spp_spp.allFinite());

    return integrate_amp([&]() {
        const auto& overlap_C2D_cfg1n_cfg2n = neutron.calc_overlap_pnp();
        const auto& onebody_C2D_cfg1p_cfg2p = proton.calc_one_body_pnp(OneBody_C2D_spp_spp);

        // ColMajor: (p₁,p₂) → rows; (n₁,n₂) → columns.
        const Eigen::Map<const Eigen::VectorXcd> overlap_C1D_cfg1ncfg2n(overlap_C2D_cfg1n_cfg2n.data(), Ncfg1ncfg2n_I);
        const Eigen::Map<const Eigen::VectorXcd> onebody_C1D_cfg1pcfg2p(onebody_C2D_cfg1p_cfg2p.data(), Ncfg1pcfg2p_I);
        Eigen::Map<Eigen::MatrixXcd> onebody_C2D_cfg1pcfg2p_cfg1ncfg2n(tmp_C4D_cfg1p_cfg2p_cfg1n_cfg2n.data(), Ncfg1pcfg2p_I, Ncfg1ncfg2n_I);
        onebody_C2D_cfg1pcfg2p_cfg1ncfg2n.noalias() = onebody_C1D_cfg1pcfg2p * overlap_C1D_cfg1ncfg2n.transpose();
    });
}

inline const Eigen::Tensor<doubleC, 6, Eigen::ColMajor>& HFBProjectionNucleus::calc_one_body_pn(const Eigen::MatrixXcd& OneBody_C2D_spp_spn) {
    auto& neutron = projection_neutron;
    auto& proton = projection_proton;
    assert(OneBody_C2D_spp_spn.rows() == proton.Nsp_I && OneBody_C2D_spp_spn.cols() == neutron.Nsp_I && OneBody_C2D_spp_spn.allFinite());

    // ColMajor: (p₁,p₂) → rows; (n₁,n₂) → columns.
    Eigen::Map<Eigen::MatrixXcd> tmp_C2D_cfg1pcfg2p_cfg1ncfg2n(tmp_C4D_cfg1p_cfg2p_cfg1n_cfg2n.data(), Ncfg1pcfg2p_I, Ncfg1ncfg2n_I);

    // Aν Oᵀ ∈ ℂ^{Ncfg1ncfg2n×Nspp}.
    Eigen::MatrixXcd tmp_C2D_cfg1ncfg2n_spp{};
    tmp_C2D_cfg1ncfg2n_spp.resize(Ncfg1ncfg2n_I, proton.Nsp_I);

    return integrate_amp([&]() {
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
    });
}

inline const Eigen::Tensor<doubleC, 6, Eigen::ColMajor>& HFBProjectionNucleus::calc_two_body_Q_nnnn(int twoLambda_I, const Eigen::Tensor<double, 3, Eigen::ColMajor>& Q_F3D_spn_spn_mu) {
    auto& neutron = projection_neutron;
    auto& proton = projection_proton;

    // Qλμ ∈ ℝ^{Nsp×Nsp}; μ = −λ,…,λ.
    assert(twoLambda_I >= 0 && twoLambda_I % 2 == 0);
    assert(Q_F3D_spn_spn_mu.dimension(0) == neutron.Nsp_I && Q_F3D_spn_spn_mu.dimension(1) == neutron.Nsp_I && Q_F3D_spn_spn_mu.dimension(2) == twoLambda_I + 1);
    assert(Eigen::Map<const Eigen::VectorXd>(Q_F3D_spn_spn_mu.data(), Q_F3D_spn_spn_mu.size()).allFinite());

    return integrate_amp([&]() {
        // Hν²(Ω) Nπ(Ω); gauge normalization included; no additional ½.
        const auto& overlap_C2D_cfg1p_cfg2p = proton.calc_overlap_pnp();
        const auto& twobody_C2D_cfg1n_cfg2n = neutron.calc_two_body_Q_pnp(twoLambda_I, Q_F3D_spn_spn_mu);

        // ColMajor: (p₁,p₂) → rows; (n₁,n₂) → columns.
        const Eigen::Map<const Eigen::VectorXcd> twobody_C1D_cfg1ncfg2n(twobody_C2D_cfg1n_cfg2n.data(), Ncfg1ncfg2n_I);
        const Eigen::Map<const Eigen::VectorXcd> overlap_C1D_cfg1pcfg2p(overlap_C2D_cfg1p_cfg2p.data(), Ncfg1pcfg2p_I);
        Eigen::Map<Eigen::MatrixXcd> twobody_C2D_cfg1pcfg2p_cfg1ncfg2n(tmp_C4D_cfg1p_cfg2p_cfg1n_cfg2n.data(), Ncfg1pcfg2p_I, Ncfg1ncfg2n_I);
        twobody_C2D_cfg1pcfg2p_cfg1ncfg2n.noalias() = overlap_C1D_cfg1pcfg2p * twobody_C1D_cfg1ncfg2n.transpose();
    });
}

inline const Eigen::Tensor<doubleC, 6, Eigen::ColMajor>& HFBProjectionNucleus::calc_two_body_Q_pppp(int twoLambda_I, const Eigen::Tensor<double, 3, Eigen::ColMajor>& Q_F3D_spp_spp_mu) {
    auto& neutron = projection_neutron;
    auto& proton = projection_proton;

    // Qλμ ∈ ℝ^{Nsp×Nsp}; μ = −λ,…,λ.
    assert(twoLambda_I >= 0 && twoLambda_I % 2 == 0);
    assert(Q_F3D_spp_spp_mu.dimension(0) == proton.Nsp_I && Q_F3D_spp_spp_mu.dimension(1) == proton.Nsp_I && Q_F3D_spp_spp_mu.dimension(2) == twoLambda_I + 1);
    assert(Eigen::Map<const Eigen::VectorXd>(Q_F3D_spp_spp_mu.data(), Q_F3D_spp_spp_mu.size()).allFinite());

    return integrate_amp([&]() {
        // Nν(Ω) Hπ²(Ω); gauge normalization included; no additional ½.
        const auto& overlap_C2D_cfg1n_cfg2n = neutron.calc_overlap_pnp();
        const auto& twobody_C2D_cfg1p_cfg2p = proton.calc_two_body_Q_pnp(twoLambda_I, Q_F3D_spp_spp_mu);

        // ColMajor: (p₁,p₂) → rows; (n₁,n₂) → columns.
        const Eigen::Map<const Eigen::VectorXcd> overlap_C1D_cfg1ncfg2n(overlap_C2D_cfg1n_cfg2n.data(), Ncfg1ncfg2n_I);
        const Eigen::Map<const Eigen::VectorXcd> twobody_C1D_cfg1pcfg2p(twobody_C2D_cfg1p_cfg2p.data(), Ncfg1pcfg2p_I);
        Eigen::Map<Eigen::MatrixXcd> twobody_C2D_cfg1pcfg2p_cfg1ncfg2n(tmp_C4D_cfg1p_cfg2p_cfg1n_cfg2n.data(), Ncfg1pcfg2p_I, Ncfg1ncfg2n_I);
        twobody_C2D_cfg1pcfg2p_cfg1ncfg2n.noalias() = twobody_C1D_cfg1pcfg2p * overlap_C1D_cfg1ncfg2n.transpose();
    });
}

inline const Eigen::Tensor<doubleC, 6, Eigen::ColMajor>& HFBProjectionNucleus::calc_two_body_Q_nppn(int twoLambda_I, const Eigen::Tensor<double, 3, Eigen::ColMajor>& Q_F3D_spn_spn_mu, const Eigen::Tensor<double, 3, Eigen::ColMajor>& Q_F3D_spp_spp_mu) {
    auto& neutron = projection_neutron;
    auto& proton = projection_proton;

    // Qν,Qπ share λ and μ; single-particle dimensions may differ.
    assert(twoLambda_I >= 0 && twoLambda_I % 2 == 0);
    assert(Q_F3D_spn_spn_mu.dimension(0) == neutron.Nsp_I && Q_F3D_spn_spn_mu.dimension(1) == neutron.Nsp_I && Q_F3D_spn_spn_mu.dimension(2) == twoLambda_I + 1);
    assert(Q_F3D_spp_spp_mu.dimension(0) == proton.Nsp_I && Q_F3D_spp_spp_mu.dimension(1) == proton.Nsp_I && Q_F3D_spp_spp_mu.dimension(2) == twoLambda_I + 1);
    assert(Eigen::Map<const Eigen::VectorXd>(Q_F3D_spn_spn_mu.data(), Q_F3D_spn_spn_mu.size()).allFinite());
    assert(Eigen::Map<const Eigen::VectorXd>(Q_F3D_spp_spp_mu.data(), Q_F3D_spp_spp_mu.size()).allFinite());

    // Qν† = Qνᵀ for real input; workspaces reused over Ω,μ.
    Eigen::MatrixXcd Qdag_C2D_spn_spn{};
    Eigen::MatrixXcd Q_C2D_spp_spp{};
    Qdag_C2D_spn_spn.resize(neutron.Nsp_I, neutron.Nsp_I);
    Q_C2D_spp_spp.resize(proton.Nsp_I, proton.Nsp_I);
    Eigen::Map<Eigen::MatrixXcd> tmp_C2D_cfg1pcfg2p_cfg1ncfg2n(tmp_C4D_cfg1p_cfg2p_cfg1n_cfg2n.data(), Ncfg1pcfg2p_I, Ncfg1ncfg2n_I);

    return integrate_amp([&]() {
        tmp_C2D_cfg1pcfg2p_cfg1ncfg2n.setZero();
        for (int mu_I = 0; mu_I <= twoLambda_I; ++mu_I) {
            // Oμ = Qνμ† Qπμ at the same Euler node.
            const Eigen::Map<const Eigen::MatrixXd> Q_F2D_spn_spn(Q_F3D_spn_spn_mu.data() + static_cast<Eigen::Index>(mu_I) * neutron.Nsp_I * neutron.Nsp_I, neutron.Nsp_I, neutron.Nsp_I);
            const Eigen::Map<const Eigen::MatrixXd> Q_F2D_spp_spp(Q_F3D_spp_spp_mu.data() + static_cast<Eigen::Index>(mu_I) * proton.Nsp_I * proton.Nsp_I, proton.Nsp_I, proton.Nsp_I);
            Qdag_C2D_spn_spn = Q_F2D_spn_spn.transpose().cast<doubleC>();
            Q_C2D_spp_spp = Q_F2D_spp_spp.cast<doubleC>();
            const auto& onebody_C2D_cfg1n_cfg2n = neutron.calc_one_body_pnp(Qdag_C2D_spn_spn);
            const auto& onebody_C2D_cfg1p_cfg2p = proton.calc_one_body_pnp(Q_C2D_spp_spp);

            // Both kernels include overlaps and gauge normalizations.
            const Eigen::Map<const Eigen::VectorXcd> onebody_C1D_cfg1ncfg2n(onebody_C2D_cfg1n_cfg2n.data(), Ncfg1ncfg2n_I);
            const Eigen::Map<const Eigen::VectorXcd> onebody_C1D_cfg1pcfg2p(onebody_C2D_cfg1p_cfg2p.data(), Ncfg1pcfg2p_I);
            tmp_C2D_cfg1pcfg2p_cfg1ncfg2n.noalias() += onebody_C1D_cfg1pcfg2p * onebody_C1D_cfg1ncfg2n.transpose();
        }
    });
}

inline const Eigen::Tensor<doubleC, 6, Eigen::ColMajor>& HFBProjectionNucleus::calc_two_body_P_nnnn(int twoLambda_I, const Eigen::Tensor<double, 3, Eigen::ColMajor>& P_F3D_spn_spn_mu) {
    auto& neutron = projection_neutron;
    auto& proton = projection_proton;

    // Pλμ ∈ ℝ^{Nsp×Nsp}; μ = −λ,…,λ.
    assert(twoLambda_I >= 0 && twoLambda_I % 2 == 0);
    assert(P_F3D_spn_spn_mu.dimension(0) == neutron.Nsp_I && P_F3D_spn_spn_mu.dimension(1) == neutron.Nsp_I && P_F3D_spn_spn_mu.dimension(2) == twoLambda_I + 1);
    assert(Eigen::Map<const Eigen::VectorXd>(P_F3D_spn_spn_mu.data(), P_F3D_spn_spn_mu.size()).allFinite());

    return integrate_amp([&]() {
        // Hν²(Ω) Nπ(Ω); gauge normalization and pairing ¼ included.
        const auto& overlap_C2D_cfg1p_cfg2p = proton.calc_overlap_pnp();
        const auto& twobody_C2D_cfg1n_cfg2n = neutron.calc_two_body_P_pnp(twoLambda_I, P_F3D_spn_spn_mu);

        // ColMajor: (p₁,p₂) → rows; (n₁,n₂) → columns.
        const Eigen::Map<const Eigen::VectorXcd> twobody_C1D_cfg1ncfg2n(twobody_C2D_cfg1n_cfg2n.data(), Ncfg1ncfg2n_I);
        const Eigen::Map<const Eigen::VectorXcd> overlap_C1D_cfg1pcfg2p(overlap_C2D_cfg1p_cfg2p.data(), Ncfg1pcfg2p_I);
        Eigen::Map<Eigen::MatrixXcd> twobody_C2D_cfg1pcfg2p_cfg1ncfg2n(tmp_C4D_cfg1p_cfg2p_cfg1n_cfg2n.data(), Ncfg1pcfg2p_I, Ncfg1ncfg2n_I);
        twobody_C2D_cfg1pcfg2p_cfg1ncfg2n.noalias() = overlap_C1D_cfg1pcfg2p * twobody_C1D_cfg1ncfg2n.transpose();
    });
}

inline const Eigen::Tensor<doubleC, 6, Eigen::ColMajor>& HFBProjectionNucleus::calc_two_body_P_pppp(int twoLambda_I, const Eigen::Tensor<double, 3, Eigen::ColMajor>& P_F3D_spp_spp_mu) {
    auto& neutron = projection_neutron;
    auto& proton = projection_proton;

    // Pλμ ∈ ℝ^{Nsp×Nsp}; μ = −λ,…,λ.
    assert(twoLambda_I >= 0 && twoLambda_I % 2 == 0);
    assert(P_F3D_spp_spp_mu.dimension(0) == proton.Nsp_I && P_F3D_spp_spp_mu.dimension(1) == proton.Nsp_I && P_F3D_spp_spp_mu.dimension(2) == twoLambda_I + 1);
    assert(Eigen::Map<const Eigen::VectorXd>(P_F3D_spp_spp_mu.data(), P_F3D_spp_spp_mu.size()).allFinite());

    return integrate_amp([&]() {
        // Nν(Ω) Hπ²(Ω); gauge normalization and pairing ¼ included.
        const auto& overlap_C2D_cfg1n_cfg2n = neutron.calc_overlap_pnp();
        const auto& twobody_C2D_cfg1p_cfg2p = proton.calc_two_body_P_pnp(twoLambda_I, P_F3D_spp_spp_mu);

        // ColMajor: (p₁,p₂) → rows; (n₁,n₂) → columns.
        const Eigen::Map<const Eigen::VectorXcd> overlap_C1D_cfg1ncfg2n(overlap_C2D_cfg1n_cfg2n.data(), Ncfg1ncfg2n_I);
        const Eigen::Map<const Eigen::VectorXcd> twobody_C1D_cfg1pcfg2p(twobody_C2D_cfg1p_cfg2p.data(), Ncfg1pcfg2p_I);
        Eigen::Map<Eigen::MatrixXcd> twobody_C2D_cfg1pcfg2p_cfg1ncfg2n(tmp_C4D_cfg1p_cfg2p_cfg1n_cfg2n.data(), Ncfg1pcfg2p_I, Ncfg1ncfg2n_I);
        twobody_C2D_cfg1pcfg2p_cfg1ncfg2n.noalias() = twobody_C1D_cfg1pcfg2p * overlap_C1D_cfg1ncfg2n.transpose();
    });
}
