/**
 * @file    hfb_projection_nucleus.hpp
 * @author  cailea
 * @date    2026-09-17
 * @brief   Joint neutron-proton projection interfaces.
 */

#pragma once

#include <cassert>
#include <cmath>
#include <functional>
#include <utility>

#include "hfb_projection.hpp"

/**
 * @brief  Project neutron-proton configurations using gauge and Euler quadrature.
 * @math   |Φκ⟩ = |Φν,κν⟩ ⊗ |Φπ,κπ⟩.
 * @note   Species-local indices start at zero; Euler meshes are shared.
 * @note   Projected tensors are allocated locally and returned by value.
 * @note   Particle-number targets apply to the right state.
 * @note   Requires even reference vacua and initialized projection data.
 * @note   Vacuum overlaps must remain nonzero at quadrature nodes.
 */
class HFBProjectionNucleus {
public:
    HFBProjection projection_neutron;
    HFBProjection projection_proton;

    Eigen::Tensor<doubleC, 7, Eigen::ColMajor> result_C7D_cfg1p_cfg2p_cfg1n_cfg2n_alpha_beta_gamma{};

private:
    Eigen::Index Ncfg1pcfg2p_I = 0;
    Eigen::Index Ncfg1ncfg2n_I = 0;
    Eigen::Index Ncfg1pcfg2pcfg1ncfg2n_I = 0;

    public:
    /**
     * @brief  Construct empty neutron-proton projection workspaces.
     * @math   Cₙ = Cₚ = ∅.
     * @output Empty workspaces.
     * @note   Assign a configured object before calculations.
     */
    HFBProjectionNucleus() = default;

    /**
     * @brief  Move species projections into owned workspaces.
     * @math   |Φκ⟩ = |Φν,κν⟩ ⊗ |Φπ,κπ⟩.
     * @output Object owning neutron and proton projections.
     */
    HFBProjectionNucleus(HFBProjection projection_neutron_, HFBProjection projection_proton_)
    : projection_neutron(std::move(projection_neutron_)), projection_proton(std::move(projection_proton_)) {
        assert(projection_neutron.Nalpha_I == projection_proton.Nalpha_I);
        assert(projection_neutron.Nbeta_I == projection_proton.Nbeta_I);
        assert(projection_neutron.Ngamma_I == projection_proton.Ngamma_I);

        Ncfg1pcfg2p_I = static_cast<Eigen::Index>(projection_proton.Ncfg1_I) * projection_proton.Ncfg2_I;
        Ncfg1ncfg2n_I = static_cast<Eigen::Index>(projection_neutron.Ncfg1_I) * projection_neutron.Ncfg2_I;
        Ncfg1pcfg2pcfg1ncfg2n_I = Ncfg1pcfg2p_I * Ncfg1ncfg2n_I;

    }

    /**
     * @brief Build nuclear overlaps using gauge quadrature.
     * @math result_p₁p₂n₁n₂(Ω) = Nπ_p₁p₂(Ω) Nν_n₁n₂(Ω).
     * @output Reference to result_C7D_cfg1p_cfg2p_cfg1n_cfg2n_alpha_beta_gamma.
     * @note Subsequent builds overwrite this buffer.
     */
    const Eigen::Tensor<doubleC, 7, Eigen::ColMajor>& build_overlap();

    /**
     * @brief Build nn kernels using Pfaffians and gauge quadrature.
     * @math result(Ω) = Nπ(Ω) Hν¹(Ω).
     * @output Reference to result_C7D_cfg1p_cfg2p_cfg1n_cfg2n_alpha_beta_gamma.
     * @note Subsequent builds overwrite this buffer.
     */
    const Eigen::Tensor<doubleC, 7, Eigen::ColMajor>& build_one_body_nn(const Eigen::MatrixXcd& OneBody_C2D_spn_spn);

    /**
     * @brief Build pp kernels using Pfaffians and gauge quadrature.
     * @math result(Ω) = Hπ¹(Ω) Nν(Ω).
     * @output Reference to result_C7D_cfg1p_cfg2p_cfg1n_cfg2n_alpha_beta_gamma.
     * @note Subsequent builds overwrite this buffer.
     */
    const Eigen::Tensor<doubleC, 7, Eigen::ColMajor>& build_one_body_pp(const Eigen::MatrixXcd& OneBody_C2D_spp_spp);

    /**
     * @brief Build pn kernels using Pfaffians and gauge quadrature.
     * @math result(Ω) = Cπ(Ω) [(-1)^r₁ Aν(Ω) Oᵀ]ᵀ.
     * @output Reference to result_C7D_cfg1p_cfg2p_cfg1n_cfg2n_alpha_beta_gamma.
     * @note Subsequent builds overwrite this buffer.
     */
    const Eigen::Tensor<doubleC, 7, Eigen::ColMajor>& build_one_body_pn(const Eigen::MatrixXcd& OneBody_C2D_spp_spn);

    /**
     * @brief Build nn kernels using Pfaffians and gauge quadrature.
     * @math result(Ω) = Nπ(Ω) Hν²(Ω).
     * @output Reference to result_C7D_cfg1p_cfg2p_cfg1n_cfg2n_alpha_beta_gamma.
     * @note Subsequent builds overwrite this buffer.
     */
    const Eigen::Tensor<doubleC, 7, Eigen::ColMajor>& build_two_body_nn(const Eigen::VectorXd& gQ_F1D_i, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& Q_C3D_spn_spn_i, const Eigen::VectorXd& gP_F1D_i, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& P_C3D_spn_spn_i);

    /**
     * @brief Build pp kernels using Pfaffians and gauge quadrature.
     * @math result(Ω) = Hπ²(Ω) Nν(Ω).
     * @output Reference to result_C7D_cfg1p_cfg2p_cfg1n_cfg2n_alpha_beta_gamma.
     * @note Subsequent builds overwrite this buffer.
     */
    const Eigen::Tensor<doubleC, 7, Eigen::ColMajor>& build_two_body_pp(const Eigen::VectorXd& gQ_F1D_i, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& Q_C3D_spp_spp_i, const Eigen::VectorXd& gP_F1D_i, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& P_C3D_spp_spp_i);

    /**
     * @brief Build pn kernels using Pfaffians and gauge quadrature.
     * @math Hpn = ½ Σi gQi (Qνi† Qπi + Qνi Qπi†).
     * @output Reference to result_C7D_cfg1p_cfg2p_cfg1n_cfg2n_alpha_beta_gamma.
     * @note Subsequent builds overwrite this buffer.
     */
    const Eigen::Tensor<doubleC, 7, Eigen::ColMajor>& build_two_body_pn(const Eigen::VectorXd& gQ_F1D_i, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& Q_C3D_spn_spn_i, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& Q_C3D_spp_spp_i);

    /**
     * @brief Integrate cached nuclear kernels using Euler quadrature.
     * @math result = (2I+1)/VΩ ΣΩ wΩ Dᴵ*(Ω) cache(Ω).
     * @output Projected tensor for one I, returned by value.
     * @note Cache includes both species' particle-number projections.
     */
    Eigen::Tensor<doubleC, 6, Eigen::ColMajor> integrate_amp(int TargetTwoI_I, const Eigen::Tensor<doubleC, 7, Eigen::ColMajor>& cache_C7D_cfg1p_cfg2p_cfg1n_cfg2n_alpha_beta_gamma);

private:
    /**
     * @brief Build nuclear kernels using gauge quadrature.
     * @math cache_ab(Ω) = ⟨Φ₁;a|O Pνᴺ Pπᶻ R(Ω)|Φ₂;b⟩.
     * @output Reference to result(cfg1p,cfg2p,cfg1n,cfg2n,α,β,γ).
     * @note Excludes K selection, Euler weights, and Wigner factors.
     */
    const Eigen::Tensor<doubleC, 7, Eigen::ColMajor>& build_amp(const std::function<void(HFBProjectionPNP&, HFBProjectionPNP&, Eigen::Tensor<doubleC, 4, Eigen::ColMajor>&)>& kernel_Func);
};

inline const Eigen::Tensor<doubleC, 7, Eigen::ColMajor>& HFBProjectionNucleus::build_amp(const std::function<void(HFBProjectionPNP&, HFBProjectionPNP&, Eigen::Tensor<doubleC, 4, Eigen::ColMajor>&)>& kernel_Func) {
    auto& neutron = projection_neutron;
    auto& proton = projection_proton;

    // Ων = Ωπ.
    assert(neutron.Nalpha_I > 0 && neutron.Nbeta_I > 0 && neutron.Ngamma_I > 0);
    assert(neutron.Nalpha_I == proton.Nalpha_I && neutron.Nbeta_I == proton.Nbeta_I && neutron.Ngamma_I == proton.Ngamma_I);
    assert(neutron.alpha_F1D_alpha.isApprox(proton.alpha_F1D_alpha));
    assert(neutron.beta_F1D_beta.isApprox(proton.beta_F1D_beta));
    assert(neutron.gamma_F1D_gamma.isApprox(proton.gamma_F1D_gamma));

    // cache ∈ ℂ^{cfg1p×cfg2p×cfg1n×cfg2n×Nα×Nβ×Nγ}.
    result_C7D_cfg1p_cfg2p_cfg1n_cfg2n_alpha_beta_gamma.resize(proton.Ncfg1_I, proton.Ncfg2_I, neutron.Ncfg1_I, neutron.Ncfg2_I, neutron.Nalpha_I, neutron.Nbeta_I, neutron.Ngamma_I);

    // Thread-local PNP workspaces; disjoint angular cache blocks.
    #pragma omp parallel
    {
        HFBProjectionPNP neutron_thread(neutron.Nsp_I, neutron.hfb_pfaffian.config1_I2D_cfg1_cqp1, neutron.hfb_pfaffian.config2_I2D_cfg2_cqp2, neutron.TargetN_I, neutron.Nphi_I);
        HFBProjectionPNP proton_thread(proton.Nsp_I, proton.hfb_pfaffian.config1_I2D_cfg1_cqp1, proton.hfb_pfaffian.config2_I2D_cfg2_cqp2, proton.TargetN_I, proton.Nphi_I);
        neutron_thread.U1_C2D_sp_qp1 = neutron.U1_C2D_sp_qp1;
        neutron_thread.V1_C2D_sp_qp1 = neutron.V1_C2D_sp_qp1;
        neutron_thread.U2_C2D_sp_qp2 = neutron.U2_C2D_sp_qp2;
        neutron_thread.V2_C2D_sp_qp2 = neutron.V2_C2D_sp_qp2;
        proton_thread.U1_C2D_sp_qp1 = proton.U1_C2D_sp_qp1;
        proton_thread.V1_C2D_sp_qp1 = proton.V1_C2D_sp_qp1;
        proton_thread.U2_C2D_sp_qp2 = proton.U2_C2D_sp_qp2;
        proton_thread.V2_C2D_sp_qp2 = proton.V2_C2D_sp_qp2;
        Eigen::Tensor<doubleC, 4, Eigen::ColMajor> tmp_thread_C4D_cfg1p_cfg2p_cfg1n_cfg2n(proton.Ncfg1_I, proton.Ncfg2_I, neutron.Ncfg1_I, neutron.Ncfg2_I);
        auto kernel_thread_Func = kernel_Func;

        // Ω → Dν,Dπ → cache(Ω).
        #pragma omp for schedule(static)
        for (int beta_I = 0; beta_I < neutron.Nbeta_I; ++beta_I) {
            for (int alpha_I = 0; alpha_I < neutron.Nalpha_I; ++alpha_I) {
                for (int gamma_I = 0; gamma_I < neutron.Ngamma_I; ++gamma_I) {
                    {
                        // D = Rz(α) Ry(β) Rz(γ).
                        const Eigen::Map<const Eigen::MatrixXcd> RzAlpha_C2D_sp_sp(neutron.Rz_C3D_sp_sp_alpha.data() + static_cast<Eigen::Index>(alpha_I) * neutron.Nsp_I * neutron.Nsp_I, neutron.Nsp_I, neutron.Nsp_I);
                        const Eigen::Map<const Eigen::MatrixXcd> RyBeta_C2D_sp_sp(neutron.Ry_C3D_sp_sp_beta.data() + static_cast<Eigen::Index>(beta_I) * neutron.Nsp_I * neutron.Nsp_I, neutron.Nsp_I, neutron.Nsp_I);
                        const Eigen::Map<const Eigen::MatrixXcd> RzGamma_C2D_sp_sp(neutron.Rz_C3D_sp_sp_gamma.data() + static_cast<Eigen::Index>(gamma_I) * neutron.Nsp_I * neutron.Nsp_I, neutron.Nsp_I, neutron.Nsp_I);
                        neutron_thread.D_C2D_sp_sp.noalias() = RzAlpha_C2D_sp_sp * RyBeta_C2D_sp_sp * RzGamma_C2D_sp_sp;
                    }
                    {
                        // D = Rz(α) Ry(β) Rz(γ).
                        const Eigen::Map<const Eigen::MatrixXcd> RzAlpha_C2D_sp_sp(proton.Rz_C3D_sp_sp_alpha.data() + static_cast<Eigen::Index>(alpha_I) * proton.Nsp_I * proton.Nsp_I, proton.Nsp_I, proton.Nsp_I);
                        const Eigen::Map<const Eigen::MatrixXcd> RyBeta_C2D_sp_sp(proton.Ry_C3D_sp_sp_beta.data() + static_cast<Eigen::Index>(beta_I) * proton.Nsp_I * proton.Nsp_I, proton.Nsp_I, proton.Nsp_I);
                        const Eigen::Map<const Eigen::MatrixXcd> RzGamma_C2D_sp_sp(proton.Rz_C3D_sp_sp_gamma.data() + static_cast<Eigen::Index>(gamma_I) * proton.Nsp_I * proton.Nsp_I, proton.Nsp_I, proton.Nsp_I);
                        proton_thread.D_C2D_sp_sp.noalias() = RzAlpha_C2D_sp_sp * RyBeta_C2D_sp_sp * RzGamma_C2D_sp_sp;
                    }
                    kernel_thread_Func(neutron_thread, proton_thread, tmp_thread_C4D_cfg1p_cfg2p_cfg1n_cfg2n);

                    // ColMajor: cfg1p → cfg2p → cfg1n → cfg2n → α → β → γ.
                    const Eigen::Index offset_I = (alpha_I + static_cast<Eigen::Index>(neutron.Nalpha_I) * (beta_I + static_cast<Eigen::Index>(neutron.Nbeta_I) * gamma_I)) * Ncfg1pcfg2pcfg1ncfg2n_I;
                    Eigen::TensorMap<Eigen::Tensor<doubleC, 4, Eigen::ColMajor>> cache_C4D_cfg1p_cfg2p_cfg1n_cfg2n(result_C7D_cfg1p_cfg2p_cfg1n_cfg2n_alpha_beta_gamma.data() + offset_I, proton.Ncfg1_I, proton.Ncfg2_I, neutron.Ncfg1_I, neutron.Ncfg2_I);
                    cache_C4D_cfg1p_cfg2p_cfg1n_cfg2n = tmp_thread_C4D_cfg1p_cfg2p_cfg1n_cfg2n;
                }
            }
        }
    }
    return result_C7D_cfg1p_cfg2p_cfg1n_cfg2n_alpha_beta_gamma;
}

inline Eigen::Tensor<doubleC, 6, Eigen::ColMajor> HFBProjectionNucleus::integrate_amp(int TargetTwoI_I, const Eigen::Tensor<doubleC, 7, Eigen::ColMajor>& cache_C7D_cfg1p_cfg2p_cfg1n_cfg2n_alpha_beta_gamma) {
    auto& neutron = projection_neutron;
    auto& proton = projection_proton;

    // Ων = Ωπ; wν = wπ.
    assert(neutron.alpha_F1D_alpha.isApprox(proton.alpha_F1D_alpha));
    assert(neutron.beta_F1D_beta.isApprox(proton.beta_F1D_beta));
    assert(neutron.gamma_F1D_gamma.isApprox(proton.gamma_F1D_gamma));
    assert(neutron.weight_F1D_alpha.isApprox(proton.weight_F1D_alpha));
    assert(neutron.weight_F1D_beta.isApprox(proton.weight_F1D_beta));
    assert(neutron.weight_F1D_gamma.isApprox(proton.weight_F1D_gamma));

    assert(neutron.TwoK1_I1D_cfg1.size() == 0 || neutron.TwoK1_I1D_cfg1.size() == neutron.Ncfg1_I);
    assert(neutron.TwoK2_I1D_cfg2.size() == 0 || neutron.TwoK2_I1D_cfg2.size() == neutron.Ncfg2_I);
    assert(proton.TwoK1_I1D_cfg1.size() == 0 || proton.TwoK1_I1D_cfg1.size() == proton.Ncfg1_I);
    assert(proton.TwoK2_I1D_cfg2.size() == 0 || proton.TwoK2_I1D_cfg2.size() == proton.Ncfg2_I);

    // cache ∈ ℂ^{cfg1p×cfg2p×cfg1n×cfg2n×Nα×Nβ×Nγ}.
    assert(cache_C7D_cfg1p_cfg2p_cfg1n_cfg2n_alpha_beta_gamma.dimension(0) == proton.Ncfg1_I);
    assert(cache_C7D_cfg1p_cfg2p_cfg1n_cfg2n_alpha_beta_gamma.dimension(1) == proton.Ncfg2_I);
    assert(cache_C7D_cfg1p_cfg2p_cfg1n_cfg2n_alpha_beta_gamma.dimension(2) == neutron.Ncfg1_I);
    assert(cache_C7D_cfg1p_cfg2p_cfg1n_cfg2n_alpha_beta_gamma.dimension(3) == neutron.Ncfg2_I);
    assert(cache_C7D_cfg1p_cfg2p_cfg1n_cfg2n_alpha_beta_gamma.dimension(4) == neutron.Nalpha_I);
    assert(cache_C7D_cfg1p_cfg2p_cfg1n_cfg2n_alpha_beta_gamma.dimension(5) == neutron.Nbeta_I);
    assert(cache_C7D_cfg1p_cfg2p_cfg1n_cfg2n_alpha_beta_gamma.dimension(6) == neutron.Ngamma_I);

    // VΩ = Σwα Σwβ Σwγ; wβ includes sinβ.
    const double volume_F = neutron.weight_F1D_alpha.sum() * neutron.weight_F1D_beta.sum() * neutron.weight_F1D_gamma.sum();
    assert(std::isfinite(volume_F) && volume_F > 0.0);

    assert(TargetTwoI_I >= 0);
    const double normalization_F = (TargetTwoI_I + 1.0) / volume_F;
    Eigen::Tensor<doubleC, 6, Eigen::ColMajor> result_C6D_cfg1p_cfg2p_cfg1n_cfg2n_K1_K2(proton.Ncfg1_I, proton.Ncfg2_I, neutron.Ncfg1_I, neutron.Ncfg2_I, TargetTwoI_I + 1, TargetTwoI_I + 1);
    result_C6D_cfg1p_cfg2p_cfg1n_cfg2n_K1_K2.setZero();

    /**
     * @brief Select twice-K values from configuration labels.
     * @math 2K = 2Kn + 2Kp; −2I ≤ 2K ≤ 2I.
     * @output Inclusive TwoK range; invalid labels yield an empty range.
     */
    const auto select_K = [&](int TargetTwoI_I, const Eigen::VectorXi& TwoKn_I1D_cfg, const Eigen::VectorXi& TwoKp_I1D_cfg, int cfgn_I, int cfgp_I) -> std::pair<int, int> {
        if (TwoKn_I1D_cfg.size() == 0 || TwoKp_I1D_cfg.size() == 0) { return {-TargetTwoI_I, TargetTwoI_I}; }

        const int TwoK_I = TwoKn_I1D_cfg(cfgn_I) + TwoKp_I1D_cfg(cfgp_I);
        if (std::abs(TwoK_I) > TargetTwoI_I || (TwoK_I + TargetTwoI_I) % 2 != 0) { return {1, -1}; }

        return {TwoK_I, TwoK_I};
    };

    // result = (2I+1)/VΩ ΣΩ wΩ Dᴵ*(Ω) cache(Ω).
    RepresentationSpin wignerD_X(TargetTwoI_I);
    for (int beta_I = 0; beta_I < neutron.Nbeta_I; ++beta_I) {
        const auto& wignerDyBeta_C2D_K_K = wignerD_X.calc_Ry(neutron.beta_F1D_beta(beta_I));
        for (int alpha_I = 0; alpha_I < neutron.Nalpha_I; ++alpha_I) {
            for (int gamma_I = 0; gamma_I < neutron.Ngamma_I; ++gamma_I) {
                const double weight_F = neutron.weight_F1D_alpha(alpha_I) * neutron.weight_F1D_beta(beta_I) * neutron.weight_F1D_gamma(gamma_I);
                for (int cfg2n_I = 0; cfg2n_I < neutron.Ncfg2_I; ++cfg2n_I) {
                    for (int cfg2p_I = 0; cfg2p_I < proton.Ncfg2_I; ++cfg2p_I) {
                        const auto [TwoK2min_I, TwoK2max_I] = select_K(TargetTwoI_I, neutron.TwoK2_I1D_cfg2, proton.TwoK2_I1D_cfg2, cfg2n_I, cfg2p_I);
                        for (int cfg1n_I = 0; cfg1n_I < neutron.Ncfg1_I; ++cfg1n_I) {
                            for (int cfg1p_I = 0; cfg1p_I < proton.Ncfg1_I; ++cfg1p_I) {
                                const auto [TwoK1min_I, TwoK1max_I] = select_K(TargetTwoI_I, neutron.TwoK1_I1D_cfg1, proton.TwoK1_I1D_cfg1, cfg1n_I, cfg1p_I);
                                for (int TwoK2_I = TwoK2min_I; TwoK2_I <= TwoK2max_I; TwoK2_I += 2) {
                                    const int K2_I = (TwoK2_I + TargetTwoI_I) / 2;
                                    const double K2_F = 0.5 * TwoK2_I;
                                    for (int TwoK1_I = TwoK1min_I; TwoK1_I <= TwoK1max_I; TwoK1_I += 2) {
                                        // Dᴵ* = exp(iK₁α) dᴵ* exp(iK₂γ).
                                        const int K1_I = (TwoK1_I + TargetTwoI_I) / 2;
                                        const double K1_F = 0.5 * TwoK1_I;
                                        const doubleC factorOmega_C = std::conj(wignerDyBeta_C2D_K_K(K1_I, K2_I)) * std::exp(doubleC(0.0, K1_F * neutron.alpha_F1D_alpha(alpha_I) + K2_F * neutron.gamma_F1D_gamma(gamma_I)));
                                        const doubleC weight_factor_C = normalization_F * weight_F * factorOmega_C;
                                        result_C6D_cfg1p_cfg2p_cfg1n_cfg2n_K1_K2(cfg1p_I, cfg2p_I, cfg1n_I, cfg2n_I, K1_I, K2_I) += weight_factor_C * cache_C7D_cfg1p_cfg2p_cfg1n_cfg2n_alpha_beta_gamma(cfg1p_I, cfg2p_I, cfg1n_I, cfg2n_I, alpha_I, beta_I, gamma_I);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return result_C6D_cfg1p_cfg2p_cfg1n_cfg2n_K1_K2;
}

inline const Eigen::Tensor<doubleC, 7, Eigen::ColMajor>& HFBProjectionNucleus::build_overlap() {
    return build_amp( [&](HFBProjectionPNP& neutron, HFBProjectionPNP& proton, Eigen::Tensor<doubleC, 4, Eigen::ColMajor>& tmp_C4D_cfg1p_cfg2p_cfg1n_cfg2n) {
        const auto& overlap_C2D_cfg1n_cfg2n = neutron.calc_overlap_pnp();
        const auto& overlap_C2D_cfg1p_cfg2p = proton.calc_overlap_pnp();

        // ColMajor: (p₁,p₂) → rows; (n₁,n₂) → columns.
        const Eigen::Map<const Eigen::VectorXcd> overlap_C1D_cfg1ncfg2n(overlap_C2D_cfg1n_cfg2n.data(), Ncfg1ncfg2n_I);
        const Eigen::Map<const Eigen::VectorXcd> overlap_C1D_cfg1pcfg2p(overlap_C2D_cfg1p_cfg2p.data(), Ncfg1pcfg2p_I);
        Eigen::Map<Eigen::MatrixXcd> overlap_C2D_cfg1pcfg2p_cfg1ncfg2n(tmp_C4D_cfg1p_cfg2p_cfg1n_cfg2n.data(), Ncfg1pcfg2p_I, Ncfg1ncfg2n_I);
        overlap_C2D_cfg1pcfg2p_cfg1ncfg2n.noalias() = overlap_C1D_cfg1pcfg2p * overlap_C1D_cfg1ncfg2n.transpose();
    });
}

inline const Eigen::Tensor<doubleC, 7, Eigen::ColMajor>& HFBProjectionNucleus::build_one_body_nn(const Eigen::MatrixXcd& OneBody_C2D_spn_spn) {
    auto& neutron = projection_neutron;
    auto& proton = projection_proton;
    assert(OneBody_C2D_spn_spn.rows() == neutron.Nsp_I && OneBody_C2D_spn_spn.cols() == neutron.Nsp_I && OneBody_C2D_spn_spn.allFinite());

    return build_amp([&](HFBProjectionPNP& neutron, HFBProjectionPNP& proton, Eigen::Tensor<doubleC, 4, Eigen::ColMajor>& tmp_C4D_cfg1p_cfg2p_cfg1n_cfg2n) {
        const auto& overlap_C2D_cfg1p_cfg2p = proton.calc_overlap_pnp();
        const auto& onebody_C2D_cfg1n_cfg2n = neutron.calc_one_body_pnp(OneBody_C2D_spn_spn);

        // ColMajor: (p₁,p₂) → rows; (n₁,n₂) → columns.
        const Eigen::Map<const Eigen::VectorXcd> onebody_C1D_cfg1ncfg2n(onebody_C2D_cfg1n_cfg2n.data(), Ncfg1ncfg2n_I);
        const Eigen::Map<const Eigen::VectorXcd> overlap_C1D_cfg1pcfg2p(overlap_C2D_cfg1p_cfg2p.data(), Ncfg1pcfg2p_I);
        Eigen::Map<Eigen::MatrixXcd> onebody_C2D_cfg1pcfg2p_cfg1ncfg2n(tmp_C4D_cfg1p_cfg2p_cfg1n_cfg2n.data(), Ncfg1pcfg2p_I, Ncfg1ncfg2n_I);
        onebody_C2D_cfg1pcfg2p_cfg1ncfg2n.noalias() = overlap_C1D_cfg1pcfg2p * onebody_C1D_cfg1ncfg2n.transpose();
    });
}

inline const Eigen::Tensor<doubleC, 7, Eigen::ColMajor>& HFBProjectionNucleus::build_one_body_pp(const Eigen::MatrixXcd& OneBody_C2D_spp_spp) {
    auto& neutron = projection_neutron;
    auto& proton = projection_proton;
    assert(OneBody_C2D_spp_spp.rows() == proton.Nsp_I && OneBody_C2D_spp_spp.cols() == proton.Nsp_I && OneBody_C2D_spp_spp.allFinite());

    return build_amp([&](HFBProjectionPNP& neutron, HFBProjectionPNP& proton, Eigen::Tensor<doubleC, 4, Eigen::ColMajor>& tmp_C4D_cfg1p_cfg2p_cfg1n_cfg2n) {
        const auto& overlap_C2D_cfg1n_cfg2n = neutron.calc_overlap_pnp();
        const auto& onebody_C2D_cfg1p_cfg2p = proton.calc_one_body_pnp(OneBody_C2D_spp_spp);

        // ColMajor: (p₁,p₂) → rows; (n₁,n₂) → columns.
        const Eigen::Map<const Eigen::VectorXcd> overlap_C1D_cfg1ncfg2n(overlap_C2D_cfg1n_cfg2n.data(), Ncfg1ncfg2n_I);
        const Eigen::Map<const Eigen::VectorXcd> onebody_C1D_cfg1pcfg2p(onebody_C2D_cfg1p_cfg2p.data(), Ncfg1pcfg2p_I);
        Eigen::Map<Eigen::MatrixXcd> onebody_C2D_cfg1pcfg2p_cfg1ncfg2n(tmp_C4D_cfg1p_cfg2p_cfg1n_cfg2n.data(), Ncfg1pcfg2p_I, Ncfg1ncfg2n_I);
        onebody_C2D_cfg1pcfg2p_cfg1ncfg2n.noalias() = onebody_C1D_cfg1pcfg2p * overlap_C1D_cfg1ncfg2n.transpose();
    });
}

inline const Eigen::Tensor<doubleC, 7, Eigen::ColMajor>& HFBProjectionNucleus::build_one_body_pn(const Eigen::MatrixXcd& OneBody_C2D_spp_spn) {
    auto& neutron = projection_neutron;
    auto& proton = projection_proton;
    assert(OneBody_C2D_spp_spn.rows() == proton.Nsp_I && OneBody_C2D_spp_spn.cols() == neutron.Nsp_I && OneBody_C2D_spp_spn.allFinite());

    // ColMajor: (p₁,p₂) → rows; (n₁,n₂) → columns.

    // Aν Oᵀ ∈ ℂ^{Ncfg1ncfg2n×Nspp}.
    Eigen::MatrixXcd tmp_C2D_cfg1ncfg2n_spp{};
    tmp_C2D_cfg1ncfg2n_spp.resize(Ncfg1ncfg2n_I, proton.Nsp_I);

    return build_amp([&, tmp_C2D_cfg1ncfg2n_spp = std::move(tmp_C2D_cfg1ncfg2n_spp)](HFBProjectionPNP& neutron, HFBProjectionPNP& proton, Eigen::Tensor<doubleC, 4, Eigen::ColMajor>& tmp_C4D_cfg1p_cfg2p_cfg1n_cfg2n) mutable {
        Eigen::Map<Eigen::MatrixXcd> tmp_C2D_cfg1pcfg2p_cfg1ncfg2n(tmp_C4D_cfg1p_cfg2p_cfg1n_cfg2n.data(), Ncfg1pcfg2p_I, Ncfg1ncfg2n_I);

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

inline const Eigen::Tensor<doubleC, 7, Eigen::ColMajor>& HFBProjectionNucleus::build_two_body_nn(const Eigen::VectorXd& gQ_F1D_i, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& Q_C3D_spn_spn_i, const Eigen::VectorXd& gP_F1D_i, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& P_C3D_spn_spn_i) {
    return build_amp([&](HFBProjectionPNP& neutron, HFBProjectionPNP& proton, Eigen::Tensor<doubleC, 4, Eigen::ColMajor>& tmp_C4D_cfg1p_cfg2p_cfg1n_cfg2n) {
        // Hν²(Ω) Nπ(Ω); weighted Q/P kernel; pairing ¼ included.
        const auto& overlap_C2D_cfg1p_cfg2p = proton.calc_overlap_pnp();
        const auto& twobody_C2D_cfg1n_cfg2n = neutron.calc_two_body_pnp(gQ_F1D_i, Q_C3D_spn_spn_i, gP_F1D_i, P_C3D_spn_spn_i);

        // ColMajor: (p₁,p₂) → rows; (n₁,n₂) → columns.
        const Eigen::Map<const Eigen::VectorXcd> twobody_C1D_cfg1ncfg2n(twobody_C2D_cfg1n_cfg2n.data(), Ncfg1ncfg2n_I);
        const Eigen::Map<const Eigen::VectorXcd> overlap_C1D_cfg1pcfg2p(overlap_C2D_cfg1p_cfg2p.data(), Ncfg1pcfg2p_I);
        Eigen::Map<Eigen::MatrixXcd> twobody_C2D_cfg1pcfg2p_cfg1ncfg2n(tmp_C4D_cfg1p_cfg2p_cfg1n_cfg2n.data(), Ncfg1pcfg2p_I, Ncfg1ncfg2n_I);
        twobody_C2D_cfg1pcfg2p_cfg1ncfg2n.noalias() = overlap_C1D_cfg1pcfg2p * twobody_C1D_cfg1ncfg2n.transpose();
    });
}

inline const Eigen::Tensor<doubleC, 7, Eigen::ColMajor>& HFBProjectionNucleus::build_two_body_pp(const Eigen::VectorXd& gQ_F1D_i, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& Q_C3D_spp_spp_i, const Eigen::VectorXd& gP_F1D_i, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& P_C3D_spp_spp_i) {
    return build_amp([&](HFBProjectionPNP& neutron, HFBProjectionPNP& proton, Eigen::Tensor<doubleC, 4, Eigen::ColMajor>& tmp_C4D_cfg1p_cfg2p_cfg1n_cfg2n) {
        // Nν(Ω) Hπ²(Ω); weighted Q/P kernel; pairing ¼ included.
        const auto& overlap_C2D_cfg1n_cfg2n = neutron.calc_overlap_pnp();
        const auto& twobody_C2D_cfg1p_cfg2p = proton.calc_two_body_pnp(gQ_F1D_i, Q_C3D_spp_spp_i, gP_F1D_i, P_C3D_spp_spp_i);

        // ColMajor: (p₁,p₂) → rows; (n₁,n₂) → columns.
        const Eigen::Map<const Eigen::VectorXcd> overlap_C1D_cfg1ncfg2n(overlap_C2D_cfg1n_cfg2n.data(), Ncfg1ncfg2n_I);
        const Eigen::Map<const Eigen::VectorXcd> twobody_C1D_cfg1pcfg2p(twobody_C2D_cfg1p_cfg2p.data(), Ncfg1pcfg2p_I);
        Eigen::Map<Eigen::MatrixXcd> twobody_C2D_cfg1pcfg2p_cfg1ncfg2n(tmp_C4D_cfg1p_cfg2p_cfg1n_cfg2n.data(), Ncfg1pcfg2p_I, Ncfg1ncfg2n_I);
        twobody_C2D_cfg1pcfg2p_cfg1ncfg2n.noalias() = twobody_C1D_cfg1pcfg2p * overlap_C1D_cfg1ncfg2n.transpose();
    });
}

inline const Eigen::Tensor<doubleC, 7, Eigen::ColMajor>& HFBProjectionNucleus::build_two_body_pn(const Eigen::VectorXd& gQ_F1D_i, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& Q_C3D_spn_spn_i, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& Q_C3D_spp_spp_i) {
    auto& neutron = projection_neutron;
    auto& proton = projection_proton;

    // NQν = NQπ = NgQ.
    assert(gQ_F1D_i.allFinite());
    assert(Q_C3D_spn_spn_i.dimension(0) == neutron.Nsp_I && Q_C3D_spn_spn_i.dimension(1) == neutron.Nsp_I && Q_C3D_spn_spn_i.dimension(2) == gQ_F1D_i.size());
    assert(Q_C3D_spp_spp_i.dimension(0) == proton.Nsp_I && Q_C3D_spp_spp_i.dimension(1) == proton.Nsp_I && Q_C3D_spp_spp_i.dimension(2) == gQ_F1D_i.size());
    assert(Eigen::Map<const Eigen::VectorXcd>(Q_C3D_spn_spn_i.data(), Q_C3D_spn_spn_i.size()).allFinite());
    assert(Eigen::Map<const Eigen::VectorXcd>(Q_C3D_spp_spp_i.data(), Q_C3D_spp_spp_i.size()).allFinite());

    // (Qν†,Qπ) and (Qν,Qπ†) reuse species workspaces.
    Eigen::MatrixXcd Qdag_C2D_spn_spn{};
    Eigen::MatrixXcd Q_C2D_spp_spp{};
    Qdag_C2D_spn_spn.resize(neutron.Nsp_I, neutron.Nsp_I);
    Q_C2D_spp_spp.resize(proton.Nsp_I, proton.Nsp_I);

    return build_amp([&, Qdag_C2D_spn_spn = std::move(Qdag_C2D_spn_spn), Q_C2D_spp_spp = std::move(Q_C2D_spp_spp)](HFBProjectionPNP& neutron, HFBProjectionPNP& proton, Eigen::Tensor<doubleC, 4, Eigen::ColMajor>& tmp_C4D_cfg1p_cfg2p_cfg1n_cfg2n) mutable {
        Eigen::Map<Eigen::MatrixXcd> tmp_C2D_cfg1pcfg2p_cfg1ncfg2n(tmp_C4D_cfg1p_cfg2p_cfg1n_cfg2n.data(), Ncfg1pcfg2p_I, Ncfg1ncfg2n_I);

        tmp_C2D_cfg1pcfg2p_cfg1ncfg2n.setZero();
        for (Eigen::Index iQ_I = 0; iQ_I < gQ_F1D_i.size(); ++iQ_I) {
            // Hi = ½ gQi (Qνi† Qπi + Qνi Qπi†).
            const Eigen::Map<const Eigen::MatrixXcd> Q_C2D_spn_spn(Q_C3D_spn_spn_i.data() + static_cast<Eigen::Index>(iQ_I) * neutron.Nsp_I * neutron.Nsp_I, neutron.Nsp_I, neutron.Nsp_I);
            const Eigen::Map<const Eigen::MatrixXcd> Qp_C2D_spp_spp(Q_C3D_spp_spp_i.data() + static_cast<Eigen::Index>(iQ_I) * proton.Nsp_I * proton.Nsp_I, proton.Nsp_I, proton.Nsp_I);
            // ½ gQi Qνi† Qπi.
            {
                Qdag_C2D_spn_spn = Q_C2D_spn_spn.adjoint();
                Q_C2D_spp_spp = Qp_C2D_spp_spp;
                const auto& onebody_C2D_cfg1n_cfg2n = neutron.calc_one_body_pnp(Qdag_C2D_spn_spn);
                const auto& onebody_C2D_cfg1p_cfg2p = proton.calc_one_body_pnp(Q_C2D_spp_spp);

                // Both kernels include overlaps and gauge normalizations.
                const Eigen::Map<const Eigen::VectorXcd> onebody_C1D_cfg1ncfg2n(onebody_C2D_cfg1n_cfg2n.data(), Ncfg1ncfg2n_I);
                const Eigen::Map<const Eigen::VectorXcd> onebody_C1D_cfg1pcfg2p(onebody_C2D_cfg1p_cfg2p.data(), Ncfg1pcfg2p_I);
                tmp_C2D_cfg1pcfg2p_cfg1ncfg2n.noalias() += 0.5 * gQ_F1D_i(iQ_I) * onebody_C1D_cfg1pcfg2p * onebody_C1D_cfg1ncfg2n.transpose();
            }
            // ½ gQi Qνi Qπi†.
            {
                Qdag_C2D_spn_spn = Q_C2D_spn_spn;
                Q_C2D_spp_spp = Qp_C2D_spp_spp.adjoint();
                const auto& onebody_C2D_cfg1n_cfg2n = neutron.calc_one_body_pnp(Qdag_C2D_spn_spn);
                const auto& onebody_C2D_cfg1p_cfg2p = proton.calc_one_body_pnp(Q_C2D_spp_spp);

                // Both kernels include overlaps and gauge normalizations.
                const Eigen::Map<const Eigen::VectorXcd> onebody_C1D_cfg1ncfg2n(onebody_C2D_cfg1n_cfg2n.data(), Ncfg1ncfg2n_I);
                const Eigen::Map<const Eigen::VectorXcd> onebody_C1D_cfg1pcfg2p(onebody_C2D_cfg1p_cfg2p.data(), Ncfg1pcfg2p_I);
                tmp_C2D_cfg1pcfg2p_cfg1ncfg2n.noalias() += 0.5 * gQ_F1D_i(iQ_I) * onebody_C1D_cfg1pcfg2p * onebody_C1D_cfg1ncfg2n.transpose();
            }
        }
    });
}
