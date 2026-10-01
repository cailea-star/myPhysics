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
     * @brief  Integrate neutron Q/P kernels using Pfaffians and quadrature.
     * @math   K = Σi gQi KQi + Σi gPi KPi; P† = ½Σab Pab c†a c†b.
     * @output Tensor indexed by cfg1p, cfg2p, cfg1n, cfg2n, K1, K2.
     * @note   gQ,gP include signs; pairing ¼ is internal.
     * @note   Q/P lists may differ in length; zero coefficients skip terms.
     */
    const Eigen::Tensor<doubleC, 6, Eigen::ColMajor>& calc_two_body_nnnn(const Eigen::VectorXd& gQ_F1D_i, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& Q_C3D_spn_spn_i, const Eigen::VectorXd& gP_F1D_i, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& P_C3D_spn_spn_i);

    /**
     * @brief  Integrate proton Q/P kernels using Pfaffians and quadrature.
     * @math   K = Σi gQi KQi + Σi gPi KPi; P† = ½Σab Pab c†a c†b.
     * @output Tensor indexed by cfg1p, cfg2p, cfg1n, cfg2n, K1, K2.
     * @note   gQ,gP include signs; pairing ¼ is internal.
     * @note   Q/P lists may differ in length; zero coefficients skip terms.
     */
    const Eigen::Tensor<doubleC, 6, Eigen::ColMajor>& calc_two_body_pppp(const Eigen::VectorXd& gQ_F1D_i, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& Q_C3D_spp_spp_i, const Eigen::VectorXd& gP_F1D_i, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& P_C3D_spp_spp_i);

    /**
     * @brief  Integrate mixed QQ kernels using Pfaffians and quadrature.
     * @math   O = Σμ Q̂ν,λμ† Q̂π,λμ; cν†cν cπ†cπ = cν†cπ†cπcν.
     * @math   Qτ ∈ ℂ^{Nspτ × Nspτ × (2λ+1)}; twoLambda_I = 2λ ≥ 0, λ ∈ ℤ; μ = −λ,…,λ.
     * @output Tensor indexed by cfg1p, cfg2p, cfg1n, cfg2n, K1, K2.
     * @note   Consistent complete spherical tensors; coupling strength is external.
     */
    const Eigen::Tensor<doubleC, 6, Eigen::ColMajor>& calc_two_body_Q_nppn(int twoLambda_I, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& Q_C3D_spn_spn_mu, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& Q_C3D_spp_spp_mu);

private:
    /**
     * @brief  Integrate nuclear kernels using Euler quadrature.
     * @math   result = (2I+1)/VΩ ΣΩ wΩ Dᴵ*(Ω) K(Ω).
     * @output Updated result_C6D_cfg1p_cfg2p_cfg1n_cfg2n_K1_K2.
     * @note   Callback receives thread-local species and kernel workspaces.
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

    assert(neutron.TwoK1_I1D_cfg1.size() == 0 || neutron.TwoK1_I1D_cfg1.size() == neutron.Ncfg1_I);
    assert(neutron.TwoK2_I1D_cfg2.size() == 0 || neutron.TwoK2_I1D_cfg2.size() == neutron.Ncfg2_I);
    assert(proton.TwoK1_I1D_cfg1.size() == 0 || proton.TwoK1_I1D_cfg1.size() == proton.Ncfg1_I);
    assert(proton.TwoK2_I1D_cfg2.size() == 0 || proton.TwoK2_I1D_cfg2.size() == proton.Ncfg2_I);

    // VΩ = Σwα Σwβ Σwγ; wβ includes sinβ.
    const double volume_F = neutron.weight_F1D_alpha.sum() * neutron.weight_F1D_beta.sum() * neutron.weight_F1D_gamma.sum();
    assert(std::isfinite(volume_F) && volume_F > 0.0);
    const double normalization_F = (TargetTwoI_I + 1.0) / volume_F;

    result_C6D_cfg1p_cfg2p_cfg1n_cfg2n_K1_K2.setZero();

    /**
     * @brief Select twice-K values from configuration labels.
     * @math 2K = 2Kn + 2Kp; −2I ≤ 2K ≤ 2I.
     * @output Inclusive TwoK range; invalid labels yield an empty range.
     */
    const auto select_K = [&](const Eigen::VectorXi& TwoKn_I1D_cfg, const Eigen::VectorXi& TwoKp_I1D_cfg, int cfgn_I, int cfgp_I) -> std::pair<int, int> {
        if (TwoKn_I1D_cfg.size() == 0 || TwoKp_I1D_cfg.size() == 0) { return {-TargetTwoI_I, TargetTwoI_I}; }

        const int TwoK_I = TwoKn_I1D_cfg(cfgn_I) + TwoKp_I1D_cfg(cfgp_I);
        if (std::abs(TwoK_I) > TargetTwoI_I || (TwoK_I + TargetTwoI_I) % 2 != 0) { return {1, -1}; }

        return {TwoK_I, TwoK_I};
    };

    // Ω → Dν,Dπ → K(Ω) → weighted projection contribution.
    const auto calc_at_angle = [&](int alpha_I, int beta_I, int gamma_I, HFBProjectionPNP& neutron_pnp, HFBProjectionPNP& proton_pnp, Eigen::Tensor<doubleC, 4, Eigen::ColMajor>& tmp_thread_C4D_cfg1p_cfg2p_cfg1n_cfg2n, const Eigen::MatrixXcd& RyI_C2D_K_K, auto& kernel_thread_Func) {
        {
            // D = Rz(α) Ry(β) Rz(γ).
            const Eigen::Map<const Eigen::MatrixXcd> RzAlpha_C2D_sp_sp(neutron.Rz_C3D_sp_sp_alpha.data() + static_cast<Eigen::Index>(alpha_I) * neutron.Nsp_I * neutron.Nsp_I, neutron.Nsp_I, neutron.Nsp_I);
            const Eigen::Map<const Eigen::MatrixXcd> RyBeta_C2D_sp_sp(neutron.Ry_C3D_sp_sp_beta.data() + static_cast<Eigen::Index>(beta_I) * neutron.Nsp_I * neutron.Nsp_I, neutron.Nsp_I, neutron.Nsp_I);
            const Eigen::Map<const Eigen::MatrixXcd> RzGamma_C2D_sp_sp(neutron.Rz_C3D_sp_sp_gamma.data() + static_cast<Eigen::Index>(gamma_I) * neutron.Nsp_I * neutron.Nsp_I, neutron.Nsp_I, neutron.Nsp_I);
            neutron_pnp.D_C2D_sp_sp.noalias() = RzAlpha_C2D_sp_sp * RyBeta_C2D_sp_sp * RzGamma_C2D_sp_sp;
        }
        {
            // D = Rz(α) Ry(β) Rz(γ).
            const Eigen::Map<const Eigen::MatrixXcd> RzAlpha_C2D_sp_sp(proton.Rz_C3D_sp_sp_alpha.data() + static_cast<Eigen::Index>(alpha_I) * proton.Nsp_I * proton.Nsp_I, proton.Nsp_I, proton.Nsp_I);
            const Eigen::Map<const Eigen::MatrixXcd> RyBeta_C2D_sp_sp(proton.Ry_C3D_sp_sp_beta.data() + static_cast<Eigen::Index>(beta_I) * proton.Nsp_I * proton.Nsp_I, proton.Nsp_I, proton.Nsp_I);
            const Eigen::Map<const Eigen::MatrixXcd> RzGamma_C2D_sp_sp(proton.Rz_C3D_sp_sp_gamma.data() + static_cast<Eigen::Index>(gamma_I) * proton.Nsp_I * proton.Nsp_I, proton.Nsp_I, proton.Nsp_I);
            proton_pnp.D_C2D_sp_sp.noalias() = RzAlpha_C2D_sp_sp * RyBeta_C2D_sp_sp * RzGamma_C2D_sp_sp;
        }

        kernel_thread_Func(neutron_pnp, proton_pnp, tmp_thread_C4D_cfg1p_cfg2p_cfg1n_cfg2n);

        // result += wΩ Dᴵ*(Ω) K(Ω); only accumulation is serialized.
        const double weight_F = normalization_F * neutron.weight_F1D_alpha(alpha_I) * neutron.weight_F1D_beta(beta_I) * neutron.weight_F1D_gamma(gamma_I);
        #pragma omp critical(psm_amp_accumulate)
        {
            for (int cfg2n_I = 0; cfg2n_I < neutron.Ncfg2_I; ++cfg2n_I) {
                for (int cfg2p_I = 0; cfg2p_I < proton.Ncfg2_I; ++cfg2p_I) {
                    const auto [TwoK2min_I, TwoK2max_I] = select_K(neutron.TwoK2_I1D_cfg2, proton.TwoK2_I1D_cfg2, cfg2n_I, cfg2p_I);
                    for (int cfg1n_I = 0; cfg1n_I < neutron.Ncfg1_I; ++cfg1n_I) {
                        for (int cfg1p_I = 0; cfg1p_I < proton.Ncfg1_I; ++cfg1p_I) {
                            const auto [TwoK1min_I, TwoK1max_I] = select_K(neutron.TwoK1_I1D_cfg1, proton.TwoK1_I1D_cfg1, cfg1n_I, cfg1p_I);
                            for (int TwoK2_I = TwoK2min_I; TwoK2_I <= TwoK2max_I; TwoK2_I += 2) {
                                const int K2_I = (TwoK2_I + TargetTwoI_I) / 2;
                                const double K2_F = 0.5 * TwoK2_I;
                                for (int TwoK1_I = TwoK1min_I; TwoK1_I <= TwoK1max_I; TwoK1_I += 2) {
                                    // Dᴵ* = exp(iK₁α) dᴵ* exp(iK₂γ).
                                    const int K1_I = (TwoK1_I + TargetTwoI_I) / 2;
                                    const double K1_F = 0.5 * TwoK1_I;
                                    const doubleC factorOmega_C = std::conj(RyI_C2D_K_K(K1_I, K2_I)) * std::exp(doubleC(0.0, K1_F * neutron.alpha_F1D_alpha(alpha_I) + K2_F * neutron.gamma_F1D_gamma(gamma_I)));
                                    const doubleC weight_factor_C = weight_F * factorOmega_C;
                                    result_C6D_cfg1p_cfg2p_cfg1n_cfg2n_K1_K2(cfg1p_I, cfg2p_I, cfg1n_I, cfg2n_I, K1_I, K2_I) += weight_factor_C * tmp_thread_C4D_cfg1p_cfg2p_cfg1n_cfg2n(cfg1p_I, cfg2p_I, cfg1n_I, cfg2n_I);
                                }
                            }
                        }
                    }
                }
            }
        }
    };


    // Thread-local contractions and kernels; rotation grids are read-only.
    #pragma omp parallel
    {
        HFBProjectionPNP neutron_thread(neutron.TargetN_I, neutron.Nsp_I, neutron.hfb_pfaffian.config1_I2D_cfg1_cqp1, neutron.hfb_pfaffian.config2_I2D_cfg2_cqp2, neutron.Nphi_I);
        HFBProjectionPNP proton_thread(proton.TargetN_I, proton.Nsp_I, proton.hfb_pfaffian.config1_I2D_cfg1_cqp1, proton.hfb_pfaffian.config2_I2D_cfg2_cqp2, proton.Nphi_I);
        neutron_thread.U1_C2D_sp_qp1 = neutron.U1_C2D_sp_qp1;
        neutron_thread.V1_C2D_sp_qp1 = neutron.V1_C2D_sp_qp1;
        neutron_thread.U2_C2D_sp_qp2 = neutron.U2_C2D_sp_qp2;
        neutron_thread.V2_C2D_sp_qp2 = neutron.V2_C2D_sp_qp2;
        proton_thread.U1_C2D_sp_qp1 = proton.U1_C2D_sp_qp1;
        proton_thread.V1_C2D_sp_qp1 = proton.V1_C2D_sp_qp1;
        proton_thread.U2_C2D_sp_qp2 = proton.U2_C2D_sp_qp2;
        proton_thread.V2_C2D_sp_qp2 = proton.V2_C2D_sp_qp2;
        RepresentationSpin representation_thread(TargetTwoI_I);
        Eigen::Tensor<doubleC, 4, Eigen::ColMajor> tmp_thread_C4D_cfg1p_cfg2p_cfg1n_cfg2n(proton.Ncfg1_I, proton.Ncfg2_I, neutron.Ncfg1_I, neutron.Ncfg2_I);
        auto kernel_thread_Func = kernel_Func;

        // dᴵ(β) is shared by both species.
        #pragma omp for schedule(static)
        for (int beta_I = 0; beta_I < neutron.Nbeta_I; ++beta_I) {
            const auto& RyI_C2D_K_K = representation_thread.calc_Ry(neutron.beta_F1D_beta(beta_I));
            for (int alpha_I = 0; alpha_I < neutron.Nalpha_I; ++alpha_I) {
                for (int gamma_I = 0; gamma_I < neutron.Ngamma_I; ++gamma_I) {
                    calc_at_angle(alpha_I, beta_I, gamma_I, neutron_thread, proton_thread, tmp_thread_C4D_cfg1p_cfg2p_cfg1n_cfg2n, RyI_C2D_K_K, kernel_thread_Func);
                }
            }
        }
    }
    return result_C6D_cfg1p_cfg2p_cfg1n_cfg2n_K1_K2;
}

inline const Eigen::Tensor<doubleC, 6, Eigen::ColMajor>& HFBProjectionNucleus::calc_overlap() {
    auto& neutron = projection_neutron;
    auto& proton = projection_proton;

    return integrate_amp([&](HFBProjectionPNP& neutron, HFBProjectionPNP& proton, Eigen::Tensor<doubleC, 4, Eigen::ColMajor>& tmp_C4D_cfg1p_cfg2p_cfg1n_cfg2n) {
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

    return integrate_amp([&](HFBProjectionPNP& neutron, HFBProjectionPNP& proton, Eigen::Tensor<doubleC, 4, Eigen::ColMajor>& tmp_C4D_cfg1p_cfg2p_cfg1n_cfg2n) {
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

    return integrate_amp([&](HFBProjectionPNP& neutron, HFBProjectionPNP& proton, Eigen::Tensor<doubleC, 4, Eigen::ColMajor>& tmp_C4D_cfg1p_cfg2p_cfg1n_cfg2n) {
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

    // Aν Oᵀ ∈ ℂ^{Ncfg1ncfg2n×Nspp}.
    Eigen::MatrixXcd tmp_C2D_cfg1ncfg2n_spp{};
    tmp_C2D_cfg1ncfg2n_spp.resize(Ncfg1ncfg2n_I, proton.Nsp_I);

    return integrate_amp([&, tmp_C2D_cfg1ncfg2n_spp = std::move(tmp_C2D_cfg1ncfg2n_spp)](HFBProjectionPNP& neutron, HFBProjectionPNP& proton, Eigen::Tensor<doubleC, 4, Eigen::ColMajor>& tmp_C4D_cfg1p_cfg2p_cfg1n_cfg2n) mutable {
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

inline const Eigen::Tensor<doubleC, 6, Eigen::ColMajor>& HFBProjectionNucleus::calc_two_body_nnnn(const Eigen::VectorXd& gQ_F1D_i, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& Q_C3D_spn_spn_i, const Eigen::VectorXd& gP_F1D_i, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& P_C3D_spn_spn_i) {
    return integrate_amp([&](HFBProjectionPNP& neutron, HFBProjectionPNP& proton, Eigen::Tensor<doubleC, 4, Eigen::ColMajor>& tmp_C4D_cfg1p_cfg2p_cfg1n_cfg2n) {
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

inline const Eigen::Tensor<doubleC, 6, Eigen::ColMajor>& HFBProjectionNucleus::calc_two_body_pppp(const Eigen::VectorXd& gQ_F1D_i, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& Q_C3D_spp_spp_i, const Eigen::VectorXd& gP_F1D_i, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& P_C3D_spp_spp_i) {
    return integrate_amp([&](HFBProjectionPNP& neutron, HFBProjectionPNP& proton, Eigen::Tensor<doubleC, 4, Eigen::ColMajor>& tmp_C4D_cfg1p_cfg2p_cfg1n_cfg2n) {
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

inline const Eigen::Tensor<doubleC, 6, Eigen::ColMajor>& HFBProjectionNucleus::calc_two_body_Q_nppn(int twoLambda_I, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& Q_C3D_spn_spn_mu, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& Q_C3D_spp_spp_mu) {
    auto& neutron = projection_neutron;
    auto& proton = projection_proton;

    // Qν,Qπ share λ and μ; single-particle dimensions may differ.
    assert(twoLambda_I >= 0 && twoLambda_I % 2 == 0);
    assert(Q_C3D_spn_spn_mu.dimension(0) == neutron.Nsp_I && Q_C3D_spn_spn_mu.dimension(1) == neutron.Nsp_I && Q_C3D_spn_spn_mu.dimension(2) == twoLambda_I + 1);
    assert(Q_C3D_spp_spp_mu.dimension(0) == proton.Nsp_I && Q_C3D_spp_spp_mu.dimension(1) == proton.Nsp_I && Q_C3D_spp_spp_mu.dimension(2) == twoLambda_I + 1);
    assert(Eigen::Map<const Eigen::VectorXcd>(Q_C3D_spn_spn_mu.data(), Q_C3D_spn_spn_mu.size()).allFinite());
    assert(Eigen::Map<const Eigen::VectorXcd>(Q_C3D_spp_spp_mu.data(), Q_C3D_spp_spp_mu.size()).allFinite());

    // Qν† = (Qν*)ᵀ; workspaces reused over Ω,μ.
    Eigen::MatrixXcd Qdag_C2D_spn_spn{};
    Eigen::MatrixXcd Q_C2D_spp_spp{};
    Qdag_C2D_spn_spn.resize(neutron.Nsp_I, neutron.Nsp_I);
    Q_C2D_spp_spp.resize(proton.Nsp_I, proton.Nsp_I);

    return integrate_amp([&, Qdag_C2D_spn_spn = std::move(Qdag_C2D_spn_spn), Q_C2D_spp_spp = std::move(Q_C2D_spp_spp)](HFBProjectionPNP& neutron, HFBProjectionPNP& proton, Eigen::Tensor<doubleC, 4, Eigen::ColMajor>& tmp_C4D_cfg1p_cfg2p_cfg1n_cfg2n) mutable {
        Eigen::Map<Eigen::MatrixXcd> tmp_C2D_cfg1pcfg2p_cfg1ncfg2n(tmp_C4D_cfg1p_cfg2p_cfg1n_cfg2n.data(), Ncfg1pcfg2p_I, Ncfg1ncfg2n_I);

        tmp_C2D_cfg1pcfg2p_cfg1ncfg2n.setZero();
        for (int mu_I = 0; mu_I <= twoLambda_I; ++mu_I) {
            // Oμ = Qνμ† Qπμ at the same Euler node.
            const Eigen::Map<const Eigen::MatrixXcd> Q_C2D_spn_spn(Q_C3D_spn_spn_mu.data() + static_cast<Eigen::Index>(mu_I) * neutron.Nsp_I * neutron.Nsp_I, neutron.Nsp_I, neutron.Nsp_I);
            const Eigen::Map<const Eigen::MatrixXcd> Qp_C2D_spp_spp(Q_C3D_spp_spp_mu.data() + static_cast<Eigen::Index>(mu_I) * proton.Nsp_I * proton.Nsp_I, proton.Nsp_I, proton.Nsp_I);
            Qdag_C2D_spn_spn = Q_C2D_spn_spn.adjoint();
            Q_C2D_spp_spp = Qp_C2D_spp_spp;
            const auto& onebody_C2D_cfg1n_cfg2n = neutron.calc_one_body_pnp(Qdag_C2D_spn_spn);
            const auto& onebody_C2D_cfg1p_cfg2p = proton.calc_one_body_pnp(Q_C2D_spp_spp);

            // Both kernels include overlaps and gauge normalizations.
            const Eigen::Map<const Eigen::VectorXcd> onebody_C1D_cfg1ncfg2n(onebody_C2D_cfg1n_cfg2n.data(), Ncfg1ncfg2n_I);
            const Eigen::Map<const Eigen::VectorXcd> onebody_C1D_cfg1pcfg2p(onebody_C2D_cfg1p_cfg2p.data(), Ncfg1pcfg2p_I);
            tmp_C2D_cfg1pcfg2p_cfg1ncfg2n.noalias() += onebody_C1D_cfg1pcfg2p * onebody_C1D_cfg1ncfg2n.transpose();
        }
    });
}
