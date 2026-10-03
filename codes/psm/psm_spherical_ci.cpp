/**
 * @file    psm_spherical_ci.cpp
 * @author  cailea
 * @date    2026-09-28
 * @brief   Projection kernels and PSM configuration mixing.
 */

#include "psm_spherical.hpp"
#include "spherical_rotation.hpp"
#include <Eigen/Eigenvalues>
#include <cassert>
#include <cmath>
#include <utility>

/**
 * @brief Assemble neutron-proton projection workspaces.
 * @math P = Pₙ^(N−Ncore) Pₚ^(Z−Zcore) Pᴵ.
 * @output Generated configurations and initialized projection_nucleus.
 */
void PSMSpherical::build_projection(int Nbeta_I, int Nphin_I, int Nphip_I) {
    assert(Un_F2D_2spn_2qpn.rows() > 0 && Up_F2D_2spp_2qpp.rows() > 0);
    // (Eqp,2K,cutoffs) → neutron/proton configurations.
    psm_setting.configuration_neutron.build_config(Eqpn_F1D_2qpn, TwoKn_I1D_2qpn);
    psm_setting.configuration_proton.build_config(Eqpp_F1D_2qpp, TwoKp_I1D_2qpp);
    assert(!psm_setting.configuration_neutron.config_I2D_cfg_cqp.empty() && !psm_setting.configuration_proton.config_I2D_cfg_cqp.empty());

    // 2K_cfg = Σqp 2K_qp; Φ₁ = Φ₂.
    TwoKn_I1D_cfgn.resize(psm_setting.configuration_neutron.config_I2D_cfg_cqp.size());
    TwoKn_I1D_cfgn.setZero();
    TwoKp_I1D_cfgp.resize(psm_setting.configuration_proton.config_I2D_cfg_cqp.size());
    TwoKp_I1D_cfgp.setZero();
    for (Eigen::Index cfgn_I = 0; cfgn_I < TwoKn_I1D_cfgn.size(); ++cfgn_I) {
        for (const int qp_I : psm_setting.configuration_neutron.config_I2D_cfg_cqp[cfgn_I]) { TwoKn_I1D_cfgn(cfgn_I) += TwoKn_I1D_2qpn(qp_I); }
    }
    for (Eigen::Index cfgp_I = 0; cfgp_I < TwoKp_I1D_cfgp.size(); ++cfgp_I) {
        for (const int qp_I : psm_setting.configuration_proton.config_I2D_cfg_cqp[cfgp_I]) { TwoKp_I1D_cfgp(cfgp_I) += TwoKp_I1D_2qpp(qp_I); }
    }

    // Axial projection: Nα = Nγ = 1.
    SphericalRotation rotation_neutron(psm_setting.sphericalsetting_neutron.labels_S1D_sp, 1, Nbeta_I, 1);
    SphericalRotation rotation_proton(psm_setting.sphericalsetting_proton.labels_S1D_sp, 1, Nbeta_I, 1);
    rotation_neutron.build();
    rotation_proton.build();

    // (N−Ncore,Z−Zcore,Cₙ,Cₚ) → nuclear projection.
    projection_nucleus = HFBProjectionNucleus(
        HFBProjection(static_cast<int>(Un_F2D_2spn_2qpn.rows()), psm_setting.configuration_neutron.config_I2D_cfg_cqp, psm_setting.configuration_neutron.config_I2D_cfg_cqp, psm_setting.N_I - psm_setting.Ncore_I, Nphin_I, 1, Nbeta_I, 1),
        HFBProjection(static_cast<int>(Up_F2D_2spp_2qpp.rows()), psm_setting.configuration_proton.config_I2D_cfg_cqp, psm_setting.configuration_proton.config_I2D_cfg_cqp, psm_setting.Z_I - psm_setting.Zcore_I, Nphip_I, 1, Nbeta_I, 1));

    // Φ₁ = Φ₂: U₁ = U₂, V₁ = V₂.
    projection_nucleus.projection_neutron.update_UV(Un_F2D_2spn_2qpn, Vn_F2D_2spn_2qpn, Un_F2D_2spn_2qpn, Vn_F2D_2spn_2qpn);
    projection_nucleus.projection_proton.update_UV(Up_F2D_2spp_2qpp, Vp_F2D_2spp_2qpp, Up_F2D_2spp_2qpp, Vp_F2D_2spp_2qpp);

    // {Rz(α),Ry(β),Rz(γ),Ω,wΩ} → neutron projection.
    projection_nucleus.projection_neutron.update_alpha(std::move(rotation_neutron.Rz_C3D_2sp_2sp_alpha), std::move(rotation_neutron.alpha_F1D_alpha), std::move(rotation_neutron.weight_F1D_alpha));
    projection_nucleus.projection_neutron.update_beta(std::move(rotation_neutron.Ry_C3D_2sp_2sp_beta), std::move(rotation_neutron.beta_F1D_beta), std::move(rotation_neutron.weight_F1D_beta));
    projection_nucleus.projection_neutron.update_gamma(std::move(rotation_neutron.Rz_C3D_2sp_2sp_gamma), std::move(rotation_neutron.gamma_F1D_gamma), std::move(rotation_neutron.weight_F1D_gamma));

    // {Rz(α),Ry(β),Rz(γ),Ω,wΩ} → proton projection.
    projection_nucleus.projection_proton.update_alpha(std::move(rotation_proton.Rz_C3D_2sp_2sp_alpha), std::move(rotation_proton.alpha_F1D_alpha), std::move(rotation_proton.weight_F1D_alpha));
    projection_nucleus.projection_proton.update_beta(std::move(rotation_proton.Ry_C3D_2sp_2sp_beta), std::move(rotation_proton.beta_F1D_beta), std::move(rotation_proton.weight_F1D_beta));
    projection_nucleus.projection_proton.update_gamma(std::move(rotation_proton.Rz_C3D_2sp_2sp_gamma), std::move(rotation_proton.gamma_F1D_gamma), std::move(rotation_proton.weight_F1D_gamma));
}

/**
 * @brief Build CI angle kernels using gauge quadrature.
 * @math G₂ = γG₀; (ρ,Q,G) → N(Ω),H(Ω).
 * @output N_C5D and H_C5D beta caches.
 * @note Requires solved HFBCS and initialized projection.
 */
void PSMSpherical::build_ci(double gamma_F) {
    // Axial cache: α = γ = 0.
    assert(projection_nucleus.projection_neutron.Nalpha_I == 1 && projection_nucleus.projection_neutron.Ngamma_I == 1);
    assert(projection_nucleus.projection_neutron.alpha_F1D_alpha(0) == 0.0 && projection_nucleus.projection_neutron.gamma_F1D_gamma(0) == 0.0);

    psm_setting.set_CI_QQ(rhon_F2D_2spn_2spn, rhop_F2D_2spp_2spp, Qn_F3D_2spn_2spn_mu, Qp_F3D_2spp_2spp_mu);
    psm_setting.set_CI_PP(gamma_F);

    // N_ab(Ω) = ⟨Φ_a|PₙPₚR(Ω)|Φ_b⟩.
    const Eigen::Index Ncfgn_I = projection_nucleus.projection_neutron.Ncfg1_I;
    const Eigen::Index Ncfgp_I = projection_nucleus.projection_proton.Ncfg1_I;
    assert(Ncfgn_I > 0 && Ncfgp_I > 0);
    assert(Ncfgn_I == projection_nucleus.projection_neutron.Ncfg2_I && Ncfgp_I == projection_nucleus.projection_proton.Ncfg2_I);

    // H_ab(Ω) = ⟨Φ_a|HPₙPₚR(Ω)|Φ_b⟩.
    assert(psm_setting.N_I >= psm_setting.Ncore_I && psm_setting.Z_I >= psm_setting.Zcore_I);
    assert(projection_nucleus.projection_neutron.Nsp_I > 0 && projection_nucleus.projection_proton.Nsp_I > 0);
    assert(psm_setting.N_I - psm_setting.Ncore_I == projection_nucleus.projection_neutron.TargetN_I && psm_setting.Z_I - psm_setting.Zcore_I == projection_nucleus.projection_proton.TargetN_I);

    // h₀ ∈ ℝ → complex projection inputs.
    const Eigen::MatrixXcd h0n_C2D_2spn_2spn = h0n_F2D_2spn_2spn.cast<doubleC>();
    const Eigen::MatrixXcd h0p_C2D_2spp_2spp = h0p_F2D_2spp_2spp.cast<doubleC>();

    // Q₂μ = δ_NN′(r/b)²√(4π/5)Y₂μ; P₂μ = Q₂μ P₀.
    assert(std::isfinite(psm_setting.chi2_nn_F) && std::isfinite(psm_setting.chi2_np_F) && std::isfinite(psm_setting.chi2_pp_F));
    assert(std::isfinite(psm_setting.G0_nn_F) && std::isfinite(psm_setting.G0_pp_F) && std::isfinite(psm_setting.G2_nn_F) && std::isfinite(psm_setting.G2_pp_F));

    // Real spherical operators → complex projection inputs.
    const Eigen::Tensor<doubleC, 3, Eigen::ColMajor> Qn_C3D_2spn_2spn_i = Qn_F3D_2spn_2spn_mu.cast<doubleC>();
    const Eigen::Tensor<doubleC, 3, Eigen::ColMajor> Qp_C3D_2spp_2spp_i = Qp_F3D_2spp_2spp_mu.cast<doubleC>();

    // Q = (Q₂₋₂,…,Q₂₂); P = (P₀₀,P₂₋₂,…,P₂₂).
    const Eigen::VectorXd gQn_F1D_i = Eigen::VectorXd::Constant(5, -0.5 * psm_setting.chi2_nn_F);
    const Eigen::VectorXd gQp_F1D_i = Eigen::VectorXd::Constant(5, -0.5 * psm_setting.chi2_pp_F);
    Eigen::VectorXd gPn_F1D_i = Eigen::VectorXd::Constant(6, -psm_setting.G2_nn_F);
    Eigen::VectorXd gPp_F1D_i = Eigen::VectorXd::Constant(6, -psm_setting.G2_pp_F);
    gPn_F1D_i(0) = -psm_setting.G0_nn_F;
    gPp_F1D_i(0) = -psm_setting.G0_pp_F;
    const Eigen::Tensor<doubleC, 3, Eigen::ColMajor> Pn_C3D_2spn_2spn_i = P0n_F3D_2spn_2spn_mu.concatenate(P2n_F3D_2spn_2spn_mu, 2).cast<doubleC>();
    const Eigen::Tensor<doubleC, 3, Eigen::ColMajor> Pp_C3D_2spp_2spp_i = P0p_F3D_2spp_2spp_mu.concatenate(P2p_F3D_2spp_2spp_mu, 2).cast<doubleC>();

    // operator = (N,H¹,H²,Q₋₂,Q₋₂†,…,Q₂,Q₂†).
    Eigen::Tensor<doubleC, 3, Eigen::ColMajor> resultn_C3D_cfgn_cfgn_operator(Ncfgn_I, Ncfgn_I, 13);
    Eigen::Tensor<doubleC, 3, Eigen::ColMajor> resultp_C3D_cfgp_cfgp_operator(Ncfgp_I, Ncfgp_I, 13);
    Eigen::MatrixXcd Qdag_C2D_2spn_2spn(h0n_C2D_2spn_2spn.rows(), h0n_C2D_2spn_2spn.cols());
    Eigen::MatrixXcd Qdag_C2D_2spp_2spp(h0p_C2D_2spp_2spp.rows(), h0p_C2D_2spp_2spp.cols());

    // Nuclear operator = (N,H); one Ω traversal.
    Eigen::Tensor<doubleC, 8, Eigen::ColMajor> result_C8D_cfgn1_cfgn2_cfgp1_cfgp2_operator_gamma_alpha_beta(Ncfgn_I, Ncfgn_I, Ncfgp_I, Ncfgp_I, 2, 1, 1, projection_nucleus.projection_neutron.Nbeta_I);
    projection_nucleus.build_amp(result_C8D_cfgn1_cfgn2_cfgp1_cfgp2_operator_gamma_alpha_beta, [&](int, int, int, auto& resultAtOmega_C5D_cfgn1_cfgn2_cfgp1_cfgp2_operator) {
        // One contraction update per φ supplies all operators.
        projection_nucleus.projection_neutron.integrate_pnp(Eigen::TensorMap<Eigen::Tensor<doubleC, 3, Eigen::ColMajor>>(resultn_C3D_cfgn_cfgn_operator.data(), resultn_C3D_cfgn_cfgn_operator.dimensions()), projection_nucleus.projection_neutron.TargetN_I, [&](int, auto& resultAtPhi_C3D_cfg1_cfg2_operator) {
            Eigen::Map<Eigen::MatrixXcd> overlap_C2D_cfg1_cfg2(resultAtPhi_C3D_cfg1_cfg2_operator.data(), Ncfgn_I, Ncfgn_I);
            Eigen::Map<Eigen::MatrixXcd> onebody_C2D_cfg1_cfg2(resultAtPhi_C3D_cfg1_cfg2_operator.data() + Ncfgn_I * Ncfgn_I, Ncfgn_I, Ncfgn_I);
            Eigen::Map<Eigen::MatrixXcd> twobody_C2D_cfg1_cfg2(resultAtPhi_C3D_cfg1_cfg2_operator.data() + 2 * Ncfgn_I * Ncfgn_I, Ncfgn_I, Ncfgn_I);
            projection_nucleus.projection_neutron.hfb_pfaffian.calc_overlap(overlap_C2D_cfg1_cfg2);
            projection_nucleus.projection_neutron.hfb_pfaffian.calc_one_body(onebody_C2D_cfg1_cfg2, h0n_C2D_2spn_2spn);
            projection_nucleus.projection_neutron.hfb_pfaffian.calc_two_body(twobody_C2D_cfg1_cfg2, gQn_F1D_i, gPn_F1D_i, Qn_C3D_2spn_2spn_i, Pn_C3D_2spn_2spn_i);

            // operator 3+2i → Qᵢ; 4+2i → Qᵢ†.
            for (Eigen::Index iQ_I = 0; iQ_I < 5; ++iQ_I) {
                const Eigen::Index offset_I = iQ_I * h0n_C2D_2spn_2spn.size();
                const Eigen::Map<const Eigen::MatrixXcd> Q_C2D_sp_sp(Qn_C3D_2spn_2spn_i.data() + offset_I, h0n_C2D_2spn_2spn.rows(), h0n_C2D_2spn_2spn.cols());
                Qdag_C2D_2spn_2spn = Q_C2D_sp_sp.adjoint();
                Eigen::Map<Eigen::MatrixXcd> Q_C2D_cfg1_cfg2(resultAtPhi_C3D_cfg1_cfg2_operator.data() + (3 + 2 * iQ_I) * Ncfgn_I * Ncfgn_I, Ncfgn_I, Ncfgn_I);
                Eigen::Map<Eigen::MatrixXcd> Qdag_C2D_cfg1_cfg2(resultAtPhi_C3D_cfg1_cfg2_operator.data() + (4 + 2 * iQ_I) * Ncfgn_I * Ncfgn_I, Ncfgn_I, Ncfgn_I);
                projection_nucleus.projection_neutron.hfb_pfaffian.calc_one_body(Q_C2D_cfg1_cfg2, Q_C2D_sp_sp);
                projection_nucleus.projection_neutron.hfb_pfaffian.calc_one_body(Qdag_C2D_cfg1_cfg2, Qdag_C2D_2spn_2spn);
            }
        });

        // One contraction update per φ supplies all operators.
        projection_nucleus.projection_proton.integrate_pnp(Eigen::TensorMap<Eigen::Tensor<doubleC, 3, Eigen::ColMajor>>(resultp_C3D_cfgp_cfgp_operator.data(), resultp_C3D_cfgp_cfgp_operator.dimensions()), projection_nucleus.projection_proton.TargetN_I, [&](int, auto& resultAtPhi_C3D_cfg1_cfg2_operator) {
            Eigen::Map<Eigen::MatrixXcd> overlap_C2D_cfg1_cfg2(resultAtPhi_C3D_cfg1_cfg2_operator.data(), Ncfgp_I, Ncfgp_I);
            Eigen::Map<Eigen::MatrixXcd> onebody_C2D_cfg1_cfg2(resultAtPhi_C3D_cfg1_cfg2_operator.data() + Ncfgp_I * Ncfgp_I, Ncfgp_I, Ncfgp_I);
            Eigen::Map<Eigen::MatrixXcd> twobody_C2D_cfg1_cfg2(resultAtPhi_C3D_cfg1_cfg2_operator.data() + 2 * Ncfgp_I * Ncfgp_I, Ncfgp_I, Ncfgp_I);
            projection_nucleus.projection_proton.hfb_pfaffian.calc_overlap(overlap_C2D_cfg1_cfg2);
            projection_nucleus.projection_proton.hfb_pfaffian.calc_one_body(onebody_C2D_cfg1_cfg2, h0p_C2D_2spp_2spp);
            projection_nucleus.projection_proton.hfb_pfaffian.calc_two_body(twobody_C2D_cfg1_cfg2, gQp_F1D_i, gPp_F1D_i, Qp_C3D_2spp_2spp_i, Pp_C3D_2spp_2spp_i);

            // operator 3+2i → Qᵢ; 4+2i → Qᵢ†.
            for (Eigen::Index iQ_I = 0; iQ_I < 5; ++iQ_I) {
                const Eigen::Index offset_I = iQ_I * h0p_C2D_2spp_2spp.size();
                const Eigen::Map<const Eigen::MatrixXcd> Q_C2D_sp_sp(Qp_C3D_2spp_2spp_i.data() + offset_I, h0p_C2D_2spp_2spp.rows(), h0p_C2D_2spp_2spp.cols());
                Qdag_C2D_2spp_2spp = Q_C2D_sp_sp.adjoint();
                Eigen::Map<Eigen::MatrixXcd> Q_C2D_cfg1_cfg2(resultAtPhi_C3D_cfg1_cfg2_operator.data() + (3 + 2 * iQ_I) * Ncfgp_I * Ncfgp_I, Ncfgp_I, Ncfgp_I);
                Eigen::Map<Eigen::MatrixXcd> Qdag_C2D_cfg1_cfg2(resultAtPhi_C3D_cfg1_cfg2_operator.data() + (4 + 2 * iQ_I) * Ncfgp_I * Ncfgp_I, Ncfgp_I, Ncfgp_I);
                projection_nucleus.projection_proton.hfb_pfaffian.calc_one_body(Q_C2D_cfg1_cfg2, Q_C2D_sp_sp);
                projection_nucleus.projection_proton.hfb_pfaffian.calc_one_body(Qdag_C2D_cfg1_cfg2, Qdag_C2D_2spp_2spp);
            }
        });

        // Columns: (N,H¹,H²,Qᵢ,Qᵢ†); rows: configuration pairs.
        const Eigen::Map<const Eigen::MatrixXcd> resultn_C2D_cfgncfgn_operator(resultn_C3D_cfgn_cfgn_operator.data(), Ncfgn_I * Ncfgn_I, 13);
        const Eigen::Map<const Eigen::MatrixXcd> resultp_C2D_cfgpcfgp_operator(resultp_C3D_cfgp_cfgp_operator.data(), Ncfgp_I * Ncfgp_I, 13);
        Eigen::Map<Eigen::MatrixXcd> N_C2D_cfgncfgn_cfgpcfgp(resultAtOmega_C5D_cfgn1_cfgn2_cfgp1_cfgp2_operator.data(), Ncfgn_I * Ncfgn_I, Ncfgp_I * Ncfgp_I);
        Eigen::Map<Eigen::MatrixXcd> H_C2D_cfgncfgn_cfgpcfgp(resultAtOmega_C5D_cfgn1_cfgn2_cfgp1_cfgp2_operator.data() + Ncfgn_I * Ncfgn_I * Ncfgp_I * Ncfgp_I, Ncfgn_I * Ncfgn_I, Ncfgp_I * Ncfgp_I);

        // N = Nν Nπᵀ; H = (Hν¹+Hν²)Nπᵀ + Nν(Hπ¹+Hπ²)ᵀ.
        N_C2D_cfgncfgn_cfgpcfgp.noalias() = resultn_C2D_cfgncfgn_operator.col(0) * resultp_C2D_cfgpcfgp_operator.col(0).transpose();
        H_C2D_cfgncfgn_cfgpcfgp.noalias() = (resultn_C2D_cfgncfgn_operator.col(1) + resultn_C2D_cfgncfgn_operator.col(2)) * resultp_C2D_cfgpcfgp_operator.col(0).transpose();
        H_C2D_cfgncfgn_cfgpcfgp.noalias() += resultn_C2D_cfgncfgn_operator.col(0) * (resultp_C2D_cfgpcfgp_operator.col(1) + resultp_C2D_cfgpcfgp_operator.col(2)).transpose();

        // Hνπ = −χnp/2 Σᵢ [N(Qᵢ†)P(Qᵢ)ᵀ + N(Qᵢ)P(Qᵢ†)ᵀ].
        for (Eigen::Index iQ_I = 0; iQ_I < 5; ++iQ_I) {
            H_C2D_cfgncfgn_cfgpcfgp.noalias() -= 0.5 * psm_setting.chi2_np_F * resultn_C2D_cfgncfgn_operator.col(4 + 2 * iQ_I) * resultp_C2D_cfgpcfgp_operator.col(3 + 2 * iQ_I).transpose();
            H_C2D_cfgncfgn_cfgpcfgp.noalias() -= 0.5 * psm_setting.chi2_np_F * resultn_C2D_cfgncfgn_operator.col(3 + 2 * iQ_I) * resultp_C2D_cfgpcfgp_operator.col(4 + 2 * iQ_I).transpose();
        }
    });
    N_C5D_cfgn_cfgn_cfgp_cfgp_beta = result_C8D_cfgn1_cfgn2_cfgp1_cfgp2_operator_gamma_alpha_beta.chip(0, 6).chip(0, 5).chip(0, 4);
    H_C5D_cfgn_cfgn_cfgp_cfgp_beta = result_C8D_cfgn1_cfgn2_cfgp1_cfgp2_operator_gamma_alpha_beta.chip(0, 6).chip(0, 5).chip(1, 4);
}

/**
 * @brief Solve configuration mixing using canonical norm orthogonalization.
 * @math N = U n U†; X = Uᵣ nᵣ^(-1/2); (X†HX)Xci = Xci E.
 * @math f = X Xci; Xci†Xci = 1; r = eigenN, ν = eigenH.
 * @output Eci_F1D_eigenH and Xci_C2D_eigenN_eigenH.
 * @note Retains nᵢ > 1e-10 nmax; energies ascend.
 */
void PSMSpherical::solve_ci(int TargetTwoI_I) {
    assert(TargetTwoI_I >= 0);
    assert(N_C5D_cfgn_cfgn_cfgp_cfgp_beta.size() > 0 && H_C5D_cfgn_cfgn_cfgp_cfgp_beta.size() > 0);

    Eigen::Tensor<doubleC, 5, Eigen::ColMajor> N_C5D_cfgn_cfgn_cfgp_cfgp_operator(N_C5D_cfgn_cfgn_cfgp_cfgp_beta.dimension(0), N_C5D_cfgn_cfgn_cfgp_cfgp_beta.dimension(1), N_C5D_cfgn_cfgn_cfgp_cfgp_beta.dimension(2), N_C5D_cfgn_cfgn_cfgp_cfgp_beta.dimension(3), 1);
    projection_nucleus.integrate_amp_with_K(N_C5D_cfgn_cfgn_cfgp_cfgp_operator, TargetTwoI_I, TwoKn_I1D_cfgn, TwoKn_I1D_cfgn, TwoKp_I1D_cfgp, TwoKp_I1D_cfgp, N_C5D_cfgn_cfgn_cfgp_cfgp_beta.reshape(Eigen::array<Eigen::Index, 6>{N_C5D_cfgn_cfgn_cfgp_cfgp_beta.dimension(0), N_C5D_cfgn_cfgn_cfgp_cfgp_beta.dimension(1), N_C5D_cfgn_cfgn_cfgp_cfgp_beta.dimension(2), N_C5D_cfgn_cfgn_cfgp_cfgp_beta.dimension(3), 1, N_C5D_cfgn_cfgn_cfgp_cfgp_beta.dimension(4)}));
    const Eigen::Tensor<doubleC, 4, Eigen::ColMajor> N_C4D_cfgn_cfgp_cfgn_cfgp = N_C5D_cfgn_cfgn_cfgp_cfgp_operator.chip(0, 4).shuffle(Eigen::array<int, 4>{0, 2, 1, 3});
    Eigen::Tensor<doubleC, 5, Eigen::ColMajor> H_C5D_cfgn_cfgn_cfgp_cfgp_operator(H_C5D_cfgn_cfgn_cfgp_cfgp_beta.dimension(0), H_C5D_cfgn_cfgn_cfgp_cfgp_beta.dimension(1), H_C5D_cfgn_cfgn_cfgp_cfgp_beta.dimension(2), H_C5D_cfgn_cfgn_cfgp_cfgp_beta.dimension(3), 1);
    projection_nucleus.integrate_amp_with_K(H_C5D_cfgn_cfgn_cfgp_cfgp_operator, TargetTwoI_I, TwoKn_I1D_cfgn, TwoKn_I1D_cfgn, TwoKp_I1D_cfgp, TwoKp_I1D_cfgp, H_C5D_cfgn_cfgn_cfgp_cfgp_beta.reshape(Eigen::array<Eigen::Index, 6>{H_C5D_cfgn_cfgn_cfgp_cfgp_beta.dimension(0), H_C5D_cfgn_cfgn_cfgp_cfgp_beta.dimension(1), H_C5D_cfgn_cfgn_cfgp_cfgp_beta.dimension(2), H_C5D_cfgn_cfgn_cfgp_cfgp_beta.dimension(3), 1, H_C5D_cfgn_cfgn_cfgp_cfgp_beta.dimension(4)}));
    const Eigen::Tensor<doubleC, 4, Eigen::ColMajor> H_C4D_cfgn_cfgp_cfgn_cfgp = H_C5D_cfgn_cfgn_cfgp_cfgp_operator.chip(0, 4).shuffle(Eigen::array<int, 4>{0, 2, 1, 3});
    // (cfgn,cfgp,cfgn′,cfgp′) → (cfgncfgp,cfgncfgp′); zero-copy.
    const Eigen::Index Ncfgncfgp_I = N_C4D_cfgn_cfgp_cfgn_cfgp.dimension(0) * N_C4D_cfgn_cfgp_cfgn_cfgp.dimension(1);
    const Eigen::Map<const Eigen::MatrixXcd> N_C2D_cfgncfgp_cfgncfgp(N_C4D_cfgn_cfgp_cfgn_cfgp.data(), Ncfgncfgp_I, Ncfgncfgp_I);
    const Eigen::Map<const Eigen::MatrixXcd> H_C2D_cfgncfgp_cfgncfgp(H_C4D_cfgn_cfgp_cfgn_cfgp.data(), Ncfgncfgp_I, Ncfgncfgp_I);

    // N → (N+N†)/2; H → (H+H†)/2.
    const Eigen::MatrixXcd Nhermitian_C2D_cfgncfgp_cfgncfgp = 0.5 * (N_C2D_cfgncfgp_cfgncfgp + N_C2D_cfgncfgp_cfgncfgp.adjoint());
    const Eigen::MatrixXcd Hhermitian_C2D_cfgncfgp_cfgncfgp = 0.5 * (H_C2D_cfgncfgp_cfgncfgp + H_C2D_cfgncfgp_cfgncfgp.adjoint());
    const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> N_solver(Nhermitian_C2D_cfgncfgp_cfgncfgp);
    assert(N_solver.info() == Eigen::Success);
    const double normCut_F = 1.0e-10 * N_solver.eigenvalues().maxCoeff();
    assert(normCut_F > 0.0 && N_solver.eigenvalues().minCoeff() >= -normCut_F);
    const Eigen::Index Nstate_I = (N_solver.eigenvalues().array() > normCut_F).count();
    assert(Nstate_I > 0);

    // X = Uᵣ nᵣ^(-1/2); ascending norm eigenvalues.
    const Eigen::MatrixXcd U_C2D_cfgncfgp_eigenN = N_solver.eigenvectors().rightCols(Nstate_I);
    const Eigen::VectorXd n_F1D_eigenN = N_solver.eigenvalues().tail(Nstate_I);
    Eigen::MatrixXcd X_C2D_cfgncfgp_eigenN = U_C2D_cfgncfgp_eigenN * n_F1D_eigenN.array().sqrt().inverse().matrix().asDiagonal();
    Eigen::MatrixXcd H_C2D_eigenN_eigenN = X_C2D_cfgncfgp_eigenN.adjoint() * Hhermitian_C2D_cfgncfgp_cfgncfgp * X_C2D_cfgncfgp_eigenN;
    H_C2D_eigenN_eigenN = (0.5 * (H_C2D_eigenN_eigenN + H_C2D_eigenN_eigenN.adjoint())).eval();
    const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> H_solver(H_C2D_eigenN_eigenN);
    assert(H_solver.info() == Eigen::Success);

    // (X†HX)Xci = Xci E; Xci†Xci = 1.
    Eci_F1D_eigenH = H_solver.eigenvalues();
    Xci_C2D_eigenN_eigenH = H_solver.eigenvectors();
}
