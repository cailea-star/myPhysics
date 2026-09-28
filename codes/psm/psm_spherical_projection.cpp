/**
 * @file    psm_spherical_projection.cpp
 * @author  cailea
 * @date    2026-09-28
 * @brief   Configuration construction from Nilsson-BCS states.
 */

#include "psm_spherical.hpp"
#include "spherical_rotation.hpp"
#include "spherical_basis.hpp"
#include <gsl/gsl_sf_coupling.h>
#include <algorithm>
#include <cmath>
#include <utility>

/**
 * @brief Set neutron cutoffs and enumerate configurations.
 * @math Nqp = size(Eqp); (Ncqp,ECut,NCut,TwoKCut) → configurations.
 * @output Updated neutron cutoffs and configuration list.
 */
void PSMSpherical::build_config_neutron(const Eigen::VectorXi& Ncqp_I1D_Ncqp, const Eigen::VectorXd& ECut_F1D_Ncqp, const Eigen::VectorXi& NCut_I1D_Ncqp, const Eigen::VectorXi& TwoKCut_I1D_Ncqp) {
    const int Nqp_I = static_cast<int>(Eqpn_F1D_2qpn.size());
    assert(Nqp_I > 0);
    configuration_neutron = HFBConfiguration(Nqp_I, Ncqp_I1D_Ncqp, ECut_F1D_Ncqp, NCut_I1D_Ncqp, TwoKCut_I1D_Ncqp);
    configuration_neutron.build_config(Eqpn_F1D_2qpn, TwoKn_I1D_2qpn);
}

/**
 * @brief Set proton cutoffs and enumerate configurations.
 * @math Nqp = size(Eqp); (Ncqp,ECut,NCut,TwoKCut) → configurations.
 * @output Updated proton cutoffs and configuration list.
 */
void PSMSpherical::build_config_proton(const Eigen::VectorXi& Ncqp_I1D_Ncqp, const Eigen::VectorXd& ECut_F1D_Ncqp, const Eigen::VectorXi& NCut_I1D_Ncqp, const Eigen::VectorXi& TwoKCut_I1D_Ncqp) {
    const int Nqp_I = static_cast<int>(Eqpp_F1D_2qpp.size());
    assert(Nqp_I > 0);
    configuration_proton = HFBConfiguration(Nqp_I, Ncqp_I1D_Ncqp, ECut_F1D_Ncqp, NCut_I1D_Ncqp, TwoKCut_I1D_Ncqp);
    configuration_proton.build_config(Eqpp_F1D_2qpp, TwoKp_I1D_2qpp);
}

/**
 * @brief Assemble neutron-proton projection workspaces.
 * @math P = Pₙ^(N−Ncore) Pₚ^(Z−Zcore) Pᴵ.
 * @output Initialized projection_nucleus with shared left/right configurations.
 */
void PSMSpherical::build_projection(int TargetN_I, int TargetZ_I, int TargetTwoI_I, int Nalpha_I, int Nbeta_I, int Ngamma_I, int Nphin_I, int Nphip_I) {
    assert(Un_F2D_2spn_2qpn.rows() > 0 && Up_F2D_2spp_2qpp.rows() > 0);
    assert(!configuration_neutron.config_I2D_cfg_cqp.empty() && !configuration_proton.config_I2D_cfg_cqp.empty());

    // {labels,Nα,Nβ,Nγ} → rotation caches.
    SphericalRotation rotation_neutron(sphericalsetting_neutron.labels_S1D_sp, Nalpha_I, Nbeta_I, Ngamma_I);
    SphericalRotation rotation_proton(sphericalsetting_proton.labels_S1D_sp, Nalpha_I, Nbeta_I, Ngamma_I);
    rotation_neutron.build();
    rotation_proton.build();

    // (N−Ncore,Z−Zcore,2I,Cₙ,Cₚ) → nuclear projection.
    projection_nucleus = HFBProjectionNucleus(
        HFBProjection(TargetN_I - Ncore_I, TargetTwoI_I, static_cast<int>(Un_F2D_2spn_2qpn.rows()), configuration_neutron.config_I2D_cfg_cqp, configuration_neutron.config_I2D_cfg_cqp, Nphin_I, Nalpha_I, Nbeta_I, Ngamma_I),
        HFBProjection(TargetZ_I - Zcore_I, TargetTwoI_I, static_cast<int>(Up_F2D_2spp_2qpp.rows()), configuration_proton.config_I2D_cfg_cqp, configuration_proton.config_I2D_cfg_cqp, Nphip_I, Nalpha_I, Nbeta_I, Ngamma_I));

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
 * @brief Build the projected configuration overlap tensor.
 * @math (p₁,p₂,n₁,n₂,K₁,K₂) → (p₁,n₁,K₁,p₂,n₂,K₂).
 * @output N_C6D_cfgp_cfgn_K_cfgp_cfgn_K; K indexes −I,…,I.
 */
void PSMSpherical::build_norm() {
    const Eigen::Index Ncfgn_I = projection_nucleus.projection_neutron.Ncfg1_I;
    const Eigen::Index Ncfgp_I = projection_nucleus.projection_proton.Ncfg1_I;
    assert(Ncfgn_I > 0 && Ncfgp_I > 0);
    assert(Ncfgn_I == projection_nucleus.projection_neutron.Ncfg2_I && Ncfgp_I == projection_nucleus.projection_proton.Ncfg2_I);

    // (p₁,p₂,n₁,n₂,K₁,K₂) → (p₁,n₁,K₁,p₂,n₂,K₂).
    const auto& N_C6D_cfg1p_cfg2p_cfg1n_cfg2n_K1_K2 = projection_nucleus.calc_overlap();
    N_C6D_cfgp_cfgn_K_cfgp_cfgn_K = N_C6D_cfg1p_cfg2p_cfg1n_cfg2n_K1_K2.shuffle(Eigen::array<int, 6>{0, 2, 4, 1, 3, 5});
}

/**
 * @brief Construct spherical single-particle fields and project them.
 * @math h₀ = ℏω₀(N+3/2) − κℏωc[2l·s+μ(l²−N(N+3)/2)].
 * @output Honebody_C6D_cfgp_cfgn_K_cfgp_cfgn_K.
 * @note Total N,Z determine A,κ,μ,ωc; b determines ω₀.
 */
void PSMSpherical::build_onebody(int TargetN_I, int TargetZ_I) {
    assert(TargetN_I >= Ncore_I && TargetZ_I >= Zcore_I);
    assert(projection_nucleus.projection_neutron.Nsp_I > 0 && projection_nucleus.projection_proton.Nsp_I > 0);
    assert(TargetN_I - Ncore_I == projection_nucleus.projection_neutron.TargetN_I && TargetZ_I - Zcore_I == projection_nucleus.projection_proton.TargetN_I);

    // A = N+Z; ℏωc,τ = 41.4678 A^(-1/3)(2Nτ/A)^(1/3).
    const double A_F = static_cast<double>(TargetN_I) + TargetZ_I;
    assert(A_F > 0.0);
    const double hbarOmegacn_F = 41.4678 * std::cbrt(2.0 * TargetN_I / A_F) / std::cbrt(A_F);
    const double hbarOmegacp_F = 41.4678 * std::cbrt(2.0 * TargetZ_I / A_F) / std::cbrt(A_F);
    Eigen::MatrixXcd h0n_C2D_2spn_2spn{};
    Eigen::MatrixXcd h0p_C2D_2spp_2spp{};

    // ⟨Nljm|h₀|N'l'j'm'⟩ = δ_NN'δ_ll'δ_jj'δ_mm' ε_Nlj.
    const auto build_h0 = [&](const SphericalSetting& setting, double hbarOmegac_F, double Akappa_F, double Bkappa_F, double Amu_F, double Bmu_F, const std::array<double, 9>& Fkappa_F1D_Nshell, const std::array<double, 9>& Fmu_F1D_Nshell, Eigen::MatrixXcd& h0_C2D_2sp_2sp) {
        const Eigen::Index Nsp_I = static_cast<Eigen::Index>(setting.labels_S1D_sp.size());
        const double hbarOmega0_F = 41.4678 / (setting.b_F * setting.b_F);
        h0_C2D_2sp_2sp.resize(2 * Nsp_I, 2 * Nsp_I);
        h0_C2D_2sp_2sp.setZero();
        for (Eigen::Index sp_I = 0; sp_I < Nsp_I; ++sp_I) {
            const auto& label = setting.labels_S1D_sp[sp_I];
            const auto shell_iter = std::find(HFBCSData::Nshell_I1D_Nshell.begin(), HFBCSData::Nshell_I1D_Nshell.end(), label.N_I);
            assert(shell_iter != HFBCSData::Nshell_I1D_Nshell.end());
            const std::size_t shell_I = static_cast<std::size_t>(shell_iter - HFBCSData::Nshell_I1D_Nshell.begin());

            // κ,μ = (Aκ,μ − A Bκ,μ)Fκ,μ(N).
            const double kappa_F = (Akappa_F - A_F * Bkappa_F) * Fkappa_F1D_Nshell[shell_I];
            const double mu_F = (Amu_F - A_F * Bmu_F) * Fmu_F1D_Nshell[shell_I];
            const double j_F = 0.5 * label.twoj_I;
            const double l2_F = label.l_I * (label.l_I + 1.0);
            const double twoLS_F = j_F * (j_F + 1.0) - l2_F - 0.75;
            const double l2Average_F = 0.5 * label.N_I * (label.N_I + 3.0);

            // ε_Nlj = ℏω₀(N+3/2) − κℏωc[2l·s+μ(l²−⟨l²⟩_N)].
            const double spherical_F = hbarOmega0_F * (label.N_I + 1.5) - kappa_F * hbarOmegac_F * (twoLS_F + mu_F * (l2_F - l2Average_F));
            h0_C2D_2sp_2sp(sp_I, sp_I) = spherical_F;
            h0_C2D_2sp_2sp(sp_I + Nsp_I, sp_I + Nsp_I) = spherical_F;
        }
    };

    // (+m,−m): h₀ₙ,h₀ₚ share the U,V basis ordering.
    build_h0(sphericalsetting_neutron, hbarOmegacn_F, HFBCSData::Akappa_n_F, HFBCSData::Bkappa_n_F, HFBCSData::Amu_n_F, HFBCSData::Bmu_n_F, HFBCSData::Fkappa_n_F1D_Nshell, HFBCSData::Fmu_n_F1D_Nshell, h0n_C2D_2spn_2spn);
    build_h0(sphericalsetting_proton, hbarOmegacp_F, HFBCSData::Akappa_p_F, HFBCSData::Bkappa_p_F, HFBCSData::Amu_p_F, HFBCSData::Bmu_p_F, HFBCSData::Fkappa_p_F1D_Nshell, HFBCSData::Fmu_p_F1D_Nshell, h0p_C2D_2spp_2spp);

    // Store neutron output before the shared projection buffer changes.
    Honebody_C6D_cfgp_cfgn_K_cfgp_cfgn_K = projection_nucleus.calc_one_body_nn(h0n_C2D_2spn_2spn).shuffle(Eigen::array<int, 6>{0, 2, 4, 1, 3, 5});
    Honebody_C6D_cfgp_cfgn_K_cfgp_cfgn_K += projection_nucleus.calc_one_body_pp(h0p_C2D_2spp_2spp).shuffle(Eigen::array<int, 6>{0, 2, 4, 1, 3, 5});
}

/**
 * @brief Build quadrupole operators using radial quadrature and GSL 3j.
 * @math H² = −χnn:Qn†Qn:/2 − χnp Qn†Qp − χpp:Qp†Qp:/2 − Στλ Gλ,τ Pλ,τ†Pλ,τ.
 * @output Htwobody_C6D_cfgp_cfgn_K_cfgp_cfgn_K.
 * @note Q₂μ = δ_NN′(r/b)²√(4π/5)Y₂μ; P₂μ = Q₂μ P₀.
 * @note Independent strengths; λ = 0,2; no BCS-strength matching.
 */
void PSMSpherical::build_twobody(double chi2_nn_F, double chi2_np_F, double chi2_pp_F, double G0_nn_F, double G0_pp_F, double G2_nn_F, double G2_pp_F) {
    assert(projection_nucleus.projection_neutron.Nsp_I > 0 && projection_nucleus.projection_proton.Nsp_I > 0);
    assert(std::isfinite(chi2_nn_F) && std::isfinite(chi2_np_F) && std::isfinite(chi2_pp_F));
    assert(std::isfinite(G0_nn_F) && std::isfinite(G0_pp_F) && std::isfinite(G2_nn_F) && std::isfinite(G2_pp_F));

    // Q₂μ,P₂μ: μ = −2,…,2; P₀: μ = 0.
    Eigen::Tensor<double, 3, Eigen::ColMajor> Q2n_F3D_2spn_2spn_mu{};
    Eigen::Tensor<double, 3, Eigen::ColMajor> Q2p_F3D_2spp_2spp_mu{};
    Eigen::Tensor<double, 3, Eigen::ColMajor> P0n_F3D_2spn_2spn_mu{};
    Eigen::Tensor<double, 3, Eigen::ColMajor> P0p_F3D_2spp_2spp_mu{};
    Eigen::Tensor<double, 3, Eigen::ColMajor> P2n_F3D_2spn_2spn_mu{};
    Eigen::Tensor<double, 3, Eigen::ColMajor> P2p_F3D_2spp_2spp_mu{};

    // (+m,−m) → Q₂μ,P₀,P₂μ in the U,V basis.
    const auto build_QP = [](const SphericalSetting& setting, Eigen::Tensor<double, 3, Eigen::ColMajor>& Q2_F3D_2sp_2sp_mu, Eigen::Tensor<double, 3, Eigen::ColMajor>& P0_F3D_2sp_2sp_mu, Eigen::Tensor<double, 3, Eigen::ColMajor>& P2_F3D_2sp_2sp_mu) {
        const Eigen::Index Nsp_I = static_cast<Eigen::Index>(setting.labels_S1D_sp.size());
        std::vector<SphericalSPLabel> labels_S1D_2sp = setting.labels_S1D_sp;
        labels_S1D_2sp.reserve(2 * Nsp_I);
        for (const auto& label : setting.labels_S1D_sp) {
            labels_S1D_2sp.emplace_back(label.N_I, label.twoj_I, -label.twom_I);
        }

        // ∫r²dr R₁R₂(r/b)²; retain ΔN = 0 as PSM.
        const SphericalLaguerreBasis radial_basis(setting.b_F, setting.Nr_I, setting.labels_S1D_sp);
        const Eigen::VectorXd r2Weights_F1D_r = radial_basis.w_F1D_r.array() * radial_basis.r_F1D_r.array().square() / (setting.b_F * setting.b_F);
        Eigen::MatrixXd r2_F2D_sp_sp{};
        r2_F2D_sp_sp.resize(Nsp_I, Nsp_I);
        for (Eigen::Index sp2_I = 0; sp2_I < Nsp_I; ++sp2_I) {
            for (Eigen::Index sp1_I = 0; sp1_I < Nsp_I; ++sp1_I) {
                const double deltaN_F = static_cast<double>(setting.labels_S1D_sp[sp1_I].N_I == setting.labels_S1D_sp[sp2_I].N_I);
                r2_F2D_sp_sp(sp1_I, sp2_I) = deltaN_F * (radial_basis.phi_F2D_sp_r.row(sp1_I).array() * radial_basis.phi_F2D_sp_r.row(sp2_I).array() * r2Weights_F1D_r.transpose().array()).sum();
            }
        }

        // ⟨l₁m₁|√(4π/5)Y₂μ|l₂m₂⟩ = (−1)^m₁√[(2l₁+1)(2l₂+1)](l₁ 2 l₂;0 0 0)(l₁ 2 l₂;−m₁ μ m₂).
        Q2_F3D_2sp_2sp_mu.resize(2 * Nsp_I, 2 * Nsp_I, 5);
        for (int mu_I = 0; mu_I < 5; ++mu_I) {
            for (Eigen::Index sp2_I = 0; sp2_I < 2 * Nsp_I; ++sp2_I) {
                for (Eigen::Index sp1_I = 0; sp1_I < 2 * Nsp_I; ++sp1_I) {
                    const auto& label1 = labels_S1D_2sp[sp1_I];
                    const auto& label2 = labels_S1D_2sp[sp2_I];
                    const int mlUp_I = (label1.twom_I - 1) / 2;
                    const int mlDn_I = (label1.twom_I + 1) / 2;
                    const double phaseUp_F = 1.0 - 2.0 * static_cast<double>(mlUp_I % 2 != 0);
                    const double phaseDn_F = 1.0 - 2.0 * static_cast<double>(mlDn_I % 2 != 0);
                    const double norm_F = std::sqrt((2.0 * label1.l_I + 1.0) * (2.0 * label2.l_I + 1.0));
                    const double coupling0_F = gsl_sf_coupling_3j(2 * label1.l_I, 4, 2 * label2.l_I, 0, 0, 0);
                    const double couplingUp_F = gsl_sf_coupling_3j(2 * label1.l_I, 4, 2 * label2.l_I, -label1.twom_I + 1, 2 * (mu_I - 2), label2.twom_I - 1);
                    const double couplingDn_F = gsl_sf_coupling_3j(2 * label1.l_I, 4, 2 * label2.l_I, -label1.twom_I - 1, 2 * (mu_I - 2), label2.twom_I + 1);
                    const double angularUp_F = label1.CGSpinUp_F * label2.CGSpinUp_F * phaseUp_F * couplingUp_F;
                    const double angularDn_F = label1.CGSpinDn_F * label2.CGSpinDn_F * phaseDn_F * couplingDn_F;
                    Q2_F3D_2sp_2sp_mu(sp1_I, sp2_I, mu_I) = r2_F2D_sp_sp(sp1_I % Nsp_I, sp2_I % Nsp_I) * norm_F * coupling0_F * (angularUp_F + angularDn_F);
                }
            }
        }

        // P₀ = [[0,η],[−η,0]]; P† = ½Σab Pab c†a c†b.
        P0_F3D_2sp_2sp_mu.resize(2 * Nsp_I, 2 * Nsp_I, 1);
        P0_F3D_2sp_2sp_mu.setZero();
        for (Eigen::Index sp_I = 0; sp_I < Nsp_I; ++sp_I) {
            P0_F3D_2sp_2sp_mu(sp_I, sp_I + Nsp_I, 0) = setting.eta_F1D_sp(sp_I);
            P0_F3D_2sp_2sp_mu(sp_I + Nsp_I, sp_I, 0) = -setting.eta_F1D_sp(sp_I);
        }

        // P₂μ = Q₂μ P₀; P₂μᵀ = −P₂μ.
        P2_F3D_2sp_2sp_mu.resize(2 * Nsp_I, 2 * Nsp_I, 5);
        for (int mu_I = 0; mu_I < 5; ++mu_I) {
            for (Eigen::Index sp2_I = 0; sp2_I < Nsp_I; ++sp2_I) {
                for (Eigen::Index sp1_I = 0; sp1_I < 2 * Nsp_I; ++sp1_I) {
                    P2_F3D_2sp_2sp_mu(sp1_I, sp2_I, mu_I) = -Q2_F3D_2sp_2sp_mu(sp1_I, sp2_I + Nsp_I, mu_I) * setting.eta_F1D_sp(sp2_I);
                    P2_F3D_2sp_2sp_mu(sp1_I, sp2_I + Nsp_I, mu_I) = Q2_F3D_2sp_2sp_mu(sp1_I, sp2_I, mu_I) * setting.eta_F1D_sp(sp2_I);
                }
            }
        }
    };

    // Species-local operators; strengths remain external.
    build_QP(sphericalsetting_neutron, Q2n_F3D_2spn_2spn_mu, P0n_F3D_2spn_2spn_mu, P2n_F3D_2spn_2spn_mu);
    build_QP(sphericalsetting_proton, Q2p_F3D_2spp_2spp_mu, P0p_F3D_2spp_2spp_mu, P2p_F3D_2spp_2spp_mu);

    // QQ: same-species contraction removed; mixed species counted once.
    Htwobody_C6D_cfgp_cfgn_K_cfgp_cfgn_K = projection_nucleus.calc_two_body_Q_nnnn(4, Q2n_F3D_2spn_2spn_mu).shuffle(Eigen::array<int, 6>{0, 2, 4, 1, 3, 5}) * doubleC(-0.5 * chi2_nn_F, 0.0);
    Htwobody_C6D_cfgp_cfgn_K_cfgp_cfgn_K += projection_nucleus.calc_two_body_Q_nppn(4, Q2n_F3D_2spn_2spn_mu, Q2p_F3D_2spp_2spp_mu).shuffle(Eigen::array<int, 6>{0, 2, 4, 1, 3, 5}) * doubleC(-chi2_np_F, 0.0);
    Htwobody_C6D_cfgp_cfgn_K_cfgp_cfgn_K += projection_nucleus.calc_two_body_Q_pppp(4, Q2p_F3D_2spp_2spp_mu).shuffle(Eigen::array<int, 6>{0, 2, 4, 1, 3, 5}) * doubleC(-0.5 * chi2_pp_F, 0.0);

    // Pairing kernels include ¼; apply −G₀ and −G₂.
    Htwobody_C6D_cfgp_cfgn_K_cfgp_cfgn_K += projection_nucleus.calc_two_body_P_nnnn(0, P0n_F3D_2spn_2spn_mu).shuffle(Eigen::array<int, 6>{0, 2, 4, 1, 3, 5}) * doubleC(-G0_nn_F, 0.0);
    Htwobody_C6D_cfgp_cfgn_K_cfgp_cfgn_K += projection_nucleus.calc_two_body_P_pppp(0, P0p_F3D_2spp_2spp_mu).shuffle(Eigen::array<int, 6>{0, 2, 4, 1, 3, 5}) * doubleC(-G0_pp_F, 0.0);
    Htwobody_C6D_cfgp_cfgn_K_cfgp_cfgn_K += projection_nucleus.calc_two_body_P_nnnn(4, P2n_F3D_2spn_2spn_mu).shuffle(Eigen::array<int, 6>{0, 2, 4, 1, 3, 5}) * doubleC(-G2_nn_F, 0.0);
    Htwobody_C6D_cfgp_cfgn_K_cfgp_cfgn_K += projection_nucleus.calc_two_body_P_pppp(4, P2p_F3D_2spp_2spp_mu).shuffle(Eigen::array<int, 6>{0, 2, 4, 1, 3, 5}) * doubleC(-G2_pp_F, 0.0);
}
