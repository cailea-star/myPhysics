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

    // Axial projection: Nα = Nγ = 1.
    SphericalRotation rotation_neutron(psm_setting.sphericalsetting_neutron.labels_S1D_sp, 1, Nbeta_I, 1);
    SphericalRotation rotation_proton(psm_setting.sphericalsetting_proton.labels_S1D_sp, 1, Nbeta_I, 1);
    rotation_neutron.build();
    rotation_proton.build();

    // (N−Ncore,Z−Zcore,Cₙ,Cₚ) → nuclear projection.
    projection_nucleus = HFBProjectionNucleus(
        HFBProjection(static_cast<int>(Un_F2D_2spn_2qpn.rows()), psm_setting.configuration_neutron.config_I2D_cfg_cqp, psm_setting.configuration_neutron.config_I2D_cfg_cqp, psm_setting.N_I - psm_setting.Ncore_I, Nphin_I, 1, Nbeta_I, 1),
        HFBProjection(static_cast<int>(Up_F2D_2spp_2qpp.rows()), psm_setting.configuration_proton.config_I2D_cfg_cqp, psm_setting.configuration_proton.config_I2D_cfg_cqp, psm_setting.Z_I - psm_setting.Zcore_I, Nphip_I, 1, Nbeta_I, 1));

    // 2K_cfg = Σqp 2K_qp; Φ₁ = Φ₂.
    const auto set_config_K = [](HFBProjection& projection, const auto& config_I2D_cfg_cqp, const Eigen::VectorXi& TwoK_I1D_qp) {
        projection.TwoK1_I1D_cfg1.resize(config_I2D_cfg_cqp.size());
        projection.TwoK1_I1D_cfg1.setZero();
        for (int cfg_I = 0; cfg_I < static_cast<int>(config_I2D_cfg_cqp.size()); ++cfg_I) {
            for (const int qp_I : config_I2D_cfg_cqp[cfg_I]) {
                projection.TwoK1_I1D_cfg1(cfg_I) += TwoK_I1D_qp(qp_I);
            }
        }
        projection.TwoK2_I1D_cfg2 = projection.TwoK1_I1D_cfg1;
    };
    set_config_K(projection_nucleus.projection_neutron, psm_setting.configuration_neutron.config_I2D_cfg_cqp, TwoKn_I1D_2qpn);
    set_config_K(projection_nucleus.projection_proton, psm_setting.configuration_proton.config_I2D_cfg_cqp, TwoKp_I1D_2qpp);

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
 * @output N_C7D and H_C7D angle caches.
 * @note Requires solved HFBCS and initialized projection.
 */
void PSMSpherical::build_ci(double gamma_F) {
    psm_setting.set_CI_QQ(rhon_F2D_2spn_2spn, rhop_F2D_2spp_2spp, Qn_F3D_2spn_2spn_mu, Qp_F3D_2spp_2spp_mu);
    psm_setting.set_CI_PP(gamma_F);

    // N_ab(Ω) = ⟨Φ_a|PₙPₚR(Ω)|Φ_b⟩.
    const Eigen::Index Ncfgn_I = projection_nucleus.projection_neutron.Ncfg1_I;
    const Eigen::Index Ncfgp_I = projection_nucleus.projection_proton.Ncfg1_I;
    assert(Ncfgn_I > 0 && Ncfgp_I > 0);
    assert(Ncfgn_I == projection_nucleus.projection_neutron.Ncfg2_I && Ncfgp_I == projection_nucleus.projection_proton.Ncfg2_I);

    N_C7D_cfgn_cfgn_cfgp_cfgp_alpha_beta_gamma = projection_nucleus.build_overlap();

    // H_ab(Ω) = ⟨Φ_a|HPₙPₚR(Ω)|Φ_b⟩.
    assert(psm_setting.N_I >= psm_setting.Ncore_I && psm_setting.Z_I >= psm_setting.Zcore_I);
    assert(projection_nucleus.projection_neutron.Nsp_I > 0 && projection_nucleus.projection_proton.Nsp_I > 0);
    assert(psm_setting.N_I - psm_setting.Ncore_I == projection_nucleus.projection_neutron.TargetN_I && psm_setting.Z_I - psm_setting.Zcore_I == projection_nucleus.projection_proton.TargetN_I);

    // h₀ ∈ ℝ → complex projection inputs.
    const Eigen::MatrixXcd h0n_C2D_2spn_2spn = h0n_F2D_2spn_2spn.cast<doubleC>();
    const Eigen::MatrixXcd h0p_C2D_2spp_2spp = h0p_F2D_2spp_2spp.cast<doubleC>();

    // H¹(Ω) = Hν¹(Ω) Nπ(Ω) + Nν(Ω) Hπ¹(Ω).
    H_C7D_cfgn_cfgn_cfgp_cfgp_alpha_beta_gamma = projection_nucleus.build_one_body_nn(h0n_C2D_2spn_2spn);
    H_C7D_cfgn_cfgn_cfgp_cfgp_alpha_beta_gamma += projection_nucleus.build_one_body_pp(h0p_C2D_2spp_2spp);

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

    // H²(Ω) = Hν²(Ω) Nπ(Ω) + Nν(Ω) Hπ²(Ω) + Hνπ(Ω).
    H_C7D_cfgn_cfgn_cfgp_cfgp_alpha_beta_gamma += projection_nucleus.build_two_body_nn(gQn_F1D_i, Qn_C3D_2spn_2spn_i, gPn_F1D_i, Pn_C3D_2spn_2spn_i);
    H_C7D_cfgn_cfgn_cfgp_cfgp_alpha_beta_gamma += projection_nucleus.build_two_body_pp(gQp_F1D_i, Qp_C3D_2spp_2spp_i, gPp_F1D_i, Pp_C3D_2spp_2spp_i);
    H_C7D_cfgn_cfgn_cfgp_cfgp_alpha_beta_gamma += projection_nucleus.build_two_body_pn(Eigen::VectorXd::Constant(5, -psm_setting.chi2_np_F), Qn_C3D_2spn_2spn_i, Qp_C3D_2spp_2spp_i);
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
    assert(N_C7D_cfgn_cfgn_cfgp_cfgp_alpha_beta_gamma.size() > 0 && H_C7D_cfgn_cfgn_cfgp_cfgp_alpha_beta_gamma.size() > 0);

    const Eigen::Tensor<doubleC, 6, Eigen::ColMajor> N_C6D_cfgn_cfgp_K_cfgn_cfgp_K = projection_nucleus.integrate_amp(TargetTwoI_I, N_C7D_cfgn_cfgn_cfgp_cfgp_alpha_beta_gamma).shuffle(Eigen::array<int, 6>{0, 2, 4, 1, 3, 5});
    const Eigen::Tensor<doubleC, 6, Eigen::ColMajor> H_C6D_cfgn_cfgp_K_cfgn_cfgp_K = projection_nucleus.integrate_amp(TargetTwoI_I, H_C7D_cfgn_cfgn_cfgp_cfgp_alpha_beta_gamma).shuffle(Eigen::array<int, 6>{0, 2, 4, 1, 3, 5});
    // (cfgn,cfgp,K,cfgn′,cfgp′,K′) → (cfgncfgpK,cfgncfgpK′); zero-copy.
    const Eigen::Index NcfgncfgpK_I = N_C6D_cfgn_cfgp_K_cfgn_cfgp_K.dimension(0) * N_C6D_cfgn_cfgp_K_cfgn_cfgp_K.dimension(1) * N_C6D_cfgn_cfgp_K_cfgn_cfgp_K.dimension(2);
    const Eigen::Map<const Eigen::MatrixXcd> N_C2D_cfgncfgpK_cfgncfgpK(N_C6D_cfgn_cfgp_K_cfgn_cfgp_K.data(), NcfgncfgpK_I, NcfgncfgpK_I);
    const Eigen::Map<const Eigen::MatrixXcd> H_C2D_cfgncfgpK_cfgncfgpK(H_C6D_cfgn_cfgp_K_cfgn_cfgp_K.data(), NcfgncfgpK_I, NcfgncfgpK_I);

    // N → (N+N†)/2; H → (H+H†)/2.
    const Eigen::MatrixXcd Nhermitian_C2D_cfgncfgpK_cfgncfgpK = 0.5 * (N_C2D_cfgncfgpK_cfgncfgpK + N_C2D_cfgncfgpK_cfgncfgpK.adjoint());
    const Eigen::MatrixXcd Hhermitian_C2D_cfgncfgpK_cfgncfgpK = 0.5 * (H_C2D_cfgncfgpK_cfgncfgpK + H_C2D_cfgncfgpK_cfgncfgpK.adjoint());
    const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> N_solver(Nhermitian_C2D_cfgncfgpK_cfgncfgpK);
    assert(N_solver.info() == Eigen::Success);
    const double normCut_F = 1.0e-10 * N_solver.eigenvalues().maxCoeff();
    assert(normCut_F > 0.0 && N_solver.eigenvalues().minCoeff() >= -normCut_F);
    const Eigen::Index Nstate_I = (N_solver.eigenvalues().array() > normCut_F).count();
    assert(Nstate_I > 0);

    // X = Uᵣ nᵣ^(-1/2); ascending norm eigenvalues.
    const Eigen::MatrixXcd U_C2D_cfgncfgpK_eigenN = N_solver.eigenvectors().rightCols(Nstate_I);
    const Eigen::VectorXd n_F1D_eigenN = N_solver.eigenvalues().tail(Nstate_I);
    Eigen::MatrixXcd X_C2D_cfgncfgpK_eigenN = U_C2D_cfgncfgpK_eigenN * n_F1D_eigenN.array().sqrt().inverse().matrix().asDiagonal();
    Eigen::MatrixXcd H_C2D_eigenN_eigenN = X_C2D_cfgncfgpK_eigenN.adjoint() * Hhermitian_C2D_cfgncfgpK_cfgncfgpK * X_C2D_cfgncfgpK_eigenN;
    H_C2D_eigenN_eigenN = (0.5 * (H_C2D_eigenN_eigenN + H_C2D_eigenN_eigenN.adjoint())).eval();
    const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> H_solver(H_C2D_eigenN_eigenN);
    assert(H_solver.info() == Eigen::Success);

    // (X†HX)Xci = Xci E; Xci†Xci = 1.
    Eci_F1D_eigenH = H_solver.eigenvalues();
    Xci_C2D_eigenN_eigenH = H_solver.eigenvectors();
}
