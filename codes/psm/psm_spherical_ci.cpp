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
void PSMSpherical::build_projection(int TargetTwoI_I, int Nalpha_I, int Nbeta_I, int Ngamma_I, int Nphin_I, int Nphip_I) {
    assert(Un_F2D_2spn_2qpn.rows() > 0 && Up_F2D_2spp_2qpp.rows() > 0);
    // (Eqp,2K,cutoffs) → neutron/proton configurations.
    psm_setting.configuration_neutron.build_config(Eqpn_F1D_2qpn, TwoKn_I1D_2qpn);
    psm_setting.configuration_proton.build_config(Eqpp_F1D_2qpp, TwoKp_I1D_2qpp);
    assert(!psm_setting.configuration_neutron.config_I2D_cfg_cqp.empty() && !psm_setting.configuration_proton.config_I2D_cfg_cqp.empty());

    // {labels,Nα,Nβ,Nγ} → rotation caches.
    SphericalRotation rotation_neutron(psm_setting.sphericalsetting_neutron.labels_S1D_sp, Nalpha_I, Nbeta_I, Ngamma_I);
    SphericalRotation rotation_proton(psm_setting.sphericalsetting_proton.labels_S1D_sp, Nalpha_I, Nbeta_I, Ngamma_I);
    rotation_neutron.build();
    rotation_proton.build();

    // (N−Ncore,Z−Zcore,2I,Cₙ,Cₚ) → nuclear projection.
    projection_nucleus = HFBProjectionNucleus(
        HFBProjection(psm_setting.N_I - psm_setting.Ncore_I, TargetTwoI_I, static_cast<int>(Un_F2D_2spn_2qpn.rows()), psm_setting.configuration_neutron.config_I2D_cfg_cqp, psm_setting.configuration_neutron.config_I2D_cfg_cqp, Nphin_I, Nalpha_I, Nbeta_I, Ngamma_I),
        HFBProjection(psm_setting.Z_I - psm_setting.Zcore_I, TargetTwoI_I, static_cast<int>(Up_F2D_2spp_2qpp.rows()), psm_setting.configuration_proton.config_I2D_cfg_cqp, psm_setting.configuration_proton.config_I2D_cfg_cqp, Nphip_I, Nalpha_I, Nbeta_I, Ngamma_I));

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
 * @brief Construct and project one-body and two-body Hamiltonians.
 * @math h₀ = ℏω₀(N+3/2) − κℏωc[2l·s+μ(l²−N(N+3)/2)].
 * @math H² = −χnn:Qn†Qn:/2 − χnp Qn†Qp − χpp:Qp†Qp:/2 − Στλ Gλ,τ Pλ,τ†Pλ,τ.
 * @output H_C6D_cfgp_cfgn_K_cfgp_cfgn_K.
 * @note Fields and strengths from psm_setting; λ = 0,2.
 */
void PSMSpherical::build_hamiltonian() {
    assert(psm_setting.N_I >= psm_setting.Ncore_I && psm_setting.Z_I >= psm_setting.Zcore_I);
    assert(projection_nucleus.projection_neutron.Nsp_I > 0 && projection_nucleus.projection_proton.Nsp_I > 0);
    assert(psm_setting.N_I - psm_setting.Ncore_I == projection_nucleus.projection_neutron.TargetN_I && psm_setting.Z_I - psm_setting.Zcore_I == projection_nucleus.projection_proton.TargetN_I);

    // h₀ ∈ ℝ → complex projection inputs.
    const Eigen::MatrixXcd h0n_C2D_2spn_2spn = h0n_F2D_2spn_2spn.cast<doubleC>();
    const Eigen::MatrixXcd h0p_C2D_2spp_2spp = h0p_F2D_2spp_2spp.cast<doubleC>();

    // Store neutron output before the shared projection buffer changes.
    H_C6D_cfgp_cfgn_K_cfgp_cfgn_K = projection_nucleus.calc_one_body_nn(h0n_C2D_2spn_2spn).shuffle(Eigen::array<int, 6>{0, 2, 4, 1, 3, 5});
    H_C6D_cfgp_cfgn_K_cfgp_cfgn_K += projection_nucleus.calc_one_body_pp(h0p_C2D_2spp_2spp).shuffle(Eigen::array<int, 6>{0, 2, 4, 1, 3, 5});

    // Q₂μ = δ_NN′(r/b)²√(4π/5)Y₂μ; P₂μ = Q₂μ P₀.
    assert(std::isfinite(psm_setting.chi2_nn_F) && std::isfinite(psm_setting.chi2_np_F) && std::isfinite(psm_setting.chi2_pp_F));
    assert(std::isfinite(psm_setting.G0_nn_F) && std::isfinite(psm_setting.G0_pp_F) && std::isfinite(psm_setting.G2_nn_F) && std::isfinite(psm_setting.G2_pp_F));

    // QQ: same-species contraction removed; mixed species counted once.
    H_C6D_cfgp_cfgn_K_cfgp_cfgn_K += projection_nucleus.calc_two_body_Q_nnnn(4, Qn_F3D_2spn_2spn_mu).shuffle(Eigen::array<int, 6>{0, 2, 4, 1, 3, 5}) * doubleC(-0.5 * psm_setting.chi2_nn_F, 0.0);
    H_C6D_cfgp_cfgn_K_cfgp_cfgn_K += projection_nucleus.calc_two_body_Q_nppn(4, Qn_F3D_2spn_2spn_mu, Qp_F3D_2spp_2spp_mu).shuffle(Eigen::array<int, 6>{0, 2, 4, 1, 3, 5}) * doubleC(-psm_setting.chi2_np_F, 0.0);
    H_C6D_cfgp_cfgn_K_cfgp_cfgn_K += projection_nucleus.calc_two_body_Q_pppp(4, Qp_F3D_2spp_2spp_mu).shuffle(Eigen::array<int, 6>{0, 2, 4, 1, 3, 5}) * doubleC(-0.5 * psm_setting.chi2_pp_F, 0.0);

    // Pairing kernels include ¼; apply −G₀ and −G₂.
    H_C6D_cfgp_cfgn_K_cfgp_cfgn_K += projection_nucleus.calc_two_body_P_nnnn(0, P0n_F3D_2spn_2spn_mu).shuffle(Eigen::array<int, 6>{0, 2, 4, 1, 3, 5}) * doubleC(-psm_setting.G0_nn_F, 0.0);
    H_C6D_cfgp_cfgn_K_cfgp_cfgn_K += projection_nucleus.calc_two_body_P_pppp(0, P0p_F3D_2spp_2spp_mu).shuffle(Eigen::array<int, 6>{0, 2, 4, 1, 3, 5}) * doubleC(-psm_setting.G0_pp_F, 0.0);
    H_C6D_cfgp_cfgn_K_cfgp_cfgn_K += projection_nucleus.calc_two_body_P_nnnn(4, P2n_F3D_2spn_2spn_mu).shuffle(Eigen::array<int, 6>{0, 2, 4, 1, 3, 5}) * doubleC(-psm_setting.G2_nn_F, 0.0);
    H_C6D_cfgp_cfgn_K_cfgp_cfgn_K += projection_nucleus.calc_two_body_P_pppp(4, P2p_F3D_2spp_2spp_mu).shuffle(Eigen::array<int, 6>{0, 2, 4, 1, 3, 5}) * doubleC(-psm_setting.G2_pp_F, 0.0);
}

/**
 * @brief Solve configuration mixing using canonical norm orthogonalization.
 * @math N = U n U†; X = Uᵣ nᵣ^(-1/2); (X†HX)c = Ec.
 * @math f = Xc; g = Uᵣc = N^(1/2)f; g†g = 1.
 * @output EPSM_F1D_state and gPSM_C4D_cfgp_cfgn_K_state.
 * @note Retains nᵢ > 1e-10 nmax; energies ascend.
 */
void PSMSpherical::solve_ci(double gamma_F) {
    psm_setting.set_CI_QQ(rhon_F2D_2spn_2spn, rhop_F2D_2spp_2spp, Qn_F3D_2spn_2spn_mu, Qp_F3D_2spp_2spp_mu);
    psm_setting.set_CI_PP(gamma_F);
    build_norm();
    build_hamiltonian();

    // a = cfgp + Ncfgp(cfgn + Ncfgn K); ColMajor.
    const Eigen::Index Ncfgp_I = N_C6D_cfgp_cfgn_K_cfgp_cfgn_K.dimension(0);
    const Eigen::Index Ncfgn_I = N_C6D_cfgp_cfgn_K_cfgp_cfgn_K.dimension(1);
    const Eigen::Index NK_I = N_C6D_cfgp_cfgn_K_cfgp_cfgn_K.dimension(2);
    const Eigen::Index Nbasis_I = Ncfgp_I * Ncfgn_I * NK_I;
    assert(Nbasis_I > 0);
    assert(N_C6D_cfgp_cfgn_K_cfgp_cfgn_K.dimension(3) == Ncfgp_I && N_C6D_cfgp_cfgn_K_cfgp_cfgn_K.dimension(4) == Ncfgn_I && N_C6D_cfgp_cfgn_K_cfgp_cfgn_K.dimension(5) == NK_I);
    assert(H_C6D_cfgp_cfgn_K_cfgp_cfgn_K.dimensions() == N_C6D_cfgp_cfgn_K_cfgp_cfgn_K.dimensions());
    const Eigen::Map<const Eigen::MatrixXcd> N_C2D_basis_basis(N_C6D_cfgp_cfgn_K_cfgp_cfgn_K.data(), Nbasis_I, Nbasis_I);
    const Eigen::Map<const Eigen::MatrixXcd> H_C2D_basis_basis(H_C6D_cfgp_cfgn_K_cfgp_cfgn_K.data(), Nbasis_I, Nbasis_I);
    assert(N_C2D_basis_basis.allFinite() && H_C2D_basis_basis.allFinite());

    // N → (N+N†)/2; H → (H+H†)/2.
    const Eigen::MatrixXcd Nhermitian_C2D_basis_basis = 0.5 * (N_C2D_basis_basis + N_C2D_basis_basis.adjoint());
    const Eigen::MatrixXcd Hhermitian_C2D_basis_basis = 0.5 * (H_C2D_basis_basis + H_C2D_basis_basis.adjoint());
    const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> norm_solver(Nhermitian_C2D_basis_basis);
    assert(norm_solver.info() == Eigen::Success);
    const double normCut_F = 1.0e-10 * norm_solver.eigenvalues().maxCoeff();
    assert(normCut_F > 0.0 && norm_solver.eigenvalues().minCoeff() >= -normCut_F);
    const Eigen::Index Nstate_I = (norm_solver.eigenvalues().array() > normCut_F).count();
    assert(Nstate_I > 0);

    // X = Uᵣ nᵣ^(-1/2); ascending norm eigenvalues.
    const Eigen::MatrixXcd U_C2D_basis_norm = norm_solver.eigenvectors().rightCols(Nstate_I);
    const Eigen::VectorXd inverseSqrtN_F1D_norm = norm_solver.eigenvalues().tail(Nstate_I).array().sqrt().inverse();
    Eigen::MatrixXcd X_C2D_basis_norm{};
    X_C2D_basis_norm.noalias() = U_C2D_basis_norm * inverseSqrtN_F1D_norm.asDiagonal();
    Eigen::MatrixXcd HX_C2D_basis_norm{};
    HX_C2D_basis_norm.noalias() = Hhermitian_C2D_basis_basis * X_C2D_basis_norm;
    Eigen::MatrixXcd H_C2D_norm_norm{};
    H_C2D_norm_norm.noalias() = X_C2D_basis_norm.adjoint() * HX_C2D_basis_norm;
    H_C2D_norm_norm = (0.5 * (H_C2D_norm_norm + H_C2D_norm_norm.adjoint())).eval();
    const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> energy_solver(H_C2D_norm_norm);
    assert(energy_solver.info() == Eigen::Success);

    // g = Uᵣc; gν†gν′ = δνν′ in the original basis labels.
    EPSM_F1D_state = energy_solver.eigenvalues();
    gPSM_C4D_cfgp_cfgn_K_state.resize(Ncfgp_I, Ncfgn_I, NK_I, Nstate_I);
    Eigen::Map<Eigen::MatrixXcd> gPSM_C2D_basis_state(gPSM_C4D_cfgp_cfgn_K_state.data(), Nbasis_I, Nstate_I);
    gPSM_C2D_basis_state.noalias() = U_C2D_basis_norm * energy_solver.eigenvectors();
}
