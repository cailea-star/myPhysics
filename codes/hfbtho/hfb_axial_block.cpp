/**
 * @file    hfb_axial_block.cpp
 * @author  cailea
 * @date    2026-05-06
 * @brief   Implement axial HFB block updates.
 */

#include <cmath>
#include <vector>

#include <Eigen/Eigenvalues>

#include "hfb_axial.hpp"

void AxialHFBBlock::update_Gamma_Delta_from_field(const AxialHFBField& field_, const AxialBasis& global_basis_) {
    const int Nz_I = static_cast<int>(field_.vcent_F2D_z_r.rows());
    const int Nr_I = static_cast<int>(field_.vcent_F2D_z_r.cols());
    const int Nbsp_I = static_cast<int>(labels_S1D_bsp.size());
    Gamma_F2D_bsp_bsp.setZero();
    Delta_F2D_bsp_bsp.setZero();

    // (r,z) → (Γ,Δ).
    const int spBegin_I = indices_I1D_bsp.front();
    for (int r_I = 0; r_I < Nr_I; ++r_I) {
        for (int z_I = 0; z_I < Nz_I; ++z_I) {
            const double w_F = global_basis_.w_F2D_z_r(z_I, r_I);
            const double rinv_F = 1.0 / global_basis_.r_F1D_r(r_I);
            const double rinv2_F = rinv_F * rinv_F;

            // (sp,z,r) → {φ,∂φ,∂²φ}.
            const Eigen::Map<const Eigen::VectorXd> phi_F1D_bsp(&global_basis_.phi_F3D_sp_z_r(spBegin_I, z_I, r_I), Nbsp_I);
            const Eigen::Map<const Eigen::VectorXd> dphidr_F1D_bsp(&global_basis_.dphidr_F3D_sp_z_r(spBegin_I, z_I, r_I), Nbsp_I);
            const Eigen::Map<const Eigen::VectorXd> dphidz_F1D_bsp(&global_basis_.dphidz_F3D_sp_z_r(spBegin_I, z_I, r_I), Nbsp_I);
            const Eigen::Map<const Eigen::VectorXd> ddphidr_F1D_bsp(&global_basis_.ddphidr_F3D_sp_z_r(spBegin_I, z_I, r_I), Nbsp_I);
            const Eigen::Map<const Eigen::VectorXd> ddphidz_F1D_bsp(&global_basis_.ddphidz_F3D_sp_z_r(spBegin_I, z_I, r_I), Nbsp_I);

            // (z,r) → {v}.
            const double vmass_F = field_.vmass_F2D_z_r(z_I, r_I);
            const double vcent_F = field_.vcent_F2D_z_r(z_I, r_I);
            const double vpair_F = field_.vpair_F2D_z_r(z_I, r_I);
            const double vD2_F = field_.vD2_F2D_z_r(z_I, r_I);
            const double vDr_F = field_.vDr_F2D_z_r(z_I, r_I);
            const double vDz_F = field_.vDz_F2D_z_r(z_I, r_I);
            const double vdJ_F = field_.vdJ_F2D_z_r(z_I, r_I);
            const double vJphiz_F = field_.vJphiz_F2D_z_r(z_I, r_I);
            const double vJzphi_F = field_.vJzphi_F2D_z_r(z_I, r_I);
            const double vJphir_F = field_.vJphir_F2D_z_r(z_I, r_I);
            const double vJrphi_F = field_.vJrphi_F2D_z_r(z_I, r_I);

            // Σ_1=Σ_2 → (Γ_12,Δ_12).
            const auto add_same_spin_Func = [&](const std::vector<int>& indices_I1D_bsp_, double spinSign_F) {
                const int Nspin_I = static_cast<int>(indices_I1D_bsp_.size());
                for (int position2_I = 0; position2_I < Nspin_I; ++position2_I) {
                    for (int position1_I = position2_I; position1_I < Nspin_I; ++position1_I) {
                        const int bsp1_I = indices_I1D_bsp_[position1_I];
                        const int bsp2_I = indices_I1D_bsp_[position2_I];
                        const int Lambda_I = labels_S1D_bsp[bsp1_I].Lambda_I;

                        // (φ_1,φ_2,Λ) → transition densities.
                        const double phi1_F = phi_F1D_bsp(bsp1_I);
                        const double phi2_F = phi_F1D_bsp(bsp2_I);
                        const double dphidr1_F = dphidr_F1D_bsp(bsp1_I);
                        const double dphidr2_F = dphidr_F1D_bsp(bsp2_I);
                        const double dphidz1_F = dphidz_F1D_bsp(bsp1_I);
                        const double dphidz2_F = dphidz_F1D_bsp(bsp2_I);
                        const double Lambda2_F = Lambda_I * Lambda_I;
                        const double laplacianPhi1_F = ddphidr_F1D_bsp(bsp1_I) + dphidr1_F * rinv_F + ddphidz_F1D_bsp(bsp1_I) - rinv2_F * Lambda2_F * phi1_F;
                        const double laplacianPhi2_F = ddphidr_F1D_bsp(bsp2_I) + dphidr2_F * rinv_F + ddphidz_F1D_bsp(bsp2_I) - rinv2_F * Lambda2_F * phi2_F;
                        const double rho12_F = phi1_F * phi2_F;
                        const double tau12_F = dphidr1_F * dphidr2_F + dphidz1_F * dphidz2_F + rinv2_F * Lambda2_F * rho12_F;
                        const double rhoD2_12_F = laplacianPhi1_F * phi2_F + phi1_F * laplacianPhi2_F + 2.0 * tau12_F;
                        const double rhoDr_12_F = dphidr1_F * phi2_F + phi1_F * dphidr2_F;
                        const double rhoDz_12_F = dphidz1_F * phi2_F + phi1_F * dphidz2_F;
                        const double Jphiz12_F = spinSign_F * 2.0 * Lambda_I * rinv_F * rho12_F;
                        const double dJ12_F = spinSign_F * Lambda_I * rinv_F * rhoDr_12_F;

                        // ({v},{ρ_12,τ_12,J_12}) → (Γ_12,Δ_12).
                        const double Delta12_F = vpair_F * rho12_F;
                        const double Gamma12_F = vcent_F * rho12_F + vmass_F * tau12_F + vD2_F * rhoD2_12_F + vDr_F * rhoDr_12_F + vDz_F * rhoDz_12_F + vdJ_F * dJ12_F + vJphiz_F * Jphiz12_F;
                        Gamma_F2D_bsp_bsp(bsp1_I, bsp2_I) += Gamma12_F * w_F;
                        Delta_F2D_bsp_bsp(bsp1_I, bsp2_I) += Delta12_F * w_F;
                        const int offDiagonal_I = static_cast<int>(bsp1_I != bsp2_I);
                        Gamma_F2D_bsp_bsp(bsp2_I, bsp1_I) += Gamma12_F * w_F * offDiagonal_I;
                        Delta_F2D_bsp_bsp(bsp2_I, bsp1_I) += Delta12_F * w_F * offDiagonal_I;
                    }
                }
            };
            add_same_spin_Func(indices_I1D_bup, 1.0);
            add_same_spin_Func(indices_I1D_bdn, -1.0);

            // Σ_↑≠Σ_↓ → Γ_↑↓.
            for (int up_I = 0; up_I < static_cast<int>(indices_I1D_bup.size()); ++up_I) {
                for (int down_I = 0; down_I < static_cast<int>(indices_I1D_bdn.size()); ++down_I) {
                    const int bspUp_I = indices_I1D_bup[up_I];
                    const int bspDown_I = indices_I1D_bdn[down_I];
                    const int LambdaUp_I = labels_S1D_bsp[bspUp_I].Lambda_I;
                    const int LambdaDown_I = labels_S1D_bsp[bspDown_I].Lambda_I;

                    // (φ_↑,φ_↓,Λ_↑,Λ_↓) → J_↑↓.
                    const double phiUp_F = phi_F1D_bsp(bspUp_I);
                    const double phiDown_F = phi_F1D_bsp(bspDown_I);
                    const double dphidrUp_F = dphidr_F1D_bsp(bspUp_I);
                    const double dphidrDown_F = dphidr_F1D_bsp(bspDown_I);
                    const double dphidzUp_F = dphidz_F1D_bsp(bspUp_I);
                    const double dphidzDown_F = dphidz_F1D_bsp(bspDown_I);
                    const double JrphiUpDown_F = dphidrUp_F * phiDown_F - phiUp_F * dphidrDown_F;
                    const double JzphiUpDown_F = dphidzUp_F * phiDown_F - phiUp_F * dphidzDown_F;
                    const double JphirUpDown_F = rinv_F * (LambdaUp_I * phiUp_F * phiDown_F + LambdaDown_I * phiUp_F * phiDown_F);
                    const double dJUpDown_F = dphidrUp_F * dphidzDown_F - dphidzUp_F * dphidrDown_F - rinv_F * (LambdaUp_I * phiUp_F * dphidzDown_F + LambdaDown_I * phiDown_F * dphidzUp_F);

                    // ({v},{J_↑↓}) → Γ_↑↓.
                    const double GammaUpDown_F = vdJ_F * dJUpDown_F + vJzphi_F * JzphiUpDown_F + vJphir_F * JphirUpDown_F + vJrphi_F * JrphiUpDown_F;
                    Gamma_F2D_bsp_bsp(bspUp_I, bspDown_I) += GammaUpDown_F * w_F;
                    Gamma_F2D_bsp_bsp(bspDown_I, bspUp_I) += GammaUpDown_F * w_F;
                }
            }
        }
    }
}

void AxialHFBBlock::add_lipkin_nogami(double lambda2_F) {
    // Γ → Γ + 2λ₂(2ρ-I).
    Gamma_F2D_bsp_bsp.noalias() += 4.0 * lambda2_F * rho_F2D_bsp_bsp;
    Gamma_F2D_bsp_bsp.diagonal().array() -= 2.0 * lambda2_F;
}

void AxialHFBBlock::update_UV_E_rho_kappa(const HFBSettings& hfbsettings_, double lambda_F) {
    const int Nbsp_I = static_cast<int>(labels_S1D_bsp.size());

    // Γ-λI → h.
    Eigen::MatrixXd h_F2D_bsp_bsp = Gamma_F2D_bsp_bsp;
    h_F2D_bsp_bsp.diagonal().array() -= lambda_F;

    // (h,Δ) → ℋ.
    Eigen::MatrixXd Hhfb_F2D_2bsp_2bsp(2 * Nbsp_I, 2 * Nbsp_I);
    Hhfb_F2D_2bsp_2bsp.topLeftCorner(Nbsp_I, Nbsp_I) = h_F2D_bsp_bsp;
    Hhfb_F2D_2bsp_2bsp.topRightCorner(Nbsp_I, Nbsp_I) = Delta_F2D_bsp_bsp;
    Hhfb_F2D_2bsp_2bsp.bottomLeftCorner(Nbsp_I, Nbsp_I) = Delta_F2D_bsp_bsp;
    Hhfb_F2D_2bsp_2bsp.bottomRightCorner(Nbsp_I, Nbsp_I) = -h_F2D_bsp_bsp;

    // ℋ → (E,X).
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver_(Hhfb_F2D_2bsp_2bsp);
    const Eigen::VectorXd& eigenvalues_F1D_state = eigensolver_.eigenvalues();
    const Eigen::MatrixXd& eigenvectors_F2D_2bsp_state = eigensolver_.eigenvectors();
    constexpr double EspCutTolerance_F = 1.0e-6;
    static const double EspCutTail_F = std::log(1.0 / EspCutTolerance_F - 1.0) / 100.0;

    // E>0 → (E_bqp,U,V).
    Eqp_F1D_bqp.setZero();
    U_F2D_bsp_bqp.setZero();
    V_F2D_bsp_bqp.setZero();
    int Nbqp_I = 0;
    for (int state_I = 0; state_I < 2 * Nbsp_I; ++state_I) {
        if (eigenvalues_F1D_state(state_I) > 0.0) {
            Eqp_F1D_bqp(Nbqp_I) = eigenvalues_F1D_state(state_I);
            U_F2D_bsp_bqp.col(Nbqp_I) = eigenvectors_F2D_2bsp_state.col(state_I).head(Nbsp_I);
            V_F2D_bsp_bqp.col(Nbqp_I) = eigenvectors_F2D_2bsp_state.col(state_I).tail(Nbsp_I);
            ++Nbqp_I;
        }
    }

    // (E_bqp,V,λ,T,E_cut) → (f_U,f_V,f_κ).
    Eigen::VectorXd factorU_F1D_bqp = Eigen::VectorXd::Zero(Nbqp_I);
    Eigen::VectorXd factorV_F1D_bqp = Eigen::VectorXd::Zero(Nbqp_I);
    Eigen::VectorXd factorK_F1D_bqp = Eigen::VectorXd::Zero(Nbqp_I);
    for (int bqp_I = 0; bqp_I < Nbqp_I; ++bqp_I) {
        const double Eqp_F = Eqp_F1D_bqp(bqp_I);
        const double Vnorm2_F = V_F2D_bsp_bqp.col(bqp_I).squaredNorm();
        const double Esp_F = Eqp_F * (1.0 - 2.0 * Vnorm2_F) + lambda_F;
        if (hfbsettings_.useEspCut_B && Esp_F > hfbsettings_.EspCut_F + EspCutTail_F) {continue;}
        const double occupation_F = (hfbsettings_.temperature_F > 0.0) ? 0.5 * (1.0 - std::tanh(0.5 * Eqp_F / hfbsettings_.temperature_F)) : 0.0;
        factorU_F1D_bqp(bqp_I) = occupation_F;
        factorV_F1D_bqp(bqp_I) = 1.0 - occupation_F;
        factorK_F1D_bqp(bqp_I) = 1.0 - 2.0 * occupation_F;
    }

    // (U,V,f_U,f_V,f_κ) → (ρ,κ).
    rho_F2D_bsp_bsp.noalias() = V_F2D_bsp_bqp * factorV_F1D_bqp.asDiagonal() * V_F2D_bsp_bqp.transpose();
    rho_F2D_bsp_bsp.noalias() += U_F2D_bsp_bqp * factorU_F1D_bqp.asDiagonal() * U_F2D_bsp_bqp.transpose();
    kappa_F2D_bsp_bsp.noalias() = V_F2D_bsp_bqp * factorK_F1D_bqp.asDiagonal() * U_F2D_bsp_bqp.transpose();

    // (ρ,κ) → ([ρ+ρᵀ]/2,-[κ+κᵀ]/2).
    rho_F2D_bsp_bsp = 0.5 * (rho_F2D_bsp_bsp + rho_F2D_bsp_bsp.transpose()).eval();
    kappa_F2D_bsp_bsp = -0.5 * (kappa_F2D_bsp_bsp + kappa_F2D_bsp_bsp.transpose()).eval();
}
