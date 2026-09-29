/**
 * @file    psm_spherical_operators.cpp
 * @author  cailea
 * @date    2026-09-29
 * @brief   Spherical PSM operators and quadrupole strengths.
 */

#include "psm_spherical.hpp"
#include "spherical_basis.hpp"
#include <gsl/gsl_sf_coupling.h>
#include <cassert>
#include <cmath>
#include <vector>

void PSMSpherical::build_Q(const SphericalSetting& sphericalsetting, Eigen::Tensor<double, 3, Eigen::ColMajor>& Q_F3D_2sp_2sp_mu) {
    const Eigen::Index Nsp_I = static_cast<Eigen::Index>(sphericalsetting.labels_S1D_sp.size());
    assert(Nsp_I > 0);
    std::vector<SphericalSPLabel> labels_S1D_2sp = sphericalsetting.labels_S1D_sp;
    labels_S1D_2sp.reserve(2 * Nsp_I);
    for (const auto& label : sphericalsetting.labels_S1D_sp) {
        labels_S1D_2sp.emplace_back(label.N_I, label.twoj_I, -label.twom_I);
    }

    // ∫r²dr R₁R₂(r/b)²; retain ΔN = 0 as PSM.
    const SphericalLaguerreBasis radial_basis(sphericalsetting.b_F, sphericalsetting.Nr_I, sphericalsetting.labels_S1D_sp);
    const Eigen::VectorXd r2Weights_F1D_r = radial_basis.w_F1D_r.array() * radial_basis.r_F1D_r.array().square() / (sphericalsetting.b_F * sphericalsetting.b_F);
    Eigen::MatrixXd r2_F2D_sp_sp{};
    r2_F2D_sp_sp.resize(Nsp_I, Nsp_I);
    for (Eigen::Index sp2_I = 0; sp2_I < Nsp_I; ++sp2_I) {
        for (Eigen::Index sp1_I = 0; sp1_I < Nsp_I; ++sp1_I) {
            const double deltaN_F = static_cast<double>(sphericalsetting.labels_S1D_sp[sp1_I].N_I == sphericalsetting.labels_S1D_sp[sp2_I].N_I);
            r2_F2D_sp_sp(sp1_I, sp2_I) = deltaN_F * (radial_basis.phi_F2D_sp_r.row(sp1_I).array() * radial_basis.phi_F2D_sp_r.row(sp2_I).array() * r2Weights_F1D_r.transpose().array()).sum();
        }
    }

    // ⟨l₁m₁|√(4π/5)Y₂μ|l₂m₂⟩ = (−1)^m₁√[(2l₁+1)(2l₂+1)](l₁ 2 l₂;0 0 0)(l₁ 2 l₂;−m₁ μ m₂).
    Q_F3D_2sp_2sp_mu.resize(2 * Nsp_I, 2 * Nsp_I, 5);
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
                Q_F3D_2sp_2sp_mu(sp1_I, sp2_I, mu_I) = r2_F2D_sp_sp(sp1_I % Nsp_I, sp2_I % Nsp_I) * norm_F * coupling0_F * (angularUp_F + angularDn_F);
            }
        }
    }

}

void PSMSpherical::build_P(const SphericalSetting& setting, const Eigen::Tensor<double, 3, Eigen::ColMajor>& Q2_F3D_2sp_2sp_mu, Eigen::Tensor<double, 3, Eigen::ColMajor>& P0_F3D_2sp_2sp_mu, Eigen::Tensor<double, 3, Eigen::ColMajor>& P2_F3D_2sp_2sp_mu) {
    const Eigen::Index Nsp_I = static_cast<Eigen::Index>(setting.labels_S1D_sp.size());
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
}

void PSMSpherical::build_h0(const SphericalSetting& setting, double hbarOmega0_F, double hbarOmega00_F, double kappa_F, double mu_F, Eigen::MatrixXd& h0_F2D_2sp_2sp) {
    const Eigen::Index Nsp_I = static_cast<Eigen::Index>(setting.labels_S1D_sp.size());
    h0_F2D_2sp_2sp.resize(2 * Nsp_I, 2 * Nsp_I);
    h0_F2D_2sp_2sp.setZero();
    for (Eigen::Index sp_I = 0; sp_I < Nsp_I; ++sp_I) {
        const auto& label = setting.labels_S1D_sp[sp_I];
        const double j_F = 0.5 * label.twoj_I;
        const double l2_F = label.l_I * (label.l_I + 1.0);
        const double twoLS_F = j_F * (j_F + 1.0) - l2_F - 0.75;
        const double l2Average_F = 0.5 * label.N_I * (label.N_I + 3.0);

        // ε_Nlj = ℏω₀(N+3/2) − κℏωc[2l·s+μ(l²−⟨l²⟩_N)].
        const double spherical_F = hbarOmega0_F * (label.N_I + 1.5) - kappa_F * hbarOmega00_F * (twoLS_F + mu_F * (l2_F - l2Average_F));
        h0_F2D_2sp_2sp(sp_I, sp_I) = spherical_F;
        h0_F2D_2sp_2sp(sp_I + Nsp_I, sp_I + Nsp_I) = spherical_F;
    }
}

void PSMSpherical::build_h0_deformation(const SphericalSetting& setting, double hbarOmega0_F, double epsilon2_F, double epsilon4_F, Eigen::MatrixXd& hDeform_F2D_2sp_2sp) {
    assert(setting.useAxialSym_B && setting.useParity_B);
    const int Nsp_I = static_cast<int>(setting.labels_S1D_sp.size());
    const double epsilon20_F = 2.0 * epsilon2_F / 3.0;
    hDeform_F2D_2sp_2sp.resize(2 * Nsp_I, 2 * Nsp_I);
    hDeform_F2D_2sp_2sp.setZero();

    // ⟨l₁m|P_L|l₂m⟩ = (−1)^m√[(2l₁+1)(2l₂+1)](l₁ L l₂;0 0 0)(l₁ L l₂;−m 0 m).
    const auto calc_legendre = [](const SphericalSPLabel& label1, const SphericalSPLabel& label2, int L_I) {
        const double norm_F = std::sqrt((2.0 * label1.l_I + 1.0) * (2.0 * label2.l_I + 1.0));
        const double coupling0_F = gsl_sf_coupling_3j(2 * label1.l_I, 2 * L_I, 2 * label2.l_I, 0, 0, 0);
        const int mlUp_I = (label1.twom_I - 1) / 2;
        const int mlDn_I = (label1.twom_I + 1) / 2;
        const double phaseUp_F = 1.0 - 2.0 * static_cast<double>(std::abs(mlUp_I) % 2);
        const double phaseDn_F = 1.0 - 2.0 * static_cast<double>(std::abs(mlDn_I) % 2);
        const double couplingUp_F = gsl_sf_coupling_3j(2 * label1.l_I, 2 * L_I, 2 * label2.l_I, -2 * mlUp_I, 0, 2 * mlUp_I);
        const double couplingDn_F = gsl_sf_coupling_3j(2 * label1.l_I, 2 * L_I, 2 * label2.l_I, -2 * mlDn_I, 0, 2 * mlDn_I);
        const double angularUp_F = label1.CGSpinUp_F * label2.CGSpinUp_F * phaseUp_F * couplingUp_F;
        const double angularDn_F = label1.CGSpinDn_F * label2.CGSpinDn_F * phaseDn_F * couplingDn_F;
        return norm_F * coupling0_F * (angularUp_F + angularDn_F);
    };

    // ∫r²dr R₁R₂(r/b)²; existing shell and axial blocks.
    const SphericalLaguerreBasis radial_basis(setting.b_F, setting.Nr_I, setting.labels_S1D_sp);
    const Eigen::VectorXd r2Weights_F1D_r = radial_basis.w_F1D_r.array() * radial_basis.r_F1D_r.array().square() / (setting.b_F * setting.b_F);
    for (const auto& indices : setting.indices_I2D_block_bsp) {
        for (int bsp1_I = 0; bsp1_I < static_cast<int>(indices.size()); ++bsp1_I) {
            const int sp1_I = indices[bsp1_I];
            const auto& label1 = setting.labels_S1D_sp[sp1_I];
            for (int bsp2_I = 0; bsp2_I <= bsp1_I; ++bsp2_I) {
                const int sp2_I = indices[bsp2_I];
                const auto& label2 = setting.labels_S1D_sp[sp2_I];
                assert(label1.N_I == label2.N_I && label1.twom_I == label2.twom_I);
                const double r2_F = (radial_basis.phi_F2D_sp_r.row(sp1_I).array() * radial_basis.phi_F2D_sp_r.row(sp2_I).array() * r2Weights_F1D_r.transpose().array()).sum();
                const double P2_F = calc_legendre(label1, label2, 2);
                const double P4_F = calc_legendre(label1, label2, 4);
                const double deformation_F = hbarOmega0_F * r2_F * (-epsilon20_F * P2_F + epsilon4_F * P4_F);

                // h⁻⁻ = ηh⁺⁺η; h⁺⁻ = h⁻⁺ = 0.
                hDeform_F2D_2sp_2sp(sp1_I, sp2_I) = deformation_F;
                hDeform_F2D_2sp_2sp(sp2_I, sp1_I) = deformation_F;
                hDeform_F2D_2sp_2sp(sp1_I + Nsp_I, sp2_I + Nsp_I) = setting.eta_F1D_sp(sp1_I) * setting.eta_F1D_sp(sp2_I) * deformation_F;
                hDeform_F2D_2sp_2sp(sp2_I + Nsp_I, sp1_I + Nsp_I) = hDeform_F2D_2sp_2sp(sp1_I + Nsp_I, sp2_I + Nsp_I);
            }
        }
    }
}
