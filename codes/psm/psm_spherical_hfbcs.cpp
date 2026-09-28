/**
 * @file    psm_spherical_hfbcs.cpp
 * @author  cailea
 * @date    2026-09-27
 * @brief   Nilsson single-particle fields and BCS pairing.
 */

#include "psm_spherical.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_sf_legendre.h>

#include "integration_gauss.hpp"
#include "spherical_basis.hpp"

/**
 * @brief Build Nilsson matrices using quadrature; solve BCS using Brent.
 * @math h = ℏω₀(N+3/2) − κℏωc[2l·s+μ(l²−N(N+3)/2)]
 * @math     + ℏω₀(r/b)²[−(2ε₂/3)P₂+ε₄P₄]; Γ = 0.
 * @output Updated species fields, eigenstates, and BCS solutions.
 * @note b fixes ℏω₀; volume correction gives ωc=ω₀/f.
 * @note BCS uses all selected shells; accuracy is 1e-10.
 */
void PSMSpherical::build_hfbcs(int TargetN_I, int TargetZ_I, double Gn_F, double Gp_F, double epsilon2_F, double epsilon4_F) {
    assert(TargetN_I >= Ncore_I && TargetZ_I >= Zcore_I);
    assert(static_cast<double>(TargetN_I) + TargetZ_I > 0.0);
    assert(std::isfinite(Gn_F) && Gn_F >= 0.0 && std::isfinite(Gp_F) && Gp_F >= 0.0);
    assert(std::isfinite(epsilon2_F) && std::isfinite(epsilon4_F));

    // A = N+Z; e = 2ε₂/3.
    const double A_F = static_cast<double>(TargetN_I) + TargetZ_I;
    const double epsilon20_F = 2.0 * epsilon2_F / 3.0;
    assert(1.0 + 0.5 * epsilon20_F > 0.0 && 1.0 - epsilon20_F > 0.0);

    // f³ = (1/2)∫[-1,1](1−eP₂+2ε₄P₄)^(-3/2)dx / [(1+e/2)√(1−e)].
    const GaussLegendreMeshes volume_meshes(40);
    double volumeIntegral_F = 0.0;
    for (int x_I = 0; x_I < volume_meshes.x_F1D_x.size(); ++x_I) {
        const double x_F = volume_meshes.x_F1D_x(x_I);
        const double shape_F = 1.0 - epsilon20_F * gsl_sf_legendre_Pl(2, x_F) + 2.0 * epsilon4_F * gsl_sf_legendre_Pl(4, x_F);
        assert(shape_F > 0.0);
        volumeIntegral_F += 0.5 * volume_meshes.w_F1D_x(x_I) / std::pow(shape_F, 1.5);
    }
    const double volumeFactor_F = std::cbrt(volumeIntegral_F / ((1.0 + 0.5 * epsilon20_F) * std::sqrt(1.0 - epsilon20_F)));

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

    // (b,{N},κ,μ,G,Nactive) → h_Nilsson → BCS.
    const auto build_species = [&](const SphericalSetting& setting, HFBCS& hfb, int Nactive_I, double G_F, double Akappa_F, double Bkappa_F, double Amu_F, double Bmu_F, const std::array<double, 9>& Fkappa_F1D_Nshell, const std::array<double, 9>& Fmu_F1D_Nshell) {
        assert(setting.useAxialSym_B && setting.useParity_B);
        assert(hfb.Nblock_I == static_cast<int>(setting.labels_S2D_block_bsp.size()));
        assert(Nactive_I >= 0 && Nactive_I <= 2 * static_cast<int>(setting.labels_S1D_sp.size()));

        // PSM: ℏ²/M = 41.4678 MeV fm²; ℏω₀ = ℏ²/(Mb²).
        const double hbarOmega0_F = 41.4678 / (setting.b_F * setting.b_F);
        const double hbarOmegac_F = hbarOmega0_F / volumeFactor_F;
        const SphericalLaguerreBasis radial_basis(setting.b_F, setting.Nr_I, setting.labels_S1D_sp);
        const Eigen::VectorXd r2Weights_F1D_r = radial_basis.w_F1D_r.array() * radial_basis.r_F1D_r.array().square() / (setting.b_F * setting.b_F);

        for (int block_I = 0; block_I < hfb.Nblock_I; ++block_I) {
            const auto& labels = setting.labels_S2D_block_bsp[block_I];
            const auto& indices = setting.indices_I2D_block_bsp[block_I];
            auto& field = hfb.fields[block_I];
            const int Nshell_I = labels.front().N_I;
            const auto shell_iter = std::find(HFBCSData::Nshell_I1D_Nshell.begin(), HFBCSData::Nshell_I1D_Nshell.end(), Nshell_I);
            assert(shell_iter != HFBCSData::Nshell_I1D_Nshell.end());
            const std::size_t shell_I = static_cast<std::size_t>(shell_iter - HFBCSData::Nshell_I1D_Nshell.begin());

            // κ,μ = (Aκ,μ − A Bκ,μ)Fκ,μ(N); ⟨l²⟩_N = N(N+3)/2.
            const double kappa_F = (Akappa_F - A_F * Bkappa_F) * Fkappa_F1D_Nshell[shell_I];
            const double mu_F = (Amu_F - A_F * Bmu_F) * Fmu_F1D_Nshell[shell_I];
            const double l2Average_F = 0.5 * Nshell_I * (Nshell_I + 3.0);
            field.GammaPosPos_F2D_bsp_bsp.setZero();
            field.h0PosPos_F2D_bsp_bsp.setZero();

            for (int bsp1_I = 0; bsp1_I < hfb.Nbsp_I1D_block[block_I]; ++bsp1_I) {
                const auto& label1 = labels[bsp1_I];
                const int sp1_I = indices[bsp1_I];
                assert(label1.N_I == Nshell_I && label1.twom_I == labels.front().twom_I);

                // ⟨1|(r/b)²P_L|2⟩ = ∫r²dr R₁R₂(r/b)² ⟨1|P_L|2⟩.
                for (int bsp2_I = 0; bsp2_I <= bsp1_I; ++bsp2_I) {
                    const auto& label2 = labels[bsp2_I];
                    const int sp2_I = indices[bsp2_I];
                    const double r2_F = (radial_basis.phi_F2D_sp_r.row(sp1_I).array() * radial_basis.phi_F2D_sp_r.row(sp2_I).array() * r2Weights_F1D_r.transpose().array()).sum();
                    const double P2_F = calc_legendre(label1, label2, 2);
                    const double P4_F = calc_legendre(label1, label2, 4);
                    const double deformation_F = hbarOmega0_F * r2_F * (-epsilon20_F * P2_F + epsilon4_F * P4_F);
                    field.h0PosPos_F2D_bsp_bsp(bsp1_I, bsp2_I) = deformation_F;
                    field.h0PosPos_F2D_bsp_bsp(bsp2_I, bsp1_I) = deformation_F;
                }

                // 2l·s = j(j+1)−l(l+1)−3/4.
                const double j_F = 0.5 * label1.twoj_I;
                const double l2_F = label1.l_I * (label1.l_I + 1.0);
                const double twoLS_F = j_F * (j_F + 1.0) - l2_F - 0.75;
                const double spherical_F = hbarOmega0_F * (Nshell_I + 1.5) - kappa_F * hbarOmegac_F * (twoLS_F + mu_F * (l2_F - l2Average_F));
                field.h0PosPos_F2D_bsp_bsp(bsp1_I, bsp1_I) += spherical_F;
            }
        }

        // (h,G,Nactive) → λ,Δ,E,u,v,ρ,κ; no pairing-window cutoff.
        hfb.G_F = G_F;
        hfb.search_lambda(Nactive_I, std::numeric_limits<double>::infinity(), 1.0e-10);
    };

    // Nactive = N−Ncore; Zactive = Z−Zcore.
    build_species(sphericalsetting_neutron, hfb_neutron, TargetN_I - Ncore_I, Gn_F, HFBCSData::Akappa_n_F, HFBCSData::Bkappa_n_F, HFBCSData::Amu_n_F, HFBCSData::Bmu_n_F, HFBCSData::Fkappa_n_F1D_Nshell, HFBCSData::Fmu_n_F1D_Nshell);
    build_species(sphericalsetting_proton, hfb_proton, TargetZ_I - Zcore_I, Gp_F, HFBCSData::Akappa_p_F, HFBCSData::Bkappa_p_F, HFBCSData::Amu_p_F, HFBCSData::Bmu_p_F, HFBCSData::Fkappa_p_F1D_Nshell, HFBCSData::Fmu_p_F1D_Nshell);
}

/**
 * @brief Expand BCS blocks and their time-reversed partners.
 * @math U = diag(fu,ηfu); V = [[0,−fv],[ηfv,0]].
 * @math E = (E₊,E₊); 2K = (2K₊,−2K₊).
 * @output Full neutron and proton U,V,Eqp,TwoK arrays.
 * @note Rows follow spherical labels; columns follow block → eigenstate.
 */
void PSMSpherical::build_UVEK() {
    // Positive-K states first; negative-K partners second.
    const auto expand_UVEK = [](const SphericalSetting& setting, const HFBCS& hfb, Eigen::MatrixXd& U_F2D_2sp_2qp, Eigen::MatrixXd& V_F2D_2sp_2qp, Eigen::VectorXd& Eqp_F1D_2qp, Eigen::VectorXi& TwoK_I1D_2qp) {
        const int Nsp_I = static_cast<int>(setting.labels_S1D_sp.size());
        const int Nqp_I = Nsp_I;
        assert(hfb.Nblock_I == static_cast<int>(setting.indices_I2D_block_bsp.size()));

        // U,V ∈ ℝ^(2Nsp×2Nqp); E,2K ∈ ℝ^(2Nqp).
        U_F2D_2sp_2qp.resize(2 * Nsp_I, 2 * Nqp_I);
        V_F2D_2sp_2qp.resize(2 * Nsp_I, 2 * Nqp_I);
        U_F2D_2sp_2qp.setZero();
        V_F2D_2sp_2qp.setZero();
        Eqp_F1D_2qp.resize(2 * Nqp_I);
        TwoK_I1D_2qp.resize(2 * Nqp_I);

        int qpBegin_I = 0;
        for (int block_I = 0; block_I < hfb.Nblock_I; ++block_I) {
            const auto& indices = setting.indices_I2D_block_bsp[block_I];
            const auto& solution = hfb.solutions[block_I];
            const int Nbqp_I = static_cast<int>(indices.size());
            assert(Nbqp_I > 0 && solution.Eqp_F1D_bhf.size() == Nbqp_I);
            assert(solution.fhf_F2D_sp_bhf.rows() == Nbqp_I && solution.fhf_F2D_sp_bhf.cols() == Nbqp_I);
            const int TwoK_I = setting.labels_S1D_sp[indices.front()].twom_I;

            // E_qp̄ = E_qp; K_qp̄ = −K_qp.
            Eqp_F1D_2qp.segment(qpBegin_I, Nbqp_I) = solution.Eqp_F1D_bhf;
            Eqp_F1D_2qp.segment(qpBegin_I + Nqp_I, Nbqp_I) = solution.Eqp_F1D_bhf;
            TwoK_I1D_2qp.segment(qpBegin_I, Nbqp_I).setConstant(TwoK_I);
            TwoK_I1D_2qp.segment(qpBegin_I + Nqp_I, Nbqp_I).setConstant(-TwoK_I);

            for (int bqp_I = 0; bqp_I < Nbqp_I; ++bqp_I) {
                const int qp_I = qpBegin_I + bqp_I;
                for (int bsp_I = 0; bsp_I < Nbqp_I; ++bsp_I) {
                    const int sp_I = indices[bsp_I];
                    const double eta_F = solution.eta_F1D_bsp(bsp_I);
                    const double fu_F = solution.fhf_F2D_sp_bhf(bsp_I, bqp_I) * solution.uPos_F1D_bhf(bqp_I);
                    const double fv_F = solution.fhf_F2D_sp_bhf(bsp_I, bqp_I) * solution.vNeg_F1D_bhf(bqp_I);

                    // U⁺⁺ = fu; U⁻⁻ = ηfu; V⁺⁻ = −fv; V⁻⁺ = ηfv.
                    U_F2D_2sp_2qp(sp_I, qp_I) = fu_F;
                    U_F2D_2sp_2qp(sp_I + Nsp_I, qp_I + Nqp_I) = eta_F * fu_F;
                    V_F2D_2sp_2qp(sp_I, qp_I + Nqp_I) = -fv_F;
                    V_F2D_2sp_2qp(sp_I + Nsp_I, qp_I) = eta_F * fv_F;
                }
            }
            qpBegin_I += Nbqp_I;
        }
        assert(qpBegin_I == Nqp_I);
    };

    expand_UVEK(sphericalsetting_neutron, hfb_neutron, Un_F2D_2spn_2qpn, Vn_F2D_2spn_2qpn, Eqpn_F1D_2qpn, TwoKn_I1D_2qpn);
    expand_UVEK(sphericalsetting_proton, hfb_proton, Up_F2D_2spp_2qpp, Vp_F2D_2spp_2qpp, Eqpp_F1D_2qpp, TwoKp_I1D_2qpp);
}
