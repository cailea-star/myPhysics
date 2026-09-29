/**
 * @file    psm_spherical_hfbcs.cpp
 * @author  cailea
 * @date    2026-09-27
 * @brief   Nilsson single-particle fields and BCS pairing.
 */

#include "psm_spherical.hpp"

#include <cassert>
#include <cmath>
#include <limits>


/**
 * @brief Assemble Nilsson blocks; solve BCS using Brent.
 * @math h = ℏω₀(N+3/2) − κℏω₀₀[2l·s+μ(l²−N(N+3)/2)]
 * @math     + ℏω₀(r/b)²[−(2ε₂/3)P₂+ε₄P₄]; Γ = 0.
 * @math U = diag(fu,ηfu); V = [[0,−fv],[ηfv,0]].
 * @math E = (E₊,E₊); 2K = (2K₊,−2K₊).
 * @output BCS solutions and full U,V,Eqp,TwoK,ρ arrays.
 * @note Uses stored spherical and deformation fields.
 * @note BCS uses all selected shells; accuracy is 1e-10.
 * @note Rows follow spherical labels; columns follow block → eigenstate.
 */
void PSMSpherical::solve_hfbcs() {
    assert(psm_setting.N_I >= psm_setting.Ncore_I && psm_setting.Z_I >= psm_setting.Zcore_I);
    assert(static_cast<double>(psm_setting.N_I) + psm_setting.Z_I > 0.0);
    assert(std::isfinite(psm_setting.G0_nn_F) && psm_setting.G0_nn_F >= 0.0 && std::isfinite(psm_setting.G0_pp_F) && psm_setting.G0_pp_F >= 0.0);
    assert(std::isfinite(psm_setting.epsilon2_F) && std::isfinite(psm_setting.epsilon4_F));

    assert(psm_setting.hbarOmega0_n_F > 0.0 && psm_setting.hbarOmega0_p_F > 0.0);

    // (h₀,h_deform,G,Nactive) → h_Nilsson → BCS.
    const auto build_species = [](const SphericalSetting& sphericalsetting, HFBCS& hfb, int Nactive_I, double G_F, const Eigen::MatrixXd& h0_F2D_2sp_2sp, const Eigen::MatrixXd& hDeform_F2D_2sp_2sp) {
        assert(sphericalsetting.useAxialSym_B && sphericalsetting.useParity_B);
        assert(hfb.Nblock_I == static_cast<int>(sphericalsetting.labels_S2D_block_bsp.size()));
        assert(Nactive_I >= 0 && Nactive_I <= 2 * static_cast<int>(sphericalsetting.labels_S1D_sp.size()));

        // h_Nilsson = h₀ + h_deform; Γ = 0.
        for (int block_I = 0; block_I < hfb.Nblock_I; ++block_I) {
            const auto& indices = sphericalsetting.indices_I2D_block_bsp[block_I];
            auto& field = hfb.fields[block_I];
            field.GammaPosPos_F2D_bsp_bsp.setZero();
            for (int bsp1_I = 0; bsp1_I < hfb.Nbsp_I1D_block[block_I]; ++bsp1_I) {
                const int sp1_I = indices[bsp1_I];
                for (int bsp2_I = 0; bsp2_I < hfb.Nbsp_I1D_block[block_I]; ++bsp2_I) {
                    const int sp2_I = indices[bsp2_I];
                    field.h0PosPos_F2D_bsp_bsp(bsp1_I, bsp2_I) = h0_F2D_2sp_2sp(sp1_I, sp2_I) + hDeform_F2D_2sp_2sp(sp1_I, sp2_I);
                }
            }
        }

        // (h,G,Nactive) → λ,Δ,E,u,v,ρ,κ; no pairing-window cutoff.
        hfb.G_F = G_F;
        hfb.search_lambda(Nactive_I, std::numeric_limits<double>::infinity(), 1.0e-10);
    };

    // Nactive = N−Ncore; Zactive = Z−Zcore.
    build_species(psm_setting.sphericalsetting_neutron, hfb_neutron, psm_setting.N_I - psm_setting.Ncore_I, psm_setting.G0_nn_F, h0n_F2D_2spn_2spn, h0Deform_n_F2D_2spn_2spn);
    build_species(psm_setting.sphericalsetting_proton, hfb_proton, psm_setting.Z_I - psm_setting.Zcore_I, psm_setting.G0_pp_F, h0p_F2D_2spp_2spp, h0Deform_p_F2D_2spp_2spp);

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

    expand_UVEK(psm_setting.sphericalsetting_neutron, hfb_neutron, Un_F2D_2spn_2qpn, Vn_F2D_2spn_2qpn, Eqpn_F1D_2qpn, TwoKn_I1D_2qpn);
    expand_UVEK(psm_setting.sphericalsetting_proton, hfb_proton, Up_F2D_2spp_2qpp, Vp_F2D_2spp_2qpp, Eqpp_F1D_2qpp, TwoKp_I1D_2qpp);

    // ρ = VVᵀ.
    rhon_F2D_2spn_2spn.noalias() = Vn_F2D_2spn_2qpn * Vn_F2D_2spn_2qpn.transpose();
    rhop_F2D_2spp_2spp.noalias() = Vp_F2D_2spp_2qpp * Vp_F2D_2spp_2qpp.transpose();
}
