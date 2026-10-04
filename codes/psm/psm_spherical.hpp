/**
 * @file    psm_spherical.hpp
 * @author  cailea
 * @date    2026-09-27
 * @brief   Projected shell model in a spherical basis.
 */

#pragma once

#include <utility>
#include <unsupported/Eigen/CXX11/Tensor>

#include "hfb_bcs.hpp"
#include "psm_spherical_setting.hpp"
#include "hfb_projection_nucleus.hpp"

class PSMSpherical {
public:
    PSMSphericalSetting psm_setting;

    // h₀ in complete (+m,−m) space.
    Eigen::MatrixXd h0n_F2D_2spn_2spn{};
    Eigen::MatrixXd h0p_F2D_2spp_2spp{};

    // h_deform = ℏω₀(r/b)²[−(2ε₂/3)P₂+ε₄P₄].
    Eigen::MatrixXd h0Deform_n_F2D_2spn_2spn{};
    Eigen::MatrixXd h0Deform_p_F2D_2spp_2spp{};

    // Q₂μ = δ_NN′(r/b)²√(4π/5)Y₂μ; μ = −2,…,2; (+m,−m).
    Eigen::Tensor<double, 3, Eigen::ColMajor> Qn_F3D_2spn_2spn_mu{};
    Eigen::Tensor<double, 3, Eigen::ColMajor> Qp_F3D_2spp_2spp_mu{};

    // P₀: μ = 0; P₂μ: μ = −2,…,2.
    Eigen::Tensor<double, 3, Eigen::ColMajor> P0n_F3D_2spn_2spn_mu{};
    Eigen::Tensor<double, 3, Eigen::ColMajor> P0p_F3D_2spp_2spp_mu{};
    Eigen::Tensor<double, 3, Eigen::ColMajor> P2n_F3D_2spn_2spn_mu{};
    Eigen::Tensor<double, 3, Eigen::ColMajor> P2p_F3D_2spp_2spp_mu{};

    HFBCS hfb_neutron;
    HFBCS hfb_proton;

    // Complete spaces: (+K,−K).
    Eigen::MatrixXd Un_F2D_2spn_2qpn{};
    Eigen::MatrixXd Vn_F2D_2spn_2qpn{};
    Eigen::MatrixXd Up_F2D_2spp_2qpp{};
    Eigen::MatrixXd Vp_F2D_2spp_2qpp{};

    // ρ = VVᵀ; complete (+K,−K) space.
    Eigen::MatrixXd rhon_F2D_2spn_2spn{};
    Eigen::MatrixXd rhop_F2D_2spp_2spp{};

    Eigen::VectorXd Eqpn_F1D_2qpn{};
    Eigen::VectorXd Eqpp_F1D_2qpp{};
    Eigen::VectorXi TwoKn_I1D_2qpn{};
    Eigen::VectorXi TwoKp_I1D_2qpp{};

    // 2K_cfg = Σqp 2K_qp.
    Eigen::VectorXi TwoKn_I1D_cfgn{};
    Eigen::VectorXi TwoKp_I1D_cfgp{};

    HFBProjectionNucleus projection_nucleus;

    // (cfgn,cfgn,cfgp,cfgp,β).
    Eigen::Tensor<doubleC, 5, Eigen::ColMajor> H_C5D_cfgn_cfgn_cfgp_cfgp_beta{};
    Eigen::Tensor<doubleC, 5, Eigen::ColMajor> N_C5D_cfgn_cfgn_cfgp_cfgp_beta{};

    // Xci(r,ν): natural orthonormal basis; Xci†Xci = 1.
    Eigen::VectorXd Eci_F1D_eigenH{};
    Eigen::MatrixXcd Xci_C2D_eigenN_eigenH{};

public:
    /**
     * @brief Store settings, construct operators, and initialize BCS workspaces.
     * @math psm_setting → h₀,h_deform,Q₂μ,P₀,P₂μ,ηₙ,ηₚ → HFBCS.
     * @output Owned settings, ten operators, and neutron/proton BCS workspaces.
     * @note Requires configured single-particle spaces and HFBCS parameters.
     * @note Settings are passed by value, then moved.
     */
    PSMSpherical(PSMSphericalSetting psm_setting_) : psm_setting(std::move(psm_setting_)),
    hfb_neutron(psm_setting.sphericalsetting_neutron.eta_F2D_block_bsp),
    hfb_proton(psm_setting.sphericalsetting_proton.eta_F2D_block_bsp) {
        // psm_setting → h₀,h_deform,Q₂μ,P₀,P₂μ in complete (+m,−m) spaces.
        build_h0(psm_setting.sphericalsetting_neutron, psm_setting.hbarOmega0_n_F, psm_setting.hbarOmega00_n_F, psm_setting.kappa_n_F, psm_setting.mu_n_F, h0n_F2D_2spn_2spn);
        build_h0(psm_setting.sphericalsetting_proton, psm_setting.hbarOmega0_p_F, psm_setting.hbarOmega00_p_F, psm_setting.kappa_p_F, psm_setting.mu_p_F, h0p_F2D_2spp_2spp);
        build_h0_deformation(psm_setting.sphericalsetting_neutron, psm_setting.hbarOmega0_n_F, psm_setting.epsilon2_F, psm_setting.epsilon4_F, h0Deform_n_F2D_2spn_2spn);
        build_h0_deformation(psm_setting.sphericalsetting_proton, psm_setting.hbarOmega0_p_F, psm_setting.epsilon2_F, psm_setting.epsilon4_F, h0Deform_p_F2D_2spp_2spp);
        build_Q(psm_setting.sphericalsetting_neutron, Qn_F3D_2spn_2spn_mu);
        build_Q(psm_setting.sphericalsetting_proton, Qp_F3D_2spp_2spp_mu);
        build_P(psm_setting.sphericalsetting_neutron, Qn_F3D_2spn_2spn_mu, P0n_F3D_2spn_2spn_mu, P2n_F3D_2spn_2spn_mu);
        build_P(psm_setting.sphericalsetting_proton, Qp_F3D_2spp_2spp_mu, P0p_F3D_2spp_2spp_mu, P2p_F3D_2spp_2spp_mu);
    }

    /**
     * @brief Build Nilsson states and solve BCS pairing.
     * @math h_Nilsson(ε₂,ε₄), Gₙ,Gₚ → λ,Δ,E,u,v.
     * @math (f,u,v,η) → U,V; E = (E₊,E₊); K = (K₊,−K₊).
     * @output HFBCS solutions and full U,V,Eqp,TwoK,ρ arrays.
     * @note Reads psm_setting; BCS uses G0_nn,G0_pp in MeV.
     * @note Full-space ordering is (+K,−K).
     */
    void solve_hfbcs();

    /**
     * @brief Build CI angle kernels using gauge quadrature.
     * @math G₂ = γG₀; H = h₀ − λₙN̂ − λₚẐ + H_QQ + H_PP.
     * @output N_C5D and H_C5D beta caches.
     * @note Requires solved HFBCS and initialized projection.
     */
    void build_ci(double gamma_F);

    /**
     * @brief Solve configuration mixing using canonical norm orthogonalization.
     * @math Hᴵf = ENᴵf; Xci†Xci = 1.
     * @output Eci_F1D_eigenH and Xci_C2D_eigenN_eigenH.
     * @note Requires build_ci(); retains nᵢ > 1e-10 nmax.
     */
    void solve_ci(int TargetTwoI_I);

    /**
     * @brief Assemble neutron-proton projection workspaces.
     * @math Nactive = N − Ncore; Zactive = Z − Zcore.
     * @output Generated configurations and initialized projection_nucleus.
     * @note Requires solved U,V,Eqp,TwoK and configured cutoffs.
     * @note Axial projection: Nalpha = Ngamma = 1.
     */
    void build_projection(int Nbeta_I, int Nphin_I, int Nphip_I);

private:
    /**
     * @brief Build spherical single-particle fields.
     * @math h₀ = ℏω₀(N+3/2) − κℏω₀₀[2l·s+μ(l²−N(N+3)/2)].
     * @output Diagonal h₀ in complete (+m,−m) space.
     */
    static void build_h0(const SphericalSetting& setting, double hbarOmega0_F, double hbarOmega00_F, double kappa_F, double mu_F, Eigen::MatrixXd& h0_F2D_2sp_2sp);

    /**
     * @brief Build deformation fields using radial quadrature.
     * @math h_deform = ℏω₀(r/b)²[−(2ε₂/3)P₂+ε₄P₄].
     * @output Deformation contribution in complete (+m,−m) space.
     * @note Retains existing shell and axial blocks.
     */
    static void build_h0_deformation(const SphericalSetting& setting, double hbarOmega0_F, double epsilon2_F, double epsilon4_F, Eigen::MatrixXd& hDeform_F2D_2sp_2sp);

    /**
     * @brief Build quadrupole matrix using radial quadrature.
     * @math Q₂μ = δ_NN′(r/b)²√(4π/5)Y₂μ; μ = −2,…,2.
     * @output Five Q₂μ matrices in complete (+m,−m) space.
     */
    static void build_Q(const SphericalSetting& sphericalsetting, Eigen::Tensor<double, 3, Eigen::ColMajor>& Q_F3D_2sp_2sp_mu);

    /**
     * @brief Build pairing matrices from quadrupole matrices.
     * @math P₀ = [[0,η],[−η,0]]; P₂μ = Q₂μ P₀.
     * @output P₀ and five P₂μ matrices.
     */
    static void build_P(const SphericalSetting& setting, const Eigen::Tensor<double, 3, Eigen::ColMajor>& Q2_F3D_2sp_2sp_mu, Eigen::Tensor<double, 3, Eigen::ColMajor>& P0_F3D_2sp_2sp_mu, Eigen::Tensor<double, 3, Eigen::ColMajor>& P2_F3D_2sp_2sp_mu);
};
