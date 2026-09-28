/**
 * @file    psm_spherical.hpp
 * @author  cailea
 * @date    2026-09-27
 * @brief   Projected shell model in a spherical basis.
 */

#pragma once

#include <array>
#include <vector>

#include "hfb_bcs.hpp"
#include "hfb_configuration.hpp"
#include "hfb_projection_nucleus.hpp"
#include "spherical_setting.hpp"

// κ_q(N,A) = (Aκ_q − A Bκ_q) Fκ_q(N).
// μ_q(N,A) = (Aμ_q − A Bμ_q) Fμ_q(N).
// q ∈ {n,p}; A: mass number; N: oscillator shell.
struct HFBCSData {
    // κ_q(A) = Aκ_q − A Bκ_q; μ_q(A) = Aμ_q − A Bμ_q.
    inline static constexpr double Akappa_n_F = 0.06433;
    inline static constexpr double Bkappa_n_F = 0.004e-3;
    inline static constexpr double Akappa_p_F = 0.07441;
    inline static constexpr double Bkappa_p_F = 0.068e-3;
    inline static constexpr double Amu_n_F = 0.59010;
    inline static constexpr double Bmu_n_F = 1.080e-3;
    inline static constexpr double Amu_p_F = 0.51180;
    inline static constexpr double Bmu_p_F = -0.560e-3;

    inline static constexpr std::array<int, 9> Nshell_I1D_Nshell{0, 1, 2, 3, 4, 5, 6, 7, 8};

    // Fκ(N) = Fμ(N) = 1.
    inline static constexpr std::array<double, 9> Fkappa_n_F1D_Nshell{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    inline static constexpr std::array<double, 9> Fkappa_p_F1D_Nshell{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    inline static constexpr std::array<double, 9> Fmu_n_F1D_Nshell{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    inline static constexpr std::array<double, 9> Fmu_p_F1D_Nshell{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
};

class PSMSpherical {
public:
    int Ncore_I = 0;
    int Zcore_I = 0;

    SphericalSetting sphericalsetting_neutron;
    SphericalSetting sphericalsetting_proton;

    HFBCS hfb_neutron;
    HFBCS hfb_proton;

    // Complete spaces: (+K,−K).
    Eigen::MatrixXd Un_F2D_2spn_2qpn{};
    Eigen::MatrixXd Vn_F2D_2spn_2qpn{};
    Eigen::MatrixXd Up_F2D_2spp_2qpp{};
    Eigen::MatrixXd Vp_F2D_2spp_2qpp{};

    Eigen::VectorXd Eqpn_F1D_2qpn{};
    Eigen::VectorXd Eqpp_F1D_2qpp{};
    Eigen::VectorXi TwoKn_I1D_2qpn{};
    Eigen::VectorXi TwoKp_I1D_2qpp{};

    HFBConfiguration configuration_neutron;
    HFBConfiguration configuration_proton;

    HFBProjectionNucleus projection_nucleus;

    // (cfgp,cfgn,K,cfgp,cfgn,K).
    Eigen::Tensor<doubleC, 6, Eigen::ColMajor> Honebody_C6D_cfgp_cfgn_K_cfgp_cfgn_K{};
    Eigen::Tensor<doubleC, 6, Eigen::ColMajor> Htwobody_C6D_cfgp_cfgn_K_cfgp_cfgn_K{};
    Eigen::Tensor<doubleC, 6, Eigen::ColMajor> N_C6D_cfgp_cfgn_K_cfgp_cfgn_K{};

public:
    /**
     * @brief Initialize species bases and BCS workspaces.
     * @math (bₙ,{Nₙ}; bₚ,{Nₚ}) → blocks (m,N), η → HFBCS.
     * @math Ncore = Nmin(Nmin+1)(Nmin+2)/3; similarly Zcore.
     * @output Species settings, core counts, and allocated BCS workspaces.
     * @note Shells below each active space are filled.
     */
    PSMSpherical(double bn_F_, double bp_F_, const std::vector<int>& Nshelln_I1D_Nshell_, const std::vector<int>& Nshellp_I1D_Nshell_) : 
    sphericalsetting_neutron(bn_F_, Nshelln_I1D_Nshell_, true, true, true), 
    sphericalsetting_proton(bp_F_, Nshellp_I1D_Nshell_, true, true, true), 
    hfb_neutron(sphericalsetting_neutron.eta_F2D_block_bsp), 
    hfb_proton(sphericalsetting_proton.eta_F2D_block_bsp) {
        const int Nminn_I = sphericalsetting_neutron.Nshell_I1D_Nshell.front();
        const int Nminp_I = sphericalsetting_proton.Nshell_I1D_Nshell.front();

        // Ncore = Σ_{N=0}^{Nmin−1}(N+1)(N+2).
        Ncore_I = Nminn_I * (Nminn_I + 1) * (Nminn_I + 2) / 3;
        Zcore_I = Nminp_I * (Nminp_I + 1) * (Nminp_I + 2) / 3;
    }

    /**
     * @brief Build Nilsson states and solve BCS pairing.
     * @math h_Nilsson(ε₂,ε₄), Gₙ,Gₚ → λ,Δ,E,u,v.
     * @output Updated neutron and proton HFBCS solutions.
     * @note TargetN,TargetZ include cores; G is in MeV.
     */
    void build_hfbcs(int TargetN_I, int TargetZ_I, double Gn_F, double Gp_F, double epsilon2_F, double epsilon4_F);

    /**
     * @brief Expand BCS blocks into complete quasiparticle spaces.
     * @math (f,u,v,η) → U,V; E = (E₊,E₊); K = (K₊,−K₊).
     * @output Full neutron and proton U,V,Eqp,TwoK arrays.
     * @note Requires solved BCS states; ordering is (+K,−K).
     */
    void build_UVEK();

    /**
     * @brief Set neutron cutoffs and enumerate configurations.
     * @math Nqp = 2Σ_b Nbsp_b; E = (E₊,E₊); K = (K₊,−K₊).
     * @output Updated neutron cutoffs and configuration list.
     * @note Requires build_UVEK(); QP order is (+K,−K).
     */
    void build_config_neutron(const Eigen::VectorXi& Ncqp_I1D_Ncqp, const Eigen::VectorXd& ECut_F1D_Ncqp, const Eigen::VectorXi& NCut_I1D_Ncqp, const Eigen::VectorXi& TwoKCut_I1D_Ncqp);

    /**
     * @brief Set proton cutoffs and enumerate configurations.
     * @math Nqp = 2Σ_b Nbsp_b; E = (E₊,E₊); K = (K₊,−K₊).
     * @output Updated proton cutoffs and configuration list.
     * @note Requires build_UVEK(); QP order is (+K,−K).
     */
    void build_config_proton(const Eigen::VectorXi& Ncqp_I1D_Ncqp, const Eigen::VectorXd& ECut_F1D_Ncqp, const Eigen::VectorXi& NCut_I1D_Ncqp, const Eigen::VectorXi& TwoKCut_I1D_Ncqp);

    /**
     * @brief Assemble neutron-proton projection workspaces.
     * @math Nactive = TargetN − Ncore; Zactive = TargetZ − Zcore.
     * @output Initialized projection_nucleus with shared left/right configurations.
     * @note Requires U,V and neutron/proton configurations.
     */
    void build_projection(int TargetN_I, int TargetZ_I, int TargetTwoI_I, int Nalpha_I, int Nbeta_I, int Ngamma_I, int Nphin_I, int Nphip_I);

    /**
     * @brief Build the projected configuration overlap tensor.
     * @math N_ab = ⟨Φ_a|Pₙ Pₚ Pᴵ|Φ_b⟩.
     * @output N_C6D_cfgp_cfgn_K_cfgp_cfgn_K; cfgp varies fastest.
     * @note Requires build_projection().
     */
    void build_norm();

    /**
     * @brief Construct spherical single-particle fields and project them.
     * @math H¹ = ⟨Φ₁|(h₀ₙ+h₀ₚ)PₙPₚPᴵ|Φ₂⟩.
     * @output Honebody_C6D_cfgp_cfgn_K_cfgp_cfgn_K.
     * @note TargetN,TargetZ include cores; requires build_projection().
     */
    void build_onebody(int TargetN_I, int TargetZ_I);

    /**
     * @brief Construct quadrupole and pairing operators and project them.
     * @math (χ₂,nn,χ₂,np,χ₂,pp,G₀,nn,G₀,pp,G₂,nn,G₂,pp) → H².
     * @output Htwobody_C6D_cfgp_cfgn_K_cfgp_cfgn_K.
     * @note Independent input strengths; requires build_projection().
     */
    void build_twobody(double chi2_nn_F, double chi2_np_F, double chi2_pp_F, double G0_nn_F, double G0_pp_F, double G2_nn_F, double G2_pp_F);
};
