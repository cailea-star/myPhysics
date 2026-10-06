/**
 * @file    psm_spherical_setting.hpp
 * @author  cailea
 * @date    2026-09-29
 * @brief   Spherical model spaces and quasiparticle configurations.
 */

#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>
#include <unsupported/Eigen/CXX11/Tensor>

#include "hfb_configuration_nucleus.hpp"
#include "spherical_setting.hpp"

// 
inline constexpr double hbar2_m_F = 41.4678; // ℏ²/m [MeV fm²]; PSM convention.

// κ_q(A) = Aκ_q − A Bκ_q.
// μ_q(A) = Aμ_q − A Bμ_q.
// q ∈ {n,p}; A = N+Z.
struct PSMParametersOneBody {
    // κ_q(A) = Aκ_q − A Bκ_q; μ_q(A) = Aμ_q − A Bμ_q.
    inline static constexpr double Akappa_n_F = 0.06433;
    inline static constexpr double Bkappa_n_F = 0.004e-3;
    inline static constexpr double Akappa_p_F = 0.07441;
    inline static constexpr double Bkappa_p_F = 0.068e-3;
    inline static constexpr double Amu_n_F = 0.59010;
    inline static constexpr double Bmu_n_F = 1.080e-3;
    inline static constexpr double Amu_p_F = 0.51180;
    inline static constexpr double Bmu_p_F = -0.560e-3;
};

struct PSMParametersPair {
    // G₀,nn/pp = [g₁ ∓ g₂(N−Z)/A]/A.
    inline static constexpr double g1_F = 20.12;
    inline static constexpr double g2_F = 13.13;
};

class PSMSphericalSetting {
public:
    // N,Z include cores; A = N+Z.
    int N_I = 0;
    int Z_I = 0;

    int Ncore_I = 0;
    int Zcore_I = 0;

    // h_HO = ℏω₀(N_shell+3/2).
    double hbarOmega0_n_F = 0.0;
    double hbarOmega0_p_F = 0.0;

    // h₀ = h_HO − κℏω₀₀[2l·s+μ(l²−⟨l²⟩)], ⟨l²⟩ = N_shell(N_shell+3)/2; l,s in ℏ units.
    double hbarOmega00_n_F = 0.0;
    double hbarOmega00_p_F = 0.0;
    double kappa_n_F = 0.0;
    double kappa_p_F = 0.0;
    double mu_n_F = 0.0;
    double mu_p_F = 0.0;

    // h_deform = ℏω₀(r/b)²[−(2ε₂/3)P₂(cosθ)+ε₄P₄(cosθ)].
    double epsilon2_F = 0.0;
    double epsilon4_F = 0.0;

    // H_QQ = −χnn:Qn†·Qn:/2 − χnp Qn†·Qp − χpp:Qp†·Qp:/2.
    double chi2_nn_F = 0.0;
    double chi2_np_F = 0.0;
    double chi2_pp_F = 0.0;

    // H_pair = −Σ_{λ=0,2} Σ_{q=n,p} Gλ,qq Pλ,q†·Pλ,q.
    double G0_nn_F = 0.0;
    double G0_pp_F = 0.0;
    double G2_nn_F = 0.0;
    double G2_pp_F = 0.0;

    SphericalSetting sphericalsetting_neutron;
    SphericalSetting sphericalsetting_proton;

    // {N_shell}_active ⊆ {N_shell}_basis.
    std::vector<int> Nshelln_I1D_Nactive{};
    std::vector<int> Nshellp_I1D_Nactive{};

    HFBConfigurationNucleus configuration_nucleus;

public:
    /**
     * @brief Construct empty PSM settings.
     * @math {Nₙ} = {Nₚ} = ∅.
     * @output Empty species spaces and zero parameters.
     */
    PSMSphericalSetting() = default;

    /**
     * @brief Set neutron space and reset nuclear configurations.
     * @math Ncore = Nmin(Nmin+1)(Nmin+2)/3; ℏω₀ = 41.4678/b².
     * @output Neutron basis, active shells, core count, frequency; empty configurations.
     * @note Empty active shells select all basis shells.
     */
    void set_sp_neutron(double b_F_, const std::vector<int>& Nshell_I1D_Nshell_, const std::vector<int>& Nshell_I1D_Nactive_ = {});

    /**
     * @brief Set proton space and reset nuclear configurations.
     * @math Zcore = Nmin(Nmin+1)(Nmin+2)/3; ℏω₀ = 41.4678/b².
     * @output Proton basis, active shells, core count, frequency; empty configurations.
     * @note Empty active shells select all basis shells.
     */
    void set_sp_proton(double b_F_, const std::vector<int>& Nshell_I1D_Nshell_, const std::vector<int>& Nshell_I1D_Nactive_ = {});

    /**
     * @brief Set nuclear composition, deformation, and single-particle parameters.
     * @math A = N+Z; ℏω₀ = 41.4678/b²; κ,μ = Aκ,μ − A Bκ,μ.
     * @math ℏω₀₀,n/p = 41.4678 A^(-1/3)(2N/A,2Z/A)^(1/3).
     * @output N,Z,ε₂,ε₄,ℏω₀,ℏω₀₀,κ,μ and monopole pairing strengths.
     * @note N,Z include cores.
     */
    void set_hfbcs(int N_I_, int Z_I_, double epsilon2_F_, double epsilon4_F_);

    /**
     * @brief Set quadrupole strengths from density contractions.
     * @math q = Tr(Q₂₀ρ); D = ℏω₀ₙqₙ+ℏω₀ₚqₚ; χ₂,ab = (2ε₂/3)ℏω₀aℏω₀b/D.
     * @output Updated χ₂,nn,χ₂,np,χ₂,pp.
     * @note Requires set_hfbcs(), full (+m,−m) densities and Q; D ≠ 0.
     */
    void set_ci_QQ(const Eigen::MatrixXd& rhon_F2D_2spn_2spn_, const Eigen::MatrixXd& rhop_F2D_2spp_2spp_, const Eigen::Tensor<double, 3, Eigen::ColMajor>& Qn_F3D_2spn_2spn_mu_, const Eigen::Tensor<double, 3, Eigen::ColMajor>& Qp_F3D_2spp_2spp_mu_);

    /**
     * @brief Set quadrupole pairing strengths.
     * @math G₂,nn = γG₀,nn; G₂,pp = γG₀,pp.
     * @output Updated G₂,nn,G₂,pp.
     * @note Requires set_hfbcs().
     */
    void set_ci_PP(double gamma_F);


    /**
     * @brief Store mixed nuclear configuration cutoffs.
     * @math (Ncqp_n,Ncqp_p,ECut,NCut,TwoKCut) → configuration_nucleus.
     * @output Stored cutoffs; empty mixed and tensor configurations.
     */
    void set_ci_cut(const Eigen::MatrixXi& Ncqp_I2D_Ncqp_np, const Eigen::VectorXd& ECut_F1D_Ncqp, const Eigen::VectorXi& NCut_I1D_Ncqp, const Eigen::VectorXi& TwoKCut_I1D_Ncqp);

    /**
     * @brief Select Kramers configurations using Nilsson-energy ordering.
     * @math Ehf → representative order; Ecfg = ΣEqp.
     * @output Mixed nuclear configurations and unique tensor factors.
     * @note Energy and K arrays follow full U,V column order.
     */
    void build_ci_config(const Eigen::VectorXd& Ehfn_F1D_2qpn, const Eigen::VectorXd& Eqpn_F1D_2qpn, const Eigen::VectorXi& TwoKn_I1D_2qpn, const Eigen::VectorXd& Ehfp_F1D_2qpp, const Eigen::VectorXd& Eqpp_F1D_2qpp, const Eigen::VectorXi& TwoKp_I1D_2qpp);

};

inline void PSMSphericalSetting::set_sp_neutron(double b_F_, const std::vector<int>& Nshell_I1D_Nshell_, const std::vector<int>& Nshell_I1D_Nactive_) {
    sphericalsetting_neutron = SphericalSetting(b_F_, Nshell_I1D_Nshell_, true, true, true);
    Nshelln_I1D_Nactive = Nshell_I1D_Nactive_;
    if (Nshelln_I1D_Nactive.empty()) { Nshelln_I1D_Nactive = Nshell_I1D_Nshell_; }
    assert(std::is_sorted(Nshelln_I1D_Nactive.begin(), Nshelln_I1D_Nactive.end()) && std::adjacent_find(Nshelln_I1D_Nactive.begin(), Nshelln_I1D_Nactive.end()) == Nshelln_I1D_Nactive.end() && std::includes(Nshell_I1D_Nshell_.begin(), Nshell_I1D_Nshell_.end(), Nshelln_I1D_Nactive.begin(), Nshelln_I1D_Nactive.end()));
    const int Nmin_I = sphericalsetting_neutron.Nshell_I1D_Nshell.front();

    // Ncore = Σ_{N=0}^{Nmin−1}(N+1)(N+2); ℏω₀ = ℏ²/(mb²).
    Ncore_I = Nmin_I * (Nmin_I + 1) * (Nmin_I + 2) / 3;
    hbarOmega0_n_F = hbar2_m_F / (b_F_ * b_F_);
    configuration_nucleus.confign_I2D_cfgnp_cqpn.clear();
    configuration_nucleus.configp_I2D_cfgnp_cqpp.clear();
    configuration_nucleus.confign_I2D_tensorcfgn_cqpn.clear();
    configuration_nucleus.configp_I2D_tensorcfgp_cqpp.clear();
}

inline void PSMSphericalSetting::set_sp_proton(double b_F_, const std::vector<int>& Nshell_I1D_Nshell_, const std::vector<int>& Nshell_I1D_Nactive_) {
    sphericalsetting_proton = SphericalSetting(b_F_, Nshell_I1D_Nshell_, true, true, true);
    Nshellp_I1D_Nactive = Nshell_I1D_Nactive_;
    if (Nshellp_I1D_Nactive.empty()) { Nshellp_I1D_Nactive = Nshell_I1D_Nshell_; }
    assert(std::is_sorted(Nshellp_I1D_Nactive.begin(), Nshellp_I1D_Nactive.end()) && std::adjacent_find(Nshellp_I1D_Nactive.begin(), Nshellp_I1D_Nactive.end()) == Nshellp_I1D_Nactive.end() && std::includes(Nshell_I1D_Nshell_.begin(), Nshell_I1D_Nshell_.end(), Nshellp_I1D_Nactive.begin(), Nshellp_I1D_Nactive.end()));
    const int Nmin_I = sphericalsetting_proton.Nshell_I1D_Nshell.front();

    // Zcore = Σ_{N=0}^{Nmin−1}(N+1)(N+2); ℏω₀ = ℏ²/(mb²).
    Zcore_I = Nmin_I * (Nmin_I + 1) * (Nmin_I + 2) / 3;
    hbarOmega0_p_F = hbar2_m_F / (b_F_ * b_F_);
    configuration_nucleus.confign_I2D_cfgnp_cqpn.clear();
    configuration_nucleus.configp_I2D_cfgnp_cqpp.clear();
    configuration_nucleus.confign_I2D_tensorcfgn_cqpn.clear();
    configuration_nucleus.configp_I2D_tensorcfgp_cqpp.clear();
}

inline void PSMSphericalSetting::set_hfbcs(int N_I_, int Z_I_, double epsilon2_F_, double epsilon4_F_) {
    assert(N_I_ >= Ncore_I && Z_I_ >= Zcore_I);
    assert(std::isfinite(epsilon2_F_) && std::isfinite(epsilon4_F_));
    const double A_F = static_cast<double>(N_I_) + Z_I_;
    assert(A_F > 0.0);

    // (N,Z,ε₂,ε₄) → stored inputs.
    N_I = N_I_;
    Z_I = Z_I_;
    epsilon2_F = epsilon2_F_;
    epsilon4_F = epsilon4_F_;

    // ℏω₀ = ℏ²/(m b²).
    hbarOmega0_n_F = hbar2_m_F / (sphericalsetting_neutron.b_F * sphericalsetting_neutron.b_F);
    hbarOmega0_p_F = hbar2_m_F / (sphericalsetting_proton.b_F * sphericalsetting_proton.b_F);

    // ℏω₀₀,n/p = (ℏ²/m) A^(-1/3)(2N/A,2Z/A)^(1/3).
    hbarOmega00_n_F = hbar2_m_F * std::cbrt(2.0 * N_I / A_F) / std::cbrt(A_F);
    hbarOmega00_p_F = hbar2_m_F * std::cbrt(2.0 * Z_I / A_F) / std::cbrt(A_F);

    // κ = Aκ − A Bκ;
    kappa_n_F = PSMParametersOneBody::Akappa_n_F - A_F * PSMParametersOneBody::Bkappa_n_F;
    kappa_p_F = PSMParametersOneBody::Akappa_p_F - A_F * PSMParametersOneBody::Bkappa_p_F;

    // μ = Aμ − A Bμ;
    mu_n_F = PSMParametersOneBody::Amu_n_F - A_F * PSMParametersOneBody::Bmu_n_F;
    mu_p_F = PSMParametersOneBody::Amu_p_F - A_F * PSMParametersOneBody::Bmu_p_F;

    // δ = (N−Z)/A; G₀,nn/pp = (g₁ ∓ g₂δ)/A.
    const double delta_F = (static_cast<double>(N_I) - Z_I) / A_F;
    G0_nn_F = (PSMParametersPair::g1_F - PSMParametersPair::g2_F * delta_F) / A_F;
    G0_pp_F = (PSMParametersPair::g1_F + PSMParametersPair::g2_F * delta_F) / A_F;
}

inline void PSMSphericalSetting::set_ci_QQ(const Eigen::MatrixXd& rhon_F2D_2spn_2spn_, const Eigen::MatrixXd& rhop_F2D_2spp_2spp_, const Eigen::Tensor<double, 3, Eigen::ColMajor>& Qn_F3D_2spn_2spn_mu_, const Eigen::Tensor<double, 3, Eigen::ColMajor>& Qp_F3D_2spp_2spp_mu_) {
    assert(hbarOmega0_n_F > 0.0 && hbarOmega0_p_F > 0.0);

    assert(Qn_F3D_2spn_2spn_mu_.dimension(0) > 0 && Qn_F3D_2spn_2spn_mu_.dimension(2) == 5);
    assert(Qp_F3D_2spp_2spp_mu_.dimension(0) > 0 && Qp_F3D_2spp_2spp_mu_.dimension(2) == 5);

    // ColMajor: μ = 0 occupies the third contiguous matrix.
    const Eigen::Map<const Eigen::MatrixXd> Q20n_F2D_2spn_2spn(Qn_F3D_2spn_2spn_mu_.data() + 2 * Qn_F3D_2spn_2spn_mu_.dimension(0) * Qn_F3D_2spn_2spn_mu_.dimension(1), Qn_F3D_2spn_2spn_mu_.dimension(0), Qn_F3D_2spn_2spn_mu_.dimension(1));
    const Eigen::Map<const Eigen::MatrixXd> Q20p_F2D_2spp_2spp(Qp_F3D_2spp_2spp_mu_.data() + 2 * Qp_F3D_2spp_2spp_mu_.dimension(0) * Qp_F3D_2spp_2spp_mu_.dimension(1), Qp_F3D_2spp_2spp_mu_.dimension(0), Qp_F3D_2spp_2spp_mu_.dimension(1));
    assert(rhon_F2D_2spn_2spn_.rows() == Q20n_F2D_2spn_2spn.rows() && rhon_F2D_2spn_2spn_.cols() == Q20n_F2D_2spn_2spn.cols());
    assert(rhop_F2D_2spp_2spp_.rows() == Q20p_F2D_2spp_2spp.rows() && rhop_F2D_2spp_2spp_.cols() == Q20p_F2D_2spp_2spp.cols());
    assert(rhon_F2D_2spn_2spn_.allFinite() && rhop_F2D_2spp_2spp_.allFinite());

    // q = Tr(Q₂₀ρ) = Σab Q₂₀,ab ρba.
    const double Q20n_F = (Q20n_F2D_2spn_2spn.array() * rhon_F2D_2spn_2spn_.transpose().array()).sum();
    const double Q20p_F = (Q20p_F2D_2spp_2spp.array() * rhop_F2D_2spp_2spp_.transpose().array()).sum();
    // D = ℏω₀ₙ⟨Q₂₀,n⟩ + ℏω₀ₚ⟨Q₂₀,p⟩; common to nn,np,pp.
    const double denominator_F = hbarOmega0_n_F * Q20n_F + hbarOmega0_p_F * Q20p_F;
    assert(std::isfinite(denominator_F) && denominator_F != 0.0);
    const double factor_F = (2.0 * epsilon2_F / 3.0) / denominator_F;

    // χ₂,ab = (2ε₂/3)ℏω₀aℏω₀b/D.
    chi2_nn_F = factor_F * hbarOmega0_n_F * hbarOmega0_n_F;
    chi2_np_F = factor_F * hbarOmega0_n_F * hbarOmega0_p_F;
    chi2_pp_F = factor_F * hbarOmega0_p_F * hbarOmega0_p_F;
}

inline void PSMSphericalSetting::set_ci_PP(double gamma_F) {
    assert(std::isfinite(gamma_F));
    // G₂,nn/pp = γG₀,nn/pp.
    G2_nn_F = gamma_F * G0_nn_F;
    G2_pp_F = gamma_F * G0_pp_F;
}

inline void PSMSphericalSetting::set_ci_cut(const Eigen::MatrixXi& Ncqp_I2D_Ncqp_np, const Eigen::VectorXd& ECut_F1D_Ncqp, const Eigen::VectorXi& NCut_I1D_Ncqp, const Eigen::VectorXi& TwoKCut_I1D_Ncqp) {
    configuration_nucleus = HFBConfigurationNucleus(Ncqp_I2D_Ncqp_np, ECut_F1D_Ncqp, NCut_I1D_Ncqp, TwoKCut_I1D_Ncqp);
}

inline void PSMSphericalSetting::build_ci_config(const Eigen::VectorXd& Ehfn_F1D_2qpn, const Eigen::VectorXd& Eqpn_F1D_2qpn, const Eigen::VectorXi& TwoKn_I1D_2qpn, const Eigen::VectorXd& Ehfp_F1D_2qpp, const Eigen::VectorXd& Eqpp_F1D_2qpp, const Eigen::VectorXi& TwoKp_I1D_2qpp) {
    assert(configuration_nucleus.Ncqp_I2D_Ncqp_np.rows() > 0);
    // Nilsson-ordered representatives: +1/2,−3/2,+5/2,… .
    const auto calc_representative = [](const SphericalSetting& setting, const Eigen::VectorXd& Ehf_F1D_2qp, const Eigen::VectorXi& TwoK_I1D_2qp, const std::vector<int>& Nshell_I1D_Nactive) {
        const int NqpHalf_I = static_cast<int>(setting.labels_S1D_sp.size());
        assert(Ehf_F1D_2qp.size() == 2 * NqpHalf_I && TwoK_I1D_2qp.size() == Ehf_F1D_2qp.size());
        // K = +1/2,−3/2,+5/2,…; qp̄ = qp + NqpHalf.
        std::vector<int> qp_I1D_representative{};
        for (int qp_I = 0; qp_I < NqpHalf_I; ++qp_I) {
            const int N_I = setting.labels_S1D_sp[qp_I].N_I;
            if (std::find(Nshell_I1D_Nactive.begin(), Nshell_I1D_Nactive.end(), N_I) == Nshell_I1D_Nactive.end()) {continue;}
            const int partner_I = (TwoK_I1D_2qp(qp_I) % 4 == 3);
            qp_I1D_representative.push_back(qp_I + partner_I * NqpHalf_I);
        }
        std::stable_sort(qp_I1D_representative.begin(), qp_I1D_representative.end(), [&](int qp1_I, int qp2_I) {return Ehf_F1D_2qp(qp1_I) < Ehf_F1D_2qp(qp2_I);});
        return qp_I1D_representative;
    };
    const auto qpn_I1D_representative = calc_representative(sphericalsetting_neutron, Ehfn_F1D_2qpn, TwoKn_I1D_2qpn, Nshelln_I1D_Nactive);
    const auto qpp_I1D_representative = calc_representative(sphericalsetting_proton, Ehfp_F1D_2qpp, TwoKp_I1D_2qpp, Nshellp_I1D_Nactive);
    auto& configuration = configuration_nucleus;
    // (N mod 2,Z mod 2) → EE,EO,OE,OO.
    if (N_I % 2 == 0 && Z_I % 2 == 0) {configuration.build_config_kramers_ee(Eqpn_F1D_2qpn, TwoKn_I1D_2qpn, qpn_I1D_representative, Eqpp_F1D_2qpp, TwoKp_I1D_2qpp, qpp_I1D_representative);}
    if (N_I % 2 == 0 && Z_I % 2 == 1) {configuration.build_config_kramers_eo(Eqpn_F1D_2qpn, TwoKn_I1D_2qpn, qpn_I1D_representative, Eqpp_F1D_2qpp, TwoKp_I1D_2qpp, qpp_I1D_representative);}
    if (N_I % 2 == 1 && Z_I % 2 == 0) {configuration.build_config_kramers_oe(Eqpn_F1D_2qpn, TwoKn_I1D_2qpn, qpn_I1D_representative, Eqpp_F1D_2qpp, TwoKp_I1D_2qpp, qpp_I1D_representative);}
    if (N_I % 2 == 1 && Z_I % 2 == 1) {configuration.build_config_kramers_oo(Eqpn_F1D_2qpn, TwoKn_I1D_2qpn, qpn_I1D_representative, Eqpp_F1D_2qpp, TwoKp_I1D_2qpp, qpp_I1D_representative);}
    assert(!configuration.confign_I2D_cfgnp_cqpn.empty());
}
