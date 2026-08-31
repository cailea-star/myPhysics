/**
 * @file    cc_fusion.hpp
 * @author  cailea
 * @date    2026-08-25
 * @brief   Coupled-channel fusion wave and probability tools.
 */

#pragma once

#include <Eigen/Core>
#include <unsupported/Eigen/CXX11/Tensor>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <functional>
#include <numbers>
#include <tuple>
#include <utility>

#include "basis_bspline.hpp"
#include "bvp_rmatrix.hpp"
#include "cc_couplings.hpp"
#include "cc_params.hpp"
#include "ivp_numerov.hpp"
#include "ivp_rungekutta.hpp"
#include "spherical_radial_boundary.hpp"

/**
 * @brief  Generate the incoming-wave fusion boundary analytically.
 * @math   (u,u')=(e^{-ikr},-ike^{-ikr}) or (e^{κr},κe^{κr})
 * @output Inner-boundary pair (u,u').
 * @note   l_I and Ze2_F are callback-only.
 */
inline SphericalRadialBoundary fusion_inner_boundary(double r_F, int l_I, double hmass_F, double Ze2_F, double E_F) {
    assert(std::isfinite(r_F) && r_F > 0.0);
    assert(l_I >= 0 && std::isfinite(hmass_F) && hmass_F > 0.0);
    assert(std::isfinite(Ze2_F) && std::isfinite(E_F));
    (void) l_I;
    (void) Ze2_F;
    double k_F = std::sqrt(std::abs(E_F) / hmass_F);
    bool isClosed_B = E_F < 0.0;
    std::complex<double> sign_C = std::complex<double>(0.0, -1.0) * static_cast<double>(!isClosed_B) + std::complex<double>(1.0, 0.0) * static_cast<double>(isClosed_B);
    std::complex<double> u_C = std::exp(sign_C * k_F * r_F);
    std::complex<double> du_C = sign_C * k_F * u_C;
    return {u_C, du_C};
}

/**
 * @brief  Integrate coupled radial equations by classical RK4.
 * @math   d[u,u']/dr=[u',F(r)u]
 * @output Channel-by-solution-by-radius wave tensor.
 */
inline Eigen::Tensor<std::complex<double>, 3, Eigen::ColMajor> fusion_wavefunction_rk4(const Eigen::Ref<const Eigen::VectorXd>& r_F1D_r, const CCParams& params, int L_I, const SphericalRadialBoundaryFunc& boundary_Func = fusion_inner_boundary) {
    int Nch_I = static_cast<int>(params.channel_1D_ch.size());
    int Nsol_I = Nch_I;
    int Nr_I = static_cast<int>(r_F1D_r.size());
    assert(Nch_I > 0 && Nr_I >= 2 && L_I >= 0);
    assert(r_F1D_r.allFinite() && boundary_Func);

    // {E_ch-V_ch(r_0)-V_L(r_0)} → diag(u_0,u'_0).
    Eigen::MatrixXcd y0_C2D_ch_2sol = Eigen::MatrixXcd::Zero(Nch_I, 2 * Nsol_I);
    Eigen::VectorXd VRmin_F1D_ch = params.calc_Vall_matrix(r_F1D_r(0)).diagonal().real();
    for (int channel_I = 0; channel_I < Nch_I; ++channel_I) {
        int channel_l_I = L_I;
        double ERmin_F = params.channel_1D_ch[channel_I].Ech_F - VRmin_F1D_ch(channel_I) - calc_Vcent(r_F1D_r(0), params.Varg_params, channel_l_I);
        auto [u_C, du_C] = boundary_Func(r_F1D_r(0), channel_l_I, params.Varg_params.hmass_F, params.Varg_params.Ze2_F, ERmin_F);
        y0_C2D_ch_2sol(channel_I, channel_I) = u_C;
        y0_C2D_ch_2sol(channel_I, Nsol_I + channel_I) = du_C;
    }

    // vec([u,u']) → vec([u',Fu]).
    Eigen::MatrixXcd F_C2D_ch_ch(Nch_I, Nch_I);
    auto dydr_Func = [&](double r_F, const Eigen::Ref<const Eigen::VectorXcd>& y_C1D_ch2sol, Eigen::Ref<Eigen::VectorXcd> dydr_C1D_ch2sol) {
        Eigen::Map<const Eigen::MatrixXcd> y_C2D_ch_2sol(y_C1D_ch2sol.data(), Nch_I, 2 * Nsol_I);
        Eigen::Map<Eigen::MatrixXcd> dydr_C2D_ch_2sol(dydr_C1D_ch2sol.data(), Nch_I, 2 * Nsol_I);
        F_C2D_ch_ch = calc_coupled_F_matrix(r_F, params, L_I);
        dydr_C2D_ch_2sol.leftCols(Nsol_I) = y_C2D_ch_2sol.rightCols(Nsol_I);
        dydr_C2D_ch_2sol.rightCols(Nsol_I).noalias() = F_C2D_ch_ch * y_C2D_ch_2sol.leftCols(Nsol_I);
    };
    Eigen::VectorXcd y0_C1D_ch2sol = Eigen::Map<const Eigen::VectorXcd>(y0_C2D_ch_2sol.data(), y0_C2D_ch_2sol.size());
    IVP_RK4State<std::complex<double>> rk4_State(dydr_Func, r_F1D_r(0), y0_C1D_ch2sol);
    Eigen::Tensor<std::complex<double>, 3, Eigen::ColMajor> u_C3D_ch_sol_r(Nch_I, Nsol_I, Nr_I);
    Eigen::TensorMap<Eigen::Tensor<const std::complex<double>, 2, Eigen::ColMajor>> u0_C2D_ch_sol(y0_C2D_ch_2sol.data(), Nch_I, Nsol_I);
    u_C3D_ch_sol_r.chip(0, 2) = u0_C2D_ch_sol;

    // r_i → RK4(r_{i+1}).
    for (int r_I = 1; r_I < Nr_I; ++r_I) {
        Eigen::VectorXcd& ycurr_C1D_ch2sol = rk4_State.step(r_F1D_r(r_I));
        Eigen::Map<const Eigen::MatrixXcd> ycurr_C2D_ch_2sol(ycurr_C1D_ch2sol.data(), Nch_I, 2 * Nsol_I);
        Eigen::TensorMap<Eigen::Tensor<const std::complex<double>, 2, Eigen::ColMajor>> ucurr_C2D_ch_sol(ycurr_C2D_ch_2sol.data(), Nch_I, Nsol_I);
        u_C3D_ch_sol_r.chip(r_I, 2) = ucurr_C2D_ch_sol;
    }
    return u_C3D_ch_sol_r;
}

/**
 * @brief  Integrate coupled radial equations by matrix Numerov.
 * @math   u''(r)=F(r)u(r)
 * @output Channel-by-solution-by-radius wave tensor.
 * @note   RK4 supplies the first two mesh values.
 */
inline Eigen::Tensor<std::complex<double>, 3, Eigen::ColMajor> fusion_wavefunction_numerov(const Eigen::Ref<const Eigen::VectorXd>& r_F1D_r, const CCParams& params, int L_I, const SphericalRadialBoundaryFunc& boundary_Func = fusion_inner_boundary) {
    int Nch_I = static_cast<int>(params.channel_1D_ch.size());
    int Nr_I = static_cast<int>(r_F1D_r.size());
    assert(Nch_I > 0 && Nr_I >= 2 && L_I >= 0);
    Eigen::Tensor<std::complex<double>, 3, Eigen::ColMajor> u01_C3D_ch_sol_r = fusion_wavefunction_rk4(r_F1D_r.head(2), params, L_I, boundary_Func);
    Eigen::Map<const Eigen::MatrixXcd> u0_C2D_ch_sol(u01_C3D_ch_sol_r.data(), Nch_I, Nch_I);
    Eigen::Map<const Eigen::MatrixXcd> u1_C2D_ch_sol(u01_C3D_ch_sol_r.data() + Nch_I * Nch_I, Nch_I, Nch_I);
    auto F_Func = [&](double r_F, Eigen::Ref<Eigen::MatrixXcd> F_C2D_ch_ch) {F_C2D_ch_ch = calc_coupled_F_matrix(r_F, params, L_I);};
    return ivp_numerov_mat<std::complex<double>>(F_Func, u0_C2D_ch_sol, u1_C2D_ch_sol, r_F1D_r);
}

/**
 * @brief  Match numerical waves to Coulomb-Hankel solutions.
 * @math   u(r)=H⁺(r)C⁺+H⁻(r)C⁻
 * @output Coefficient pair (C⁺,C⁻).
 * @note   Uses the final and seventh-final mesh points.
 */
inline std::pair<Eigen::MatrixXcd, Eigen::MatrixXcd> align_Cpm(const Eigen::Ref<const Eigen::VectorXd>& r_F1D_r, const CCParams& params, int L_I, const Eigen::Tensor<std::complex<double>, 3, Eigen::ColMajor>& u_C3D_ch_sol_r) {
    int Nch_I = static_cast<int>(params.channel_1D_ch.size());
    int Nsol_I = static_cast<int>(u_C3D_ch_sol_r.dimension(1));
    int Nr_I = static_cast<int>(r_F1D_r.size());
    assert(Nch_I > 0 && Nsol_I > 0 && Nr_I >= 7 && L_I >= 0);
    assert(u_C3D_ch_sol_r.dimension(0) == Nch_I && u_C3D_ch_sol_r.dimension(2) == Nr_I);

    // {u(r_max),u(r_fit)} → [H⁺,H⁻][C⁺;C⁻].
    Eigen::Map<const Eigen::MatrixXcd> uRmax_C2D_ch_sol(u_C3D_ch_sol_r.data() + (Nr_I - 1) * Nch_I * Nsol_I, Nch_I, Nsol_I);
    Eigen::Map<const Eigen::MatrixXcd> uRfit_C2D_ch_sol(u_C3D_ch_sol_r.data() + (Nr_I - 7) * Nch_I * Nsol_I, Nch_I, Nsol_I);
    Eigen::MatrixXcd HpRmax_C2D_ch_ch = Eigen::MatrixXcd::Zero(Nch_I, Nch_I);
    Eigen::MatrixXcd HpRfit_C2D_ch_ch = Eigen::MatrixXcd::Zero(Nch_I, Nch_I);
    Eigen::MatrixXcd HmRmax_C2D_ch_ch = Eigen::MatrixXcd::Zero(Nch_I, Nch_I);
    Eigen::MatrixXcd HmRfit_C2D_ch_ch = Eigen::MatrixXcd::Zero(Nch_I, Nch_I);
    double rmax_F = r_F1D_r(Nr_I - 1);
    double rfit_F = r_F1D_r(Nr_I - 7);
    for (int channel_I = 0; channel_I < Nch_I; ++channel_I) {
        int channel_l_I = L_I;
        double Ech_F = params.channel_1D_ch[channel_I].Ech_F;
        HpRmax_C2D_ch_ch(channel_I, channel_I) = spherical_radial_boundary_coulomb_hplus(rmax_F, channel_l_I, params.Varg_params.hmass_F, params.Varg_params.Ze2_F, Ech_F).first;
        HpRfit_C2D_ch_ch(channel_I, channel_I) = spherical_radial_boundary_coulomb_hplus(rfit_F, channel_l_I, params.Varg_params.hmass_F, params.Varg_params.Ze2_F, Ech_F).first;
        HmRmax_C2D_ch_ch(channel_I, channel_I) = spherical_radial_boundary_coulomb_hminus(rmax_F, channel_l_I, params.Varg_params.hmass_F, params.Varg_params.Ze2_F, Ech_F).first;
        HmRfit_C2D_ch_ch(channel_I, channel_I) = spherical_radial_boundary_coulomb_hminus(rfit_F, channel_l_I, params.Varg_params.hmass_F, params.Varg_params.Ze2_F, Ech_F).first;
    }
    Eigen::MatrixXcd Hpm_C2D_2ch_2ch(2 * Nch_I, 2 * Nch_I);
    Hpm_C2D_2ch_2ch.topLeftCorner(Nch_I, Nch_I) = HpRmax_C2D_ch_ch;
    Hpm_C2D_2ch_2ch.topRightCorner(Nch_I, Nch_I) = HmRmax_C2D_ch_ch;
    Hpm_C2D_2ch_2ch.bottomLeftCorner(Nch_I, Nch_I) = HpRfit_C2D_ch_ch;
    Hpm_C2D_2ch_2ch.bottomRightCorner(Nch_I, Nch_I) = HmRfit_C2D_ch_ch;
    Eigen::MatrixXcd uBoundary_C2D_2ch_sol(2 * Nch_I, Nsol_I);
    uBoundary_C2D_2ch_sol.topRows(Nch_I) = uRmax_C2D_ch_sol;
    uBoundary_C2D_2ch_sol.bottomRows(Nch_I) = uRfit_C2D_ch_sol;
    Eigen::MatrixXcd Cpm_C2D_2ch_sol = Hpm_C2D_2ch_2ch.fullPivLu().solve(uBoundary_C2D_2ch_sol);
    Eigen::MatrixXcd Cp_C2D_ch_sol = Cpm_C2D_2ch_sol.topRows(Nch_I);
    Eigen::MatrixXcd Cm_C2D_ch_sol = Cpm_C2D_2ch_sol.bottomRows(Nch_I);
    return {Cp_C2D_ch_sol, Cm_C2D_ch_sol};
}

/**
 * @brief  Compute fusion probabilities by current matching.
 * @math   P_s=Σ_ch k_ch|f_ch,s|²/[-Im(H_s^-*H_s^-′)]
 * @output Fusion probability for each incident solution.
 * @note   Requires one numerical solution per channel.
 */
inline Eigen::VectorXd fusion_probability_u(const Eigen::Ref<const Eigen::VectorXd>& r_F1D_r, const CCParams& params, int L_I, const Eigen::Tensor<std::complex<double>, 3, Eigen::ColMajor>& u_C3D_ch_sol_r) {
    int Nch_I = static_cast<int>(params.channel_1D_ch.size());
    int Nsol_I = static_cast<int>(u_C3D_ch_sol_r.dimension(1));
    assert(Nch_I > 0 && Nsol_I == Nch_I && r_F1D_r.size() >= 7 && L_I >= 0);
    Eigen::Map<const Eigen::MatrixXcd> u0_C2D_ch_sol(u_C3D_ch_sol_r.data(), Nch_I, Nsol_I);
    auto [Cp_C2D_ch_sol, Cm_C2D_ch_sol] = align_Cpm(r_F1D_r, params, L_I, u_C3D_ch_sol_r);
    (void) Cp_C2D_ch_sol;
    Eigen::MatrixXcd f0_C2D_ch_sol = Cm_C2D_ch_sol.transpose().fullPivLu().solve(u0_C2D_ch_sol.transpose()).transpose();
    Eigen::VectorXd jRmin_F1D_sol = Eigen::VectorXd::Zero(Nsol_I);
    Eigen::VectorXd jRmax_F1D_sol(Nsol_I);
    Eigen::VectorXd VRmin_F1D_ch = params.calc_Vall_matrix(r_F1D_r(0)).diagonal().real();
    double rmax_F = r_F1D_r(r_F1D_r.size() - 1);

    // (f_0,H_s^-) → (j_min,j_max).
    for (int solution_I = 0; solution_I < Nsol_I; ++solution_I) {
        int solution_l_I = L_I;
        double Esol_F = params.channel_1D_ch[solution_I].Ech_F;
        auto [HmRmax_C, dHmRmax_C] = spherical_radial_boundary_coulomb_hminus(rmax_F, solution_l_I, params.Varg_params.hmass_F, params.Varg_params.Ze2_F, Esol_F);
        jRmax_F1D_sol(solution_I) = -std::imag(std::conj(HmRmax_C) * dHmRmax_C);
        for (int channel_I = 0; channel_I < Nch_I; ++channel_I) {
            int channel_l_I = L_I;
            double ERmin_F = params.channel_1D_ch[channel_I].Ech_F - VRmin_F1D_ch(channel_I) - calc_Vcent(r_F1D_r(0), params.Varg_params, channel_l_I);
            if (ERmin_F <= 0.0) {continue;}
            double kRmin_F = std::sqrt(ERmin_F / params.Varg_params.hmass_F);
            jRmin_F1D_sol(solution_I) += kRmin_F * std::norm(f0_C2D_ch_sol(channel_I, solution_I));
        }
    }
    return jRmin_F1D_sol.array() / jRmax_F1D_sol.array();
}

/**
 * @brief  Compute fusion probabilities from the scattering matrix.
 * @math   P_s=max(0,1-Σ_ch(k_ch/k_s)|S_ch,s|²)
 * @output Fusion probability for each incident channel.
 */
inline Eigen::VectorXd fusion_probability_S(const Eigen::Ref<const Eigen::VectorXd>& r_F1D_r, const CCParams& params, int L_I, const Eigen::Ref<const Eigen::MatrixXcd>& S_C2D_ch_sol) {
    int Nch_I = static_cast<int>(params.channel_1D_ch.size());
    assert(Nch_I > 0 && r_F1D_r.size() >= 2 && L_I >= 0);
    assert(S_C2D_ch_sol.rows() == Nch_I && S_C2D_ch_sol.cols() == Nch_I);
    Eigen::VectorXd kmax_F1D_ch(Nch_I);
    for (int channel_I = 0; channel_I < Nch_I; ++channel_I) {
        assert(params.channel_1D_ch[channel_I].Ech_F > 0.0);
        kmax_F1D_ch(channel_I) = std::sqrt(params.channel_1D_ch[channel_I].Ech_F / params.Varg_params.hmass_F);
    }
    Eigen::VectorXd P_F1D_sol = Eigen::VectorXd::Ones(Nch_I);
    for (int solution_I = 0; solution_I < Nch_I; ++solution_I) {
        for (int channel_I = 0; channel_I < Nch_I; ++channel_I) {P_F1D_sol(solution_I) -= kmax_F1D_ch(channel_I) / kmax_F1D_ch(solution_I) * std::norm(S_C2D_ch_sol(channel_I, solution_I));}
        P_F1D_sol(solution_I) = std::max(0.0, P_F1D_sol(solution_I));
    }
    return P_F1D_sol;
}

/**
 * @brief  Compute fusion probabilities by B-spline R-matrix matching.
 * @math   S=(H⁺-RH⁺′)⁻¹(H⁻-RH⁻′); P=1-K^{-1}S†KS
 * @output Fusion probability for each incident channel.
 * @note   Requires endpoint-retaining B-splines.
 */
inline Eigen::VectorXd fusion_probability_R(const CCParams& params, int L_I, RMatrix& r_matrix) {
    int Nch_I = static_cast<int>(params.channel_1D_ch.size());
    const Eigen::VectorXd& r_F1D_r = r_matrix.b_basis.b_funcs.x_F1D_x;
    assert(Nch_I > 0 && r_F1D_r.size() >= 2 && L_I >= 0);
    double rmin_F = r_F1D_r(0);
    double rmax_F = r_F1D_r(r_F1D_r.size() - 1);
    Eigen::VectorXcd lambdaIn_C1D_ch(Nch_I);
    Eigen::VectorXd VRmin_F1D_ch = params.calc_Vall_matrix(rmin_F).diagonal().real();
    for (int channel_I = 0; channel_I < Nch_I; ++channel_I) {
        double ERmin_F = params.channel_1D_ch[channel_I].Ech_F - VRmin_F1D_ch(channel_I) - calc_Vcent(rmin_F, params.Varg_params, L_I);
        double k_F = std::sqrt(std::abs(ERmin_F) / params.Varg_params.hmass_F);
        lambdaIn_C1D_ch(channel_I) = ERmin_F >= 0.0 ? std::complex<double>(0.0, -k_F) : std::complex<double>(k_F, 0.0);
    }
    r_matrix.update_Ginv([&](double r_F) {return calc_coupled_F_matrix(r_F, params, L_I);}, lambdaIn_C1D_ch);
    Eigen::VectorXcd Hp_C1D_ch(Nch_I);
    Eigen::VectorXcd dHp_C1D_ch(Nch_I);
    Eigen::VectorXcd Hm_C1D_ch(Nch_I);
    Eigen::VectorXcd dHm_C1D_ch(Nch_I);
    for (int channel_I = 0; channel_I < Nch_I; ++channel_I) {
        double Ech_F = params.channel_1D_ch[channel_I].Ech_F;
        assert(Ech_F > 0.0);
        auto [Hp_C, dHp_C] = spherical_radial_boundary_coulomb_hplus(rmax_F, L_I, params.Varg_params.hmass_F, params.Varg_params.Ze2_F, Ech_F);
        auto [Hm_C, dHm_C] = spherical_radial_boundary_coulomb_hminus(rmax_F, L_I, params.Varg_params.hmass_F, params.Varg_params.Ze2_F, Ech_F);
        Hp_C1D_ch(channel_I) = Hp_C;
        dHp_C1D_ch(channel_I) = dHp_C;
        Hm_C1D_ch(channel_I) = Hm_C;
        dHm_C1D_ch(channel_I) = dHm_C;
    }
    Eigen::MatrixXcd S_C2D_ch_sol = r_matrix.calc_S_matrix(Hp_C1D_ch, dHp_C1D_ch, Hm_C1D_ch, dHm_C1D_ch);
    return fusion_probability_S(r_F1D_r, params, L_I, S_C2D_ch_sol);
}

/**
 * @brief  Sum partial waves into fusion observables.
 * @math   σ=πh_μE⁻¹Σ_L(2L+1)P_L; ⟨L⟩=Σ_LLw_L/Σ_Lw_L
 * @output Tuple (σ,⟨L⟩,P₀).
 * @note   σ uses millibarns.
 */
inline std::tuple<double, double, double> fusion_cross_section(const CCParams& params, const std::function<double(int)>& PL_Func) {
    assert(!params.channel_1D_ch.empty() && PL_Func);
    double E_F = params.channel_1D_ch[0].Ech_F;
    assert(std::isfinite(E_F) && E_F > 0.0);
    double sigmaFactor_F = std::numbers::pi * params.Varg_params.hmass_F / E_F * 10.0;
    double P0_F = PL_Func(0);
    double sumPL_F = P0_F;
    double sumLPL_F = 0.0;
    for (int L_I = 1; L_I <= 200; ++L_I) {
        double PL_F = PL_Func(L_I);
        double weight_F = (2.0 * L_I + 1.0) * PL_F;
        sumPL_F += weight_F;
        sumLPL_F += L_I * weight_F;
        if (weight_F < sumPL_F * 1.0e-4) {break;}
    }
    double Lavg_F = sumPL_F > 0.0 ? sumLPL_F / sumPL_F : 0.0;
    return {sigmaFactor_F * sumPL_F, Lavg_F, P0_F};
}
