/**
 * @file    cc_decay.hpp
 * @author  cailea
 * @date    2026-08-25
 * @brief   Coupled-channel decay matching, widths, and node counts.
 */

#pragma once

#include <Eigen/Core>
#include <Eigen/QR>
#include <Eigen/SVD>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstddef>
#include <tuple>
#include <utility>
#include <vector>

#include "cc_params.hpp"
#include "cc_couplings.hpp"
#include "ivp_logderiv.hpp"
#include "spherical_radial_boundary.hpp"

namespace {

/**
 * @brief  Build diagonal boundaries by analytic power laws.
 * @math   B(rₘᵢₙ)=[diag(uₗ),diag(uₗ′)]
 * @output Regular coupled-channel boundary matrix.
 */
inline Eigen::MatrixXcd decay_regular_boundary(double rmin_F, const CCParams& params)
{
    int Nch_I = static_cast<int>(params.channel_1D_ch.size());
    Eigen::MatrixXcd boundary_C2D_ch_2sol = Eigen::MatrixXcd::Zero(Nch_I, 2 * Nch_I);

    // l_ch → (uₗ,uₗ′) → B_ch,ch.
    for (int channel_I = 0; channel_I < Nch_I; ++channel_I) {
        auto [u_C, dudr_C] = spherical_radial_boundary_regular(rmin_F, params.channel_1D_ch[channel_I].twol_I / 2, params.Varg_params.hmass_F, params.Varg_params.Ze2_F, params.channel_1D_ch[channel_I].Ech_F);
        boundary_C2D_ch_2sol(channel_I, channel_I) = u_C;
        boundary_C2D_ch_2sol(channel_I, Nch_I + channel_I) = dudr_C;
    }
    return boundary_C2D_ch_2sol;
}

/**
 * @brief  Build real Coulomb-Hankel boundaries through GSL.
 * @math   B(rₘₐₓ)=[diag(Re Hₗ⁺),diag(Re Hₗ⁺′)]
 * @output Outgoing coupled-channel boundary matrix.
 */
inline Eigen::MatrixXcd decay_coulomb_boundary(double rmax_F, const CCParams& params)
{
    int Nch_I = static_cast<int>(params.channel_1D_ch.size());
    Eigen::MatrixXcd boundary_C2D_ch_2sol = Eigen::MatrixXcd::Zero(Nch_I, 2 * Nch_I);

    // (l_ch,E_ch,rₘₐₓ) → Re(Hₗ⁺,Hₗ⁺′) → B_ch,ch.
    for (int channel_I = 0; channel_I < Nch_I; ++channel_I) {
        auto [Hplus_C, dHplusdr_C] = spherical_radial_boundary_coulomb_hplus(rmax_F, params.channel_1D_ch[channel_I].twol_I / 2, params.Varg_params.hmass_F, params.Varg_params.Ze2_F, params.channel_1D_ch[channel_I].Ech_F);
        boundary_C2D_ch_2sol(channel_I, channel_I) = Hplus_C.real();
        boundary_C2D_ch_2sol(channel_I, Nch_I + channel_I) = dHplusdr_C.real();
    }
    return boundary_C2D_ch_2sol;
}

}

/**
 * @brief  Match log derivatives by QR-RK4 and SVD.
 * @math   Y=dQ Q⁻¹; (Yᵢₙ-Yₒᵤₜ)uₘ=0
 * @output Minimum mismatch and inner/outer coefficients.
 * @note   Requires endpoint-clamped adaptive steps.
 */
inline std::tuple<double, Eigen::VectorXcd, Eigen::VectorXcd> decay_match(double rmin_F, double rmatch_F, double rmax_F, double dr_F, const CCParams& params)
{
    int Nch_I = static_cast<int>(params.channel_1D_ch.size());
    assert(Nch_I > 0 && 0.0 < rmin_F && rmin_F < rmatch_F && rmatch_F < rmax_F && dr_F > 0.0);

    // r → F(r).
    Real2TMatFunc<std::complex<double>> F_Func = [&](double rstage_F, Eigen::Ref<Eigen::MatrixXcd> F_C2D_ch_ch) {
        F_C2D_ch_ch = calc_coupled_F_matrix(rstage_F, params);
    };

    // (Q₀,dQ₀,R₀,r₀) → (Qₘ,dQₘ,Rₘ,rₘ).
    auto propagate = [&](IVP_RK4QRState<std::complex<double>> state_QR, double r_F, double rend_F, double direction_F) {
        assert(std::abs(direction_F) == 1.0);
        while (direction_F * (rend_F - r_F) > 1.0e-12) {
            double drstep_F = direction_F * std::min({dr_F, std::abs(rend_F - r_F), 0.1 * r_F});
            state_QR.step(r_F + drstep_F);
            r_F += drstep_F;
        }
        return state_QR;
    };

    Eigen::MatrixXcd inBoundary_C2D_ch_2sol = decay_regular_boundary(rmin_F, params);
    Eigen::MatrixXcd outBoundary_C2D_ch_2sol = decay_coulomb_boundary(rmax_F, params);

    // (Bᵢₙ,Bₒᵤₜ) → (Qᵢₙ,dQᵢₙ,Rᵢₙ,Qₒᵤₜ,dQₒᵤₜ,Rₒᵤₜ).
    IVP_RK4QRState<std::complex<double>> inState_QR = propagate(IVP_RK4QRState<std::complex<double>>(F_Func, rmin_F, inBoundary_C2D_ch_2sol.leftCols(Nch_I), inBoundary_C2D_ch_2sol.rightCols(Nch_I)), rmin_F, rmatch_F, 1.0);
    IVP_RK4QRState<std::complex<double>> outState_QR = propagate(IVP_RK4QRState<std::complex<double>>(F_Func, rmax_F, outBoundary_C2D_ch_2sol.leftCols(Nch_I), outBoundary_C2D_ch_2sol.rightCols(Nch_I)), rmax_F, rmatch_F, -1.0);

    // Y=(∂ᵣu)u⁻¹=dQ Q⁻¹.
    const Eigen::MatrixXcd& QIn_C2D_ch_ch = inState_QR.Qcurr_T2D_ch_sol;
    const Eigen::MatrixXcd& dQdrIn_C2D_ch_ch = inState_QR.dQcurr_T2D_ch_sol;
    const Eigen::MatrixXcd& QOut_C2D_ch_ch = outState_QR.Qcurr_T2D_ch_sol;
    const Eigen::MatrixXcd& dQdrOut_C2D_ch_ch = outState_QR.dQcurr_T2D_ch_sol;
    const Eigen::MatrixXcd Yin_C2D_ch_ch = QIn_C2D_ch_ch.transpose().colPivHouseholderQr().solve(dQdrIn_C2D_ch_ch.transpose()).transpose();
    const Eigen::MatrixXcd Yout_C2D_ch_ch = QOut_C2D_ch_ch.transpose().colPivHouseholderQr().solve(dQdrOut_C2D_ch_ch.transpose()).transpose();

    // εₘ=min ||(Yᵢₙ-Yₒᵤₜ)uₘ||₂.
    const Eigen::JacobiSVD<Eigen::MatrixXcd> matchSolver_SVD(Yin_C2D_ch_ch - Yout_C2D_ch_ch, Eigen::ComputeFullV);
    Eigen::Index minSingular_I = Nch_I - 1;
    double matchError_F = matchSolver_SVD.singularValues()(minSingular_I);
    Eigen::VectorXcd umatch_C1D_ch = matchSolver_SVD.matrixV().col(minSingular_I);

    // uₘ=uc=(QR)c=Q(Rc)=Qa.
    Eigen::VectorXcd aIn_C1D_ch = QIn_C2D_ch_ch.colPivHouseholderQr().solve(umatch_C1D_ch);
    Eigen::VectorXcd aOut_C1D_ch = QOut_C2D_ch_ch.colPivHouseholderQr().solve(umatch_C1D_ch);
    Eigen::VectorXcd cIn_C1D_ch = inState_QR.Rcurr_T2D_sol_sol.colPivHouseholderQr().solve(aIn_C1D_ch);
    Eigen::VectorXcd cOut_C1D_ch = outState_QR.Rcurr_T2D_sol_sol.colPivHouseholderQr().solve(aOut_C1D_ch);

    // arg(cₒᵤₜ,max)=0.
    Eigen::Index maxCoefficient_I;
    cOut_C1D_ch.cwiseAbs().maxCoeff(&maxCoefficient_I);
    std::complex<double> phase_C = cOut_C1D_ch(maxCoefficient_I) / std::abs(cOut_C1D_ch(maxCoefficient_I));
    cOut_C1D_ch /= phase_C;
    cIn_C1D_ch /= phase_C;
    return {matchError_F, cIn_C1D_ch, cOut_C1D_ch};
}

/**
 * @brief  Integrate norms by QR-RK4 and trapezoids.
 * @math   u=QR; ||u||²=|R₀₀|²; N=∫||u||²dr
 * @math   Γ_ch=2P_αh_μk_ch|cₒᵤₜ,ch|²/N
 * @output Partial decay widths by channel.
 * @note   Uses preformation factor P_α=0.39.
 * @note   Requires endpoint-clamped adaptive steps.
 */
inline Eigen::VectorXd decay_width(double rmin_F, double rmatch_F, double rmax_F, double dr_F, const CCParams& params)
{
    int Nch_I = static_cast<int>(params.channel_1D_ch.size());
    assert(Nch_I > 0 && 0.0 < rmin_F && rmin_F < rmatch_F && rmatch_F < rmax_F && dr_F > 0.0);

    // (rₘᵢₙ,rₘ,rₘₐₓ,dr) → (εₘ,cᵢₙ,cₒᵤₜ).
    auto [matchError_F, cIn_C1D_ch, cOut_C1D_ch] = decay_match(rmin_F, rmatch_F, rmax_F, dr_F, params);
    (void) matchError_F;

    // r → F(r).
    Real2TMatFunc<std::complex<double>> F_Func = [&](double rstage_F, Eigen::Ref<Eigen::MatrixXcd> F_C2D_ch_ch) {
        F_C2D_ch_ch = calc_coupled_F_matrix(rstage_F, params);
    };

    // Bᵢₙcᵢₙ → (uᵢₙ,uᵢₙ′) → (Qᵢₙ,dQᵢₙ,Rᵢₙ).
    Eigen::MatrixXcd boundary_C2D_ch_2sol = decay_regular_boundary(rmin_F, params);
    Eigen::MatrixXcd y_C2D_ch_2(Nch_I, 2);
    y_C2D_ch_2.leftCols(1).noalias() = boundary_C2D_ch_2sol.leftCols(Nch_I) * cIn_C1D_ch;
    y_C2D_ch_2.rightCols(1).noalias() = boundary_C2D_ch_2sol.rightCols(Nch_I) * cIn_C1D_ch;
    IVP_RK4QRState<std::complex<double>> inState_QR(F_Func, rmin_F, y_C2D_ch_2.leftCols(1), y_C2D_ch_2.rightCols(1));

    // (N,ρ,r)=(0,|R₀₀|²,rₘᵢₙ).
    double norm_F = 0.0;
    double density_F = std::norm(inState_QR.Rcurr_T2D_sol_sol(0, 0));
    double r_F = rmin_F;

    // N ← N + dr(ρᵢ+ρᵢ₊₁)/2.
    while (rmatch_F - r_F > 1.0e-12) {
        double drstep_F = std::min({dr_F, rmatch_F - r_F, 0.1 * r_F});
        inState_QR.step(r_F + drstep_F);
        double nextDensity_F = std::norm(inState_QR.Rcurr_T2D_sol_sol(0, 0));
        norm_F += 0.5 * drstep_F * (density_F + nextDensity_F);
        density_F = nextDensity_F;
        r_F += drstep_F;
    }

    // Bₒᵤₜcₒᵤₜ → (uₒᵤₜ,uₒᵤₜ′) → (Qₒᵤₜ,dQₒᵤₜ,Rₒᵤₜ).
    boundary_C2D_ch_2sol = decay_coulomb_boundary(rmax_F, params);
    y_C2D_ch_2.leftCols(1).noalias() = boundary_C2D_ch_2sol.leftCols(Nch_I) * cOut_C1D_ch;
    y_C2D_ch_2.rightCols(1).noalias() = boundary_C2D_ch_2sol.rightCols(Nch_I) * cOut_C1D_ch;
    IVP_RK4QRState<std::complex<double>> outState_QR(F_Func, rmax_F, y_C2D_ch_2.leftCols(1), y_C2D_ch_2.rightCols(1));
    density_F = std::norm(outState_QR.Rcurr_T2D_sol_sol(0, 0));
    r_F = rmax_F;

    // N ← N - dr(ρᵢ+ρᵢ₋₁)/2; dr<0.
    while (r_F - rmatch_F > 1.0e-12) {
        double drstep_F = -std::min({dr_F, r_F - rmatch_F, 0.1 * r_F});
        outState_QR.step(r_F + drstep_F);
        double nextDensity_F = std::norm(outState_QR.Rcurr_T2D_sol_sol(0, 0));
        norm_F -= 0.5 * drstep_F * (density_F + nextDensity_F);
        density_F = nextDensity_F;
        r_F += drstep_F;
    }
    assert(std::isfinite(norm_F) && norm_F > 0.0);

    constexpr double Palpha_F = 0.39;
    Eigen::VectorXd Gamma_F1D_ch(Nch_I);

    // (E_ch,cₒᵤₜ,ch,N) → Γ_ch.
    for (int channel_I = 0; channel_I < Nch_I; ++channel_I) {
        assert(params.channel_1D_ch[channel_I].Ech_F > 0.0);
        double k_F = std::sqrt(params.channel_1D_ch[channel_I].Ech_F / params.Varg_params.hmass_F);
        Gamma_F1D_ch(channel_I) = 2.0 * Palpha_F * params.Varg_params.hmass_F * k_F * std::norm(cOut_C1D_ch(channel_I)) / norm_F;
    }
    return Gamma_F1D_ch;
}

/**
 * @brief  Count nodes after QR-RK4 propagation.
 * @math   uᵢ=QᵢRᵢ; N_ch=Σᵢ𝟙[uᵢ₋₁,ch uᵢ,ch<0]
 * @output Internal radial-node count by channel.
 * @note   Ignores amplitudes below 10⁻⁸ of maxima.
 * @note   Requires endpoint-clamped adaptive steps.
 */
inline Eigen::VectorXi decay_nodes(double rmin_F, double rmatch_F, double rmax_F, double dr_F, const CCParams& params)
{
    int Nch_I = static_cast<int>(params.channel_1D_ch.size());
    assert(Nch_I > 0 && 0.0 < rmin_F && rmin_F < rmatch_F && rmatch_F < rmax_F && dr_F > 0.0);

    // (rₘᵢₙ,rₘ,rₘₐₓ,dr) → (εₘ,cᵢₙ,cₒᵤₜ).
    auto [matchError_F, cIn_C1D_ch, cOut_C1D_ch] = decay_match(rmin_F, rmatch_F, rmax_F, dr_F, params);
    (void) matchError_F;
    (void) cOut_C1D_ch;

    // r → F(r).
    Real2TMatFunc<std::complex<double>> F_Func = [&](double rstage_F, Eigen::Ref<Eigen::MatrixXcd> F_C2D_ch_ch) {
        F_C2D_ch_ch = calc_coupled_F_matrix(rstage_F, params);
    };

    // cᵢₙ → cᵢₙ exp[-i arg(cᵢₙ,max)].
    Eigen::VectorXcd cin_C1D_ch = cIn_C1D_ch;
    Eigen::Index maxCoefficient_I;
    cin_C1D_ch.cwiseAbs().maxCoeff(&maxCoefficient_I);
    assert(std::abs(cin_C1D_ch(maxCoefficient_I)) > 0.0);
    cin_C1D_ch /= cin_C1D_ch(maxCoefficient_I) / std::abs(cin_C1D_ch(maxCoefficient_I));

    // Bᵢₙcᵢₙ → (u₀,u₀′) → (Q₀,dQ₀,R₀).
    Eigen::MatrixXcd boundary_C2D_ch_2sol = decay_regular_boundary(rmin_F, params);
    Eigen::MatrixXcd y_C2D_ch_2sol(Nch_I, 2);
    y_C2D_ch_2sol.leftCols(1).noalias() = boundary_C2D_ch_2sol.leftCols(Nch_I) * cin_C1D_ch;
    y_C2D_ch_2sol.rightCols(1).noalias() = boundary_C2D_ch_2sol.rightCols(Nch_I) * cin_C1D_ch;
    IVP_RK4QRState<std::complex<double>> state_QR(F_Func, rmin_F, y_C2D_ch_2sol.leftCols(1), y_C2D_ch_2sol.rightCols(1));
    Eigen::VectorXcd u_C1D_ch = y_C2D_ch_2sol.col(0);
    std::vector<Eigen::VectorXd> u_F2D_r_ch;

    // {u(rᵢ)} ← u₀.
    u_F2D_r_ch.reserve(static_cast<std::size_t>(std::ceil((rmatch_F - rmin_F) / std::min(dr_F, 0.1 * rmin_F))) + 1);
    u_F2D_r_ch.push_back(u_C1D_ch.real());
    double r_F = rmin_F;

    // (Qᵢ,Rᵢ) → uᵢ=QᵢRᵢ → {u(rᵢ)}.
    while (rmatch_F - r_F > 1.0e-12) {
        double drstep_F = std::min({dr_F, rmatch_F - r_F, 0.1 * r_F});
        state_QR.step(r_F + drstep_F);
        u_C1D_ch = state_QR.Qcurr_T2D_ch_sol.col(0) * state_QR.Rcurr_T2D_sol_sol(0, 0);
        u_F2D_r_ch.push_back(u_C1D_ch.real());
        r_F += drstep_F;
    }

    // uₘₐₓ,ch=maxᵢ|uᵢ,ch|.
    Eigen::VectorXd umax_F1D_ch = Eigen::VectorXd::Zero(Nch_I);
    for (const Eigen::VectorXd& u_F1D_ch : u_F2D_r_ch) {
        umax_F1D_ch = umax_F1D_ch.cwiseMax(u_F1D_ch.cwiseAbs());
    }

    Eigen::VectorXi nodes_I1D_ch = Eigen::VectorXi::Zero(Nch_I);
    constexpr double nodeTolerance_F = 1.0e-8;

    // N_ch=Σᵢ𝟙[uᵢ₋₁,ch uᵢ,ch<0]𝟙[|u|>εuₘₐₓ,ch].
    for (std::size_t r_I = 1; r_I < u_F2D_r_ch.size(); ++r_I) {
        for (int channel_I = 0; channel_I < Nch_I; ++channel_I) {
            double uprev_F = u_F2D_r_ch[r_I - 1](channel_I);
            double u_F = u_F2D_r_ch[r_I](channel_I);
            nodes_I1D_ch(channel_I) += static_cast<int>((uprev_F * u_F < 0.0) & (std::max(std::abs(uprev_F), std::abs(u_F)) > nodeTolerance_F * umax_F1D_ch(channel_I)));
        }
    }
    return nodes_I1D_ch;
}
