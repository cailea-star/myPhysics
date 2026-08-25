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
#include "ivp_rungekutta.hpp"
#include "spherical_radial_boundary.hpp"

namespace {

/**
 * @brief  Advance coupled radial solutions by one RK4 step.
 * @math   (u,u')_{n+1}=RK4[(u',Fu),r_n,dr]
 * @output Coupled radial solutions at r+dr.
 */
inline Eigen::MatrixXcd decay_rk4_step(const Eigen::Ref<const Eigen::MatrixXcd>& y_C2D_ch_2sol, double r_F, double dr_F, const CCParams& params)
{
    int Nch_I = static_cast<int>(params.channel_1D_ch.size());
    int Nsol_I = static_cast<int>(y_C2D_ch_2sol.cols() / 2);
    assert(y_C2D_ch_2sol.rows() == Nch_I && y_C2D_ch_2sol.cols() == 2 * Nsol_I && Nsol_I > 0);

    // (u,u') → (u',Fu).
    auto rhs_Func = [&](double rstage_F, const Eigen::Ref<const Eigen::VectorXcd>& y_C1D_ch2sol, Eigen::Ref<Eigen::VectorXcd> dydr_C1D_ch2sol) {
        Eigen::Map<const Eigen::MatrixXcd> y_C2D_ch_2sol(y_C1D_ch2sol.data(), Nch_I, 2 * Nsol_I);
        Eigen::Map<Eigen::MatrixXcd> dydr_C2D_ch_2sol(dydr_C1D_ch2sol.data(), Nch_I, 2 * Nsol_I);
        auto u_C2D_ch_sol = y_C2D_ch_2sol.leftCols(Nsol_I);
        auto dudr_C2D_ch_sol = y_C2D_ch_2sol.rightCols(Nsol_I);
        dydr_C2D_ch_2sol.leftCols(Nsol_I) = dudr_C2D_ch_sol;
        dydr_C2D_ch_2sol.rightCols(Nsol_I).noalias() = calc_coupled_F_matrix(rstage_F, params) * u_C2D_ch_sol;
    };

    // y_{ch,2sol} ↔ y_i.
    Eigen::VectorXcd y_C1D_ch2sol(Eigen::Map<const Eigen::VectorXcd>(y_C2D_ch_2sol.data(), y_C2D_ch_2sol.size()));
    IVP_RK4State<std::complex<double>> rk4_State(rhs_Func, r_F, y_C1D_ch2sol);
    const Eigen::VectorXcd& ynext_C1D_ch2sol = rk4_State.step(r_F + dr_F);
    return Eigen::MatrixXcd(Eigen::Map<const Eigen::MatrixXcd>(ynext_C1D_ch2sol.data(), Nch_I, 2 * Nsol_I));
}

/**
 * @brief  Build diagonal regular power-law boundaries.
 * @math   B(r_{min})=[diag(u_l),diag(u_l')]
 * @output Regular coupled-channel boundary matrix.
 */
inline Eigen::MatrixXcd decay_regular_boundary(double rmin_F, const CCParams& params)
{
    int Nch_I = static_cast<int>(params.channel_1D_ch.size());
    Eigen::MatrixXcd boundary_C2D_ch_2sol = Eigen::MatrixXcd::Zero(Nch_I, 2 * Nch_I);
    for (int channel_I = 0; channel_I < Nch_I; ++channel_I) {
        auto [u_C, dudr_C] = spherical_radial_boundary_regular(rmin_F, params.channel_1D_ch[channel_I].twol_I / 2, params.Varg_params.hmass_F, params.Varg_params.Ze2_F, params.channel_1D_ch[channel_I].Ech_F);
        boundary_C2D_ch_2sol(channel_I, channel_I) = u_C;
        boundary_C2D_ch_2sol(channel_I, Nch_I + channel_I) = dudr_C;
    }
    return boundary_C2D_ch_2sol;
}

/**
 * @brief  Build diagonal real Coulomb-Hankel boundaries.
 * @math   B(r_{max})=[diag(Re H_l^+),diag(Re ∂_rH_l^+)]
 * @output Outgoing coupled-channel boundary matrix.
 */
inline Eigen::MatrixXcd decay_coulomb_boundary(double rmax_F, const CCParams& params)
{
    int Nch_I = static_cast<int>(params.channel_1D_ch.size());
    Eigen::MatrixXcd boundary_C2D_ch_2sol = Eigen::MatrixXcd::Zero(Nch_I, 2 * Nch_I);
    for (int channel_I = 0; channel_I < Nch_I; ++channel_I) {
        auto [Hplus_C, dHplusdr_C] = spherical_radial_boundary_coulomb_hplus(rmax_F, params.channel_1D_ch[channel_I].twol_I / 2, params.Varg_params.hmass_F, params.Varg_params.Ze2_F, params.channel_1D_ch[channel_I].Ech_F);
        boundary_C2D_ch_2sol(channel_I, channel_I) = Hplus_C.real();
        boundary_C2D_ch_2sol(channel_I, Nch_I + channel_I) = dHplusdr_C.real();
    }
    return boundary_C2D_ch_2sol;
}

}

/**
 * @brief  Match regular and Coulomb solutions by QR-SVD.
 * @math   (Y_{in}-Y_{out})u_m=0
 * @output Minimum mismatch and inner/outer coefficients.
 * @note   Adaptive endpoint-controlled RK4 loops are approved.
 */
inline std::tuple<double, Eigen::VectorXcd, Eigen::VectorXcd> decay_match(double rmin_F, double rmatch_F, double rmax_F, double dr_F, const CCParams& params)
{
    int Nch_I = static_cast<int>(params.channel_1D_ch.size());
    assert(Nch_I > 0 && 0.0 < rmin_F && rmin_F < rmatch_F && rmatch_F < rmax_F && dr_F > 0.0);

    Real2TMatFunc<std::complex<double>> F_Func = [&](double rstage_F, Eigen::Ref<Eigen::MatrixXcd> F_C2D_ch_ch) {
        F_C2D_ch_ch = calc_coupled_F_matrix(rstage_F, params);
    };

    // (Q_0,dQ_0,R_0,r_0) → (Q_m,dQ_m,R_m,r_m).
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
    IVP_RK4QRState<std::complex<double>> inState_QR = propagate(IVP_RK4QRState<std::complex<double>>(F_Func, rmin_F, inBoundary_C2D_ch_2sol.leftCols(Nch_I), inBoundary_C2D_ch_2sol.rightCols(Nch_I)), rmin_F, rmatch_F, 1.0);
    IVP_RK4QRState<std::complex<double>> outState_QR = propagate(IVP_RK4QRState<std::complex<double>>(F_Func, rmax_F, outBoundary_C2D_ch_2sol.leftCols(Nch_I), outBoundary_C2D_ch_2sol.rightCols(Nch_I)), rmax_F, rmatch_F, -1.0);

    // Y=(∂_ru)u^{-1}=dQ Q^{-1}.
    const Eigen::MatrixXcd& QIn_C2D_ch_ch = inState_QR.Qcurr_T2D_ch_sol;
    const Eigen::MatrixXcd& dQdrIn_C2D_ch_ch = inState_QR.dQcurr_T2D_ch_sol;
    const Eigen::MatrixXcd& QOut_C2D_ch_ch = outState_QR.Qcurr_T2D_ch_sol;
    const Eigen::MatrixXcd& dQdrOut_C2D_ch_ch = outState_QR.dQcurr_T2D_ch_sol;
    const Eigen::MatrixXcd Yin_C2D_ch_ch = QIn_C2D_ch_ch.transpose().colPivHouseholderQr().solve(dQdrIn_C2D_ch_ch.transpose()).transpose();
    const Eigen::MatrixXcd Yout_C2D_ch_ch = QOut_C2D_ch_ch.transpose().colPivHouseholderQr().solve(dQdrOut_C2D_ch_ch.transpose()).transpose();

    // min ||(Y_in-Y_out)u_m||_2.
    const Eigen::JacobiSVD<Eigen::MatrixXcd> matchSolver_SVD(Yin_C2D_ch_ch - Yout_C2D_ch_ch, Eigen::ComputeFullV);
    Eigen::Index minSingular_I = Nch_I - 1;
    double matchError_F = matchSolver_SVD.singularValues()(minSingular_I);
    Eigen::VectorXcd umatch_C1D_ch = matchSolver_SVD.matrixV().col(minSingular_I);

    // u_m = u c = (QR) c = Q (R c) = Q a.
    Eigen::VectorXcd aIn_C1D_ch = QIn_C2D_ch_ch.colPivHouseholderQr().solve(umatch_C1D_ch);
    Eigen::VectorXcd aOut_C1D_ch = QOut_C2D_ch_ch.colPivHouseholderQr().solve(umatch_C1D_ch);
    Eigen::VectorXcd cIn_C1D_ch = inState_QR.Rcurr_T2D_sol_sol.colPivHouseholderQr().solve(aIn_C1D_ch);
    Eigen::VectorXcd cOut_C1D_ch = outState_QR.Rcurr_T2D_sol_sol.colPivHouseholderQr().solve(aOut_C1D_ch);

    // arg(c_out,max)=0.
    Eigen::Index maxCoefficient_I;
    cOut_C1D_ch.cwiseAbs().maxCoeff(&maxCoefficient_I);
    std::complex<double> phase_C = cOut_C1D_ch(maxCoefficient_I) / std::abs(cOut_C1D_ch(maxCoefficient_I));
    cOut_C1D_ch /= phase_C;
    cIn_C1D_ch /= phase_C;
    return {matchError_F, cIn_C1D_ch, cOut_C1D_ch};
}

/**
 * @brief  Integrate norm by adaptive RK4 and trapezoids.
 * @math   Γ_ch=2P_αh_μk_ch|c_{out,ch}|²/||u||²
 * @output Partial decay widths by channel.
 * @note   Uses P_α=0.39.
 * @note   Adaptive endpoint-controlled RK4 loops are approved.
 */
inline Eigen::VectorXd decay_width(double rmin_F, double rmatch_F, double rmax_F, double dr_F, const CCParams& params)
{
    int Nch_I = static_cast<int>(params.channel_1D_ch.size());
    assert(Nch_I > 0 && 0.0 < rmin_F && rmin_F < rmatch_F && rmatch_F < rmax_F && dr_F > 0.0);

    auto [matchError_F, cIn_C1D_ch, cOut_C1D_ch] = decay_match(rmin_F, rmatch_F, rmax_F, dr_F, params);
    (void) matchError_F;
    Eigen::MatrixXcd boundary_C2D_ch_2sol = decay_regular_boundary(rmin_F, params);
    Eigen::MatrixXcd y_C2D_ch_2sol(Nch_I, 2);
    y_C2D_ch_2sol.leftCols(1).noalias() = boundary_C2D_ch_2sol.leftCols(Nch_I) * cIn_C1D_ch;
    y_C2D_ch_2sol.rightCols(1).noalias() = boundary_C2D_ch_2sol.rightCols(Nch_I) * cIn_C1D_ch;
    double norm_F = 0.0;
    double density_F = y_C2D_ch_2sol.col(0).squaredNorm();
    double r_F = rmin_F;
    while (rmatch_F - r_F > 1.0e-12) {
        double drstep_F = std::min({dr_F, rmatch_F - r_F, 0.1 * r_F});
        y_C2D_ch_2sol = decay_rk4_step(y_C2D_ch_2sol, r_F, drstep_F, params);
        double nextDensity_F = y_C2D_ch_2sol.col(0).squaredNorm();
        norm_F += 0.5 * drstep_F * (density_F + nextDensity_F);
        density_F = nextDensity_F;
        r_F += drstep_F;
    }

    boundary_C2D_ch_2sol = decay_coulomb_boundary(rmax_F, params);
    y_C2D_ch_2sol.leftCols(1).noalias() = boundary_C2D_ch_2sol.leftCols(Nch_I) * cOut_C1D_ch;
    y_C2D_ch_2sol.rightCols(1).noalias() = boundary_C2D_ch_2sol.rightCols(Nch_I) * cOut_C1D_ch;
    density_F = y_C2D_ch_2sol.col(0).squaredNorm();
    r_F = rmax_F;
    while (r_F - rmatch_F > 1.0e-12) {
        double drstep_F = -std::min({dr_F, r_F - rmatch_F, 0.1 * r_F});
        y_C2D_ch_2sol = decay_rk4_step(y_C2D_ch_2sol, r_F, drstep_F, params);
        double nextDensity_F = y_C2D_ch_2sol.col(0).squaredNorm();
        norm_F -= 0.5 * drstep_F * (density_F + nextDensity_F);
        density_F = nextDensity_F;
        r_F += drstep_F;
    }
    assert(std::isfinite(norm_F) && norm_F > 0.0);

    constexpr double Palpha_F = 0.39;
    Eigen::VectorXd Gamma_F1D_ch(Nch_I);
    for (int channel_I = 0; channel_I < Nch_I; ++channel_I) {
        assert(params.channel_1D_ch[channel_I].Ech_F > 0.0);
        double k_F = std::sqrt(params.channel_1D_ch[channel_I].Ech_F / params.Varg_params.hmass_F);
        Gamma_F1D_ch(channel_I) = 2.0 * Palpha_F * params.Varg_params.hmass_F * k_F * std::norm(cOut_C1D_ch(channel_I)) / norm_F;
    }
    return Gamma_F1D_ch;
}

/**
 * @brief  Count significant sign changes after adaptive RK4.
 * @math   N_{ch}=Σ_r 1[u_{r-1,ch}u_{r,ch}<0]
 * @output Internal radial-node count by channel.
 * @note   Ignores crossings below 10⁻⁸ of channel maxima.
 * @note   Adaptive endpoint-controlled RK4 loops are approved.
 */
inline Eigen::VectorXi decay_nodes(double rmin_F, double rmatch_F, double rmax_F, double dr_F, const CCParams& params)
{
    int Nch_I = static_cast<int>(params.channel_1D_ch.size());
    assert(Nch_I > 0 && 0.0 < rmin_F && rmin_F < rmatch_F && rmatch_F < rmax_F && dr_F > 0.0);

    auto [matchError_F, cIn_C1D_ch, cOut_C1D_ch] = decay_match(rmin_F, rmatch_F, rmax_F, dr_F, params);
    (void) matchError_F;
    (void) cOut_C1D_ch;
    Eigen::VectorXcd cin_C1D_ch = cIn_C1D_ch;
    Eigen::Index maxCoefficient_I;
    cin_C1D_ch.cwiseAbs().maxCoeff(&maxCoefficient_I);
    assert(std::abs(cin_C1D_ch(maxCoefficient_I)) > 0.0);
    cin_C1D_ch /= cin_C1D_ch(maxCoefficient_I) / std::abs(cin_C1D_ch(maxCoefficient_I));

    Eigen::MatrixXcd boundary_C2D_ch_2sol = decay_regular_boundary(rmin_F, params);
    Eigen::MatrixXcd y_C2D_ch_2sol(Nch_I, 2);
    y_C2D_ch_2sol.leftCols(1).noalias() = boundary_C2D_ch_2sol.leftCols(Nch_I) * cin_C1D_ch;
    y_C2D_ch_2sol.rightCols(1).noalias() = boundary_C2D_ch_2sol.rightCols(Nch_I) * cin_C1D_ch;
    std::vector<Eigen::VectorXd> u_F2D_r_ch;
    u_F2D_r_ch.reserve(static_cast<std::size_t>(std::ceil((rmatch_F - rmin_F) / std::min(dr_F, 0.1 * rmin_F))) + 1);
    u_F2D_r_ch.push_back(y_C2D_ch_2sol.col(0).real());
    double r_F = rmin_F;
    while (rmatch_F - r_F > 1.0e-12) {
        double drstep_F = std::min({dr_F, rmatch_F - r_F, 0.1 * r_F});
        y_C2D_ch_2sol = decay_rk4_step(y_C2D_ch_2sol, r_F, drstep_F, params);
        u_F2D_r_ch.push_back(y_C2D_ch_2sol.col(0).real());
        r_F += drstep_F;
    }

    Eigen::VectorXd umax_F1D_ch = Eigen::VectorXd::Zero(Nch_I);
    for (const Eigen::VectorXd& u_F1D_ch : u_F2D_r_ch) {
        umax_F1D_ch = umax_F1D_ch.cwiseMax(u_F1D_ch.cwiseAbs());
    }

    Eigen::VectorXi nodes_I1D_ch = Eigen::VectorXi::Zero(Nch_I);
    constexpr double nodeTolerance_F = 1.0e-8;
    for (std::size_t r_I = 1; r_I < u_F2D_r_ch.size(); ++r_I) {
        for (int channel_I = 0; channel_I < Nch_I; ++channel_I) {
            double uprev_F = u_F2D_r_ch[r_I - 1](channel_I);
            double u_F = u_F2D_r_ch[r_I](channel_I);
            nodes_I1D_ch(channel_I) += static_cast<int>((uprev_F * u_F < 0.0) & (std::max(std::abs(uprev_F), std::abs(u_F)) > nodeTolerance_F * umax_F1D_ch(channel_I)));
        }
    }
    return nodes_I1D_ch;
}
