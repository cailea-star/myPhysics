/**
 * @file    cc_couplings.hpp
 * @author  cailea
 * @date    2026-08-23
 * @brief   Coupled-channel operators, potentials, and radial equations.
 */

#pragma once

#include <Eigen/Core>
#include <Eigen/Eigenvalues>

#include <cassert>
#include <cmath>
#include <complex>
#include <numbers>
#include <vector>

#include <gsl/gsl_sf_coupling.h>

#include "cc_params.hpp"

/**
 * @brief  Evaluate reduced spherical harmonics via Wigner 3j.
 * @math   ⟨l₁||Y_λ||l₂⟩ = (-1)^l₁ √[(2l₁+1)(2λ+1)(2l₂+1)/(4π)] (l₁ λ l₂; 0 0 0)
 * @output Reduced spherical-harmonic matrix element.
 * @note   GSL arguments use doubled angular momenta.
 */
inline double calc_reduced_Y_matrix_element(const Channel& channel1, const Channel& channel2, int twolambda_I)
{
    assert(twolambda_I >= 0 && twolambda_I % 2 == 0);
    int twol1_I = channel1.twol_I;
    int twol2_I = channel2.twol_I;
    double phase_l1_F = std::pow(-1.0, twol1_I / 2);

    // √[(2l₁+1)(2λ+1)(2l₂+1)/(4π)].
    double Ynorm_F = std::sqrt((twol1_I + 1.0) * (twolambda_I + 1.0) * (twol2_I + 1.0) / (4.0 * std::numbers::pi));
    double wigner3j_F = gsl_sf_coupling_3j(twol1_I, twolambda_I, twol2_I, 0, 0, 0);
    double reducedY_F = phase_l1_F * Ynorm_F * wigner3j_F;
    return reducedY_F;
}

/**
 * @brief  Evaluate reduced Wigner-D elements via Wigner 3j.
 * @math   ⟨I₁K₁||D^λ_ν||I₂K₂⟩ = √[(2I₁+1)(2I₂+1)] (-1)^(I₁-K₁) (I₁ λ I₂; -K₁ ν K₂)
 * @output Reduced Wigner-D matrix element.
 * @note   GSL arguments use doubled angular momenta.
 */
inline double calc_reduced_D_matrix_element(const Channel& channel1, const Channel& channel2, int twolambda_I, int twonu_I = 0)
{
    assert(twolambda_I >= 0 && twolambda_I % 2 == 0);
    assert(std::abs(twonu_I) <= twolambda_I && (twolambda_I - twonu_I) % 2 == 0);
    int twoIt1_I = channel1.twoIt_I;
    int twoIt2_I = channel2.twoIt_I;
    int twoKt1_I = channel1.twoKt_I;
    int twoKt2_I = channel2.twoKt_I;
    bool isKAllowed_B = twoKt2_I + twonu_I == twoKt1_I;
    double phase_It1Kt1_F = std::pow(-1.0, (twoIt1_I - twoKt1_I) / 2);
    double Dnorm_F = std::sqrt((twoIt1_I + 1.0) * (twoIt2_I + 1.0));

    // (I₁ λ I₂; -K₁ ν K₂).
    double wigner3j_F = gsl_sf_coupling_3j(twoIt1_I, twolambda_I, twoIt2_I, -twoKt1_I, twonu_I, twoKt2_I);
    double reducedD_F = static_cast<double>(isKAllowed_B) * phase_It1Kt1_F * Dnorm_F * wigner3j_F;
    return reducedD_F;
}

/**
 * @brief  Evaluate one coupled-channel scalar matrix element.
 * @math   ⟨ch₁|O_λ|ch₂⟩ = δ_J O_λ ⟨l₁||Y_λ||l₂⟩⟨I₁K₁||D^λ_ν||I₂K₂⟩{l₁ I₁ J; I₂ l₂ λ}/√(2λ+1)
 * @output Coupled-channel matrix element.
 * @note   GSL arguments use doubled angular momenta.
 */
inline std::complex<double> calc_coupling_matrix_element(std::complex<double> Olambda_C, const Channel& channel1, const Channel& channel2, int twolambda_I, int twonu_I = 0)
{
    assert(twolambda_I >= 0 && twolambda_I % 2 == 0);
    assert(std::abs(twonu_I) <= twolambda_I && (twolambda_I - twonu_I) % 2 == 0);
    bool isJConserved_B = channel1.twoJ_I == channel2.twoJ_I;
    int twol1_I = channel1.twol_I;
    int twol2_I = channel2.twol_I;
    int twoIt1_I = channel1.twoIt_I;
    int twoIt2_I = channel2.twoIt_I;
    int twoJ1_I = channel1.twoJ_I;
    double lambdaNorm_F = 1.0 / std::sqrt(twolambda_I + 1.0);

    // {l₁ I₁ J; I₂ l₂ λ}.
    double wigner6j_F = gsl_sf_coupling_6j(twol1_I, twoIt1_I, twoJ1_I, twoIt2_I, twol2_I, twolambda_I);
    double reducedY_F = calc_reduced_Y_matrix_element(channel1, channel2, twolambda_I);
    double reducedD_F = calc_reduced_D_matrix_element(channel1, channel2, twolambda_I, twonu_I);
    std::complex<double> coupling_C = static_cast<double>(isJConserved_B) * lambdaNorm_F * wigner6j_F * reducedY_F * reducedD_F * Olambda_C;
    return coupling_C;
}

/**
 * @brief  Build a symmetric channel-space operator matrix.
 * @math   O_{ch₁ch₂} = Σ_{λ=0,2,4} ⟨ch₁|O_λ|ch₂⟩
 * @output Symmetric channel-space operator matrix.
 */
inline Eigen::MatrixXcd calc_coupling_matrix(std::complex<double> O0_C, std::complex<double> O2_C, std::complex<double> O4_C, const std::vector<Channel>& channels, int twonu_I = 0)
{
    int Nch_I = static_cast<int>(channels.size());
    Eigen::MatrixXcd O_C2D_ch_ch = Eigen::MatrixXcd::Zero(Nch_I, Nch_I);

    // O₁₂ = O₂₁.
    for (int channel1_I = 0; channel1_I < Nch_I; ++channel1_I) {
        for (int channel2_I = channel1_I; channel2_I < Nch_I; ++channel2_I) {
            std::complex<double> O12_C = calc_coupling_matrix_element(O0_C, channels[channel1_I], channels[channel2_I], 0);
            O12_C += calc_coupling_matrix_element(O2_C, channels[channel1_I], channels[channel2_I], 4, twonu_I);
            O12_C += calc_coupling_matrix_element(O4_C, channels[channel1_I], channels[channel2_I], 8, twonu_I);
            O_C2D_ch_ch(channel1_I, channel2_I) = O12_C;
            O_C2D_ch_ch(channel2_I, channel1_I) = O12_C;
        }
    }
    return O_C2D_ch_ch;
}

inline Eigen::MatrixXcd CCParams::calc_Vall_matrix(double r_F) const
{
    assert(std::isfinite(r_F) && r_F > 0.0);
    int Nch_I = static_cast<int>(channel_1D_ch.size());
    double Rcoul_F = Varg_params.Rcoul_p_F + Varg_params.Rcoul_t_F;
    assert(std::isfinite(Rcoul_F) && Rcoul_F > 0.0);
    double rlarge_F = std::fmax(r_F, Rcoul_F);
    double rsmall_F = std::fmin(r_F, Rcoul_F);

    // V_C,λ(r), λ ∈ {0, 2, 4}.
    double Vcoul0_F = calc_Vcoul(r_F, Varg_params) * std::sqrt(4.0 * std::numbers::pi);
    double Vcoul2_F = 3.0 / std::sqrt(5.0) * Varg_params.Ze2_F * std::pow(rsmall_F / rlarge_F, 2) / rlarge_F * Varg_params.beta2tcoul_F;
    double Vcoul4_F = 3.0 / std::sqrt(9.0) * Varg_params.Ze2_F * std::pow(rsmall_F / rlarge_F, 4) / rlarge_F * Varg_params.beta4tcoul_F;
    Eigen::MatrixXcd Vall_C2D_ch_ch = calc_coupling_matrix(Vcoul0_F, Vcoul2_F, Vcoul4_F, channel_1D_ch);

    // R_t,λ, λ ∈ {0, 2, 4}.
    double Rt0_F = Varg_params.Rt0_F * std::sqrt(4.0 * std::numbers::pi);
    double Rt2_F = Varg_params.Rt0_F * Varg_params.beta2t_F * std::sqrt(5.0);
    double Rt4_F = Varg_params.Rt0_F * Varg_params.beta4t_F * std::sqrt(9.0);
    Eigen::MatrixXd Rt_F2D_ch_ch = calc_coupling_matrix(Rt0_F, Rt2_F, Rt4_F, channel_1D_ch).real();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> RtSolver_Eig(Rt_F2D_ch_ch);
    assert(RtSolver_Eig.info() == Eigen::Success);
    Eigen::VectorXd Rt_F1D_eigen = RtSolver_Eig.eigenvalues();
    Eigen::MatrixXd U_F2D_ch_eigen = RtSolver_Eig.eigenvectors();
    Eigen::VectorXd Vnucl_F1D_eigen(Nch_I);
    for (int eigen_I = 0; eigen_I < Nch_I; ++eigen_I) {
        Vnucl_F1D_eigen(eigen_I) = calc_Vnucl(r_F, Varg_params, Rt_F1D_eigen(eigen_I));
    }

    // V_N = U diag[V_N(R_t)] U^T.
    Eigen::MatrixXd Vnucl_F2D_ch_ch(Nch_I, Nch_I);
    Vnucl_F2D_ch_ch.noalias() = U_F2D_ch_eigen * Vnucl_F1D_eigen.asDiagonal() * U_F2D_ch_eigen.transpose();
    Vall_C2D_ch_ch += Vnucl_F2D_ch_ch.cast<std::complex<double>>();
    return Vall_C2D_ch_ch;
}

/**
 * @brief  Build the fixed-partial-wave radial-equation matrix.
 * @math   F_{ch,ch'}(r) = [V_{ch,ch'}(r) + δ_{ch,ch'}(V_L(r)-E_ch)]/h_μ
 * @output Coupled radial-equation matrix in fm^-2.
 */
inline Eigen::MatrixXcd calc_coupled_F_matrix(double r_F, const CCParams& params, int L_I)
{
    assert(std::isfinite(r_F) && r_F > 0.0 && L_I >= 0);
    int Nch_I = static_cast<int>(params.channel_1D_ch.size());
    Eigen::MatrixXcd F_C2D_ch_ch = params.calc_Vall_matrix(r_F);

    // F_ch,ch += [V_L(r) - E_ch]/h_μ.
    for (int channel_I = 0; channel_I < Nch_I; ++channel_I) {
        F_C2D_ch_ch(channel_I, channel_I) += calc_Vcent(r_F, params.Varg_params, L_I) - params.channel_1D_ch[channel_I].Ech_F;
    }
    F_C2D_ch_ch /= params.Varg_params.hmass_F;
    return F_C2D_ch_ch;
}

/**
 * @brief  Build the channel-partial-wave radial-equation matrix.
 * @math   F_{ch,ch'}(r) = [V_{ch,ch'}(r) + δ_{ch,ch'}(V_{l_ch}(r)-E_ch)]/h_μ
 * @output Coupled radial-equation matrix in fm^-2.
 */
inline Eigen::MatrixXcd calc_coupled_F_matrix(double r_F, const CCParams& params)
{
    assert(std::isfinite(r_F) && r_F > 0.0);
    int Nch_I = static_cast<int>(params.channel_1D_ch.size());
    Eigen::MatrixXcd F_C2D_ch_ch = params.calc_Vall_matrix(r_F);

    // F_ch,ch += [V_lch(r) - E_ch]/h_μ.
    for (int channel_I = 0; channel_I < Nch_I; ++channel_I) {
        int lch_I = params.channel_1D_ch[channel_I].twol_I / 2;
        F_C2D_ch_ch(channel_I, channel_I) += calc_Vcent(r_F, params.Varg_params, lch_I) - params.channel_1D_ch[channel_I].Ech_F;
    }
    F_C2D_ch_ch /= params.Varg_params.hmass_F;
    return F_C2D_ch_ch;
}

/**
 * @brief  Evaluate the coupled first-order RK4 derivative.
 * @math   d(u,u')/dr = (u',Fu)
 * @output Coupled first-order state derivative.
 */
inline Eigen::MatrixXcd calc_coupled_rk4_derivative(double r_F, const CCParams& params, int L_I, const Eigen::Ref<const Eigen::MatrixXcd>& y_C2D_ch_2sol)
{
    assert(std::isfinite(r_F) && r_F > 0.0 && L_I >= 0);
    int Nch_I = static_cast<int>(params.channel_1D_ch.size());
    int Nsol_I = static_cast<int>(y_C2D_ch_2sol.cols() / 2);
    assert(y_C2D_ch_2sol.rows() == Nch_I && y_C2D_ch_2sol.cols() == 2 * Nsol_I && Nsol_I > 0);
    auto u_C2D_ch_sol = y_C2D_ch_2sol.leftCols(Nsol_I);
    auto dudr_C2D_ch_sol = y_C2D_ch_2sol.rightCols(Nsol_I);
    Eigen::MatrixXcd F_C2D_ch_ch = calc_coupled_F_matrix(r_F, params, L_I);
    Eigen::MatrixXcd dydr_C2D_ch_2sol(Nch_I, 2 * Nsol_I);
    dydr_C2D_ch_2sol.leftCols(Nsol_I) = dudr_C2D_ch_sol;
    dydr_C2D_ch_2sol.rightCols(Nsol_I).noalias() = F_C2D_ch_ch * u_C2D_ch_sol;
    return dydr_C2D_ch_2sol;
}
