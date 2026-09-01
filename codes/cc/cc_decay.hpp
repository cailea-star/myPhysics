/**
 * @file    cc_decay.hpp
 * @author  cailea
 * @date    2026-08-25
 * @brief   Coupled-channel decay matching, widths, and node counts.
 */

#pragma once

#include <Eigen/Core>
#include <Eigen/QR>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <Eigen/SVD>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <tuple>
#include <utility>

#include "basis_bspline.hpp"
#include "bvp_rmatrix.hpp"
#include "cc_params.hpp"
#include "cc_couplings.hpp"
#include "ivp_logderiv.hpp"
#include "spherical_radial_boundary.hpp"

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

/**
 * @brief  Match log derivatives by QR-RK4 and SVD.
 * @math   Y=dQ Q⁻¹; (Yᵢₙ-Yₒᵤₜ)uₘ=0
 * @output Minimum mismatch and normalized matched wave.
 * @note   Requires a strictly increasing radial grid containing the match point.
 */
inline std::tuple<double, Eigen::MatrixXcd> decay_match(const CCParams& params, const Eigen::Ref<const Eigen::VectorXd>& r_F1D_r, int rmatch_I)
{
    int Nch_I = static_cast<int>(params.channel_1D_ch.size());
    const int Nr_I = static_cast<int>(r_F1D_r.size());
    assert(Nch_I > 0 && Nr_I > 2 && rmatch_I > 0 && rmatch_I < Nr_I - 1 && r_F1D_r(0) > 0.0);
    for (int r_I = 1; r_I < Nr_I; ++r_I) {assert(r_F1D_r(r_I) > r_F1D_r(r_I - 1));}
    const double rmin_F = r_F1D_r(0);
    const double rmax_F = r_F1D_r(Nr_I - 1);

    // r → F(r).
    Real2TMatFunc<std::complex<double>> F_Func = [&](double rstage_F, Eigen::Ref<Eigen::MatrixXcd> F_C2D_ch_ch) {
        F_C2D_ch_ch = calc_coupled_F_matrix(rstage_F, params);
    };

    Eigen::MatrixXcd inBoundary_C2D_ch_2sol = decay_regular_boundary(rmin_F, params);
    Eigen::MatrixXcd outBoundary_C2D_ch_2sol = decay_coulomb_boundary(rmax_F, params);

    // (Bᵢₙ,Bₒᵤₜ) → (Qᵢₙ,dQᵢₙ,Rᵢₙ,Qₒᵤₜ,dQₒᵤₜ,Rₒᵤₜ).
    IVP_RK4QRState<std::complex<double>> inState_QR(F_Func, rmin_F, inBoundary_C2D_ch_2sol.leftCols(Nch_I), inBoundary_C2D_ch_2sol.rightCols(Nch_I));
    for (int r_I = 1; r_I <= rmatch_I; ++r_I) {inState_QR.step(r_F1D_r(r_I));}
    IVP_RK4QRState<std::complex<double>> outState_QR(F_Func, rmax_F, outBoundary_C2D_ch_2sol.leftCols(Nch_I), outBoundary_C2D_ch_2sol.rightCols(Nch_I));
    for (int r_I = Nr_I - 2; r_I >= rmatch_I; --r_I) {outState_QR.step(r_F1D_r(r_I));}

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
    Eigen::VectorXcd cIn_C1D_sol = inState_QR.Rcurr_T2D_sol_sol.colPivHouseholderQr().solve(aIn_C1D_ch);
    Eigen::VectorXcd cOut_C1D_sol = outState_QR.Rcurr_T2D_sol_sol.colPivHouseholderQr().solve(aOut_C1D_ch);

    // arg(cₒᵤₜ,max)=0.
    Eigen::Index maxCoefficient_I;
    cOut_C1D_sol.cwiseAbs().maxCoeff(&maxCoefficient_I);
    std::complex<double> phase_C = cOut_C1D_sol(maxCoefficient_I) / std::abs(cOut_C1D_sol(maxCoefficient_I));
    cOut_C1D_sol /= phase_C;
    cIn_C1D_sol /= phase_C;

    // (cᵢₙ,cₒᵤₜ) → {u(rᵢ)}, N=∫||u||²dr.
    Eigen::MatrixXcd umatch_C2D_ch_r(Nch_I, Nr_I);
    Eigen::MatrixXcd y_C2D_ch_2(Nch_I, 2);
    y_C2D_ch_2.leftCols(1).noalias() = inBoundary_C2D_ch_2sol.leftCols(Nch_I) * cIn_C1D_sol;
    y_C2D_ch_2.rightCols(1).noalias() = inBoundary_C2D_ch_2sol.rightCols(Nch_I) * cIn_C1D_sol;
    IVP_RK4QRState<std::complex<double>> inWaveState_QR(F_Func, rmin_F, y_C2D_ch_2.leftCols(1), y_C2D_ch_2.rightCols(1));
    umatch_C2D_ch_r.col(0) = y_C2D_ch_2.col(0);
    double norm_F = 0.0;
    double density_F = y_C2D_ch_2.col(0).squaredNorm();
    for (int r_I = 1; r_I <= rmatch_I; ++r_I) {
        inWaveState_QR.step(r_F1D_r(r_I));
        const Eigen::VectorXcd u_C1D_ch = inWaveState_QR.Qcurr_T2D_ch_sol.col(0) * inWaveState_QR.Rcurr_T2D_sol_sol(0, 0);
        const double nextDensity_F = u_C1D_ch.squaredNorm();
        norm_F += 0.5 * (r_F1D_r(r_I) - r_F1D_r(r_I - 1)) * (density_F + nextDensity_F);
        density_F = nextDensity_F;
        umatch_C2D_ch_r.col(r_I) = u_C1D_ch;
    }

    y_C2D_ch_2.leftCols(1).noalias() = outBoundary_C2D_ch_2sol.leftCols(Nch_I) * cOut_C1D_sol;
    y_C2D_ch_2.rightCols(1).noalias() = outBoundary_C2D_ch_2sol.rightCols(Nch_I) * cOut_C1D_sol;
    IVP_RK4QRState<std::complex<double>> outWaveState_QR(F_Func, rmax_F, y_C2D_ch_2.leftCols(1), y_C2D_ch_2.rightCols(1));
    umatch_C2D_ch_r.col(Nr_I - 1) = y_C2D_ch_2.col(0);
    density_F = y_C2D_ch_2.col(0).squaredNorm();
    for (int r_I = Nr_I - 2; r_I >= rmatch_I; --r_I) {
        outWaveState_QR.step(r_F1D_r(r_I));
        const Eigen::VectorXcd u_C1D_ch = outWaveState_QR.Qcurr_T2D_ch_sol.col(0) * outWaveState_QR.Rcurr_T2D_sol_sol(0, 0);
        const double nextDensity_F = u_C1D_ch.squaredNorm();
        norm_F += 0.5 * (r_F1D_r(r_I + 1) - r_F1D_r(r_I)) * (density_F + nextDensity_F);
        density_F = nextDensity_F;
        if (r_I > rmatch_I) {umatch_C2D_ch_r.col(r_I) = u_C1D_ch;}
    }
    assert(std::isfinite(norm_F) && norm_F > 0.0);

    umatch_C2D_ch_r /= std::sqrt(norm_F);
    return {matchError_F, umatch_C2D_ch_r};
}

/**
 * @brief  Match the full Bloch matrix by sparse LU and inverse iteration.
 * @math   H_Bc=0; u_ch(r)=Σ_b c_ch,b B_b(r)
 * @output Signed minimum eigenvalue magnitude and normalized wave.
 */
inline std::tuple<double, Eigen::MatrixXcd> decay_match_R(const CCParams& params, RMatrix& r_matrix)
{
    const BSplineBasis& b_basis = r_matrix.b_basis;
    const int Nch_I = static_cast<int>(params.channel_1D_ch.size());
    const int Nb_I = static_cast<int>(b_basis.B_F2D_grid_b.cols());
    const int Nchb_I = Nch_I * Nb_I;
    const int Nr_I = static_cast<int>(b_basis.b_funcs.x_F1D_x.size());
    assert(Nch_I > 0 && b_basis.b_funcs.x_F1D_x.size() > 1);

    const double rmin_F = b_basis.b_funcs.x_F1D_x(0);
    const double rmax_F = b_basis.b_funcs.x_F1D_x(b_basis.b_funcs.x_F1D_x.size() - 1);
    const Eigen::MatrixXcd innerBoundary_C2D_ch_2sol = decay_regular_boundary(rmin_F, params);
    const Eigen::MatrixXcd outerBoundary_C2D_ch_2sol = decay_coulomb_boundary(rmax_F, params);
    const Eigen::VectorXcd lambdaIn_C1D_ch = innerBoundary_C2D_ch_2sol.rightCols(Nch_I).diagonal().array() / innerBoundary_C2D_ch_2sol.leftCols(Nch_I).diagonal().array();
    const Eigen::VectorXcd lambdaOut_C1D_ch = outerBoundary_C2D_ch_2sol.rightCols(Nch_I).diagonal().array() / outerBoundary_C2D_ch_2sol.leftCols(Nch_I).diagonal().array();

    // G_Σ⁻¹ → H_B=G_Σ⁻¹-Λ_Gδ_bBδ_b'B.
    const Eigen::Tensor<std::complex<double>, 4, Eigen::ColMajor>& Ginv_C4D_ch_b_ch_b = r_matrix.update_Ginv([&](double r_F) {return calc_coupled_F_matrix(r_F, params);}, lambdaIn_C1D_ch);
    const Eigen::Map<const Eigen::MatrixXcd> Ginv_C2D_chb_chb(Ginv_C4D_ch_b_ch_b.data(), Nchb_I, Nchb_I);
    const double imagScale_F = std::max({1.0, Ginv_C2D_chb_chb.real().cwiseAbs().maxCoeff(), lambdaOut_C1D_ch.real().cwiseAbs().maxCoeff()});
    assert(Ginv_C2D_chb_chb.imag().cwiseAbs().maxCoeff() <= 1.0e-12 * imagScale_F);
    assert(lambdaOut_C1D_ch.imag().cwiseAbs().maxCoeff() <= 1.0e-12 * imagScale_F);
    Eigen::MatrixXd HB_F2D_chb_chb = Ginv_C2D_chb_chb.real();
    HB_F2D_chb_chb.bottomRightCorner(Nch_I, Nch_I).diagonal() -= lambdaOut_C1D_ch.real();

    // H_B → sign(det H_B), min|λ_i|, c_ch,b.
    const Eigen::SparseMatrix<double> HB_F2D_chb_chb_sparse = HB_F2D_chb_chb.sparseView(0.0, 0.0);
    Eigen::SparseLU<Eigen::SparseMatrix<double>> HBsolver_SparseLU;
    HBsolver_SparseLU.compute(HB_F2D_chb_chb_sparse);
    assert(HBsolver_SparseLU.info() == Eigen::Success);
    constexpr int NinverseIteration_I = 20;
    Eigen::VectorXd coeff_F1D_chb = Eigen::VectorXd::LinSpaced(Nchb_I, 1.0, 2.0).normalized();
    Eigen::VectorXd nextCoeff_F1D_chb(Nchb_I);
    for (int iteration_I = 0; iteration_I < NinverseIteration_I; ++iteration_I) {
        nextCoeff_F1D_chb = HBsolver_SparseLU.solve(coeff_F1D_chb);
        assert(HBsolver_SparseLU.info() == Eigen::Success && nextCoeff_F1D_chb.allFinite() && nextCoeff_F1D_chb.norm() > 0.0);
        coeff_F1D_chb = nextCoeff_F1D_chb.normalized();
    }
    const Eigen::VectorXd HBcoeff_F1D_chb = HB_F2D_chb_chb * coeff_F1D_chb;
    const double matchError_F = HBsolver_SparseLU.signDeterminant() * HBcoeff_F1D_chb.norm();
    Eigen::MatrixXd coeff_F2D_ch_b = Eigen::Map<const Eigen::MatrixXd>(coeff_F1D_chb.data(), Nch_I, Nb_I);

    // c_ch,b → u_ch,g; ∫||u||²dr=1.
    Eigen::MatrixXd uquad_F2D_ch_grid = coeff_F2D_ch_b * b_basis.B_F2D_grid_b.transpose();
    double norm_F = (uquad_F2D_ch_grid.cwiseAbs2().colwise().sum().array() * b_basis.w_F1D_grid.transpose().array()).sum();
    assert(std::isfinite(norm_F) && norm_F > 0.0);
    coeff_F2D_ch_b /= std::sqrt(norm_F);

    Eigen::MatrixXd B_F2D_r_b(Nr_I, Nb_I);
    for (int r_I = 0; r_I < Nr_I; ++r_I) {
        for (int b_I = 0; b_I < Nb_I; ++b_I) {B_F2D_r_b(r_I, b_I) = b_basis.b_funcs.basis_function(b_basis.b_funcs.x_F1D_x(r_I), b_I)(0);}
    }
    Eigen::MatrixXcd umatch_C2D_ch_r = (coeff_F2D_ch_b * B_F2D_r_b.transpose()).cast<std::complex<double>>();
    return {matchError_F, umatch_C2D_ch_r};
}

/**
 * @brief  Calculate widths from a normalized matched wave.
 * @math   Γ_ch=2P_αh_μk_ch|u_ch(r_max)/G_ch(r_max)|²
 * @output Partial decay widths by channel.
 * @note   Uses preformation factor P_α=0.39.
 * @note   The last radial-grid point is the outer boundary.
 */
inline Eigen::VectorXd decay_width(const CCParams& params, const Eigen::Ref<const Eigen::VectorXd>& r_F1D_r, const Eigen::Ref<const Eigen::MatrixXcd>& umatch_C2D_ch_r)
{
    int Nch_I = static_cast<int>(params.channel_1D_ch.size());
    assert(Nch_I > 0 && r_F1D_r.size() > 1 && r_F1D_r.size() == umatch_C2D_ch_r.cols() && umatch_C2D_ch_r.rows() == Nch_I);
    const double rmax_F = r_F1D_r(r_F1D_r.size() - 1);
    assert(rmax_F > 0.0);
    const Eigen::MatrixXcd boundary_C2D_ch_2sol = decay_coulomb_boundary(rmax_F, params);
    const Eigen::VectorXcd G_C1D_ch = boundary_C2D_ch_2sol.leftCols(Nch_I).diagonal();
    const Eigen::VectorXcd cOut_C1D_ch = umatch_C2D_ch_r.rightCols(1).array() / G_C1D_ch.array();
    constexpr double Palpha_F = 0.39;
    Eigen::VectorXd Gamma_F1D_ch(Nch_I);

    // (E_ch,cₒᵤₜ,ch) → Γ_ch; ||u||=1.
    for (int channel_I = 0; channel_I < Nch_I; ++channel_I) {
        assert(params.channel_1D_ch[channel_I].Ech_F > 0.0);
        double k_F = std::sqrt(params.channel_1D_ch[channel_I].Ech_F / params.Varg_params.hmass_F);
        Gamma_F1D_ch(channel_I) = 2.0 * Palpha_F * params.Varg_params.hmass_F * k_F * std::norm(cOut_C1D_ch(channel_I));
    }
    return Gamma_F1D_ch;
}

/**
 * @brief  Count nodes from a normalized matched wave.
 * @math   N_ch=Σᵢ𝟙[uᵢ₋₁,ch uᵢ,ch<0]
 * @output Internal radial-node count by channel.
 * @note   Ignores amplitudes below 10⁻⁸ of maxima.
 */
inline Eigen::VectorXi decay_nodes(const CCParams& params, const Eigen::Ref<const Eigen::MatrixXcd>& umatch_C2D_ch_r)
{
    int Nch_I = static_cast<int>(params.channel_1D_ch.size());
    assert(Nch_I > 0 && umatch_C2D_ch_r.rows() == Nch_I && umatch_C2D_ch_r.cols() > 1);
    Eigen::MatrixXcd uNode_C2D_ch_r = umatch_C2D_ch_r;
    Eigen::Index maxChannel_I;
    Eigen::Index maxR_I;
    uNode_C2D_ch_r.cwiseAbs().maxCoeff(&maxChannel_I, &maxR_I);
    assert(std::abs(uNode_C2D_ch_r(maxChannel_I, maxR_I)) > 0.0);
    const std::complex<double> phase_C = uNode_C2D_ch_r(maxChannel_I, maxR_I) / std::abs(uNode_C2D_ch_r(maxChannel_I, maxR_I));
    const Eigen::MatrixXd u_F2D_ch_r = (uNode_C2D_ch_r / phase_C).real();

    // uₘₐₓ,ch=maxᵢ|uᵢ,ch|.
    const Eigen::VectorXd umax_F1D_ch = u_F2D_ch_r.cwiseAbs().rowwise().maxCoeff();

    Eigen::VectorXi nodes_I1D_ch = Eigen::VectorXi::Zero(Nch_I);
    constexpr double nodeTolerance_F = 1.0e-8;

    // N_ch=Σᵢ𝟙[uᵢ₋₁,ch uᵢ,ch<0]𝟙[|u|>εuₘₐₓ,ch].
    for (int r_I = 1; r_I < u_F2D_ch_r.cols(); ++r_I) {
        for (int channel_I = 0; channel_I < Nch_I; ++channel_I) {
            double uprev_F = u_F2D_ch_r(channel_I, r_I - 1);
            double u_F = u_F2D_ch_r(channel_I, r_I);
            nodes_I1D_ch(channel_I) += static_cast<int>((uprev_F * u_F < 0.0) & (std::max(std::abs(uprev_F), std::abs(u_F)) > nodeTolerance_F * umax_F1D_ch(channel_I)));
        }
    }
    return nodes_I1D_ch;
}
