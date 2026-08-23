/**
 * @file    ivp_logderiv.hpp
 * @author  cailea
 * @date    2026-08-19
 * @brief   QR-stabilized matrix log-derivative propagation.
 */

#pragma once

#include <cassert>
#include <cmath>
#include <functional>

#include <Eigen/Dense>
#include <Eigen/QR>
#include <unsupported/Eigen/CXX11/Tensor>

#include "ivp_rungekutta.hpp"

template<typename T> using Real2TMatFunc = std::function<void(double, Eigen::Ref<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>)>;

/**
 * @brief  Propagate Y by RK4 with QR reorthogonalization.
 * @math   Ψ''=FΨ, Y=Ψ'Ψ⁻¹
 * @output Log-derivative tensor Y_T3D_ch_ch_x on the input grid.
 * @note   Uses Ψ(x₀)=I and Ψ'(x₀)=Y₀.
 * @note   Requires nonsingular Ψ at every output point.
 */
template<typename T>
Eigen::Tensor<T, 3, Eigen::ColMajor> ivp_logderiv_rk4(const Real2TMatFunc<T>& F_Func, const Eigen::Ref<const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>& Y0_T2D_ch_ch, const Eigen::Ref<const Eigen::VectorXd>& x_F1D_x) {
    // (N_x,N_ch) → reusable RK4/QR buffers.
    int Nx_I = static_cast<int>(x_F1D_x.size());
    int Nch_I = static_cast<int>(Y0_T2D_ch_ch.rows());
    assert(Nx_I >= 2 && x_F1D_x.allFinite());
    assert(Nch_I > 0 && Y0_T2D_ch_ch.cols() == Nch_I && Y0_T2D_ch_ch.allFinite());
    int Nchsol_I = Nch_I * Nch_I;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> F_T2D_ch_ch(Nch_I, Nch_I);
    Eigen::HouseholderQR<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> qr_(Nch_I, Nch_I);
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Q_T2D_ch_sol(Nch_I, Nch_I);
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> R_T2D_sol_sol(Nch_I, Nch_I);
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> dpsiRinv_T2D_ch_sol(Nch_I, Nch_I);

    // u = {vec(Ψ), vec(Ψ')}.
    Eigen::Vector<T, Eigen::Dynamic> u0_T1D_2chsol(2 * Nchsol_I);
    {
        Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> u0_T2D_2ch_sol(u0_T1D_2chsol.data(), 2 * Nch_I, Nch_I);
        u0_T2D_2ch_sol.topRows(Nch_I).setIdentity();
        u0_T2D_2ch_sol.bottomRows(Nch_I) = Y0_T2D_ch_ch;
    }

    // u' = {vec(Ψ'), vec(FΨ)}.
    RealTVec2TVecFunc<T> dudx_Func = [&](double x_F, const Eigen::Ref<const Eigen::Vector<T, Eigen::Dynamic>>& u_T1D_2chsol, Eigen::Ref<Eigen::Vector<T, Eigen::Dynamic>> dudx_T1D_2chsol) {
        Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> u_T2D_2ch_sol(u_T1D_2chsol.data(), 2 * Nch_I, Nch_I);
        F_Func(x_F, F_T2D_ch_ch);
        assert(F_T2D_ch_ch.allFinite());
        Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> dudx_T2D_2ch_sol(dudx_T1D_2chsol.data(), 2 * Nch_I, Nch_I);
        dudx_T2D_2ch_sol.topRows(Nch_I) = u_T2D_2ch_sol.bottomRows(Nch_I);
        dudx_T2D_2ch_sol.bottomRows(Nch_I).noalias() = F_T2D_ch_ch * u_T2D_2ch_sol.topRows(Nch_I);
    };
    IVP_RK4State<T> rk4_State(dudx_Func, x_F1D_x(0), u0_T1D_2chsol);

    // Y(x₀) = Y₀.
    Eigen::Tensor<T, 3, Eigen::ColMajor> Y_T3D_ch_ch_x(Nch_I, Nch_I, Nx_I);
    Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> Yx_T2D_ch_ch(Y_T3D_ch_ch_x.data(), Nch_I, Nch_I);
    Yx_T2D_ch_ch = Y0_T2D_ch_ch;

    // (Ψ_i,Ψ'_i) → (Ψ_{i+1},Ψ'_{i+1}) → Y_{i+1}.
    for (int x_I = 0; x_I < Nx_I - 1; ++x_I) {
        Eigen::Vector<T, Eigen::Dynamic>& u_T1D_2chsol = rk4_State.step(x_F1D_x(x_I + 1));

        // Ψ = QR, Ψ' → Ψ'R⁻¹.
        Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> u_T2D_2ch_sol(u_T1D_2chsol.data(), 2 * Nch_I, Nch_I);
        qr_.compute(u_T2D_2ch_sol.topRows(Nch_I));
        Q_T2D_ch_sol = qr_.householderQ();
        R_T2D_sol_sol = qr_.matrixQR().topLeftCorner(Nch_I, Nch_I).template triangularView<Eigen::Upper>();
        assert(R_T2D_sol_sol.diagonal().array().abs().minCoeff() > 0.0);
        dpsiRinv_T2D_ch_sol = R_T2D_sol_sol.transpose().template triangularView<Eigen::Lower>().solve(u_T2D_2ch_sol.bottomRows(Nch_I).transpose()).transpose();
        u_T2D_2ch_sol.topRows(Nch_I) = Q_T2D_ch_sol;
        u_T2D_2ch_sol.bottomRows(Nch_I) = dpsiRinv_T2D_ch_sol;

        // Y = Ψ'R⁻¹Q†.
        Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> Yx_T2D_ch_ch(Y_T3D_ch_ch_x.data() + (x_I + 1) * Nchsol_I, Nch_I, Nch_I);
        Yx_T2D_ch_ch.noalias() = dpsiRinv_T2D_ch_sol * Q_T2D_ch_sol.adjoint();
    }
    return Y_T3D_ch_ch_x;
}
