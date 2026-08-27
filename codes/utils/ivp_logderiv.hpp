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

template<typename T>
class IVP_RK4QRState {
public:
    double xcurr_F;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Qcurr_T2D_ch_sol;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> dQcurr_T2D_ch_sol;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Rcurr_T2D_sol_sol;
private:
    IVP_RK4State<T> rk4_State;
    Eigen::HouseholderQR<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> qrSolver_QR;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Rtmp_T2D_sol_sol;

public:
    /**
     * @brief  Initialize a QR-stabilized matrix RK4 recurrence.
     * @math   u₀=Q₀R₀, du₀=dQ₀R₀
     * @output Initialized QR-stabilized RK4 state.
     * @note   Requires full-column-rank u₀.
     */
    IVP_RK4QRState(const Real2TMatFunc<T>& F_Func_, double x0_F, const Eigen::Ref<const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>& u0_T2D_ch_sol, const Eigen::Ref<const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>& du0_T2D_ch_sol)
    : xcurr_F(x0_F), Qcurr_T2D_ch_sol(u0_T2D_ch_sol.rows(), u0_T2D_ch_sol.cols()), dQcurr_T2D_ch_sol(du0_T2D_ch_sol.rows(), du0_T2D_ch_sol.cols()), Rcurr_T2D_sol_sol(u0_T2D_ch_sol.cols(), u0_T2D_ch_sol.cols()), rk4_State(), qrSolver_QR(u0_T2D_ch_sol.rows(), u0_T2D_ch_sol.cols()), Rtmp_T2D_sol_sol(u0_T2D_ch_sol.cols(), u0_T2D_ch_sol.cols()) {
        int Nch_I = static_cast<int>(Qcurr_T2D_ch_sol.rows());
        int Nsol_I = static_cast<int>(Qcurr_T2D_ch_sol.cols());
        assert(F_Func_ && Nch_I >= Nsol_I && Nsol_I > 0 && std::isfinite(x0_F));
        assert(du0_T2D_ch_sol.rows() == Nch_I && du0_T2D_ch_sol.cols() == Nsol_I);
        assert(u0_T2D_ch_sol.allFinite() && du0_T2D_ch_sol.allFinite());

        // u₀=Q₀R₀, du₀=dQ₀R₀.
        qrSolver_QR.compute(u0_T2D_ch_sol);
        Qcurr_T2D_ch_sol = qrSolver_QR.householderQ() * Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Identity(Nch_I, Nsol_I);
        Rcurr_T2D_sol_sol = qrSolver_QR.matrixQR().topLeftCorner(Nsol_I, Nsol_I).template triangularView<Eigen::Upper>();
        assert((Rcurr_T2D_sol_sol.diagonal().array().abs() > 0.0).all());
        dQcurr_T2D_ch_sol = Rcurr_T2D_sol_sol.transpose().template triangularView<Eigen::Lower>().solve(du0_T2D_ch_sol.transpose()).transpose();

        // (Q₀,dQ₀) → y₀=vec([Q₀,dQ₀]).
        Eigen::Vector<T, Eigen::Dynamic> y0_T1D_ch2sol(2 * Nch_I * Nsol_I);
        Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> y0_T2D_ch_2sol(y0_T1D_ch2sol.data(), Nch_I, 2 * Nsol_I);
        y0_T2D_ch_2sol.leftCols(Nsol_I) = Qcurr_T2D_ch_sol;
        y0_T2D_ch_2sol.rightCols(Nsol_I) = dQcurr_T2D_ch_sol;

        // [Q,dQ]' = [dQ,FQ].
        RealTVec2TVecFunc<T> dydx_Func = [F_Func_, Nch_I, Nsol_I, F_T2D_ch_ch = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(Nch_I, Nch_I)](double xstage_F, const Eigen::Ref<const Eigen::Vector<T, Eigen::Dynamic>>& y_T1D_ch2sol, Eigen::Ref<Eigen::Vector<T, Eigen::Dynamic>> dydx_T1D_ch2sol) mutable {
            Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> y_T2D_ch_2sol(y_T1D_ch2sol.data(), Nch_I, 2 * Nsol_I);
            Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> dydx_T2D_ch_2sol(dydx_T1D_ch2sol.data(), Nch_I, 2 * Nsol_I);
            F_Func_(xstage_F, F_T2D_ch_ch);
            assert(F_T2D_ch_ch.allFinite());
            dydx_T2D_ch_2sol.leftCols(Nsol_I) = y_T2D_ch_2sol.rightCols(Nsol_I);
            dydx_T2D_ch_2sol.rightCols(Nsol_I).noalias() = F_T2D_ch_ch * y_T2D_ch_2sol.leftCols(Nsol_I);
        };
        rk4_State = IVP_RK4State<T>(dydx_Func, x0_F, y0_T1D_ch2sol);
    }

    /**
     * @brief  Advance the matrix state by QR-stabilized RK4.
     * @math   (Q,dQ,R) → (Q⁺,dQ⁺,R⁺)
     * @output Mutable reference to the updated QR state.
     */
    IVP_RK4QRState<T>& step(double xnext_F);
};

template<typename T>
IVP_RK4QRState<T>& IVP_RK4QRState<T>::step(double xnext_F) {
    int Nch_I = static_cast<int>(Qcurr_T2D_ch_sol.rows());
    int Nsol_I = static_cast<int>(Qcurr_T2D_ch_sol.cols());
    Eigen::Vector<T, Eigen::Dynamic>& y_T1D_ch2sol = rk4_State.step(xnext_F);
    Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> y_T2D_ch_2sol(y_T1D_ch2sol.data(), Nch_I, 2 * Nsol_I);
    auto u_T2D_ch_sol = y_T2D_ch_2sol.leftCols(Nsol_I);
    auto du_T2D_ch_sol = y_T2D_ch_2sol.rightCols(Nsol_I);

    // ũ=Q⁺R_step, dQ⁺=dũR_step⁻¹, R⁺=R_stepR.
    qrSolver_QR.compute(u_T2D_ch_sol);
    Rtmp_T2D_sol_sol = qrSolver_QR.matrixQR().topLeftCorner(Nsol_I, Nsol_I).template triangularView<Eigen::Upper>();
    Qcurr_T2D_ch_sol.noalias() = qrSolver_QR.householderQ() * Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Identity(Nch_I, Nsol_I);
    dQcurr_T2D_ch_sol = Rtmp_T2D_sol_sol.transpose().template triangularView<Eigen::Lower>().solve(du_T2D_ch_sol.transpose()).transpose();
    Rcurr_T2D_sol_sol = Rtmp_T2D_sol_sol * Rcurr_T2D_sol_sol;
    u_T2D_ch_sol = Qcurr_T2D_ch_sol;
    du_T2D_ch_sol = dQcurr_T2D_ch_sol;
    xcurr_F = xnext_F;
    return *this;
}

/**
 * @brief  Propagate Y by RK4 with QR reorthogonalization.
 * @math   Ψ''=FΨ, Y=Ψ'Ψ⁻¹
 * @output Log-derivative tensor Y_T3D_ch_ch_x on the input grid.
 * @note   Uses Ψ(x₀)=I and Ψ'(x₀)=Y₀.
 * @note   Requires nonsingular Ψ at every output point.
 */
template<typename T>
Eigen::Tensor<T, 3, Eigen::ColMajor> ivp_logderiv_rk4(const Real2TMatFunc<T>& F_Func, const Eigen::Ref<const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>& Y0_T2D_ch_ch, const Eigen::Ref<const Eigen::VectorXd>& x_F1D_x) {
    int Nx_I = static_cast<int>(x_F1D_x.size());
    int Nch_I = static_cast<int>(Y0_T2D_ch_ch.rows());
    assert(Nx_I >= 2 && x_F1D_x.allFinite());
    assert(Nch_I > 0 && Y0_T2D_ch_ch.cols() == Nch_I && Y0_T2D_ch_ch.allFinite());
    int Nchch_I = Nch_I * Nch_I;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Q0_T2D_ch_sol = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Identity(Nch_I, Nch_I);
    IVP_RK4QRState<T> qr_State(F_Func, x_F1D_x(0), Q0_T2D_ch_sol, Y0_T2D_ch_ch);

    // Y(x₀) = Y₀.
    Eigen::Tensor<T, 3, Eigen::ColMajor> Y_T3D_ch_ch_x(Nch_I, Nch_I, Nx_I);
    Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> Yx_T2D_ch_ch(Y_T3D_ch_ch_x.data(), Nch_I, Nch_I);
    Yx_T2D_ch_ch = Y0_T2D_ch_ch;

    // (Q_i,dQ_i,R_i) → (Q_{i+1},dQ_{i+1},R_{i+1}) → Y_{i+1}.
    for (int x_I = 0; x_I < Nx_I - 1; ++x_I) {
        qr_State.step(x_F1D_x(x_I + 1));
        Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> Yx_T2D_ch_ch(Y_T3D_ch_ch_x.data() + (x_I + 1) * Nchch_I, Nch_I, Nch_I);
        Yx_T2D_ch_ch = qr_State.Qcurr_T2D_ch_sol.transpose().colPivHouseholderQr().solve(qr_State.dQcurr_T2D_ch_sol.transpose()).transpose();
    }
    return Y_T3D_ch_ch_x;
}
