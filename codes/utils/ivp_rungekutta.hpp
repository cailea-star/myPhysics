/**
 * @file    ivp_rungekutta.hpp
 * @author  cailea
 * @date    2026-01-10
 * @brief   Fourth-order Runge-Kutta solvers for scalar and vector initial-value problems.
 */

#pragma once

#include <cassert>
#include <cmath>
#include <functional>

#include <Eigen/Dense>

template<typename T> using RealT2TFunc = std::function<T(double, T)>;
template<typename T> using TVec2TVecFunc = std::function<void(double, const Eigen::Ref<const Eigen::Matrix<T, Eigen::Dynamic, 1>>&, Eigen::Ref<Eigen::Matrix<T, Eigen::Dynamic, 1>>)>;

// ==================== Internal vector RK4 interface ====================

namespace rk4_detail {

/**
 * @brief  Advance one vector initial-value problem with reusable RK4 buffers.
 * @math   yₙ₊₁=yₙ+h(k₁+2k₂+2k₃+k₄)/6
 * @output Writes yₙ₊₁ to ynext_T1D_ch.
 * @note   rkwork_T2D_ch_1234tmp stores {k₁,k₂,k₃,k₄,y_tmp}.
 */
template<typename T>
void ivp_rk4_vec_step(const TVec2TVecFunc<T>& f_Func, double t_F, double dt_F, const Eigen::Ref<const Eigen::Matrix<T, Eigen::Dynamic, 1>>& y_T1D_ch, Eigen::Ref<Eigen::Matrix<T, Eigen::Dynamic, 1>> ynext_T1D_ch, Eigen::Ref<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> rkwork_T2D_ch_1234tmp) {
    int Nch_I = static_cast<int>(y_T1D_ch.size());
    assert(Nch_I > 0 && ynext_T1D_ch.size() == Nch_I);
    assert(rkwork_T2D_ch_1234tmp.rows() == Nch_I && rkwork_T2D_ch_1234tmp.cols() == 5);
    assert(std::isfinite(t_F) && std::isfinite(dt_F) && dt_F != 0.0);

    // k₁ = f(tₙ, yₙ).
    f_Func(t_F, y_T1D_ch, rkwork_T2D_ch_1234tmp.col(0));

    // k₂ = f(tₙ + h/2, yₙ + hk₁/2).
    rkwork_T2D_ch_1234tmp.col(4) = y_T1D_ch + 0.5 * dt_F * rkwork_T2D_ch_1234tmp.col(0);
    f_Func(t_F + 0.5 * dt_F, rkwork_T2D_ch_1234tmp.col(4), rkwork_T2D_ch_1234tmp.col(1));

    // k₃ = f(tₙ + h/2, yₙ + hk₂/2).
    rkwork_T2D_ch_1234tmp.col(4) = y_T1D_ch + 0.5 * dt_F * rkwork_T2D_ch_1234tmp.col(1);
    f_Func(t_F + 0.5 * dt_F, rkwork_T2D_ch_1234tmp.col(4), rkwork_T2D_ch_1234tmp.col(2));

    // k₄ = f(tₙ + h, yₙ + hk₃).
    rkwork_T2D_ch_1234tmp.col(4) = y_T1D_ch + dt_F * rkwork_T2D_ch_1234tmp.col(2);
    f_Func(t_F + dt_F, rkwork_T2D_ch_1234tmp.col(4), rkwork_T2D_ch_1234tmp.col(3));

    // yₙ₊₁ = yₙ + h(k₁ + 2k₂ + 2k₃ + k₄)/6.
    ynext_T1D_ch = y_T1D_ch + (dt_F / 6.0) * (rkwork_T2D_ch_1234tmp.col(0) + T(2.0) * rkwork_T2D_ch_1234tmp.col(1) + T(2.0) * rkwork_T2D_ch_1234tmp.col(2) + rkwork_T2D_ch_1234tmp.col(3));
}

}

// ==================== Scalar RK4 interface ====================

/**
 * @brief  Advance one scalar initial-value problem by one RK4 step.
 * @math   yₙ₊₁=yₙ+h(k₁+2k₂+2k₃+k₄)/6
 * @output Scalar value yₙ₊₁.
 */
template<typename T>
T ivp_rk4_step(const RealT2TFunc<T>& f_Func, double t_F, T y_T, double h_F) {
    assert(std::isfinite(t_F) && std::isfinite(h_F) && h_F != 0.0);
    T k1_T = f_Func(t_F, y_T);
    T k2_T = f_Func(t_F + 0.5 * h_F, y_T + 0.5 * h_F * k1_T);
    T k3_T = f_Func(t_F + 0.5 * h_F, y_T + 0.5 * h_F * k2_T);
    T k4_T = f_Func(t_F + h_F, y_T + h_F * k3_T);
    return y_T + (h_F / 6.0) * (k1_T + T(2.0) * k2_T + T(2.0) * k3_T + k4_T);
}

/**
 * @brief  Solve a scalar initial-value problem on a prescribed grid with RK4.
 * @math   y(t₀)=y₀, y(tᵢ₊₁)=RK4[f,tᵢ,y(tᵢ),tᵢ₊₁-tᵢ]
 * @output Scalar solution y_T1D_t on t_F1D_t.
 */
template<typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> ivp_rk4(const RealT2TFunc<T>& f_Func, T y0_T, const Eigen::Ref<const Eigen::VectorXd>& t_F1D_t) {
    int Nt_I = static_cast<int>(t_F1D_t.size());
    assert(Nt_I >= 2 && t_F1D_t.allFinite());
    Eigen::Matrix<T, Eigen::Dynamic, 1> y_T1D_t(Nt_I);
    y_T1D_t(0) = y0_T;
    for (int t_I = 0; t_I < Nt_I - 1; ++t_I) {
        double h_F = t_F1D_t(t_I + 1) - t_F1D_t(t_I);
        y_T1D_t(t_I + 1) = ivp_rk4_step(f_Func, t_F1D_t(t_I), y_T1D_t(t_I), h_F);
    }
    return y_T1D_t;
}

// ==================== Vector RK4 interface ====================

/**
 * @brief  Advance one vector initial-value problem by one RK4 step.
 * @math   yₙ₊₁=yₙ+h(k₁+2k₂+2k₃+k₄)/6
 * @output Channel vector yₙ₊₁.
 */
template<typename T>
void ivp_rk4_vec_step(const TVec2TVecFunc<T>& f_Func, double t_F, double dt_F, const Eigen::Ref<const Eigen::Matrix<T, Eigen::Dynamic, 1>>& y_T1D_ch, Eigen::Ref<Eigen::Matrix<T, Eigen::Dynamic, 1>> ynext_T1D_ch) {
    int Nch_I = static_cast<int>(y_T1D_ch.size());
    assert(ynext_T1D_ch.size() == Nch_I);
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> rkwork_T2D_ch_1234tmp(Nch_I, 5);
    rk4_detail::ivp_rk4_vec_step<T>(f_Func, t_F, dt_F, y_T1D_ch, ynext_T1D_ch, rkwork_T2D_ch_1234tmp);
}

/**
 * @brief  Solve a vector initial-value problem on a prescribed grid with RK4.
 * @math   y(t₀)=y₀, y(tᵢ₊₁)=RK4[f,tᵢ,y(tᵢ),tᵢ₊₁-tᵢ]
 * @output Channel-by-grid solution y_T2D_ch_t.
 */
template<typename T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> ivp_rk4_vec(const TVec2TVecFunc<T>& f_Func, const Eigen::Ref<const Eigen::VectorXd>& t_F1D_t, const Eigen::Ref<const Eigen::Matrix<T, Eigen::Dynamic, 1>>& y0_T1D_ch) {
    int Nt_I = static_cast<int>(t_F1D_t.size());
    int Nch_I = static_cast<int>(y0_T1D_ch.size());
    assert(Nt_I >= 2 && t_F1D_t.allFinite());
    assert(Nch_I > 0 && y0_T1D_ch.allFinite());
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> y_T2D_ch_t(Nch_I, Nt_I);
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> rkwork_T2D_ch_1234tmp(Nch_I, 5);
    y_T2D_ch_t.col(0) = y0_T1D_ch;
    for (int t_I = 0; t_I < Nt_I - 1; ++t_I) {
        double h_F = t_F1D_t(t_I + 1) - t_F1D_t(t_I);
        rk4_detail::ivp_rk4_vec_step<T>(f_Func, t_F1D_t(t_I), h_F, y_T2D_ch_t.col(t_I), y_T2D_ch_t.col(t_I + 1), rkwork_T2D_ch_1234tmp);
    }
    return y_T2D_ch_t;
}
