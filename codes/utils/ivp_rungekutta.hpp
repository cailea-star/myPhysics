/**
 * @file    ivp_rungekutta.hpp
 * @author  cailea
 * @date    2026-01-10
 * @brief   Fourth-order Runge-Kutta solvers for scalar and vector initial-value problems.
 */

#pragma once

#include <cmath>
#include <functional>
#include <Eigen/Dense>

template<typename T> using RealT2TFunc = std::function<T(double, T)>;
template<typename T> using TVec2TVecFunc = std::function<Eigen::Matrix<T, Eigen::Dynamic, 1>(double, const Eigen::Ref<const Eigen::Matrix<T, Eigen::Dynamic, 1>>&)>;

// ==================== Scalar fourth-order Runge-Kutta solver ====================

/**
 * @brief  Advance one scalar initial-value problem by one RK4 step.
 * @math   k₁ = f(tₙ, yₙ), k₂ = f(tₙ + h/2, yₙ + hk₁/2), k₃ = f(tₙ + h/2, yₙ + hk₂/2), k₄ = f(tₙ + h, yₙ + hk₃); yₙ₊₁ = yₙ + h(k₁ + 2k₂ + 2k₃ + k₄)/6
 * @output Returns yₙ₊₁.
 */
template<typename T>
T ivp_rk4_step(const RealT2TFunc<T>& f_Func, double t_F, T y_T, double h_F) {
    T k1_T = f_Func(t_F, y_T);
    T k2_T = f_Func(t_F + 0.5 * h_F, y_T + 0.5 * h_F * k1_T);
    T k3_T = f_Func(t_F + 0.5 * h_F, y_T + 0.5 * h_F * k2_T);
    T k4_T = f_Func(t_F + h_F, y_T + h_F * k3_T);
    return y_T + (h_F / 6.0) * (k1_T + T(2.0) * k2_T + T(2.0) * k3_T + k4_T);
}

/**
 * @brief  Solve a scalar initial-value problem on a prescribed time grid with RK4.
 * @math   y(t₀) = y₀; y(tᵢ₊₁) = RK4[f, tᵢ, y(tᵢ), tᵢ₊₁ - tᵢ]
 * @output Returns y_T1D_t evaluated on t_F1D_t.
 */
template<typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> ivp_rk4(const RealT2TFunc<T>& f_Func, T y0_T, const Eigen::Ref<const Eigen::VectorXd>& t_F1D_t) {
    int Nt_I = static_cast<int>(t_F1D_t.size());
    Eigen::Matrix<T, Eigen::Dynamic, 1> y_T1D_t(Nt_I);
    y_T1D_t(0) = y0_T;
    for (int t_I = 0; t_I <= Nt_I - 2; ++t_I) {
        double h_F = t_F1D_t(t_I + 1) - t_F1D_t(t_I);
        y_T1D_t(t_I + 1) = ivp_rk4_step(f_Func, t_F1D_t(t_I), y_T1D_t(t_I), h_F);
    }
    return y_T1D_t;
}

// ==================== Vector fourth-order Runge-Kutta solver ====================

/**
 * @brief  Advance one vector initial-value problem by one RK4 step.
 * @math   k₁ = f(tₙ, yₙ), k₂ = f(tₙ + h/2, yₙ + hk₁/2), k₃ = f(tₙ + h/2, yₙ + hk₂/2), k₄ = f(tₙ + h, yₙ + hk₃); yₙ₊₁ = yₙ + h(k₁ + 2k₂ + 2k₃ + k₄)/6
 * @output Returns y_T1D_ch at the next time point.
 */
template<typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> ivp_rk4_vec_step(const TVec2TVecFunc<T>& f_Func, double t_F, const Eigen::Ref<const Eigen::Matrix<T, Eigen::Dynamic, 1>>& y_T1D_ch, double h_F) {
    Eigen::Matrix<T, Eigen::Dynamic, 1> k1_T1D_ch = f_Func(t_F, y_T1D_ch);
    Eigen::Matrix<T, Eigen::Dynamic, 1> k2_T1D_ch = f_Func(t_F + 0.5 * h_F, y_T1D_ch + 0.5 * h_F * k1_T1D_ch);
    Eigen::Matrix<T, Eigen::Dynamic, 1> k3_T1D_ch = f_Func(t_F + 0.5 * h_F, y_T1D_ch + 0.5 * h_F * k2_T1D_ch);
    Eigen::Matrix<T, Eigen::Dynamic, 1> k4_T1D_ch = f_Func(t_F + h_F, y_T1D_ch + h_F * k3_T1D_ch);
    return y_T1D_ch + (h_F / 6.0) * (k1_T1D_ch + T(2.0) * k2_T1D_ch + T(2.0) * k3_T1D_ch + k4_T1D_ch);
}

/**
 * @brief  Solve a vector initial-value problem on a prescribed time grid with RK4.
 * @math   y(t₀) = y₀; y(tᵢ₊₁) = RK4[f, tᵢ, y(tᵢ), tᵢ₊₁ - tᵢ]
 * @output Returns y_T2D_ch_t evaluated on t_F1D_t.
 */
template<typename T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> ivp_rk4_vec(const TVec2TVecFunc<T>& f_Func, const Eigen::Ref<const Eigen::Matrix<T, Eigen::Dynamic, 1>>& y0_T1D_ch, const Eigen::Ref<const Eigen::VectorXd>& t_F1D_t) {
    int Nt_I = static_cast<int>(t_F1D_t.size());
    int Nch_I = static_cast<int>(y0_T1D_ch.size());
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> y_T2D_ch_t(Nch_I, Nt_I);
    y_T2D_ch_t.col(0) = y0_T1D_ch;
    for (int t_I = 0; t_I <= Nt_I - 2; ++t_I) {
        double h_F = t_F1D_t(t_I + 1) - t_F1D_t(t_I);
        Eigen::Matrix<T, Eigen::Dynamic, 1> y_T1D_ch = y_T2D_ch_t.col(t_I);
        y_T2D_ch_t.col(t_I + 1) = ivp_rk4_vec_step<T>(f_Func, t_F1D_t(t_I), y_T1D_ch, h_F);
    }
    return y_T2D_ch_t;
}
