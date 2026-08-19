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
template<typename T> using RealTVec2TVecFunc = std::function<void(double, const Eigen::Ref<const Eigen::Vector<T, Eigen::Dynamic>>&, Eigen::Ref<Eigen::Vector<T, Eigen::Dynamic>>)>;

// ==================== Vector RK4 state ====================

template<typename T>
class IVP_RK4State {
public:
    double tcurr_F;
private:
    RealTVec2TVecFunc<T> f_Func;
    Eigen::Vector<T, Eigen::Dynamic> ycurr_T1D_i;
    Eigen::Vector<T, Eigen::Dynamic> k1_T1D_i;
    Eigen::Vector<T, Eigen::Dynamic> k2_T1D_i;
    Eigen::Vector<T, Eigen::Dynamic> k3_T1D_i;
    Eigen::Vector<T, Eigen::Dynamic> k4_T1D_i;
    Eigen::Vector<T, Eigen::Dynamic> ytmp_T1D_i;

public:
    /**
     * @brief  Initialize a vector RK4 recurrence.
     * @math   {t₀,y₀}
     * @output Initialized RK4 state.
     */
    IVP_RK4State(const RealTVec2TVecFunc<T>& f_Func, double t0_F, const Eigen::Ref<const Eigen::Vector<T, Eigen::Dynamic>>& y0_T1D_i)
    : f_Func(f_Func), tcurr_F(t0_F), ycurr_T1D_i(y0_T1D_i), k1_T1D_i(y0_T1D_i.size()), k2_T1D_i(y0_T1D_i.size()), k3_T1D_i(y0_T1D_i.size()), k4_T1D_i(y0_T1D_i.size()), ytmp_T1D_i(y0_T1D_i.size()) {
        assert(f_Func);
        assert(std::isfinite(t0_F));
        assert(y0_T1D_i.size() > 0 && y0_T1D_i.allFinite());
    }

    /**
     * @brief  Advance the vector state by one RK4 step.
     * @math   yₙ₊₁=yₙ+dt(k₁+2k₂+2k₃+k₄)/6
     * @output Mutable reference to the updated state yₙ₊₁.
     * @note   The returned reference is owned by this state.
     */
    Eigen::Vector<T, Eigen::Dynamic>& step(double tnext_F);
};

template<typename T>
Eigen::Vector<T, Eigen::Dynamic>& IVP_RK4State<T>::step(double tnext_F) {
    double dt_F = tnext_F - tcurr_F;
    assert(std::isfinite(tnext_F) && std::isfinite(dt_F) && dt_F != 0.0);

    // k₁ = f(tₙ, yₙ).
    f_Func(tcurr_F, ycurr_T1D_i, k1_T1D_i);

    // k₂ = f(tₙ + dt/2, yₙ + dt k₁/2).
    ytmp_T1D_i = ycurr_T1D_i + 0.5 * dt_F * k1_T1D_i;
    f_Func(tcurr_F + 0.5 * dt_F, ytmp_T1D_i, k2_T1D_i);

    // k₃ = f(tₙ + dt/2, yₙ + dt k₂/2).
    ytmp_T1D_i = ycurr_T1D_i + 0.5 * dt_F * k2_T1D_i;
    f_Func(tcurr_F + 0.5 * dt_F, ytmp_T1D_i, k3_T1D_i);

    // k₄ = f(tₙ + dt, yₙ + dt k₃).
    ytmp_T1D_i = ycurr_T1D_i + dt_F * k3_T1D_i;
    f_Func(tnext_F, ytmp_T1D_i, k4_T1D_i);

    // yₙ₊₁ = yₙ + dt(k₁ + 2k₂ + 2k₃ + k₄)/6.
    ycurr_T1D_i += (dt_F / 6.0) * (k1_T1D_i + T(2.0) * k2_T1D_i + T(2.0) * k3_T1D_i + k4_T1D_i);
    tcurr_F = tnext_F;
    return ycurr_T1D_i;
}

// ==================== Scalar RK4 interface ====================

/**
 * @brief  Advance one scalar initial-value problem by one RK4 step.
 * @math   yₙ₊₁=yₙ+dt(k₁+2k₂+2k₃+k₄)/6
 * @output Scalar value yₙ₊₁.
 */
template<typename T>
T ivp_rk4_step(const RealT2TFunc<T>& f_Func, double t_F, T y_T, double dt_F) {
    assert(std::isfinite(t_F) && std::isfinite(dt_F) && dt_F != 0.0);
    T k1_T = f_Func(t_F, y_T);
    T k2_T = f_Func(t_F + 0.5 * dt_F, y_T + 0.5 * dt_F * k1_T);
    T k3_T = f_Func(t_F + 0.5 * dt_F, y_T + 0.5 * dt_F * k2_T);
    T k4_T = f_Func(t_F + dt_F, y_T + dt_F * k3_T);
    return y_T + (dt_F / 6.0) * (k1_T + T(2.0) * k2_T + T(2.0) * k3_T + k4_T);
}

/**
 * @brief  Solve a scalar initial-value problem on a prescribed grid with RK4.
 * @math   y(t₀)=y₀, y(tᵢ₊₁)=RK4[f,tᵢ,y(tᵢ),tᵢ₊₁-tᵢ]
 * @output Scalar solution y_T1D_t on t_F1D_t.
 */
template<typename T>
Eigen::Vector<T, Eigen::Dynamic> ivp_rk4(const RealT2TFunc<T>& f_Func, T y0_T, const Eigen::Ref<const Eigen::VectorXd>& t_F1D_t) {
    int Nt_I = static_cast<int>(t_F1D_t.size());
    assert(Nt_I >= 2 && t_F1D_t.allFinite());
    Eigen::Vector<T, Eigen::Dynamic> y_T1D_t(Nt_I);
    y_T1D_t(0) = y0_T;
    for (int t_I = 0; t_I < Nt_I - 1; ++t_I) {
        double dt_F = t_F1D_t(t_I + 1) - t_F1D_t(t_I);
        y_T1D_t(t_I + 1) = ivp_rk4_step(f_Func, t_F1D_t(t_I), y_T1D_t(t_I), dt_F);
    }
    return y_T1D_t;
}

// ==================== Vector RK4 interface ====================

/**
 * @brief  Solve a vector initial-value problem on a prescribed grid with RK4.
 * @math   y(t₀)=y₀, y(tᵢ₊₁)=RK4[f,tᵢ,y(tᵢ),tᵢ₊₁-tᵢ]
 * @output State-by-grid solution y_T2D_i_t.
 */
template<typename T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> ivp_rk4_vec(const RealTVec2TVecFunc<T>& f_Func, const Eigen::Ref<const Eigen::VectorXd>& t_F1D_t, const Eigen::Ref<const Eigen::Vector<T, Eigen::Dynamic>>& y0_T1D_i) {
    int Nt_I = static_cast<int>(t_F1D_t.size());
    int Ni_I = static_cast<int>(y0_T1D_i.size());
    assert(Nt_I >= 2 && t_F1D_t.allFinite());
    assert(Ni_I > 0 && y0_T1D_i.allFinite());
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> y_T2D_i_t(Ni_I, Nt_I);
    IVP_RK4State<T> rk4_State(f_Func, t_F1D_t(0), y0_T1D_i);
    y_T2D_i_t.col(0) = y0_T1D_i;
    for (int t_I = 0; t_I < Nt_I - 1; ++t_I) {y_T2D_i_t.col(t_I + 1) = rk4_State.step(t_F1D_t(t_I + 1));}
    return y_T2D_i_t;
}
