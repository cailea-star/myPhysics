/**
 * @file    ivp_numerov.hpp
 * @author  cailea
 * @date    2026-01-06
 * @brief   Scalar and matrix Numerov solvers.
 */

#pragma once

#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <functional>
#include <unsupported/Eigen/CXX11/Tensor>

template<typename T> using Real2TMatFunc = std::function<void(double, Eigen::Ref<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>)>;

// ==================== Matrix Numerov state ====================

template<typename T>
class IVP_NumerovState {
private:
    Real2TMatFunc<T> F_Func;
    double xcurr_F;
    double dx_F;
    double dx2_F;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Fprev_T2D_ch_ch;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Fcurr_T2D_ch_ch;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Fnext_T2D_ch_ch;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> yprev_T2D_ch_sol;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> ycurr_T2D_ch_sol;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> ynext_T2D_ch_sol;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> lhs_T2D_ch_ch;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> rhs_T2D_ch_sol;
    Eigen::PartialPivLU<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> lhs_LU;

public:
    /**
     * @brief  Initialize a matrix Numerov recurrence.
     * @math   {F,x₀,x₁,Y₀,Y₁}
     * @output Initialized Numerov state at x₁.
     */
    IVP_NumerovState(const Real2TMatFunc<T>& F_Func, double x0_F, double x1_F, const Eigen::Ref<const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>& y0_T2D_ch_sol, const Eigen::Ref<const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>& y1_T2D_ch_sol)
    : F_Func(F_Func), xcurr_F(x1_F), dx_F(x1_F - x0_F), dx2_F(dx_F * dx_F), Fprev_T2D_ch_ch(y0_T2D_ch_sol.rows(), y0_T2D_ch_sol.rows()), Fcurr_T2D_ch_ch(y0_T2D_ch_sol.rows(), y0_T2D_ch_sol.rows()), Fnext_T2D_ch_ch(y0_T2D_ch_sol.rows(), y0_T2D_ch_sol.rows()), yprev_T2D_ch_sol(y0_T2D_ch_sol), ycurr_T2D_ch_sol(y1_T2D_ch_sol), ynext_T2D_ch_sol(y0_T2D_ch_sol.rows(), y0_T2D_ch_sol.cols()), lhs_T2D_ch_ch(y0_T2D_ch_sol.rows(), y0_T2D_ch_sol.rows()), rhs_T2D_ch_sol(y0_T2D_ch_sol.rows(), y0_T2D_ch_sol.cols()), lhs_LU(y0_T2D_ch_sol.rows()) {
        assert(this->F_Func);
        assert(std::isfinite(x0_F) && std::isfinite(x1_F) && std::isfinite(dx_F) && dx_F != 0.0);
        assert(y0_T2D_ch_sol.rows() > 0 && y0_T2D_ch_sol.cols() > 0);
        assert(y1_T2D_ch_sol.rows() == y0_T2D_ch_sol.rows() && y1_T2D_ch_sol.cols() == y0_T2D_ch_sol.cols());
        this->F_Func(x0_F, Fprev_T2D_ch_ch);
        this->F_Func(x1_F, Fcurr_T2D_ch_ch);
        assert(Fprev_T2D_ch_ch.allFinite() && Fcurr_T2D_ch_ch.allFinite());
    }

    /**
     * @brief  Advance the matrix Numerov recurrence by one step.
     * @math   (I-dx²Fₙ₊₁/12)Yₙ₊₁=(2I+5dx²Fₙ/6)Yₙ-(I-dx²Fₙ₋₁/12)Yₙ₋₁
     * @output Constant reference to the updated solution Yₙ₊₁.
     * @note   The returned reference is owned by this state.
     */
    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& step(double xnext_F);
};

template<typename T>
const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& IVP_NumerovState<T>::step(double xnext_F) {
    assert(std::isfinite(xnext_F) && std::abs(xnext_F - xcurr_F - dx_F) <= 1.0e-10 * std::abs(dx_F));
    F_Func(xnext_F, Fnext_T2D_ch_ch);
    assert(Fnext_T2D_ch_ch.allFinite());

    // lhs = I - dx² Fₙ₊₁ / 12.
    lhs_T2D_ch_ch = -dx2_F * Fnext_T2D_ch_ch / T(12.0);
    lhs_T2D_ch_ch.diagonal().array() += T(1.0);

    // rhs = 5 dx² Fₙ Yₙ / 6 + dx² Fₙ₋₁ Yₙ₋₁ / 12 + 2 Yₙ - Yₙ₋₁.
    rhs_T2D_ch_sol.noalias() = T(5.0) * dx2_F * Fcurr_T2D_ch_ch * ycurr_T2D_ch_sol / T(6.0);
    rhs_T2D_ch_sol.noalias() += dx2_F * Fprev_T2D_ch_ch * yprev_T2D_ch_sol / T(12.0);
    rhs_T2D_ch_sol += T(2.0) * ycurr_T2D_ch_sol - yprev_T2D_ch_sol;

    // lhs Yₙ₊₁ = rhs.
    lhs_LU.compute(lhs_T2D_ch_ch);
    ynext_T2D_ch_sol = lhs_LU.solve(rhs_T2D_ch_sol);

    // {n-1,n,n+1} → {n,n+1,*}.
    Fprev_T2D_ch_ch.swap(Fcurr_T2D_ch_ch);
    Fcurr_T2D_ch_ch.swap(Fnext_T2D_ch_ch);
    yprev_T2D_ch_sol.swap(ycurr_T2D_ch_sol);
    ycurr_T2D_ch_sol.swap(ynext_T2D_ch_sol);
    xcurr_F = xnext_F;
    return ycurr_T2D_ch_sol;
}

/**
 * @brief  Check a uniform Numerov grid.
 * @math   h=x_{i+1}-x_i≠0
 * @output Uniform grid spacing.
 */
inline double ivp_numerov_assert_uniform_grid(const Eigen::Ref<const Eigen::VectorXd>& x_F1D_x) {
    int Nx_I = static_cast<int>(x_F1D_x.size());
    assert(Nx_I >= 2);
    double dx_F = x_F1D_x(1) - x_F1D_x(0);
    assert(std::isfinite(dx_F) && dx_F != 0.0);
    for (int x_I = 1; x_I < Nx_I - 1; ++x_I) {assert(std::abs(x_F1D_x(x_I + 1) - x_F1D_x(x_I) - dx_F) <= 1.0e-10 * std::abs(dx_F));}
    return dx_F;
}

// ==================== Scalar Numerov solver ====================

/**
 * @brief  Perform one scalar Numerov step.
 * @math   (1-h²F_{n+1}/12)y_{n+1}=2(1+5h²F_n/12)y_n-(1-h²F_{n-1}/12)y_{n-1}
 * @output Scalar value y_{n+1}.
 */
template<typename T>
T ivp_numerov_step(T Fprev_T, T Fcurr_T, T Fnext_T, T yprev_T, T ycurr_T, double dx_F) {
    assert(std::isfinite(dx_F) && dx_F != 0.0);
    double dx2_F = dx_F * dx_F;
    T lhs_T = T(1.0) - dx2_F * Fnext_T / T(12.0);
    T rhs_T = T(2.0) * (T(1.0) + T(5.0) * dx2_F * Fcurr_T / T(12.0)) * ycurr_T - (T(1.0) - dx2_F * Fprev_T / T(12.0)) * yprev_T;
    assert(std::abs(lhs_T) > 0.0);
    return rhs_T / lhs_T;
}

/**
 * @brief  Solve a scalar Numerov equation from sampled F.
 * @math   y''(x)=F(x)y(x)
 * @output Scalar solution y(x) on the input mesh.
 */
template<typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> ivp_numerov(const Eigen::Ref<const Eigen::Matrix<T, Eigen::Dynamic, 1>>& F_F1D_x, T y0_T, T y1_T, const Eigen::Ref<const Eigen::VectorXd>& x_F1D_x) {
    int Nx_I = static_cast<int>(x_F1D_x.size());
    double dx_F = ivp_numerov_assert_uniform_grid(x_F1D_x);
    assert(F_F1D_x.size() == Nx_I);

    Eigen::Matrix<T, Eigen::Dynamic, 1> y_T1D_x(Nx_I);
    y_T1D_x(0) = y0_T;
    y_T1D_x(1) = y1_T;
    T Fprev_T = F_F1D_x(0);
    T Fcurr_T = F_F1D_x(1);
    for (int x_I = 1; x_I <= Nx_I - 2; ++x_I) {
        T Fnext_T = F_F1D_x(x_I + 1);
        y_T1D_x(x_I + 1) = ivp_numerov_step(Fprev_T, Fcurr_T, Fnext_T, y_T1D_x(x_I - 1), y_T1D_x(x_I), dx_F);
        Fprev_T = Fcurr_T;
        Fcurr_T = Fnext_T;
    }
    return y_T1D_x;
}

/**
 * @brief  Solve a scalar Numerov equation from F(x).
 * @math   y''(x)=F(x)y(x)
 * @output Scalar solution y(x) on the input mesh.
 */
template<typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> ivp_numerov(const std::function<T(double)>& F_Func, T y0_T, T y1_T, const Eigen::Ref<const Eigen::VectorXd>& x_F1D_x) {
    int Nx_I = static_cast<int>(x_F1D_x.size());
    Eigen::Matrix<T, Eigen::Dynamic, 1> F_F1D_x(Nx_I);
    for (int x_I = 0; x_I < Nx_I; ++x_I) {F_F1D_x(x_I) = F_Func(x_F1D_x(x_I));}
    return ivp_numerov<T>(F_F1D_x, y0_T, y1_T, x_F1D_x);
}

// ==================== Matrix Numerov solver ====================

/**
 * @brief  Solve matrix Numerov equations from F(x).
 * @math   Y''(x)=F(x)Y(x)
 * @output Channel-by-solution-by-mesh solution tensor.
 */
template<typename T>
Eigen::Tensor<T, 3, Eigen::ColMajor> ivp_numerov_mat(const Real2TMatFunc<T>& F_Func, const Eigen::Ref<const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>& y0_T2D_ch_sol, const Eigen::Ref<const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>& y1_T2D_ch_sol, const Eigen::Ref<const Eigen::VectorXd>& x_F1D_x) {
    int Nx_I = static_cast<int>(x_F1D_x.size());
    int Nch_I = static_cast<int>(y0_T2D_ch_sol.rows());
    int Nsol_I = static_cast<int>(y0_T2D_ch_sol.cols());
    ivp_numerov_assert_uniform_grid(x_F1D_x);
    assert(Nch_I > 0 && Nsol_I > 0);
    assert(y1_T2D_ch_sol.rows() == Nch_I && y1_T2D_ch_sol.cols() == Nsol_I);

    Eigen::Tensor<T, 3, Eigen::ColMajor> y_T3D_ch_sol_x(Nch_I, Nsol_I, Nx_I);
    Eigen::TensorMap<Eigen::Tensor<const T, 2, Eigen::ColMajor>> y0map_T2D_ch_sol(y0_T2D_ch_sol.data(), Nch_I, Nsol_I);
    Eigen::TensorMap<Eigen::Tensor<const T, 2, Eigen::ColMajor>> y1map_T2D_ch_sol(y1_T2D_ch_sol.data(), Nch_I, Nsol_I);
    y_T3D_ch_sol_x.chip(0, 2) = y0map_T2D_ch_sol;
    y_T3D_ch_sol_x.chip(1, 2) = y1map_T2D_ch_sol;
    using YMap_T2D_ch_sol = Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>;
    IVP_NumerovState<T> numerov_State(F_Func, x_F1D_x(0), x_F1D_x(1), y0_T2D_ch_sol, y1_T2D_ch_sol);
    for (int x_I = 1; x_I <= Nx_I - 2; ++x_I) {
        YMap_T2D_ch_sol ycurrMap_T2D_ch_sol(&y_T3D_ch_sol_x(0, 0, x_I + 1), Nch_I, Nsol_I);
        ycurrMap_T2D_ch_sol = numerov_State.step(x_F1D_x(x_I + 1));
    }
    return y_T3D_ch_sol_x;
}
