/**
 * @file    ivp_numerov.hpp
 * @author  cailea
 * @date    2026-01-06
 * @brief   Scalar, vector, and matrix Numerov solvers.
 */

#pragma once

#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <functional>
#include <unsupported/Eigen/CXX11/Tensor>

template<typename T> using Real2TMatFunc = std::function<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(double)>;

/**
 * @brief  Check a uniform Numerov grid.
 * @math   h=x_{i+1}-x_i≠0
 * @output Uniform grid spacing.
 */
inline double ivp_numerov_assert_uniform_grid(const Eigen::Ref<const Eigen::VectorXd>& x_F1D_x) {
    int Nx_I = static_cast<int>(x_F1D_x.size());
    assert(Nx_I >= 2);
    double h_F = x_F1D_x(1) - x_F1D_x(0);
    assert(std::isfinite(h_F) && h_F != 0.0);
    for (int x_I = 1; x_I < Nx_I - 1; ++x_I) {assert(std::abs(x_F1D_x(x_I + 1) - x_F1D_x(x_I) - h_F) <= 1.0e-10 * std::abs(h_F));}
    return h_F;
}

// ==================== Scalar Numerov solver ====================

/**
 * @brief  Perform one scalar Numerov step.
 * @math   (1-h²F_{n+1}/12)y_{n+1}=2(1+5h²F_n/12)y_n-(1-h²F_{n-1}/12)y_{n-1}
 * @output Scalar value y_{n+1}.
 */
template<typename T>
T ivp_numerov_step(T Fprev_T, T Fcurr_T, T Fnext_T, T yprev_T, T ycurr_T, double h_F) {
    assert(std::isfinite(h_F) && h_F != 0.0);
    double h2_F = h_F * h_F;
    T lhs_T = T(1.0) - h2_F * Fnext_T / T(12.0);
    T rhs_T = T(2.0) * (T(1.0) + T(5.0) * h2_F * Fcurr_T / T(12.0)) * ycurr_T - (T(1.0) - h2_F * Fprev_T / T(12.0)) * yprev_T;
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
    double h_F = ivp_numerov_assert_uniform_grid(x_F1D_x);
    assert(F_F1D_x.size() == Nx_I);

    Eigen::Matrix<T, Eigen::Dynamic, 1> y_T1D_x(Nx_I);
    y_T1D_x(0) = y0_T;
    y_T1D_x(1) = y1_T;
    T Fprev_T = F_F1D_x(0);
    T Fcurr_T = F_F1D_x(1);
    for (int x_I = 1; x_I <= Nx_I - 2; ++x_I) {
        T Fnext_T = F_F1D_x(x_I + 1);
        y_T1D_x(x_I + 1) = ivp_numerov_step(Fprev_T, Fcurr_T, Fnext_T, y_T1D_x(x_I - 1), y_T1D_x(x_I), h_F);
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

// ==================== Vector Numerov solver ====================

/**
 * @brief  Perform one vector Numerov step.
 * @math   (I-h²F_{n+1}/12)y_{n+1}=(2I+5h²F_n/6)y_n-(I-h²F_{n-1}/12)y_{n-1}
 * @output Channel vector y_{n+1}.
 */
template<typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> ivp_numerov_vec_step(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& Fprev_T2D_ch_ch, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& Fcurr_T2D_ch_ch, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& Fnext_T2D_ch_ch, const Eigen::Ref<const Eigen::Matrix<T, Eigen::Dynamic, 1>>& yprev_T1D_ch, const Eigen::Ref<const Eigen::Matrix<T, Eigen::Dynamic, 1>>& ycurr_T1D_ch, double h_F) {
    int Nch_I = static_cast<int>(ycurr_T1D_ch.size());
    assert(std::isfinite(h_F) && h_F != 0.0);
    assert(Nch_I > 0 && yprev_T1D_ch.size() == Nch_I);
    assert(Fprev_T2D_ch_ch.rows() == Nch_I && Fprev_T2D_ch_ch.cols() == Nch_I);
    assert(Fcurr_T2D_ch_ch.rows() == Nch_I && Fcurr_T2D_ch_ch.cols() == Nch_I);
    assert(Fnext_T2D_ch_ch.rows() == Nch_I && Fnext_T2D_ch_ch.cols() == Nch_I);
    double h2_F = h_F * h_F;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> I_T2D_ch_ch = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Identity(Nch_I, Nch_I);
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> lhs_T2D_ch_ch = I_T2D_ch_ch - h2_F * Fnext_T2D_ch_ch / T(12.0);
    Eigen::Matrix<T, Eigen::Dynamic, 1> rhs_T1D_ch = (T(2.0) * I_T2D_ch_ch + T(5.0) * h2_F * Fcurr_T2D_ch_ch / T(6.0)) * ycurr_T1D_ch - (I_T2D_ch_ch - h2_F * Fprev_T2D_ch_ch / T(12.0)) * yprev_T1D_ch;
    auto lhs_LU = lhs_T2D_ch_ch.lu();
    return lhs_LU.solve(rhs_T1D_ch);
}

/**
 * @brief  Solve a vector Numerov equation.
 * @math   y''(x)=F(x)y(x)
 * @output Channel-by-mesh solution matrix.
 */
template<typename T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> ivp_numerov_vec(const Real2TMatFunc<T>& F_Func, const Eigen::Ref<const Eigen::Matrix<T, Eigen::Dynamic, 1>>& y0_T1D_ch, const Eigen::Ref<const Eigen::Matrix<T, Eigen::Dynamic, 1>>& y1_T1D_ch, const Eigen::Ref<const Eigen::VectorXd>& x_F1D_x) {
    int Nx_I = static_cast<int>(x_F1D_x.size());
    int Nch_I = static_cast<int>(y0_T1D_ch.size());
    double h_F = ivp_numerov_assert_uniform_grid(x_F1D_x);
    assert(Nch_I > 0 && y1_T1D_ch.size() == Nch_I);

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> y_T2D_ch_x(Nch_I, Nx_I);
    y_T2D_ch_x.col(0) = y0_T1D_ch;
    y_T2D_ch_x.col(1) = y1_T1D_ch;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Fprev_T2D_ch_ch = F_Func(x_F1D_x(0));
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Fcurr_T2D_ch_ch = F_Func(x_F1D_x(1));
    assert(Fprev_T2D_ch_ch.rows() == Nch_I && Fprev_T2D_ch_ch.cols() == Nch_I);
    assert(Fcurr_T2D_ch_ch.rows() == Nch_I && Fcurr_T2D_ch_ch.cols() == Nch_I);
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Fnext_T2D_ch_ch(Nch_I, Nch_I);
    Eigen::Matrix<T, Eigen::Dynamic, 1> yprev_T1D_ch(Nch_I);
    Eigen::Matrix<T, Eigen::Dynamic, 1> ycurr_T1D_ch(Nch_I);
    for (int x_I = 1; x_I <= Nx_I - 2; ++x_I) {
        Fnext_T2D_ch_ch.noalias() = F_Func(x_F1D_x(x_I + 1));
        assert(Fnext_T2D_ch_ch.rows() == Nch_I && Fnext_T2D_ch_ch.cols() == Nch_I);
        yprev_T1D_ch.noalias() = y_T2D_ch_x.col(x_I - 1);
        ycurr_T1D_ch.noalias() = y_T2D_ch_x.col(x_I);
        y_T2D_ch_x.col(x_I + 1) = ivp_numerov_vec_step<T>(Fprev_T2D_ch_ch, Fcurr_T2D_ch_ch, Fnext_T2D_ch_ch, yprev_T1D_ch, ycurr_T1D_ch, h_F);
        Fprev_T2D_ch_ch.swap(Fcurr_T2D_ch_ch);
        Fcurr_T2D_ch_ch.swap(Fnext_T2D_ch_ch);
    }
    return y_T2D_ch_x;
}

// ==================== Matrix Numerov solver ====================

/**
 * @brief  Perform one matrix Numerov step.
 * @math   (I-h²F_{n+1}/12)Y_{n+1}=(2I+5h²F_n/6)Y_n-(I-h²F_{n-1}/12)Y_{n-1}
 * @output Channel-by-solution matrix Y_{n+1}.
 */
template<typename T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> ivp_numerov_mat_step(const Eigen::Ref<const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>& Fprev_T2D_ch_ch, const Eigen::Ref<const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>& Fcurr_T2D_ch_ch, const Eigen::Ref<const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>& Fnext_T2D_ch_ch, const Eigen::Ref<const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>& yprev_T2D_ch_sol, const Eigen::Ref<const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>& ycurr_T2D_ch_sol, double h_F) {
    int Nch_I = static_cast<int>(ycurr_T2D_ch_sol.rows());
    int Nsol_I = static_cast<int>(ycurr_T2D_ch_sol.cols());
    assert(std::isfinite(h_F) && h_F != 0.0);
    assert(Nch_I > 0 && Nsol_I > 0);
    assert(yprev_T2D_ch_sol.rows() == Nch_I && yprev_T2D_ch_sol.cols() == Nsol_I);
    assert(Fprev_T2D_ch_ch.rows() == Nch_I && Fprev_T2D_ch_ch.cols() == Nch_I);
    assert(Fcurr_T2D_ch_ch.rows() == Nch_I && Fcurr_T2D_ch_ch.cols() == Nch_I);
    assert(Fnext_T2D_ch_ch.rows() == Nch_I && Fnext_T2D_ch_ch.cols() == Nch_I);
    double h2_F = h_F * h_F;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> I_T2D_ch_ch = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Identity(Nch_I, Nch_I);
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> lhs_T2D_ch_ch = I_T2D_ch_ch - h2_F * Fnext_T2D_ch_ch / T(12.0);
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> rhs_T2D_ch_sol = (T(2.0) * I_T2D_ch_ch + T(5.0) * h2_F * Fcurr_T2D_ch_ch / T(6.0)) * ycurr_T2D_ch_sol - (I_T2D_ch_ch - h2_F * Fprev_T2D_ch_ch / T(12.0)) * yprev_T2D_ch_sol;
    auto lhs_LU = lhs_T2D_ch_ch.lu();
    return lhs_LU.solve(rhs_T2D_ch_sol);
}

/**
 * @brief  Solve matrix Numerov equations from sampled F.
 * @math   Y''(x)=F(x)Y(x)
 * @output Channel-by-solution-by-mesh solution tensor.
 */
template<typename T>
Eigen::Tensor<T, 3, Eigen::ColMajor> ivp_numerov_mat(const Eigen::Tensor<T, 3, Eigen::ColMajor>& F_T3D_ch_ch_x, const Eigen::Ref<const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>& y0_T2D_ch_sol, const Eigen::Ref<const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>& y1_T2D_ch_sol, const Eigen::Ref<const Eigen::VectorXd>& x_F1D_x) {
    int Nx_I = static_cast<int>(x_F1D_x.size());
    int Nch_I = static_cast<int>(y0_T2D_ch_sol.rows());
    int Nsol_I = static_cast<int>(y0_T2D_ch_sol.cols());
    double h_F = ivp_numerov_assert_uniform_grid(x_F1D_x);
    assert(Nch_I > 0 && Nsol_I > 0);
    assert(y1_T2D_ch_sol.rows() == Nch_I && y1_T2D_ch_sol.cols() == Nsol_I);
    assert(F_T3D_ch_ch_x.dimension(0) == Nch_I && F_T3D_ch_ch_x.dimension(1) == Nch_I && F_T3D_ch_ch_x.dimension(2) == Nx_I);

    Eigen::Tensor<T, 3, Eigen::ColMajor> y_T3D_ch_sol_x(Nch_I, Nsol_I, Nx_I);
    Eigen::TensorMap<Eigen::Tensor<const T, 2, Eigen::ColMajor>> y0map_T2D_ch_sol(y0_T2D_ch_sol.data(), Nch_I, Nsol_I);
    Eigen::TensorMap<Eigen::Tensor<const T, 2, Eigen::ColMajor>> y1map_T2D_ch_sol(y1_T2D_ch_sol.data(), Nch_I, Nsol_I);
    y_T3D_ch_sol_x.chip(0, 2) = y0map_T2D_ch_sol;
    y_T3D_ch_sol_x.chip(1, 2) = y1map_T2D_ch_sol;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> yprev_T2D_ch_sol = y0_T2D_ch_sol;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> ycurr_T2D_ch_sol = y1_T2D_ch_sol;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> ynext_T2D_ch_sol(Nch_I, Nsol_I);
    using FMap_T2D_ch_ch = Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>;
    for (int x_I = 1; x_I <= Nx_I - 2; ++x_I) {
        FMap_T2D_ch_ch Fprev_T2D_ch_ch(&F_T3D_ch_ch_x(0, 0, x_I - 1), Nch_I, Nch_I);
        FMap_T2D_ch_ch Fcurr_T2D_ch_ch(&F_T3D_ch_ch_x(0, 0, x_I), Nch_I, Nch_I);
        FMap_T2D_ch_ch Fnext_T2D_ch_ch(&F_T3D_ch_ch_x(0, 0, x_I + 1), Nch_I, Nch_I);
        ynext_T2D_ch_sol.noalias() = ivp_numerov_mat_step<T>(Fprev_T2D_ch_ch, Fcurr_T2D_ch_ch, Fnext_T2D_ch_ch, yprev_T2D_ch_sol, ycurr_T2D_ch_sol, h_F);
        Eigen::TensorMap<Eigen::Tensor<T, 2, Eigen::ColMajor>> ynextMap_T2D_ch_sol(ynext_T2D_ch_sol.data(), Nch_I, Nsol_I);
        y_T3D_ch_sol_x.chip(x_I + 1, 2) = ynextMap_T2D_ch_sol;
        yprev_T2D_ch_sol.swap(ycurr_T2D_ch_sol);
        ycurr_T2D_ch_sol.swap(ynext_T2D_ch_sol);
    }
    return y_T3D_ch_sol_x;
}

/**
 * @brief  Solve matrix Numerov equations from F(x).
 * @math   Y''(x)=F(x)Y(x)
 * @output Channel-by-solution-by-mesh solution tensor.
 */
template<typename T>
Eigen::Tensor<T, 3, Eigen::ColMajor> ivp_numerov_mat(const Real2TMatFunc<T>& F_Func, const Eigen::Ref<const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>& y0_T2D_ch_sol, const Eigen::Ref<const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>& y1_T2D_ch_sol, const Eigen::Ref<const Eigen::VectorXd>& x_F1D_x) {
    int Nx_I = static_cast<int>(x_F1D_x.size());
    int Nch_I = static_cast<int>(y0_T2D_ch_sol.rows());
    assert(y1_T2D_ch_sol.rows() == Nch_I && y1_T2D_ch_sol.cols() == y0_T2D_ch_sol.cols());
    Eigen::Tensor<T, 3, Eigen::ColMajor> F_T3D_ch_ch_x(Nch_I, Nch_I, Nx_I);
    for (int x_I = 0; x_I < Nx_I; ++x_I) {
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> F_T2D_ch_ch = F_Func(x_F1D_x(x_I));
        assert(F_T2D_ch_ch.rows() == Nch_I && F_T2D_ch_ch.cols() == Nch_I);
        for (int chi_I = 0; chi_I < Nch_I; ++chi_I) {
            for (int chj_I = 0; chj_I < Nch_I; ++chj_I) {F_T3D_ch_ch_x(chi_I, chj_I, x_I) = F_T2D_ch_ch(chi_I, chj_I);}
        }
    }
    return ivp_numerov_mat<T>(F_T3D_ch_ch_x, y0_T2D_ch_sol, y1_T2D_ch_sol, x_F1D_x);
}
