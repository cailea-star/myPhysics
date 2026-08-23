/**
 * @file    derivative.hpp
 * @author  cailea
 * @date    2026-04-09
 * @brief   Finite-difference derivative utilities.
 */

#pragma once
#include <Eigen/Core>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>

namespace {
/**
* @brief  Validate a uniform mesh and return its spacing.
* @math   h = x_{i+1} - x_i != 0
* @output Uniform mesh spacing.
*/
template<typename T>
inline double assert_uniform_mesh(const Eigen::Vector<T, Eigen::Dynamic>& f_T1D_x, const Eigen::VectorXd& x_F1D_x, Eigen::Index Nx_min_I) {
    assert(f_T1D_x.size() == x_F1D_x.size());
    assert(x_F1D_x.size() >= Nx_min_I);

    // {x_i} → h; Δx_i-h → 0.
    double dx_F = x_F1D_x(1) - x_F1D_x(0);
    assert(std::isfinite(dx_F));
    assert(dx_F != 0.0);
    double tol_F = 1.0e-10 * std::abs(dx_F);
    for (Eigen::Index x_I = 1; x_I < x_F1D_x.size() - 1; ++x_I) {
        double step_F = x_F1D_x(x_I + 1) - x_F1D_x(x_I);
        assert(std::isfinite(step_F));
        assert(step_F != 0.0);
        assert(std::abs(step_F - dx_F) <= tol_F);
    }
    return dx_F;
}

/**
* @brief  Compute first derivatives by mixed-order finite differences.
* @math   f'(x_i) ≈ D_h f_i
* @output First derivative at x_i.
*/
template<typename T>
T derivative1_impl(const Eigen::Vector<T, Eigen::Dynamic>& f_T1D_x, double dx_F, Eigen::Index x_I) {
    // d=min(i,N-1-i): d=0→D¹₂, d=1→D¹₃, d≥2→D¹₅.
    const Eigen::Index Nx_I = f_T1D_x.size();
    assert(Nx_I > 1);
    const Eigen::Index leftDistance_I = x_I;
    const Eigen::Index rightDistance_I = Nx_I - 1 - x_I;
    const Eigen::Index edgeDistance_I = std::min(leftDistance_I, rightDistance_I);
    const T fm2_T = f_T1D_x(std::max<Eigen::Index>(0, x_I - 2));
    const T fm1_T = f_T1D_x(std::max<Eigen::Index>(0, x_I - 1));
    const T fp1_T = f_T1D_x(std::min<Eigen::Index>(Nx_I - 1, x_I + 1));
    const T fp2_T = f_T1D_x(std::min<Eigen::Index>(Nx_I - 1, x_I + 2));
    const double isLeft_F = static_cast<double>(leftDistance_I == 0);
    const double isRight_F = static_cast<double>(rightDistance_I == 0);
    const double isNear_F = static_cast<double>(edgeDistance_I == 1);
    const double isInterior_F = static_cast<double>(edgeDistance_I >= 2);
    const T dfdxEdge_T = (isLeft_F * (f_T1D_x(1) - f_T1D_x(0)) + isRight_F * (f_T1D_x(Nx_I - 1) - f_T1D_x(Nx_I - 2))) / dx_F;
    const T dfdxNear_T = isNear_F * (fp1_T - fm1_T) / (2.0 * dx_F);
    const T dfdxInterior_T = isInterior_F * (T(8.0) * (fp1_T - fm1_T) - (fp2_T - fm2_T)) / (12.0 * dx_F);
    return dfdxEdge_T + dfdxNear_T + dfdxInterior_T;
}

/**
* @brief  Compute second derivatives by mixed-order finite differences.
* @math   f''(x_i) ≈ D_h^2 f_i
* @output Second derivative at x_i.
*/
template<typename T>
T derivative2_impl(const Eigen::Vector<T, Eigen::Dynamic>& f_T1D_x, double dx_F, Eigen::Index x_I) {
    // d=min(i,N-1-i): d=0→D²₄, d=1→D²₃, d≥2→D²₅.
    const Eigen::Index Nx_I = f_T1D_x.size();
    assert(Nx_I > 3);
    const Eigen::Index leftDistance_I = x_I;
    const Eigen::Index rightDistance_I = Nx_I - 1 - x_I;
    const Eigen::Index edgeDistance_I = std::min(leftDistance_I, rightDistance_I);
    const T fm2_T = f_T1D_x(std::max<Eigen::Index>(0, x_I - 2));
    const T fm1_T = f_T1D_x(std::max<Eigen::Index>(0, x_I - 1));
    const T fp1_T = f_T1D_x(std::min<Eigen::Index>(Nx_I - 1, x_I + 1));
    const T fp2_T = f_T1D_x(std::min<Eigen::Index>(Nx_I - 1, x_I + 2));
    const double isLeft_F = static_cast<double>(leftDistance_I == 0);
    const double isRight_F = static_cast<double>(rightDistance_I == 0);
    const double isNear_F = static_cast<double>(edgeDistance_I == 1);
    const double isInterior_F = static_cast<double>(edgeDistance_I >= 2);
    const T d2fdx2Edge_T = (isLeft_F * (T(2.0) * f_T1D_x(0) - T(5.0) * f_T1D_x(1) + T(4.0) * f_T1D_x(2) - f_T1D_x(3)) + isRight_F * (T(2.0) * f_T1D_x(Nx_I - 1) - T(5.0) * f_T1D_x(Nx_I - 2) + T(4.0) * f_T1D_x(Nx_I - 3) - f_T1D_x(Nx_I - 4))) / (dx_F * dx_F);
    const T d2fdx2Near_T = isNear_F * (fp1_T - T(2.0) * f_T1D_x(x_I) + fm1_T) / (dx_F * dx_F);
    const T d2fdx2Interior_T = isInterior_F * (T(16.0) * (fp1_T + fm1_T) - (fp2_T + fm2_T) - T(30.0) * f_T1D_x(x_I)) / (12.0 * dx_F * dx_F);
    return d2fdx2Edge_T + d2fdx2Near_T + d2fdx2Interior_T;
}
}

/**
 * @brief  Compute one first derivative by mixed-order finite differences.
 * @math   f'(x_i) ≈ D_h f_i
 * @output First derivative at x_i.
 */
template<typename T>
T derivative1(const Eigen::Vector<T, Eigen::Dynamic>& f_T1D_x, const Eigen::VectorXd& x_F1D_x, Eigen::Index x_I) {
    double dx_F = assert_uniform_mesh(f_T1D_x, x_F1D_x, 2);
    assert(x_I >= 0 && x_I < f_T1D_x.size());
    return derivative1_impl(f_T1D_x, dx_F, x_I);
}

/**
 * @brief  Compute one second derivative by mixed-order finite differences.
 * @math   f''(x_i) ≈ D_h^2 f_i
 * @output Second derivative at x_i.
 */
template<typename T>
T derivative2(const Eigen::Vector<T, Eigen::Dynamic>& f_T1D_x, const Eigen::VectorXd& x_F1D_x, Eigen::Index x_I) {
    double dx_F = assert_uniform_mesh(f_T1D_x, x_F1D_x, 4);
    assert(x_I >= 0 && x_I < f_T1D_x.size());
    return derivative2_impl(f_T1D_x, dx_F, x_I);
}

/**
 * @brief  Compute first-derivative vector by mixed-order finite differences.
 * @math   f'_i ≈ D_h f_i
 * @output First-derivative vector.
 */
template<typename T>
Eigen::Vector<T, Eigen::Dynamic> derivative1(const Eigen::Vector<T, Eigen::Dynamic>& f_T1D_x, const Eigen::VectorXd& x_F1D_x) {
    // {f_i} → {f'_i}.
    double dx_F = assert_uniform_mesh(f_T1D_x, x_F1D_x, 2);
    Eigen::Vector<T, Eigen::Dynamic> dfdx_T1D_x(f_T1D_x.size());
    for (Eigen::Index x_I = 0; x_I < f_T1D_x.size(); ++x_I) {
        dfdx_T1D_x(x_I) = derivative1_impl(f_T1D_x, dx_F, x_I);
    }
    return dfdx_T1D_x;
}

/**
 * @brief  Compute second-derivative vector by mixed-order finite differences.
 * @math   f''_i ≈ D_h^2 f_i
 * @output Second-derivative vector.
 */
template<typename T>
Eigen::Vector<T, Eigen::Dynamic> derivative2(const Eigen::Vector<T, Eigen::Dynamic>& f_T1D_x, const Eigen::VectorXd& x_F1D_x) {
    // {f_i} → {f''_i}.
    double dx_F = assert_uniform_mesh(f_T1D_x, x_F1D_x, 4);
    Eigen::Vector<T, Eigen::Dynamic> d2fdx2_T1D_x(f_T1D_x.size());
    for (Eigen::Index x_I = 0; x_I < f_T1D_x.size(); ++x_I) {
        d2fdx2_T1D_x(x_I) = derivative2_impl(f_T1D_x, dx_F, x_I);
    }
    return d2fdx2_T1D_x;
}


/**
 * @brief  Compute callable derivative by fourth-order central differences.
 * @math   f'(x) ≈ [8(f(x+h)-f(x-h))-(f(x+2h)-f(x-2h))]/(12h)
 * @output First derivative at x.
 */
template<typename T>
T derivative1(const std::function<T(double)>& f_Func, double x_F, double dx_F = 1e-3) {
    // f(x±h),f(x±2h) → f'(x).
    assert(std::isfinite(dx_F));
    assert(dx_F != 0.0);
    T fm2_T = f_Func(x_F - 2 * dx_F);
    T fm1_T = f_Func(x_F - dx_F);
    T fp1_T = f_Func(x_F + dx_F);
    T fp2_T = f_Func(x_F + 2 * dx_F);
    return (T(8.0) * (fp1_T - fm1_T) - (fp2_T - fm2_T)) / (12.0 * dx_F);
}

/**
 * @brief  Compute callable second derivative by fourth-order central differences.
 * @math   f''(x) ≈ [16(f(x+h)+f(x-h))-(f(x+2h)+f(x-2h))-30f(x)]/(12h^2)
 * @output Second derivative at x.
 */
template<typename T>
T derivative2(const std::function<T(double)>& f_Func, double x_F, double dx_F = 1e-3) {
    // f(x),f(x±h),f(x±2h) → f''(x).
    assert(std::isfinite(dx_F));
    assert(dx_F != 0.0);
    T fm2_T = f_Func(x_F - 2 * dx_F);
    T fm1_T = f_Func(x_F - dx_F);
    T f0_T = f_Func(x_F);
    T fp1_T = f_Func(x_F + dx_F);
    T fp2_T = f_Func(x_F + 2 * dx_F);
    return (T(16.0) * (fp1_T + fm1_T) - (fp2_T + fm2_T) - T(30.0) * f0_T) / (12.0 * dx_F * dx_F);
}

/**
 * @brief  Compute real callable derivative by fourth-order central differences.
 * @math   f'(x) ≈ [8(f(x+h)-f(x-h))-(f(x+2h)-f(x-2h))]/(12h)
 * @output First derivative at x.
 */
inline double derivative1(const std::function<double(double)>& f_Func, double x_F, double dx_F = 1e-3) {
    return derivative1<double>(f_Func, x_F, dx_F);
}

/**
 * @brief  Compute real callable second derivative by fourth-order central differences.
 * @math   f''(x) ≈ [16(f(x+h)+f(x-h))-(f(x+2h)+f(x-2h))-30f(x)]/(12h^2)
 * @output Second derivative at x.
 */
inline double derivative2(const std::function<double(double)>& f_Func, double x_F, double dx_F = 1e-3) {
    return derivative2<double>(f_Func, x_F, dx_F);
}
