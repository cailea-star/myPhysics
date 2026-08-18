/**
 * @file    derivative.hpp
 * @author  cailea
 * @date    2026-04-09
 * @brief   数值导数计算函数
 */

#pragma once
#include <Eigen/Core>
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
* @brief  Compute the first derivative on a validated uniform mesh.
* @math   f'(x_i) ≈ D_h f_i
* @output First derivative at x_i.
*/
template<typename T>
T derivative1_impl(const Eigen::Vector<T, Eigen::Dynamic>& f_T1D_x, double dx_F, Eigen::Index x_I) {
    Eigen::Index Nx_I = f_T1D_x.size();
    if (x_I == 0) {
        return (f_T1D_x(1) - f_T1D_x(0)) / dx_F;
    } else if (x_I == Nx_I - 1) {
        return (f_T1D_x(Nx_I - 1) - f_T1D_x(Nx_I - 2)) / dx_F;
    } else if (x_I == 1 || x_I == Nx_I - 2) {
        return (f_T1D_x(x_I + 1) - f_T1D_x(x_I - 1)) / (2.0 * dx_F);
    } else {
        return (T(8.0) * (f_T1D_x(x_I + 1) - f_T1D_x(x_I - 1)) - (f_T1D_x(x_I + 2) - f_T1D_x(x_I - 2))) / (12.0 * dx_F);
    }
}

/**
* @brief  Compute the second derivative on a validated uniform mesh.
* @math   f''(x_i) ≈ D_h^2 f_i
* @output Second derivative at x_i.
*/
template<typename T>
T derivative2_impl(const Eigen::Vector<T, Eigen::Dynamic>& f_T1D_x, double dx_F, Eigen::Index x_I) {
    Eigen::Index Nx_I = f_T1D_x.size();
    if (x_I == 0) {
        return (T(2.0) * f_T1D_x(0) - T(5.0) * f_T1D_x(1) + T(4.0) * f_T1D_x(2) - f_T1D_x(3)) / (dx_F * dx_F);
    } else if (x_I == Nx_I - 1) {
        return (T(2.0) * f_T1D_x(Nx_I - 1) - T(5.0) * f_T1D_x(Nx_I - 2) + T(4.0) * f_T1D_x(Nx_I - 3) - f_T1D_x(Nx_I - 4)) / (dx_F * dx_F);
    } else if (x_I == 1 || x_I == Nx_I - 2) {
        return (f_T1D_x(x_I + 1) - T(2.0) * f_T1D_x(x_I) + f_T1D_x(x_I - 1)) / (dx_F * dx_F);
    } else {
        return (T(16.0) * (f_T1D_x(x_I + 1) + f_T1D_x(x_I - 1)) - (f_T1D_x(x_I + 2) + f_T1D_x(x_I - 2)) - T(30.0) * f_T1D_x(x_I)) / (12.0 * dx_F * dx_F);
    }
}
}

/**
 * @brief  Compute the first derivative at one mesh point.
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
 * @brief  Compute the second derivative at one mesh point.
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
 * @brief  Compute the first derivative on a uniform mesh.
 * @math   f'_i ≈ D_h f_i
 * @output First-derivative vector.
 */
template<typename T>
Eigen::Vector<T, Eigen::Dynamic> derivative1(const Eigen::Vector<T, Eigen::Dynamic>& f_T1D_x, const Eigen::VectorXd& x_F1D_x) {
    double dx_F = assert_uniform_mesh(f_T1D_x, x_F1D_x, 2);
    Eigen::Vector<T, Eigen::Dynamic> dfdx_T1D_x(f_T1D_x.size());
    for (Eigen::Index x_I = 0; x_I < f_T1D_x.size(); ++x_I) {
        dfdx_T1D_x(x_I) = derivative1_impl(f_T1D_x, dx_F, x_I);
    }
    return dfdx_T1D_x;
}

/**
 * @brief  Compute the second derivative on a uniform mesh.
 * @math   f''_i ≈ D_h^2 f_i
 * @output Second-derivative vector.
 */
template<typename T>
Eigen::Vector<T, Eigen::Dynamic> derivative2(const Eigen::Vector<T, Eigen::Dynamic>& f_T1D_x, const Eigen::VectorXd& x_F1D_x) {
    double dx_F = assert_uniform_mesh(f_T1D_x, x_F1D_x, 4);
    Eigen::Vector<T, Eigen::Dynamic> d2fdx2_T1D_x(f_T1D_x.size());
    for (Eigen::Index x_I = 0; x_I < f_T1D_x.size(); ++x_I) {
        d2fdx2_T1D_x(x_I) = derivative2_impl(f_T1D_x, dx_F, x_I);
    }
    return d2fdx2_T1D_x;
}


/**
 * @brief  Compute the first derivative of a callable function.
 * @math   f'(x) ≈ [8(f(x+h)-f(x-h))-(f(x+2h)-f(x-2h))]/(12h)
 * @output First derivative at x.
 */
template<typename T>
T derivative1(const std::function<T(double)>& f_Func, double x_F, double dx_F = 1e-3) {
    assert(std::isfinite(dx_F));
    assert(dx_F != 0.0);
    T fm2_T = f_Func(x_F - 2 * dx_F);
    T fm1_T = f_Func(x_F - dx_F);
    T fp1_T = f_Func(x_F + dx_F);
    T fp2_T = f_Func(x_F + 2 * dx_F);
    return (T(8.0) * (fp1_T - fm1_T) - (fp2_T - fm2_T)) / (12.0 * dx_F);
}

/**
 * @brief  Compute the second derivative of a callable function.
 * @math   f''(x) ≈ [16(f(x+h)+f(x-h))-(f(x+2h)+f(x-2h))-30f(x)]/(12h^2)
 * @output Second derivative at x.
 */
template<typename T>
T derivative2(const std::function<T(double)>& f_Func, double x_F, double dx_F = 1e-3) {
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
 * @brief  Compute the first derivative of a real callable function.
 * @math   f'(x) ≈ [8(f(x+h)-f(x-h))-(f(x+2h)-f(x-2h))]/(12h)
 * @output First derivative at x.
 */
inline double derivative1(const std::function<double(double)>& f_Func, double x_F, double dx_F = 1e-3) {
    return derivative1<double>(f_Func, x_F, dx_F);
}

/**
 * @brief  Compute the second derivative of a real callable function.
 * @math   f''(x) ≈ [16(f(x+h)+f(x-h))-(f(x+2h)+f(x-2h))-30f(x)]/(12h^2)
 * @output Second derivative at x.
 */
inline double derivative2(const std::function<double(double)>& f_Func, double x_F, double dx_F = 1e-3) {
    return derivative2<double>(f_Func, x_F, dx_F);
}
