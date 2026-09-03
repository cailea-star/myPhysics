/**
 * @file    bvp_shooting.hpp
 * @author  cailea
 * @date    2026-01-29
 * @brief   Bound-state shooting methods for eigenenergies and wave functions.
 */

#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <complex>
#include <functional>
#include <limits>
#include <utility>
#include <Eigen/Dense>
#include "derivative.hpp"
#include "integration.hpp"
#include "ivp_numerov.hpp"
#include "ivp_rungekutta.hpp"
#include "minimize.hpp"
#include "root.hpp"

using doubleC = std::complex<double>;
using ShootingBoundary = std::pair<doubleC, doubleC>;
using ShootingPotentialFunc = std::function<double(double)>;
using ShootingBoundaryFunc = std::function<ShootingBoundary(double, double)>;

class ShootingProblem {
public:
    ShootingPotentialFunc V_Func;
    ShootingBoundaryFunc Bin_Func;
    ShootingBoundaryFunc Bout_Func;
    Eigen::VectorXd x_F1D_x{};
    Eigen::VectorXd V_F1D_x{};
    int node_I;
    double hmass_F;
    int xmatch_I;
    double tol_F;

    /**
     * @brief  Construct a bound-state shooting problem on a prescribed grid.
     * @math   h_μ y''(x) = [V(x) - E] y(x)
     * @output Initialized potential mesh and shooting configuration.
     * @note   Boundary functions have the form B(x,E) = (y,y').
     */
    ShootingProblem(const ShootingPotentialFunc& V_Func_, const ShootingBoundaryFunc& Bin_Func_, const ShootingBoundaryFunc& Bout_Func_, const Eigen::Ref<const Eigen::VectorXd>& x_F1D_x_, int node_I_, double hmass_F_, int xmatch_I_, double tol_F_ = 1.0e-8) {
        V_Func = V_Func_;
        Bin_Func = Bin_Func_;
        Bout_Func = Bout_Func_;
        x_F1D_x = x_F1D_x_;
        V_F1D_x.resize(x_F1D_x_.size());
        node_I = node_I_;
        hmass_F = hmass_F_;
        xmatch_I = xmatch_I_;
        tol_F = tol_F_;
        // (V,{x_i}) → {V_i}.
        assert(x_F1D_x.size() >= 3);
        assert(node_I >= 0);
        assert(std::isfinite(hmass_F) && hmass_F > 0.0);
        assert(xmatch_I > 0 && xmatch_I < x_F1D_x.size() - 1);
        assert(std::isfinite(tol_F) && tol_F > 0.0);
        for (int x_I = 0; x_I < x_F1D_x.size(); ++x_I) {V_F1D_x(x_I) = V_Func(x_F1D_x(x_I));}
    }
};

struct ShootingSolution {
    int node_I = 0;
    double E_F = 0.0;
    doubleC R_C = 0.0; // R(E) = y'_in / y_in - y'_out / y_out
    Eigen::VectorXcd y_C1D_x{};
};

using ShootingMatchFunc = std::function<ShootingSolution(double)>;

namespace shooting_detail {

/**
 * @brief  Bracket the requested node count by bisection.
 * @math   N(E_n) = n
 * @output Energy interval containing the requested state.
 */
inline std::pair<double, double> search_node_interval(const ShootingProblem& problem, const ShootingMatchFunc& match_Func, double Elo_F, double Eup_F) {
    // (V_min,V_max,L,n) → [E_lo,E_up].
    double Vmin_F = problem.V_F1D_x.minCoeff();
    double Vmax_F = problem.V_F1D_x.maxCoeff();
    double L_F = problem.x_F1D_x(problem.x_F1D_x.size() - 1) - problem.x_F1D_x(0);
    double k_F = std::acos(-1.0) * (problem.node_I + 1.0) / L_F;
    double margin_F = std::max(10.0, Vmax_F - Vmin_F);
    double Ekin_F = problem.hmass_F * k_F * k_F;
    Elo_F = std::isfinite(Elo_F) ? Elo_F : Vmin_F - margin_F;
    Eup_F = std::isfinite(Eup_F) ? Eup_F : Vmax_F + margin_F + 4.0 * Ekin_F;

    // N(E) → [E_n^lo,E_n^up].
    Real2RealFunc Nlo_Func = [&](double E_F) {
        return match_Func(E_F).node_I - problem.node_I + 0.5;
    };
    Real2RealFunc Nup_Func = [&](double E_F) {
        return match_Func(E_F).node_I - problem.node_I - 0.5;
    };
    double Enlo_F = problem.node_I == 0 ? Elo_F : root_bisection(Nlo_Func, Elo_F, Eup_F, problem.tol_F);
    double Enup_F = root_bisection(Nup_Func, Enlo_F, Eup_F, problem.tol_F);
    return {Enlo_F, Enup_F};
}

/**
 * @brief  Count nodes from real-part sign changes.
 * @math   N = # sign changes in {Re y_i : Re y_i ≠ 0}
 * @output Number of nodes.
 * @note   Restricted to real Sturm--Liouville bound states.
 */
inline int count_nodes(const Eigen::Ref<const Eigen::VectorXcd>& y_C1D_x) {
    // {Re y_i ≠ 0} → N_sign.
    assert(y_C1D_x.size() >= 2);
    int node_I = 0;
    double yprev_F = y_C1D_x(0).real();
    for (int x_I = 1; x_I < y_C1D_x.size(); ++x_I) {
        double ycurr_F = y_C1D_x(x_I).real();
        int isSignChange_I = static_cast<int>(yprev_F * ycurr_F < 0.0);
        double isNonzero_F = static_cast<double>(ycurr_F != 0.0);
        node_I += isSignChange_I;
        yprev_F = std::array<double, 2>{yprev_F, ycurr_F}[static_cast<std::size_t>(isNonzero_F)];
    }
    return node_I;
}

/**
 * @brief  Evaluate the logarithmic-derivative matching error.
 * @math   R(E) = y'_in/y_in - y'_out/y_out
 * @output Complex matching error.
 */
inline doubleC match_error(const ShootingBoundary& yin_C1D_ydy, const ShootingBoundary& yout_C1D_ydy) {
    // (y_in,y'_in,y_out,y'_out) → W/(y_in y_out).
    doubleC W_C = yin_C1D_ydy.second * yout_C1D_ydy.first - yin_C1D_ydy.first * yout_C1D_ydy.second;
    doubleC yinout_C = yin_C1D_ydy.first * yout_C1D_ydy.first;
    if (std::abs(yinout_C) <= 1.0e-300) {return doubleC(std::numeric_limits<double>::max(), 0.0);}
    return W_C / yinout_C;
}

/**
 * @brief  Normalize a wave function by trapezoidal quadrature.
 * @math   Norm² = ∫|y(x)|²dx, y → y/Norm
 * @output Normalized wave function.
 */
inline void normalize(Eigen::VectorXcd& y_C1D_x, const Eigen::Ref<const Eigen::VectorXd>& x_F1D_x) {
    // {y_i} → Norm → {y_i/Norm}.
    Eigen::VectorXd density_F1D_x = y_C1D_x.cwiseAbs2();
    double norm_F = std::sqrt(integrate_trapezoidal(density_F1D_x, x_F1D_x));
    assert(std::isfinite(norm_F) && norm_F > 0.0);
    y_C1D_x /= norm_F;
}

}

/**
 * @brief  Match Numerov solutions inward and outward.
 * @math   h_μ y''(x) = [V(x) - E]y(x), R(E) = y'_in/y_in - y'_out/y_out
 * @output Matched and normalized solution at E.
 * @note   Requires a uniform grid.
 */
inline ShootingSolution shooting_match_numerov(const ShootingProblem& problem, double E_F) {
    // (B_in(E),B_out(E)) → (y₀,y₁,y_{N-1},y_{N-2}).
    int Nx_I = static_cast<int>(problem.x_F1D_x.size());
    ShootingBoundary Bin0_C1D_ydy = problem.Bin_Func(problem.x_F1D_x(0), E_F);
    ShootingBoundary Bin1_C1D_ydy = problem.Bin_Func(problem.x_F1D_x(1), E_F);
    ShootingBoundary Bout1_C1D_ydy = problem.Bout_Func(problem.x_F1D_x(Nx_I - 1), E_F);
    ShootingBoundary Bout2_C1D_ydy = problem.Bout_Func(problem.x_F1D_x(Nx_I - 2), E_F);

    // x_min → {y_in(x_i)}.
    Eigen::VectorXcd F_C1D_x(Nx_I);
    for (int x_I = 0; x_I < Nx_I; ++x_I) {F_C1D_x(x_I) = doubleC((problem.V_F1D_x(x_I) - E_F) / problem.hmass_F, 0.0);}
    Eigen::VectorXcd yin_C1D_x = ivp_numerov<doubleC>(F_C1D_x, Bin0_C1D_ydy.first, Bin1_C1D_ydy.first, problem.x_F1D_x);

    // x_max → {y_out(x_i)}.
    Eigen::VectorXd xrev_F1D_x = problem.x_F1D_x.reverse().eval();
    Eigen::VectorXcd Frev_C1D_x = F_C1D_x.reverse().eval();
    Eigen::VectorXcd youtRev_C1D_x = ivp_numerov<doubleC>(Frev_C1D_x, Bout1_C1D_ydy.first, Bout2_C1D_ydy.first, xrev_F1D_x);
    Eigen::VectorXcd yout_C1D_x = youtRev_C1D_x.reverse().eval();

    // R(E; x_match).
    doubleC dyin_C = derivative1<doubleC>(yin_C1D_x, problem.x_F1D_x, problem.xmatch_I);
    doubleC dyout_C = derivative1<doubleC>(yout_C1D_x, problem.x_F1D_x, problem.xmatch_I);
    ShootingBoundary yinmatch_C1D_ydy{yin_C1D_x(problem.xmatch_I), dyin_C};
    ShootingBoundary youtmatch_C1D_ydy{yout_C1D_x(problem.xmatch_I), dyout_C};
    doubleC R_C = shooting_detail::match_error(yinmatch_C1D_ydy, youtmatch_C1D_ydy);
    if (std::abs(youtmatch_C1D_ydy.first) > 1.0e-300) {yout_C1D_x *= yinmatch_C1D_ydy.first / youtmatch_C1D_ydy.first;}

    // y = y_in ⊕ y_out → y / Norm.
    Eigen::VectorXcd y_C1D_x = yin_C1D_x;
    int Ntail_I = static_cast<int>(y_C1D_x.size()) - problem.xmatch_I - 1;
    y_C1D_x.tail(Ntail_I) = yout_C1D_x.tail(Ntail_I);
    ShootingSolution solution{shooting_detail::count_nodes(y_C1D_x), E_F, R_C, y_C1D_x};
    shooting_detail::normalize(solution.y_C1D_x, problem.x_F1D_x);
    return solution;
}

/**
 * @brief  Match RK4 solutions inward and outward.
 * @math   y' = z, h_μ z' = [V(x) - E]y
 * @output Matched and normalized solution at E.
 */
inline ShootingSolution shooting_match_rk4(const ShootingProblem& problem, double E_F) {
    // y'' → (y, y')'.
    RealTVec2TVecFunc<doubleC> ode_Func = [&](double x_F, const Eigen::Ref<const Eigen::VectorXcd>& ydy_C1D_ydy, Eigen::Ref<Eigen::VectorXcd> dydx_C1D_ydy) {
        dydx_C1D_ydy(0) = ydy_C1D_ydy(1);
        dydx_C1D_ydy(1) = (problem.V_Func(x_F) - E_F) * ydy_C1D_ydy(0) / problem.hmass_F;
    };

    // (B_in(E),B_out(E)) → (y_in,y′_in,y_out,y′_out).
    int Nx_I = static_cast<int>(problem.x_F1D_x.size());
    ShootingBoundary Bin_C1D_ydy = problem.Bin_Func(problem.x_F1D_x(0), E_F);
    ShootingBoundary Bout_C1D_ydy = problem.Bout_Func(problem.x_F1D_x(Nx_I - 1), E_F);
    Eigen::Vector2cd yin0_C1D_ydy;
    Eigen::Vector2cd yout0_C1D_ydy;
    yin0_C1D_ydy << Bin_C1D_ydy.first, Bin_C1D_ydy.second;
    yout0_C1D_ydy << Bout_C1D_ydy.first, Bout_C1D_ydy.second;

    // x_min → {(y_in(x_i),y′_in(x_i))}.
    Eigen::MatrixXcd yin_C2D_ydy_x = ivp_rk4_vec<doubleC>(ode_Func, problem.x_F1D_x, yin0_C1D_ydy);
    Eigen::VectorXcd yin_C1D_x = yin_C2D_ydy_x.row(0).transpose();

    // x_max → {(y_out(x_i),y′_out(x_i))}.
    Eigen::VectorXd xrev_F1D_x = problem.x_F1D_x.reverse().eval();
    Eigen::MatrixXcd youtRev_C2D_ydy_x = ivp_rk4_vec<doubleC>(ode_Func, xrev_F1D_x, yout0_C1D_ydy);
    Eigen::VectorXcd yout_C1D_x = youtRev_C2D_ydy_x.row(0).reverse().transpose().eval();
    Eigen::VectorXcd dyout_C1D_x = youtRev_C2D_ydy_x.row(1).reverse().transpose().eval();

    // R(E; x_match).
    doubleC dyin_C = yin_C2D_ydy_x(1, problem.xmatch_I);
    doubleC dyout_C = dyout_C1D_x(problem.xmatch_I);
    ShootingBoundary yinmatch_C1D_ydy{yin_C1D_x(problem.xmatch_I), dyin_C};
    ShootingBoundary youtmatch_C1D_ydy{yout_C1D_x(problem.xmatch_I), dyout_C};
    doubleC R_C = shooting_detail::match_error(yinmatch_C1D_ydy, youtmatch_C1D_ydy);
    if (std::abs(youtmatch_C1D_ydy.first) > 1.0e-300) {yout_C1D_x *= yinmatch_C1D_ydy.first / youtmatch_C1D_ydy.first;}

    // y = y_in ⊕ y_out → y / Norm.
    Eigen::VectorXcd y_C1D_x = yin_C1D_x;
    int Ntail_I = static_cast<int>(y_C1D_x.size()) - problem.xmatch_I - 1;
    y_C1D_x.tail(Ntail_I) = yout_C1D_x.tail(Ntail_I);
    ShootingSolution solution{shooting_detail::count_nodes(y_C1D_x), E_F, R_C, y_C1D_x};
    shooting_detail::normalize(solution.y_C1D_x, problem.x_F1D_x);
    return solution;
}

/**
 * @brief  Solve Numerov shooting by Brent minimization.
 * @math   E_n = arg min_E |R(E)|, N(E) = n
 * @output Eigenenergy and normalized matched solution.
 */
inline ShootingSolution shooting_numerov(const ShootingProblem& problem, double Elo_F = NAN, double Eup_F = NAN) {
    ShootingMatchFunc match_Func = [&](double E_F) {
        return shooting_match_numerov(problem, E_F);
    };

    // [E_lo,E_up] → [E_n^lo,E_n^up]; ¬finite(E_bound) → E=0.
    std::pair<double, double> Ebound_F1D_bound = shooting_detail::search_node_interval(problem, match_Func, Elo_F, Eup_F);
    if (!std::isfinite(Ebound_F1D_bound.first) || !std::isfinite(Ebound_F1D_bound.second)) {return match_Func(0.0);}

    // E_n = argmin_E |R(E)|.
    Real2RealFunc Rabs_Func = [&](double E_F) {
        return std::abs(match_Func(E_F).R_C);
    };
    double E_F = minimize_brent(Rabs_Func, Ebound_F1D_bound.first, Ebound_F1D_bound.second, problem.tol_F);
    return match_Func(E_F);
}

/**
 * @brief  Solve RK4 shooting by Brent minimization.
 * @math   E_n = arg min_E |R(E)|, N(E) = n
 * @output Eigenenergy and normalized matched solution.
 */
inline ShootingSolution shooting_rk4(const ShootingProblem& problem, double Elo_F = NAN, double Eup_F = NAN) {
    ShootingMatchFunc match_Func = [&](double E_F) {
        return shooting_match_rk4(problem, E_F);
    };

    // [E_lo,E_up] → [E_n^lo,E_n^up]; ¬finite(E_bound) → E=0.
    std::pair<double, double> Ebound_F1D_bound = shooting_detail::search_node_interval(problem, match_Func, Elo_F, Eup_F);
    if (!std::isfinite(Ebound_F1D_bound.first) || !std::isfinite(Ebound_F1D_bound.second)) {return match_Func(0.0);}

    // E_n = argmin_E |R(E)|.
    Real2RealFunc Rabs_Func = [&](double E_F) {
        return std::abs(match_Func(E_F).R_C);
    };
    double E_F = minimize_brent(Rabs_Func, Ebound_F1D_bound.first, Ebound_F1D_bound.second, problem.tol_F);
    return match_Func(E_F);
}
