/**
 * @file    root_broyden.hpp
 * @author  cailea
 * @date    2026-05-14
 * @brief   Modified Broyden iteration for nonlinear fixed-point equations.
 */

#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <cassert>
#include <cmath>
#include <functional>

using Vec2VecFunc = std::function<void(const Eigen::VectorXd&, Eigen::VectorXd&)>;

class BroydenIterator {
protected:
    Eigen::VectorXd xcurr_F1D_i;         // x_h: accepted current iterate.
    Eigen::VectorXd rcurr_F1D_i;         // r_h = G(x_h) - x_h.
    Eigen::VectorXd gamma_F1D_h;         // γ = (ΔR^TΔR)^(-1)ΔR^Tr_h.
    Eigen::VectorXd Sinv2_F1D_h;         // S_h^(-2): truncated inverse squared singular values.
    Eigen::JacobiSVD<Eigen::MatrixXd, Eigen::ComputeThinU | Eigen::ComputeThinV> dR_SVD; // ΔR = USV^T.

public:
    int Nh_I;                            // N_h: maximum retained history length.
    int hnew_I;                          // h_new = -1: history = ∅.
    Eigen::VectorXd xnext_F1D_i;         // x_{h+1} = x_h + Δx_h.
    Eigen::VectorXd rnext_F1D_i;         // r_{h+1} = G(x_{h+1}) - x_{h+1}.
    Eigen::MatrixXd dX_F2D_i_h;          // ΔX = [Δx_1, ..., Δx_{N_h}].
    Eigen::MatrixXd dR_F2D_i_h;          // ΔR = [Δr_1, ..., Δr_{N_h}].

    /**
     * @brief  Initialize a modified Broyden fixed-point iteration.
     * @math   r_0 = G(x_0) - x_0, x_1 = x_0 + α r_0
     * @output Initial x_1, r_1, and an empty secant history.
     */
    BroydenIterator(int Nh_I, const Vec2VecFunc& G_Func, double alpha_F, const Eigen::VectorXd& x0_F1D_i, const Eigen::VectorXd& G0_F1D_i)
    : xcurr_F1D_i(Eigen::VectorXd::Zero(x0_F1D_i.size())), rcurr_F1D_i(Eigen::VectorXd::Zero(x0_F1D_i.size())), gamma_F1D_h(Eigen::VectorXd::Zero(Nh_I)), Sinv2_F1D_h(Eigen::VectorXd::Zero(Nh_I)), dR_SVD(x0_F1D_i.size(), Nh_I), Nh_I(Nh_I), hnew_I(-1), xnext_F1D_i(Eigen::VectorXd::Zero(x0_F1D_i.size())), rnext_F1D_i(Eigen::VectorXd::Zero(x0_F1D_i.size())), dX_F2D_i_h(Eigen::MatrixXd::Zero(x0_F1D_i.size(), Nh_I)), dR_F2D_i_h(Eigen::MatrixXd::Zero(x0_F1D_i.size(), Nh_I)) {
        assert(Nh_I > 0 && Nh_I <= x0_F1D_i.size());

        // r_0 = G(x_0) - x_0, x_1 = x_0 + α r_0.
        rcurr_F1D_i.noalias() = G0_F1D_i - x0_F1D_i;
        xnext_F1D_i.noalias() = x0_F1D_i + alpha_F * rcurr_F1D_i;
        G_Func(xnext_F1D_i, rnext_F1D_i);
        rnext_F1D_i.noalias() -= xnext_F1D_i;
    }

    /**
     * @brief  Perform one modified Broyden iteration.
     * @math   γ = (ΔR^T ΔR)^(-1) ΔR^T r_h, x_{h+1} = x_h + α r_h - (ΔX + α ΔR) γ
     * @output Updates x_{h+1}, r_{h+1}, and the secant history; returns ||r_{h+1}||_infinity.
     * @note   The fallback threshold is an implementation safeguard, not a convergence criterion.
     */
    double iterate(const Vec2VecFunc& G_Func, double alpha_F) {
        // (x_h, r_h) ← (x_{h+1}, r_{h+1}).
        xcurr_F1D_i.swap(xnext_F1D_i);
        rcurr_F1D_i.swap(rnext_F1D_i);

        // history = ∅ → {(Δx_h, Δr_h)}.
        if (hnew_I < 0) {
            xnext_F1D_i.noalias() = xcurr_F1D_i + alpha_F * rcurr_F1D_i;
            G_Func(xnext_F1D_i, rnext_F1D_i);
            rnext_F1D_i.noalias() -= xnext_F1D_i;
            hnew_I = 0;
            dX_F2D_i_h.col(hnew_I) = xnext_F1D_i - xcurr_F1D_i;
            dR_F2D_i_h.col(hnew_I) = rnext_F1D_i - rcurr_F1D_i;
            return rnext_F1D_i.array().abs().maxCoeff();
        }

        // ΔR = USV^T, (ΔR^T ΔR)^(-1) = VS^(-2)V^T.
        dR_SVD.compute(dR_F2D_i_h);
        const auto& S_F1D_h = dR_SVD.singularValues();
        const auto& V_F2D_h_h = dR_SVD.matrixV();
        const double tol_F = 1.0e-12;
        for (int h_I = 0; h_I < Nh_I; ++h_I) {Sinv2_F1D_h(h_I) = std::abs(S_F1D_h(h_I)) > tol_F ? 1.0 / (S_F1D_h(h_I) * S_F1D_h(h_I)) : 0.0;}

        // γ = (ΔR^T ΔR)^(-1) ΔR^T r_h, Δx_h = α r_h - (ΔX + α ΔR) γ.
        gamma_F1D_h.noalias() = V_F2D_h_h * Sinv2_F1D_h.asDiagonal() * V_F2D_h_h.transpose() * (dR_F2D_i_h.transpose() * rcurr_F1D_i);
        xnext_F1D_i.noalias() = dX_F2D_i_h * gamma_F1D_h;
        xnext_F1D_i.noalias() += alpha_F * dR_F2D_i_h * gamma_F1D_h;
        xnext_F1D_i = alpha_F * rcurr_F1D_i - xnext_F1D_i;
        const double condition_F = rcurr_F1D_i.dot(xnext_F1D_i);
        if (condition_F <= -1.0) {xnext_F1D_i.noalias() = 0.5 * alpha_F * rcurr_F1D_i;} // Δx_h = α r_h / 2.
        xnext_F1D_i.noalias() += xcurr_F1D_i; // x_{h+1} = x_h + Δx_h.

        // r_{h+1} = G(x_{h+1}) - x_{h+1}.
        G_Func(xnext_F1D_i, rnext_F1D_i);
        rnext_F1D_i.noalias() -= xnext_F1D_i;

        // history ← history ∪ {(Δx_h, Δr_h)}.
        hnew_I = (hnew_I + 1) % Nh_I;
        dX_F2D_i_h.col(hnew_I) = xnext_F1D_i - xcurr_F1D_i;
        dR_F2D_i_h.col(hnew_I) = rnext_F1D_i - rcurr_F1D_i;
        return rnext_F1D_i.array().abs().maxCoeff();
    }
};
