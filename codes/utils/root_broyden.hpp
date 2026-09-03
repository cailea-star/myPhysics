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
public:
    int Nh_I;                            // N_h: maximum retained history length.
    int hnew_I = -1;                     // h_new = -1: history = ∅.
    Eigen::VectorXd xnext_F1D_i{};       // x_{h+1} = x_h + Δx_h.
    Eigen::VectorXd rnext_F1D_i{};       // r_{h+1} = G(x_{h+1}) - x_{h+1}.
    Eigen::MatrixXd dX_F2D_i_h{};        // ΔX = [Δx_1, ..., Δx_{N_h}].
    Eigen::MatrixXd dR_F2D_i_h{};        // ΔR = [Δr_1, ..., Δr_{N_h}].

protected:
    Eigen::VectorXd xcurr_F1D_i{};       // x_h: accepted current iterate.
    Eigen::VectorXd rcurr_F1D_i{};       // r_h = G(x_h) - x_h.
    Eigen::VectorXd gamma_F1D_h{};       // γ = (ΔRᵀΔR)⁺ΔRᵀr_h.
    Eigen::VectorXd Sinv2_F1D_h{};       // S_τ⁻²: truncated inverse squared singular values.
    Eigen::JacobiSVD<Eigen::MatrixXd, Eigen::ComputeThinU | Eigen::ComputeThinV> dR_SVD; // ΔR = USV^T.

public:
    /**
     * @brief  Initialize a modified Broyden fixed-point iteration.
     * @math   r_0 = G(x_0) - x_0, x_1 = x_0 + α r_0
     * @output Initial x_1, r_1, and an empty secant history.
     */
    BroydenIterator(int Nh_I_, const Vec2VecFunc& G_Func, double alpha_F, const Eigen::VectorXd& x0_F1D_i, const Eigen::VectorXd& G0_F1D_i) {
        Nh_I = Nh_I_;
        xnext_F1D_i.resize(x0_F1D_i.size());
        rnext_F1D_i.resize(x0_F1D_i.size());
        dX_F2D_i_h.resize(x0_F1D_i.size(), Nh_I_);
        dR_F2D_i_h.resize(x0_F1D_i.size(), Nh_I_);
        xcurr_F1D_i.resize(x0_F1D_i.size());
        rcurr_F1D_i.resize(x0_F1D_i.size());
        gamma_F1D_h.resize(Nh_I_);
        Sinv2_F1D_h.resize(Nh_I_);

        xnext_F1D_i.setZero();
        rnext_F1D_i.setZero();
        dX_F2D_i_h.setZero();
        dR_F2D_i_h.setZero();
        xcurr_F1D_i.setZero();
        rcurr_F1D_i.setZero();
        gamma_F1D_h.setZero();
        Sinv2_F1D_h.setZero();

        dR_SVD = Eigen::JacobiSVD<Eigen::MatrixXd, Eigen::ComputeThinU | Eigen::ComputeThinV>(x0_F1D_i.size(), Nh_I_);
        assert(Nh_I_ > 0 && Nh_I_ <= x0_F1D_i.size());

        // r_0 = G(x_0) - x_0, x_1 = x_0 + α r_0.
        rcurr_F1D_i.noalias() = G0_F1D_i - x0_F1D_i;
        xnext_F1D_i.noalias() = x0_F1D_i + alpha_F * rcurr_F1D_i;
        G_Func(xnext_F1D_i, rnext_F1D_i);
        rnext_F1D_i.noalias() -= xnext_F1D_i;
    }

    /**
     * @brief  Perform one modified Broyden iteration.
     * @math   γ=(ΔRᵀΔR)⁺ΔRᵀr_h, x_{h+1}=x_h+αr_h-(ΔX+αΔR)γ
     * @output Updates x_{h+1}, r_{h+1}, and the secant history; returns ||r_{h+1}||_infinity.
     * @note   Safeguard threshold is not a convergence criterion.
     * @note   Truncates singular values at τ = 10⁻¹².
     */
    double iterate(const Vec2VecFunc& G_Func, double alpha_F) {
        // (x_h,r_h) ← (x_{h+1},r_{h+1}).
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

        // ΔR=USVᵀ; (ΔRᵀΔR)⁺=VS_τ⁻²Vᵀ.
        dR_SVD.compute(dR_F2D_i_h);
        const auto& S_F1D_h = dR_SVD.singularValues();
        const auto& V_F2D_h_h = dR_SVD.matrixV();
        const double tol_F = 1.0e-12;
        Sinv2_F1D_h = (S_F1D_h.array().abs() > tol_F).cast<double>();
        Sinv2_F1D_h.array() /= S_F1D_h.array().square().max(tol_F * tol_F);

        // γ=(ΔRᵀΔR)⁺ΔRᵀr_h; Δx_h=αr_h-(ΔX+αΔR)γ.
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
