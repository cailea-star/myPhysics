/**
 * @file    hfb.hpp
 * @author  cailea
 * @date    2026-09-10
 * @brief   Real finite-temperature HFB matrix calculations.
 */

#pragma once

#include <cassert>
#include <cmath>
#include <functional>

#include <Eigen/Core>
#include <Eigen/Eigenvalues>

class HFB {
public:
    using OneBodyFunc = std::function<void(Eigen::MatrixXd& h0_F2D_sp_sp)>;

    using TwoBodyFunc = std::function<void(
        const Eigen::MatrixXd& rho_F2D_sp_sp,
        const Eigen::MatrixXd& kappa_F2D_sp_sp,
        Eigen::MatrixXd& Gamma_F2D_sp_sp,
        Eigen::MatrixXd& Delta_F2D_sp_sp)>;

    int Nsp_I = 0;

    Eigen::MatrixXd h0_F2D_sp_sp{};
    Eigen::MatrixXd Gamma_F2D_sp_sp{};
    Eigen::MatrixXd Delta_F2D_sp_sp{};

    Eigen::VectorXd Eqp_F1D_qp{};
    Eigen::VectorXd f_F1D_qp{};
    Eigen::MatrixXd U_F2D_sp_qp{};
    Eigen::MatrixXd V_F2D_sp_qp{};

    Eigen::MatrixXd rho_F2D_sp_sp{};
    Eigen::MatrixXd kappa_F2D_sp_sp{};

    Eigen::MatrixXd H_F2D_2sp_2sp{};
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> H_eigensolver{};

    /**
     * @brief  Build the one-body matrix.
     * @math   h₀ᵀ = h₀.
     * @output Overwritten Nsp × Nsp h0.
     */
    OneBodyFunc build_onebody{};
    TwoBodyFunc build_twobody{};

public:
    /**
     * @brief  Initialize workspace, callbacks, and the one-body matrix.
     * @math   h = h₀ + Γ.
     * @output Built h0; remaining matrices zeroed; eigensolver allocated.
     * @note   Requires Nsp > 0 and valid callbacks.
     */
    HFB(int Nsp_I_, const OneBodyFunc& build_onebody_, const TwoBodyFunc& build_twobody_) {
        assert(Nsp_I_ > 0);
        assert(build_onebody_);
        assert(build_twobody_);

        Nsp_I = Nsp_I_;
        build_onebody = build_onebody_;
        build_twobody = build_twobody_;

        // h₀, Γ, Δ ∈ ℝ^{Nsp×Nsp}.
        h0_F2D_sp_sp.resize(Nsp_I, Nsp_I);
        Gamma_F2D_sp_sp.resize(Nsp_I, Nsp_I);
        Delta_F2D_sp_sp.resize(Nsp_I, Nsp_I);
        h0_F2D_sp_sp.setZero();
        Gamma_F2D_sp_sp.setZero();
        Delta_F2D_sp_sp.setZero();

        // E, f ∈ ℝ^{Nsp}; U, V ∈ ℝ^{Nsp×Nsp}.
        Eqp_F1D_qp.resize(Nsp_I);
        f_F1D_qp.resize(Nsp_I);
        U_F2D_sp_qp.resize(Nsp_I, Nsp_I);
        V_F2D_sp_qp.resize(Nsp_I, Nsp_I);
        Eqp_F1D_qp.setZero();
        f_F1D_qp.setZero();
        U_F2D_sp_qp.setZero();
        V_F2D_sp_qp.setZero();

        // ρ, κ ∈ ℝ^{Nsp×Nsp}.
        rho_F2D_sp_sp.resize(Nsp_I, Nsp_I);
        kappa_F2D_sp_sp.resize(Nsp_I, Nsp_I);
        rho_F2D_sp_sp.setZero();
        kappa_F2D_sp_sp.setZero();

        // H ∈ ℝ^{2Nsp×2Nsp}.
        H_F2D_2sp_2sp.resize(2 * Nsp_I, 2 * Nsp_I);
        H_F2D_2sp_2sp.setZero();
        H_eigensolver = Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>(2 * Nsp_I);

        build_onebody(h0_F2D_sp_sp);
    }

    /**
     * @brief  Update fields using the two-body callback.
     * @math   (ρ,κ) → (Γ,Δ).
     * @output Overwritten Gamma and Delta.
     * @note   Γ excludes h₀; Γᵀ = Γ, Δᵀ = -Δ.
     */
    void update_Gamma_Delta();

    /**
     * @brief  Solve thermal HFB using self-adjoint eigendecomposition.
     * @math   H = [h₀+Γ-λI, Δ; -Δ, -h₀-Γ+λI].
     * @math   H[U;V] = [U;V]E; f = 1/(1+exp(E/T)).
     * @math   ρ = V(1-f)Vᵀ + UfUᵀ.
     * @math   κ = V(1-f)Uᵀ + UfVᵀ.
     * @output Updated H, U, V, Eqp, f, rho, and kappa.
     * @output Returned mean particle number N = Tr(ρ).
     * @note   T = k_B T_phys; same energy units as Eqp.
     * @note   Requires T ≥ 0 and no zero-energy modes.
     */
    double update_UV_E_rho_kappa(double lambda_F, double temperature_F);
};


inline void HFB::update_Gamma_Delta() {
    build_twobody(rho_F2D_sp_sp, kappa_F2D_sp_sp, Gamma_F2D_sp_sp, Delta_F2D_sp_sp);
}


inline double HFB::update_UV_E_rho_kappa(double lambda_F, double temperature_F) {
    assert(temperature_F >= 0.0);
    assert(h0_F2D_sp_sp.isApprox(h0_F2D_sp_sp.transpose(), 1.0e-12));
    assert(Gamma_F2D_sp_sp.isApprox(Gamma_F2D_sp_sp.transpose(), 1.0e-12));
    assert(Delta_F2D_sp_sp.isApprox(-Delta_F2D_sp_sp.transpose(), 1.0e-12));

    // H₁₁ = h₀ + Γ - λI.
    H_F2D_2sp_2sp.topLeftCorner(Nsp_I, Nsp_I) = h0_F2D_sp_sp + Gamma_F2D_sp_sp;
    H_F2D_2sp_2sp.topLeftCorner(Nsp_I, Nsp_I).diagonal().array() -= lambda_F;

    // H₁₂ = Δ; H₂₁ = -Δ; H₂₂ = -H₁₁.
    H_F2D_2sp_2sp.topRightCorner(Nsp_I, Nsp_I) = Delta_F2D_sp_sp;
    H_F2D_2sp_2sp.bottomLeftCorner(Nsp_I, Nsp_I) = -Delta_F2D_sp_sp;
    H_F2D_2sp_2sp.bottomRightCorner(Nsp_I, Nsp_I) = -H_F2D_2sp_2sp.topLeftCorner(Nsp_I, Nsp_I);

    // H[U;V] = [U;V]E.
    H_eigensolver.compute(H_F2D_2sp_2sp);
    assert(H_eigensolver.info() == Eigen::Success);
    const Eigen::VectorXd& eigenvalues_F1D_state = H_eigensolver.eigenvalues();
    const Eigen::MatrixXd& eigenvectors_F2D_2sp_state = H_eigensolver.eigenvectors();

    // E₁ ≤ ⋯ ≤ E₂Nsp; retain Nsp positive-energy modes.
    assert(eigenvalues_F1D_state(Nsp_I - 1) < 0.0);
    assert(eigenvalues_F1D_state(Nsp_I) > 0.0);
    Eqp_F1D_qp = eigenvalues_F1D_state.tail(Nsp_I);
    U_F2D_sp_qp = eigenvectors_F2D_2sp_state.topRightCorner(Nsp_I, Nsp_I);
    V_F2D_sp_qp = eigenvectors_F2D_2sp_state.bottomRightCorner(Nsp_I, Nsp_I);

    // T ≤ 10⁻¹² → f = 0; otherwise f = e⁻ᴱᐟᵀ/(1+e⁻ᴱᐟᵀ).
    f_F1D_qp.setZero();
    if (temperature_F > 1.0e-12) {
        for (int qp_I = 0; qp_I < Nsp_I; ++qp_I) {
            const double expMinusEOverT_F = std::exp(-Eqp_F1D_qp(qp_I) / temperature_F);
            f_F1D_qp(qp_I) = expMinusEOverT_F / (1.0 + expMinusEOverT_F);
        }
    }

    // ρ = V(1-f)Vᵀ + UfUᵀ.
    rho_F2D_sp_sp.noalias() = V_F2D_sp_qp * (1.0 - f_F1D_qp.array()).matrix().asDiagonal() * V_F2D_sp_qp.transpose();
    rho_F2D_sp_sp.noalias() += U_F2D_sp_qp * f_F1D_qp.asDiagonal() * U_F2D_sp_qp.transpose();

    // κ = V(1-f)Uᵀ + UfVᵀ.
    kappa_F2D_sp_sp.noalias() = V_F2D_sp_qp * (1.0 - f_F1D_qp.array()).matrix().asDiagonal() * U_F2D_sp_qp.transpose();
    kappa_F2D_sp_sp.noalias() += U_F2D_sp_qp * f_F1D_qp.asDiagonal() * V_F2D_sp_qp.transpose();

    return rho_F2D_sp_sp.trace();
}
