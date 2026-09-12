/**
 * @file    hfb.hpp
 * @author  cailea
 * @date    2026-09-10
 * @brief   Real finite-temperature HFB matrix calculations.
 */

#pragma once

#include <cassert>
#include <cmath>

#include <Eigen/Core>
#include <Eigen/Eigenvalues>

struct HFBSolution {
    Eigen::VectorXd Eqp_F1D_qp{};
    Eigen::VectorXd f_F1D_qp{};
    Eigen::MatrixXd U_F2D_sp_qp{};
    Eigen::MatrixXd V_F2D_sp_qp{};
    Eigen::MatrixXd rho_F2D_sp_sp{};
    Eigen::MatrixXd kappa_F2D_sp_sp{};
};

struct HFBField {
    Eigen::MatrixXd h0_F2D_sp_sp{};
    Eigen::MatrixXd Gamma_F2D_sp_sp{};
    Eigen::MatrixXd Delta_F2D_sp_sp{};
};

class HFB {
public:
    double lambda_F = -7.0; // Fermi energy [MeV].
    double lambda2_F = 0.0; // Lipkin-Nogami λ₂ [MeV].
    double temperature_F = 0.0;
    double EspCut_F = 60.0; // Equivalent single-particle energy cutoff [MeV].

    int Nsp_I = 0;

    HFBField hfb_field{};
    HFBSolution hfb_solution{};

    Eigen::MatrixXd H_F2D_2sp_2sp{};
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> H_eigensolver{};

public:
    /**
     * @brief  Initialize the single-species HFB workspace.
     * @math   h = h₀ + Γ.
     * @output Zeroed matrices and allocated eigensolver workspace.
     * @note   Requires Nsp > 0.
     */
    HFB(int Nsp_I_) {
        assert(Nsp_I_ > 0);

        Nsp_I = Nsp_I_;

        // h₀, Γ, Δ ∈ ℝ^{Nsp×Nsp}.
        hfb_field.h0_F2D_sp_sp.resize(Nsp_I, Nsp_I);
        hfb_field.h0_F2D_sp_sp.setZero();
        hfb_field.Gamma_F2D_sp_sp.resize(Nsp_I, Nsp_I);
        hfb_field.Gamma_F2D_sp_sp.setZero();
        hfb_field.Delta_F2D_sp_sp.resize(Nsp_I, Nsp_I);
        hfb_field.Delta_F2D_sp_sp.setZero();

        // E, f ∈ ℝ^{Nsp}; U, V ∈ ℝ^{Nsp×Nsp}.
        hfb_solution.Eqp_F1D_qp.resize(Nsp_I);
        hfb_solution.f_F1D_qp.resize(Nsp_I);
        hfb_solution.U_F2D_sp_qp.resize(Nsp_I, Nsp_I);
        hfb_solution.V_F2D_sp_qp.resize(Nsp_I, Nsp_I);
        hfb_solution.Eqp_F1D_qp.setZero();
        hfb_solution.f_F1D_qp.setZero();
        hfb_solution.U_F2D_sp_qp.setZero();
        hfb_solution.V_F2D_sp_qp.setZero();

        // ρ, κ ∈ ℝ^{Nsp×Nsp}.
        hfb_solution.rho_F2D_sp_sp.resize(Nsp_I, Nsp_I);
        hfb_solution.kappa_F2D_sp_sp.resize(Nsp_I, Nsp_I);
        hfb_solution.rho_F2D_sp_sp.setZero();
        hfb_solution.kappa_F2D_sp_sp.setZero();

        // H ∈ ℝ^{2Nsp×2Nsp}.
        H_F2D_2sp_2sp.resize(2 * Nsp_I, 2 * Nsp_I);
        H_F2D_2sp_2sp.setZero();
        H_eigensolver = Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>(2 * Nsp_I);

    }

    /**
     * @brief  Solve thermal HFB using self-adjoint eigendecomposition.
     * @math   H = [h_LN-λI, Δ; -Δ, -h_LN+λI].
     * @math   H[U;V] = [U;V]E; f = 1/(1+exp(E/T)).
     * @math   ρ = V(1-f)Vᵀ + UfUᵀ.
     * @math   κ = V(1-f)Uᵀ + UfVᵀ.
     * @output Updated H, U, V, Eqp, f, rho, and kappa.
     * @output Returned mean particle number N = Tr(ρ).
     * @note   T = k_B T_phys; same energy units as Eqp.
     * @note   Requires T ≥ 0 and no zero-energy modes.
     * @note   ε = λ + E(1-2‖V‖²); cutoff affects ρ and κ.
     * @note   EspCut_F = ∞ disables the cutoff.
     * @note   h_LN = h₀ + Γ + 4λ₂ρ - 2λ₂I; pre-solve ρ.
     */
    double update_UV_E_rho_kappa();
};


class HFBNucleus {
public:
    HFB hfb_neutron;
    HFB hfb_proton;

public:
    /**
     * @brief  Initialize species dimensions and workspaces.
     * @math   Nsp_n, Nsp_p.
     * @output Allocated solutions, fields, and workspaces.
     */
    HFBNucleus(int NspN_I_, int NspP_I_)
    : hfb_neutron(NspN_I_), hfb_proton(NspP_I_) {}

    /**
     * @brief  Initialize one-body fields in the derived model.
     * @math   h₀,n, h₀,p.
     * @output Initialized neutron and proton one-body fields.
     */
    virtual void initialize_h0() = 0;

    /**
     * @brief  Initialize HFB fields in the derived model.
     * @math   (N,Z) → (Γ_n,Δ_n,Γ_p,Δ_p)_initial.
     * @output Initialized Gamma and Delta for both species.
     */
    virtual void initialize_GammaDelta(int TargetN_I, int TargetZ_I) = 0;

    /**
     * @brief  Update both species using their joint densities.
     * @math   (ρ_n,κ_n,ρ_p,κ_p) → (Γ_n,Δ_n,Γ_p,Δ_p).
     * @output Overwritten neutron and proton fields.
     */
    virtual void update_Gamma_Delta() = 0;

    /**
     * @brief  Iterate the unblocked HFB equations.
     * @math   (N,Z) → HFB_converged.
     * @output Updated converged solver state.
     */
    void iterate(int TargetN_I, int TargetZ_I);

};

inline double HFB::update_UV_E_rho_kappa() {
    assert(temperature_F >= 0.0);
    assert(hfb_field.h0_F2D_sp_sp.rows() == Nsp_I && hfb_field.h0_F2D_sp_sp.cols() == Nsp_I);
    assert(hfb_field.Gamma_F2D_sp_sp.rows() == Nsp_I && hfb_field.Gamma_F2D_sp_sp.cols() == Nsp_I);
    assert(hfb_field.Delta_F2D_sp_sp.rows() == Nsp_I && hfb_field.Delta_F2D_sp_sp.cols() == Nsp_I);
    assert(hfb_field.h0_F2D_sp_sp.isApprox(hfb_field.h0_F2D_sp_sp.transpose(), 1.0e-12));
    assert(hfb_field.Gamma_F2D_sp_sp.isApprox(hfb_field.Gamma_F2D_sp_sp.transpose(), 1.0e-12));
    assert(hfb_field.Delta_F2D_sp_sp.isApprox(-hfb_field.Delta_F2D_sp_sp.transpose(), 1.0e-12));

    // H₁₁ = h₀ + Γ + 4λ₂ρ - (λ+2λ₂)I.
    H_F2D_2sp_2sp.topLeftCorner(Nsp_I, Nsp_I) = hfb_field.h0_F2D_sp_sp + hfb_field.Gamma_F2D_sp_sp + 4.0 * lambda2_F * hfb_solution.rho_F2D_sp_sp;
    H_F2D_2sp_2sp.topLeftCorner(Nsp_I, Nsp_I).diagonal().array() -= lambda_F + 2.0 * lambda2_F;

    // H₁₂ = Δ; H₂₁ = -Δ; H₂₂ = -H₁₁.
    H_F2D_2sp_2sp.topRightCorner(Nsp_I, Nsp_I) = hfb_field.Delta_F2D_sp_sp;
    H_F2D_2sp_2sp.bottomLeftCorner(Nsp_I, Nsp_I) = -hfb_field.Delta_F2D_sp_sp;
    H_F2D_2sp_2sp.bottomRightCorner(Nsp_I, Nsp_I) = -H_F2D_2sp_2sp.topLeftCorner(Nsp_I, Nsp_I);

    // H[U;V] = [U;V]E.
    H_eigensolver.compute(H_F2D_2sp_2sp);
    assert(H_eigensolver.info() == Eigen::Success);
    const Eigen::VectorXd& eigenvalues_F1D_state = H_eigensolver.eigenvalues();
    const Eigen::MatrixXd& eigenvectors_F2D_2sp_state = H_eigensolver.eigenvectors();

    // E₁ ≤ ⋯ ≤ E₂Nsp; retain Nsp positive-energy modes.
    assert(eigenvalues_F1D_state(Nsp_I - 1) < 0.0);
    assert(eigenvalues_F1D_state(Nsp_I) > 0.0);
    hfb_solution.Eqp_F1D_qp = eigenvalues_F1D_state.tail(Nsp_I);
    hfb_solution.U_F2D_sp_qp = eigenvectors_F2D_2sp_state.topRightCorner(Nsp_I, Nsp_I);
    hfb_solution.V_F2D_sp_qp = eigenvectors_F2D_2sp_state.bottomRightCorner(Nsp_I, Nsp_I);

    // ε = λ + E(1-2‖V‖²); ε > Ecut + tail → factors = 0.
    const double EspCutTail_F = std::log(1.0 / 1.0e-6 - 1.0) / 100.0;
    Eigen::VectorXd factorU_F1D_qp{};
    Eigen::VectorXd factorV_F1D_qp{};
    factorU_F1D_qp.resize(Nsp_I);
    factorV_F1D_qp.resize(Nsp_I);
    factorU_F1D_qp.setZero();
    factorV_F1D_qp.setZero();

    // f = (1-tanh(E/(2T)))/2; active factors = (f,1-f).
    for (int qp_I = 0; qp_I < Nsp_I; ++qp_I) {
        const double Eqp_F = hfb_solution.Eqp_F1D_qp(qp_I);
        hfb_solution.f_F1D_qp(qp_I) = temperature_F > 1.0e-12 ? 0.5 * (1.0 - std::tanh(0.5 * Eqp_F / temperature_F)) : 0.0;
        const double Esp_F = lambda_F + Eqp_F * (1.0 - 2.0 * hfb_solution.V_F2D_sp_qp.col(qp_I).squaredNorm());
        if (Esp_F > EspCut_F + EspCutTail_F) {continue;}
        factorU_F1D_qp(qp_I) = hfb_solution.f_F1D_qp(qp_I);
        factorV_F1D_qp(qp_I) = 1.0 - hfb_solution.f_F1D_qp(qp_I);
    }

    // ρ = V(1-f)Vᵀ + UfUᵀ.
    hfb_solution.rho_F2D_sp_sp.noalias() = hfb_solution.V_F2D_sp_qp * factorV_F1D_qp.asDiagonal() * hfb_solution.V_F2D_sp_qp.transpose();
    hfb_solution.rho_F2D_sp_sp.noalias() += hfb_solution.U_F2D_sp_qp * factorU_F1D_qp.asDiagonal() * hfb_solution.U_F2D_sp_qp.transpose();

    // κ = V(1-f)Uᵀ + UfVᵀ.
    hfb_solution.kappa_F2D_sp_sp.noalias() = hfb_solution.V_F2D_sp_qp * factorV_F1D_qp.asDiagonal() * hfb_solution.U_F2D_sp_qp.transpose();
    hfb_solution.kappa_F2D_sp_sp.noalias() += hfb_solution.U_F2D_sp_qp * factorU_F1D_qp.asDiagonal() * hfb_solution.V_F2D_sp_qp.transpose();

    return hfb_solution.rho_F2D_sp_sp.trace();
}
