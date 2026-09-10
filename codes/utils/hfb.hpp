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


class HFBNucleus {
public:
    using OneBodyFunc = std::function<void(Eigen::MatrixXd& h0_F2D_sp_sp)>;
    using TwoBodyFunc = std::function<void(const HFBSolution& solution_n, const HFBSolution& solution_p, HFBField& field_n, HFBField& field_p)>;

    HFB hfb_neutron;
    HFB hfb_proton;

    OneBodyFunc build_onebody_neutron{};
    OneBodyFunc build_onebody_proton{};
    TwoBodyFunc build_twobody{};

public:
    /**
     * @brief  Initialize species and build their one-body matrices.
     * @math   h_q = h₀,q + Γ_q; q ∈ {n,p}.
     * @output Built one-body matrices, zeroed fields, stored callbacks.
     */
    HFBNucleus(int NspN_I_, int NspP_I_, const OneBodyFunc& build_onebodyN_, const OneBodyFunc& build_onebodyP_, const TwoBodyFunc& build_twobody_)
    : hfb_neutron(NspN_I_), hfb_proton(NspP_I_) {
        assert(build_onebodyN_);
        assert(build_onebodyP_);
        assert(build_twobody_);
        build_onebody_neutron = build_onebodyN_;
        build_onebody_proton = build_onebodyP_;
        build_twobody = build_twobody_;

        // (h₀,n, h₀,p) ← one-body callbacks.
        build_onebody_neutron(hfb_neutron.hfb_field.h0_F2D_sp_sp);
        build_onebody_proton(hfb_proton.hfb_field.h0_F2D_sp_sp);
    }

    /**
     * @brief  Update both species using their joint densities.
     * @math   (ρ_n,κ_n,ρ_p,κ_p) → (Γ_n,Δ_n,Γ_p,Δ_p).
     * @output Overwritten neutron and proton fields.
     */
    void update_Gamma_Delta();

};

inline double HFB::update_UV_E_rho_kappa(double lambda_F, double temperature_F) {
    assert(temperature_F >= 0.0);
    assert(hfb_field.h0_F2D_sp_sp.rows() == Nsp_I && hfb_field.h0_F2D_sp_sp.cols() == Nsp_I);
    assert(hfb_field.Gamma_F2D_sp_sp.rows() == Nsp_I && hfb_field.Gamma_F2D_sp_sp.cols() == Nsp_I);
    assert(hfb_field.Delta_F2D_sp_sp.rows() == Nsp_I && hfb_field.Delta_F2D_sp_sp.cols() == Nsp_I);
    assert(hfb_field.h0_F2D_sp_sp.isApprox(hfb_field.h0_F2D_sp_sp.transpose(), 1.0e-12));
    assert(hfb_field.Gamma_F2D_sp_sp.isApprox(hfb_field.Gamma_F2D_sp_sp.transpose(), 1.0e-12));
    assert(hfb_field.Delta_F2D_sp_sp.isApprox(-hfb_field.Delta_F2D_sp_sp.transpose(), 1.0e-12));

    // H₁₁ = h₀ + Γ - λI.
    H_F2D_2sp_2sp.topLeftCorner(Nsp_I, Nsp_I) = hfb_field.h0_F2D_sp_sp + hfb_field.Gamma_F2D_sp_sp;
    H_F2D_2sp_2sp.topLeftCorner(Nsp_I, Nsp_I).diagonal().array() -= lambda_F;

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

    // T ≤ 10⁻¹² → f = 0; otherwise f = e⁻ᴱᐟᵀ/(1+e⁻ᴱᐟᵀ).
    hfb_solution.f_F1D_qp.setZero();
    if (temperature_F > 1.0e-12) {
        for (int qp_I = 0; qp_I < Nsp_I; ++qp_I) {
            const double expMinusEOverT_F = std::exp(-hfb_solution.Eqp_F1D_qp(qp_I) / temperature_F);
            hfb_solution.f_F1D_qp(qp_I) = expMinusEOverT_F / (1.0 + expMinusEOverT_F);
        }
    }

    // ρ = V(1-f)Vᵀ + UfUᵀ.
    hfb_solution.rho_F2D_sp_sp.noalias() = hfb_solution.V_F2D_sp_qp * (1.0 - hfb_solution.f_F1D_qp.array()).matrix().asDiagonal() * hfb_solution.V_F2D_sp_qp.transpose();
    hfb_solution.rho_F2D_sp_sp.noalias() += hfb_solution.U_F2D_sp_qp * hfb_solution.f_F1D_qp.asDiagonal() * hfb_solution.U_F2D_sp_qp.transpose();

    // κ = V(1-f)Uᵀ + UfVᵀ.
    hfb_solution.kappa_F2D_sp_sp.noalias() = hfb_solution.V_F2D_sp_qp * (1.0 - hfb_solution.f_F1D_qp.array()).matrix().asDiagonal() * hfb_solution.U_F2D_sp_qp.transpose();
    hfb_solution.kappa_F2D_sp_sp.noalias() += hfb_solution.U_F2D_sp_qp * hfb_solution.f_F1D_qp.asDiagonal() * hfb_solution.V_F2D_sp_qp.transpose();

    return hfb_solution.rho_F2D_sp_sp.trace();
}


/**
 * @brief  Update both species using their joint densities.
 * @math   (ρ_n,κ_n,ρ_p,κ_p) → (Γ_n,Δ_n,Γ_p,Δ_p).
 * @output Overwritten neutron and proton fields.
 */
inline void HFBNucleus::update_Gamma_Delta() {
    build_twobody(hfb_neutron.hfb_solution, hfb_proton.hfb_solution, hfb_neutron.hfb_field, hfb_proton.hfb_field);
}
