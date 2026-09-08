/**
 * @file    hfb_pfaffian.hpp
 * @author  cailea
 * @date    2026-09-08
 * @brief   Pfaffians and HFB matrix-element kernel interfaces.
 */

#pragma once

#include <cassert>
#include <complex>
#include <functional>
#include <vector>
#include <Eigen/Core>

using doubleC = std::complex<double>;

/**
 * @brief  HFB kernels using generalized Wick Pfaffians.
 * @math   K = ⟨Φ₁|Φ₂⟩ pf(S).
 * @note   One-body and two-body kernel implementations remain pending.
 */
class HFBPfaffian {
public:
    // Outputs: overlap = ⟨Φ₁|Φ₂(β,φ)⟩, U₂(β,φ), V₂(β,φ).
    using VacuumFunc = std::function<void(int beta_I, int phi_I, doubleC& overlap_C, Eigen::MatrixXcd& U2g_C2D_2sp_2qp, Eigen::MatrixXcd& V2g_C2D_2sp_2qp)>;

    int Nsp_I = 0;
    std::vector<std::vector<int>> Nqp1_I2D_config_qp{};
    std::vector<std::vector<int>> Nqp2_I2D_config_qp{};
    VacuumFunc vacuum_Func{};

    Eigen::MatrixXcd U1_C2D_2sp_2qp{};
    Eigen::MatrixXcd V1_C2D_2sp_2qp{};
    Eigen::MatrixXcd U2_C2D_2sp_2qp{};
    Eigen::MatrixXcd V2_C2D_2sp_2qp{};
    Eigen::MatrixXcd U2g_C2D_2sp_2qp{};
    Eigen::MatrixXcd V2g_C2D_2sp_2qp{};
    Eigen::MatrixXcd Ag_C2D_2qp_2qp{};

    doubleC overlap_C = doubleC(0.0, 0.0);
    Eigen::MatrixXcd rho_C2D_sp_sp{};
    Eigen::MatrixXcd kappa_C2D_sp_sp{};
    Eigen::MatrixXcd kappabar_C2D_sp_sp{};

public:
    /**
     * @brief  Initialize HFB vacuum data and kernel buffers.
     * @math   μ_{a,i} < μ_{a,i+1}; ρ = κ = κ̄ = 0.
     * @output Fixed U₁,V₁; Nsp × Nsp kernel buffers.
     * @note   Quasiparticle indices are zero-based and strictly increasing.
     * @note   U,V contain the complete quasiparticle basis.
     * @note   Callback supplies transformed U₂,V₂ and the vacuum overlap.
     */
    HFBPfaffian(int Nsp_I_, const std::vector<std::vector<int>>& Nqp1_I2D_config_qp_, const std::vector<std::vector<int>>& Nqp2_I2D_config_qp_, const Eigen::MatrixXcd& U1_C2D_sp_qp_, const Eigen::MatrixXcd& V1_C2D_sp_qp_, const VacuumFunc& vacuum_Func_) {
        Nsp_I = Nsp_I_;
        Nqp1_I2D_config_qp = Nqp1_I2D_config_qp_;
        Nqp2_I2D_config_qp = Nqp2_I2D_config_qp_;
        U1_C2D_2sp_2qp = U1_C2D_sp_qp_;
        V1_C2D_2sp_2qp = V1_C2D_sp_qp_;
        vacuum_Func = vacuum_Func_;
        assert(Nsp_I >= 0);

        // U₁, V₁ ∈ ℂ^{Nsp×Nsp}.
        assert(U1_C2D_2sp_2qp.rows() == Nsp_I && U1_C2D_2sp_2qp.cols() == Nsp_I);
        assert(V1_C2D_2sp_2qp.rows() == Nsp_I && V1_C2D_2sp_2qp.cols() == Nsp_I);

        // U₂, V₂: untransformed vacuum buffers.
        U2_C2D_2sp_2qp.resize(Nsp_I, Nsp_I);
        V2_C2D_2sp_2qp.resize(Nsp_I, Nsp_I);
        U2_C2D_2sp_2qp.setZero();
        V2_C2D_2sp_2qp.setZero();

        // U₂(g), V₂(g): vacuum_Func outputs.
        U2g_C2D_2sp_2qp.resize(Nsp_I, Nsp_I);
        V2g_C2D_2sp_2qp.resize(Nsp_I, Nsp_I);
        U2g_C2D_2sp_2qp.setZero();
        V2g_C2D_2sp_2qp.setZero();

        // A(g) = U₁ᵀU₂(g)* + V₁ᵀV₂(g)*.
        Ag_C2D_2qp_2qp.resize(Nsp_I, Nsp_I);
        Ag_C2D_2qp_2qp.setZero();

        // ρ, κ, κ̄ ∈ ℂ^{Nsp×Nsp}.
        rho_C2D_sp_sp.resize(Nsp_I, Nsp_I);
        kappa_C2D_sp_sp.resize(Nsp_I, Nsp_I);
        kappabar_C2D_sp_sp.resize(Nsp_I, Nsp_I);
        rho_C2D_sp_sp.setZero();
        kappa_C2D_sp_sp.setZero();
        kappabar_C2D_sp_sp.setZero();
    }

    /**
     * @brief  Calculate one-body kernels using generalized Wick Pfaffians.
     * @math   K₁₂ = ⟨Φ₁;μ₁|c†₁c₂|Φ₂;μ₂(β,φ)⟩.
     * @output Unnormalized kernel, including the vacuum overlap.
     * @note   config1, config2 select zero-based left/right configurations.
     * @note   Requires |μ₁| + |μ₂| ≤ Nsp and even.
     * @note   β, φ are mesh indices; sp indices are zero-based.
     * @note   Requires nonzero vacuum overlap and quasiparticle data.
     * @note   Declaration only; implementation pending.
     */
    doubleC calc_one_body(int beta_I, int phi_I, int config1_I, int config2_I, int sp1_I, int sp2_I);

    /**
     * @brief  Calculate two-body kernels using generalized Wick Pfaffians.
     * @math   K₁₂₃₄ = ⟨Φ₁;μ₁|c†₁c†₂c₄c₃|Φ₂;μ₂(β,φ)⟩.
     * @output Unnormalized kernel; excludes operator coefficients and 1/2.
     * @note   config1, config2 select zero-based left/right configurations.
     * @note   Requires |μ₁| + |μ₂| ≤ Nsp and even.
     * @note   β, φ are mesh indices; sp indices are zero-based.
     * @note   Requires nonzero vacuum overlap and quasiparticle data.
     * @note   Declaration only; implementation pending.
     */
    doubleC calc_two_body(int beta_I, int phi_I, int config1_I, int config2_I, int sp1_I, int sp2_I, int sp3_I, int sp4_I);

    /**
     * @brief  Calculate Pfaffians using pivoted skew-symmetric elimination.
     * @math   pf(X)² = det(X); pf(∅) = 1.
     * @output pf(X); X overwritten.
     * @note   Requires finite, even-order square X with Xᵀ = -X.
     * @note   Exact-zero pivots return zero; no magnitude cutoff.
     */
    static doubleC calc_pfaffian(Eigen::MatrixXcd& X_C2D_i_i);
};

/**
 * @brief  Calculate Pfaffians using pivoted skew-symmetric elimination.
 * @math   pf(X)² = det(X); pf(∅) = 1.
 * @output pf(X); X overwritten.
 * @note   Requires finite, even-order square X with Xᵀ = -X.
 * @note   Exact-zero pivots return zero; no magnitude cutoff.
 */
inline doubleC HFBPfaffian::calc_pfaffian(Eigen::MatrixXcd& X_C2D_i_i) {
    // X ∈ ℂ^{2p×2p}; Xᵀ = -X.
    const Eigen::Index N_I = X_C2D_i_i.rows();
    assert(X_C2D_i_i.cols() == N_I && N_I % 2 == 0);
    assert(X_C2D_i_i.allFinite());
    assert(X_C2D_i_i.isApprox(-X_C2D_i_i.transpose(), 1.0e-12));
    doubleC pf_C(1.0, 0.0); // pf(∅) = 1.

    for (Eigen::Index k_I = 0; k_I < N_I; k_I += 2) {
        // Pivot: argmax_{j>k}|X_kj|.
        Eigen::Index pivotOffset_I = 0;
        const double pivotAbs_F = X_C2D_i_i.row(k_I).segment(k_I + 1, N_I - k_I - 1).cwiseAbs().maxCoeff(&pivotOffset_I);
        const Eigen::Index pivot_I = k_I + 1 + pivotOffset_I;
        if (pivotAbs_F == 0.0) {return doubleC(0.0, 0.0);}

        // X → PXPᵀ; pf(PXPᵀ) = det(P)pf(X).
        if (pivot_I != k_I + 1) {
            X_C2D_i_i.row(k_I + 1).swap(X_C2D_i_i.row(pivot_I));
            X_C2D_i_i.col(k_I + 1).swap(X_C2D_i_i.col(pivot_I));
            pf_C = -pf_C;
        }

        const doubleC pivot_C = X_C2D_i_i(k_I, k_I + 1);
        pf_C *= pivot_C;

        // X_ij ← X_ij+(X_{k+1,i}X_kj-X_ki X_{k+1,j})/X_{k,k+1}.
        for (Eigen::Index j_I = k_I + 3; j_I < N_I; ++j_I) {
            for (Eigen::Index i_I = k_I + 2; i_I < j_I; ++i_I) {
                X_C2D_i_i(i_I, j_I) += X_C2D_i_i(k_I + 1, i_I) * (X_C2D_i_i(k_I, j_I) / pivot_C) - X_C2D_i_i(k_I, i_I) * (X_C2D_i_i(k_I + 1, j_I) / pivot_C);
                X_C2D_i_i(j_I, i_I) = -X_C2D_i_i(i_I, j_I);
            }
        }
    }
    return pf_C;
}
