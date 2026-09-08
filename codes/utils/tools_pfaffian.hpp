/**
 * @file    tools_pfaffian.hpp
 * @author  cailea
 * @date    2026-09-07
 * @brief   Complex skew-symmetric matrix Pfaffians.
 */

#pragma once

#include <cassert>
#include <complex>
#include <Eigen/Core>

using doubleC = std::complex<double>;

/**
 * @brief  Calculate Pfaffians using pivoted skew-symmetric elimination.
 * @math   pf(X)² = det(X); pf(∅) = 1.
 * @output pf(X); X overwritten.
 * @note   Requires finite, even-order square X with Xᵀ = -X.
 * @note   Exact-zero pivots return zero; no magnitude cutoff.
 */
inline doubleC calc_pfaffian(Eigen::MatrixXcd& X_C2D_i_i) {
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
