/**
 * @file    group_so3_su2.hpp
 * @author  cailea
 * @date    2026-09-11
 * @brief   Angular-momentum generators and SO(3)/SU(2) rotation representations.
 */

#pragma once

#include <cassert>
#include <cmath>
#include <complex>

#include <Eigen/Core>
#include <Eigen/Eigenvalues>

/**
 * @brief  Store angular-momentum matrices and the Jy eigensystem.
 * @math   d(β) = Q diag(exp(-iβλ)) Q†; Jy = Q diag(λ) Q†.
 * @note   Matrices represent J_a/ℏ; M = -J,…,J.
 */
class RepresentationSpin {
public:
    int twoJ_I = 0;

    Eigen::MatrixXd Jx_F2D_M_M{};
    Eigen::MatrixXcd Jy_C2D_M_M{};
    Eigen::MatrixXd Jz_F2D_M_M{};
    Eigen::MatrixXd Jp_F2D_M_M{};
    Eigen::MatrixXd Jm_F2D_M_M{};

    Eigen::VectorXd JyEigenvalue_F1D_eigen{};
    Eigen::MatrixXcd JyEigenvector_C2D_M_eigen{};

    Eigen::MatrixXcd Ry_C2D_M_M{};
    Eigen::MatrixXcd Rz_C2D_M_M{};
    Eigen::MatrixXcd RzRyRz_C2D_M_M{};

public:
    /**
     * @brief  Initialize angular momentum using Hermitian eigendecomposition.
     * @math   J = twoJ_I_/2 ≥ 0; Jy = Q diag(λ) Q†.
     * @output Jp, Jm, Jz, Jx, Jy and Jy eigensystem.
     */
    explicit RepresentationSpin(int twoJ_I_) {
        assert(twoJ_I_ >= 0);
        twoJ_I = twoJ_I_;

        // dim = 2J + 1.
        Jx_F2D_M_M.resize(twoJ_I + 1, twoJ_I + 1);
        Jy_C2D_M_M.resize(twoJ_I + 1, twoJ_I + 1);
        Jz_F2D_M_M.resize(twoJ_I + 1, twoJ_I + 1);
        Jp_F2D_M_M.resize(twoJ_I + 1, twoJ_I + 1);
        Jm_F2D_M_M.resize(twoJ_I + 1, twoJ_I + 1);
        JyEigenvalue_F1D_eigen.resize(twoJ_I + 1);
        JyEigenvector_C2D_M_eigen.resize(twoJ_I + 1, twoJ_I + 1);

        // dim(R) = (2J + 1) × (2J + 1).
        Ry_C2D_M_M.resize(twoJ_I + 1, twoJ_I + 1);
        Rz_C2D_M_M.resize(twoJ_I + 1, twoJ_I + 1);
        RzRyRz_C2D_M_M.resize(twoJ_I + 1, twoJ_I + 1);

        build_Jp_Jm_Jz();
        build_Jx_Jy_Jz();
    }

    /**
     * @brief  Build ladder and diagonal angular-momentum matrices.
     * @math   (J±)_(M±1,M) = √((J∓M)(J±M+1)); (Jz)_(M,M) = M.
     * @output Jp_F2D_M_M, Jm_F2D_M_M, Jz_F2D_M_M.
     */
    void build_Jp_Jm_Jz();

    /**
     * @brief  Build Cartesian matrices and Hermitian Jy eigendecomposition.
     * @math   Jx = (Jp+Jm)/2; Jy = (Jp-Jm)/(2i); (Jz)_(M,M) = M.
     * @output Jx, Jy, Jz and Jy eigensystem.
     * @note   Requires current Jp and Jm matrices.
     */
    void build_Jx_Jy_Jz();

    /**
     * @brief  Evaluate y rotation using the Jy eigensystem.
     * @math   Ry(β) = Q diag(exp(-iβλ)) Q†.
     * @output Const reference to Ry_C2D_M_M.
     * @note   Requires the current Jy eigensystem; β in radians.
     * @note   Reference aliases the reusable rotation buffer.
     */
    const Eigen::MatrixXcd& calc_Ry(double beta_F);

    /**
     * @brief  Evaluate the Euler rotation matrix.
     * @math   R(α,β,γ) = exp(-iαJz) Ry(β) exp(-iγJz).
     * @output Const reference to RzRyRz_C2D_M_M.
     * @note   Requires the current Jy eigensystem; angles in radians.
     * @note   Reference aliases the reusable rotation buffer.
     */
    const Eigen::MatrixXcd& calc_R(double alpha_F, double beta_F, double gamma_F);
};

inline void RepresentationSpin::build_Jp_Jm_Jz() {
    const double J_F = 0.5 * twoJ_I;
    Jp_F2D_M_M.setZero();
    Jz_F2D_M_M.setZero();

    // (Jp)_(M+1,M) = √((J-M)(J+M+1)).
    for (int M_I = 0; M_I < twoJ_I; ++M_I) {
        const double M_F = M_I - J_F;
        Jp_F2D_M_M(M_I + 1, M_I) = std::sqrt((J_F - M_F) * (J_F + M_F + 1.0));
    }
    Jm_F2D_M_M = Jp_F2D_M_M.transpose();
    for (int M_I = 0; M_I <= twoJ_I; ++M_I) { Jz_F2D_M_M(M_I, M_I) = M_I - J_F; }
}

inline void RepresentationSpin::build_Jx_Jy_Jz() {
    // Jx = (Jp+Jm)/2; Jy = (Jp-Jm)/(2i); Jz|M⟩ = M|M⟩.
    Jx_F2D_M_M = 0.5 * (Jp_F2D_M_M + Jm_F2D_M_M);
    Jy_C2D_M_M = std::complex<double>(0.0, -0.5) * (Jp_F2D_M_M - Jm_F2D_M_M);
    Jz_F2D_M_M.setZero();
    for (int M_I = 0; M_I <= twoJ_I; ++M_I) { Jz_F2D_M_M(M_I, M_I) = M_I - 0.5 * twoJ_I; }

    // Jy = Q diag(λ) Q†.
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> Jy_solver(Jy_C2D_M_M);
    assert(Jy_solver.info() == Eigen::Success);
    JyEigenvalue_F1D_eigen = Jy_solver.eigenvalues();
    JyEigenvector_C2D_M_eigen = Jy_solver.eigenvectors();
}

inline const Eigen::MatrixXcd& RepresentationSpin::calc_Ry(double beta_F) {
    assert(std::isfinite(beta_F));
    // Ry(β) = Q diag(exp(-iβλ)) Q†.
    Ry_C2D_M_M.noalias() = JyEigenvector_C2D_M_eigen * (std::complex<double>(0.0, -beta_F) * JyEigenvalue_F1D_eigen.array()).exp().matrix().asDiagonal() * JyEigenvector_C2D_M_eigen.adjoint();
    return Ry_C2D_M_M;
}

inline const Eigen::MatrixXcd& RepresentationSpin::calc_R(double alpha_F, double beta_F, double gamma_F) {
    assert(std::isfinite(alpha_F) && std::isfinite(gamma_F));
    calc_Ry(beta_F);

    // Rz(α) Ry(β).
    Rz_C2D_M_M.setZero();
    Rz_C2D_M_M.diagonal() = (std::complex<double>(0.0, -alpha_F) * Jz_F2D_M_M.diagonal().array()).exp().matrix();
    RzRyRz_C2D_M_M.noalias() = Rz_C2D_M_M.diagonal().asDiagonal() * Ry_C2D_M_M;

    // R_(M',M) *= exp(-iγM).
    Rz_C2D_M_M.diagonal() = (std::complex<double>(0.0, -gamma_F) * Jz_F2D_M_M.diagonal().array()).exp().matrix();
    for (int M_I = 0; M_I <= twoJ_I; ++M_I) { RzRyRz_C2D_M_M.col(M_I) *= Rz_C2D_M_M(M_I, M_I); }
    return RzRyRz_C2D_M_M;
}
