/**
 * @file    test_integration_gauss_legendre.cpp
 * @author  cailea
 * @date    2026-04-27
 * @brief   Orthogonality test for Gauss-Legendre meshes and polynomials.
 */

#include <Eigen/Dense>
#include <algorithm>
#include <cassert>
#include <iomanip>
#include <iostream>
#include "integration_gauss.hpp"

/**
 * @brief  Test Gauss-Legendre associated-polynomial orthogonality.
 * @math   M_ij = ∑_x N_i P_i^m(x) w_x P_j^m(x) N_j = δ_ij
 * @output Prints labeled matrices and asserts the maximum orthogonality error.
 */
int main() {
    // (m,{l_i}) → M_ij → max|M-I|.
    double errMax_F = 0.0;
    double tol_F = 1.0e-10;
    std::cout << std::scientific << std::setprecision(4) << std::right;

    // (l,m,N_x)=(0:8,0,9) → M → max|M-I|.
    {
        int lmin_I = 0;
        int lmax_I = 8;
        int m_I = 0;
        int Nx_I = lmax_I + 1;
        Eigen::VectorXi l_I1D_lm(lmax_I + 1 - lmin_I);
        for (int i = lmin_I; i <= lmax_I; ++i) {l_I1D_lm(i - lmin_I) = i;}
        Eigen::VectorXi m_I1D_lm(lmax_I + 1 - lmin_I);
        for (int i = lmin_I; i <= lmax_I; ++i) {m_I1D_lm(i - lmin_I) = m_I;}
        GaussLegendreMeshes legendre_meshes(l_I1D_lm, m_I1D_lm, Nx_I);
        Eigen::VectorXd w_F1D_x = legendre_meshes.w_F1D_x;
        Eigen::VectorXd N_F1D_lm = legendre_meshes.N_F1D_lm;
        Eigen::MatrixXd P_F2D_lm_x = legendre_meshes.P_F2D_lm_x;
        Eigen::MatrixXd M_F2D_i_j = N_F1D_lm.asDiagonal() * P_F2D_lm_x * w_F1D_x.asDiagonal() * P_F2D_lm_x.transpose() * N_F1D_lm.asDiagonal();
        Eigen::MatrixXd Mref_F2D_i_j = Eigen::MatrixXd::Identity(M_F2D_i_j.rows(), M_F2D_i_j.cols());
        errMax_F = std::max(errMax_F, (M_F2D_i_j - Mref_F2D_i_j).cwiseAbs().maxCoeff());
        std::cout << "\nGauss-Legendre\n";
        std::cout << "Input: lmin = " << lmin_I << ", lmax = " << lmax_I << ", m = " << m_I << ", Nx = " << Nx_I << "\n";
        std::cout << "Reference: M\n" << Mref_F2D_i_j << "\n";
        std::cout << "Computed: M\n" << M_F2D_i_j << "\n";
    }

    // (l,m,N_x)=(3:11,3,12) → M → max|M-I|.
    {
        int lmin_I = 3;
        int lmax_I = 11;
        int m_I = 3;
        int Nx_I = lmax_I + 1;
        Eigen::VectorXi l_I1D_lm(lmax_I + 1 - lmin_I);
        for (int i = lmin_I; i <= lmax_I; ++i) {l_I1D_lm(i - lmin_I) = i;}
        Eigen::VectorXi m_I1D_lm(lmax_I + 1 - lmin_I);
        for (int i = lmin_I; i <= lmax_I; ++i) {m_I1D_lm(i - lmin_I) = m_I;}
        GaussLegendreMeshes legendre_meshes(l_I1D_lm, m_I1D_lm, Nx_I);
        Eigen::VectorXd w_F1D_x = legendre_meshes.w_F1D_x;
        Eigen::VectorXd N_F1D_lm = legendre_meshes.N_F1D_lm;
        Eigen::MatrixXd P_F2D_lm_x = legendre_meshes.P_F2D_lm_x;
        Eigen::MatrixXd M_F2D_i_j = N_F1D_lm.asDiagonal() * P_F2D_lm_x * w_F1D_x.asDiagonal() * P_F2D_lm_x.transpose() * N_F1D_lm.asDiagonal();
        Eigen::MatrixXd Mref_F2D_i_j = Eigen::MatrixXd::Identity(M_F2D_i_j.rows(), M_F2D_i_j.cols());
        errMax_F = std::max(errMax_F, (M_F2D_i_j - Mref_F2D_i_j).cwiseAbs().maxCoeff());
        std::cout << "\nGauss-Legendre\n";
        std::cout << "Input: lmin = " << lmin_I << ", lmax = " << lmax_I << ", m = " << m_I << ", Nx = " << Nx_I << "\n";
        std::cout << "Reference: M\n" << Mref_F2D_i_j << "\n";
        std::cout << "Computed: M\n" << M_F2D_i_j << "\n";
    }

    assert(errMax_F < tol_F);
    return 0;
}
