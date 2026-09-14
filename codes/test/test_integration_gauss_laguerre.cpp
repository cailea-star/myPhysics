/**
 * @file    test_integration_gauss_laguerre.cpp
 * @author  cailea
 * @date    2026-04-27
 * @brief   Orthogonality test for Gauss-Laguerre meshes and polynomials.
 */

#include <Eigen/Dense>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include "integration_gauss.hpp"

/**
 * @brief  Test generalized Gauss-Laguerre polynomial orthogonality.
 * @math   M_ij = ∑_x N_i L_i^α(x) x^α w_x L_j^α(x) N_j = δ_ij
 * @output Prints labeled matrices and asserts the maximum orthogonality error.
 */
int main() {
    // α ∈ {0,1,5,10} → max|M(α)-I|.
    double errMax_F = 0.0;
    double tol_F = 1.0e-10;
    std::cout << std::scientific << std::setprecision(4) << std::right;

    // (α,{n_i}) → M_ij(α) → stdout.
    for (double a_F : {0.0, 1.0, 5.0, 10.0}) {
        int nmax_I = 6;
        int Nx_I = nmax_I + std::ceil(a_F / 2.0) + 1;
        Eigen::VectorXi n_I1D_na(nmax_I + 1);
        for (int i = 0; i <= nmax_I; ++i) {n_I1D_na(i) = i;}
        Eigen::VectorXd a_F1D_na(nmax_I + 1);
        for (int i = 0; i <= nmax_I; ++i) {a_F1D_na(i) = a_F;}
        GaussLaguerreMeshes laguerre_meshes(n_I1D_na, a_F1D_na, Nx_I);
        Eigen::VectorXd wxa_F1D_x = laguerre_meshes.w_F1D_x.array() * laguerre_meshes.x_F1D_x.array().pow(a_F);
        Eigen::VectorXd N_F1D_na = laguerre_meshes.N_F1D_na;
        Eigen::MatrixXd L_F2D_na_x = laguerre_meshes.L_F2D_na_x;
        Eigen::MatrixXd M_F2D_i_j = N_F1D_na.asDiagonal() * L_F2D_na_x * wxa_F1D_x.asDiagonal() * L_F2D_na_x.transpose() * N_F1D_na.asDiagonal();
        Eigen::MatrixXd Mref_F2D_i_j = Eigen::MatrixXd::Identity(M_F2D_i_j.rows(), M_F2D_i_j.cols());
        errMax_F = std::max(errMax_F, (M_F2D_i_j - Mref_F2D_i_j).cwiseAbs().maxCoeff());
        std::cout << "\nGauss-Laguerre\n";
        std::cout << "Input: nmax = " << nmax_I << ", α = " << a_F << ", Nx = " << Nx_I << "\n";
        std::cout << "Reference: M\n" << Mref_F2D_i_j << "\n";
        std::cout << "Computed: M\n" << M_F2D_i_j << "\n";
    }

    // α_w=1/2; α-α_w ∈ {0,1,2} → M(α)=I.
    for (double a_F : {0.5, 1.5, 2.5}) {
        int nmax_I = 6;
        int Nx_I = nmax_I + static_cast<int>(std::ceil((a_F - 0.5) / 2.0)) + 1;
        Eigen::VectorXi n_I1D_na(nmax_I + 1);
        for (int i = 0; i <= nmax_I; ++i) {n_I1D_na(i) = i;}
        Eigen::VectorXd a_F1D_na(nmax_I + 1);
        a_F1D_na.setConstant(a_F);

        // Both constructors → identical x^(1/2)e^-x grids.
        GaussLaguerreMeshes laguerre_grid(Nx_I, 0.5);
        GaussLaguerreMeshes laguerre_meshes(n_I1D_na, a_F1D_na, Nx_I, 0.5);
        assert((laguerre_grid.x_F1D_x - laguerre_meshes.x_F1D_x).cwiseAbs().maxCoeff() < tol_F);
        assert((laguerre_grid.w_F1D_x - laguerre_meshes.w_F1D_x).cwiseAbs().maxCoeff() < tol_F);
        assert(std::abs(laguerre_grid.w_F1D_x.sum() - std::tgamma(1.5)) < tol_F);

        // x^α_w e^-x × x^(α-α_w) → x^α e^-x.
        Eigen::VectorXd wxa_F1D_x = laguerre_meshes.w_F1D_x.array() * laguerre_meshes.x_F1D_x.array().pow(a_F - 0.5);
        Eigen::VectorXd N_F1D_na = laguerre_meshes.N_F1D_na;
        Eigen::MatrixXd L_F2D_na_x = laguerre_meshes.L_F2D_na_x;
        Eigen::MatrixXd M_F2D_i_j = N_F1D_na.asDiagonal() * L_F2D_na_x * wxa_F1D_x.asDiagonal() * L_F2D_na_x.transpose() * N_F1D_na.asDiagonal();
        Eigen::MatrixXd Mref_F2D_i_j = Eigen::MatrixXd::Identity(M_F2D_i_j.rows(), M_F2D_i_j.cols());
        errMax_F = std::max(errMax_F, (M_F2D_i_j - Mref_F2D_i_j).cwiseAbs().maxCoeff());
        std::cout << "\nGauss-Laguerre: weight α = 0.5, polynomial α = " << a_F << "\n";
        std::cout << "Computed: M\n" << M_F2D_i_j << "\n";
    }

    assert(errMax_F < tol_F);
    return 0;
}
