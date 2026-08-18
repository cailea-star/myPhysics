/**
 * @file    test_integration_gauss_hermite.cpp
 * @author  cailea
 * @date    2026-04-27
 * @brief   Orthogonality test for Gauss-Hermite meshes and polynomials.
 */

#include <Eigen/Dense>
#include <cassert>
#include <iomanip>
#include <iostream>
#include "integration_gauss.hpp"

/**
 * @brief  Test Gauss-Hermite polynomial orthogonality.
 * @math   M_ij = ∑_x N_i H_i(x) w_x H_j(x) N_j = δ_ij
 * @output Prints labeled matrix corners and asserts the maximum orthogonality error.
 */
int main() {
    int nmax_I = 10;
    int Nx_I = nmax_I + 1;
    double tol_F = 1.0e-10;
    Eigen::VectorXi n_I1D_n(nmax_I + 1);
    for (int i = 0; i <= nmax_I; ++i) {n_I1D_n(i) = i;}
    GaussHermiteMeshes hermite_meshes(n_I1D_n, Nx_I);
    Eigen::VectorXd w_F1D_x = hermite_meshes.w_F1D_x;
    Eigen::VectorXd N_F1D_n = hermite_meshes.N_F1D_n;
    Eigen::MatrixXd H_F2D_n_x = hermite_meshes.H_F2D_n_x;
    Eigen::MatrixXd M_F2D_i_j = N_F1D_n.asDiagonal() * H_F2D_n_x * w_F1D_x.asDiagonal() * H_F2D_n_x.transpose() * N_F1D_n.asDiagonal();
    Eigen::MatrixXd Mref_F2D_i_j = Eigen::MatrixXd::Identity(M_F2D_i_j.rows(), M_F2D_i_j.cols());
    double errMax_F = (M_F2D_i_j - Mref_F2D_i_j).cwiseAbs().maxCoeff();

    std::cout << std::scientific << std::setprecision(4) << std::right;
    std::cout << "\nGauss-Hermite\n";
    std::cout << "Input: nmax = " << nmax_I << ", Nx = " << Nx_I << "\n";
    std::cout << "Reference: M\n";
    for (int i_I = 0; i_I < 5; ++i_I) {std::cout << Mref_F2D_i_j.block(i_I, 0, 1, 5) << " ... " << Mref_F2D_i_j.block(i_I, Mref_F2D_i_j.cols() - 5, 1, 5) << "\n";}
    std::cout << "...\n";
    for (int i_I = Mref_F2D_i_j.rows() - 5; i_I < Mref_F2D_i_j.rows(); ++i_I) {std::cout << Mref_F2D_i_j.block(i_I, 0, 1, 5) << " ... " << Mref_F2D_i_j.block(i_I, Mref_F2D_i_j.cols() - 5, 1, 5) << "\n";}
    std::cout << "Computed: M\n";
    for (int i_I = 0; i_I < 5; ++i_I) {std::cout << M_F2D_i_j.block(i_I, 0, 1, 5) << " ... " << M_F2D_i_j.block(i_I, M_F2D_i_j.cols() - 5, 1, 5) << "\n";}
    std::cout << "...\n";
    for (int i_I = M_F2D_i_j.rows() - 5; i_I < M_F2D_i_j.rows(); ++i_I) {std::cout << M_F2D_i_j.block(i_I, 0, 1, 5) << " ... " << M_F2D_i_j.block(i_I, M_F2D_i_j.cols() - 5, 1, 5) << "\n";}

    assert(errMax_F < tol_F);
    return 0;
}
