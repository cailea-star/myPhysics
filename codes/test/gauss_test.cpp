#include <cmath>
#include <iomanip>
#include <ostream>
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include "integration_gauss.hpp"


void test_hermite(int n_max) {
    // (2 nx - 1) >= n1 + n2,  nx >= nmax + 0.5; 
    int nx_points = n_max + 1;
    Eigen::VectorXi n_I1D_n(n_max + 1); 
    for (int i = 0; i <= n_max; ++i) {n_I1D_n(i) = i;}
    GaussHermiteMeshes hermite_meshes(n_I1D_n, nx_points);
    std::cout << std::scientific << std::setprecision(4) << std::right;
    std::cout << "Gauss-Hermite Meshes:\n";
    Eigen::VectorXd w_F1D_x = hermite_meshes.w_F1D_x;
    Eigen::VectorXd N_F1D_n = hermite_meshes.N_F1D_n;
    Eigen::MatrixXd H_F2D_n_x = hermite_meshes.H_F2D_n_x;
    Eigen::MatrixXd resultH_F2D_n_n = N_F1D_n.asDiagonal() * H_F2D_n_x * w_F1D_x.asDiagonal() * H_F2D_n_x.transpose() * N_F1D_n.asDiagonal();
    std::cout << resultH_F2D_n_n << "\n\n";
}

void test_laguerre(int n_max, double a_F) {
    // (2 nx - 1) >= n1 + n2 + a, nx >= nmax + a / 2.0 + 0.5;
    int nx_points = n_max + std::ceil(a_F / 2.0) + 1;
    Eigen::VectorXi n_I1D_na(n_max + 1); 
    for (int i = 0; i <= n_max; ++i) {n_I1D_na(i) = i;}
    Eigen::VectorXd a_F1D_na(n_max + 1); 
    for (int i = 0; i <= n_max; ++i) {a_F1D_na(i) = a_F;}
    GaussLaguerreMeshes laguerre_meshes(n_I1D_na, a_F1D_na, nx_points);
    std::cout << std::scientific << std::setprecision(4) << std::right;
    std::cout << "Gauss-Laguerre Meshes with a = " << a_F << ":\n";
    Eigen::VectorXd wxa_F1D_x = laguerre_meshes.w_F1D_x.array() * laguerre_meshes.x_F1D_x.array().pow(a_F);
    Eigen::VectorXd N_F1D_na = laguerre_meshes.N_F1D_na;
    Eigen::MatrixXd L_F2D_na_x = laguerre_meshes.L_F2D_na_x;
    Eigen::MatrixXd resultL_F2D_na_na = N_F1D_na.asDiagonal() * L_F2D_na_x * wxa_F1D_x.asDiagonal() * L_F2D_na_x.transpose() * N_F1D_na.asDiagonal();
    std::cout << resultL_F2D_na_na << "\n\n";
}

void test_legendre(int lmin, int l_max, int m_F) {
    // (2 nx - 1) >= l1 + l2, nx >= lmax + 0.5;
    int nx_points = l_max + 1;
    Eigen::VectorXi l_I1D_lm(l_max + 1 - lmin); 
    for (int i = lmin; i <= l_max; ++i) {l_I1D_lm(i - lmin) = i;}
    Eigen::VectorXi m_I1D_lm(l_max + 1 - lmin); 
    for (int i = lmin; i <= l_max; ++i) {m_I1D_lm(i - lmin) = m_F;}
    GaussLegendreMeshes legendre_meshes(l_I1D_lm, m_I1D_lm, nx_points);
    std::cout << std::scientific << std::setprecision(4) << std::right;
    std::cout << "Gauss-Legendre Meshes with m = " << m_F << ":\n";
    Eigen::VectorXd w_F1D_x = legendre_meshes.w_F1D_x;
    Eigen::VectorXd N_F1D_lm = legendre_meshes.N_F1D_lm;
    Eigen::MatrixXd P_F2D_lm_x = legendre_meshes.P_F2D_lm_x;
    Eigen::MatrixXd resultP_F2D_lm_lm = N_F1D_lm.asDiagonal() * P_F2D_lm_x * w_F1D_x.asDiagonal() * P_F2D_lm_x.transpose() * N_F1D_lm.asDiagonal();
    std::cout << resultP_F2D_lm_lm << "\n\n";
}

int main() {
    test_hermite(10);
    test_legendre(0, 8, 0);
    test_legendre(3, 11, 3);
    test_laguerre(6, 0.0);
    test_laguerre(6, 1.0);
    test_laguerre(6, 5.0);
    test_laguerre(6, 10.0);
}
