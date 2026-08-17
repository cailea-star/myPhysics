
#include <cmath>
#include <iomanip>
#include <ostream>
#include <iostream>
#include <Eigen/Core>
#include "ivp_numerov.hpp"

/**
 * y'' = F(x)·y 的测试函数，F(x) = -a, a > 0, 则解析解为 y(x) = A·cos(√a·x) + B·sin(√a·x);
 */
double func_scal(double x) {
    return -4.0; // F(x) = -a, a=1
}


/**
 * y'' = F(x)·y 的测试函数，F(x) = diag(1, -1), 则解析解为 
 * y1(x) = A2·cos(√a·x) + B2·sin(√a·x);
 * y2(x) = A1·exp(√a·x) + B1·exp(-√a·x);
 */
Eigen::MatrixXd Func_Mat(double x) {
    Eigen::MatrixXd F(2, 2);
    F << -4.0, 0.0,
          0.0, 4.0; 
    return F;
}


int main() {
    int nxI = 100;
    double xmin_F = 0.0;
    double xmax_F = 10.0;
    Eigen::VectorXd x_F1D_x = Eigen::VectorXd::LinSpaced(nxI, xmin_F, xmax_F);
    double ycos0 = std::cos(2.0 * x_F1D_x(0)); // A=1, B=0
    double ycos1 = std::cos(2.0 * x_F1D_x(1)); // A=1, B=0
    double ysin0 = std::sin(2.0 * x_F1D_x(0)); // A=0, B=1
    double ysin1 = std::sin(2.0 * x_F1D_x(1)); // A=0, B=1
    Eigen::VectorXd ycos_F1D_x = (x_F1D_x * 2.0).array().cos();
    Eigen::VectorXd ysin_F1D_x = (x_F1D_x * 2.0).array().sin();
    double yexpn0 = std::exp(-2.0 * x_F1D_x(0)); // A=1, B=0
    double yexpn1 = std::exp(-2.0 * x_F1D_x(1)); // A=1, B=0
    double yexpp0 = std::exp(2.0 * x_F1D_x(0)); // A=1, B=0
    double yexpp1 = std::exp(2.0 * x_F1D_x(1)); // A=0, B=1
    Eigen::VectorXd yexpp_F1D_x = (x_F1D_x * 2.0).array().exp();
    Eigen::VectorXd yexpn_F1D_x = (-x_F1D_x * 2.0).array().exp();


    std::cout << std::scientific << std::setprecision(4) << std::right;

    // 测试Numerov标量方法
    std::cout << "\n Testing ivp_numerov with F(x) = -4.0, y(x) = cos(2x):\n";
    Eigen::VectorXd ynum_F1D_x = ivp_numerov<double>(func_scal, ycos0, ycos1, x_F1D_x);
    Eigen::VectorXd err_cos_F1D_x = ynum_F1D_x - ycos_F1D_x;

    std::cout << std::setw(15) << "x" << std::setw(15) << "cos(2x)"  << std::setw(15) << "err[cos(2x)]" << "\n";
    for (int i = 0; i < nxI; ++i) {
        std::cout << std::setw(15) << x_F1D_x(i) 
        << std::setw(15) << ycos_F1D_x(i) << std::setw(15) << err_cos_F1D_x(i) << "\n";
    }


    // 测试Numerov矢量方法
    std::cout << "\n Testing ivp_numerov_vec with F(x) = diag(-4.0, 4.0), y(x) = [cos(2x), exp(2x)]:\n";
    Eigen::VectorXd yvec0_F1D_ch(2), yvec1_F1D_ch(2);
    yvec0_F1D_ch << ycos0, yexpp0;
    yvec1_F1D_ch << ycos1, yexpp1;
    Eigen::MatrixXd yvecnum_F2D_ch_x = ivp_numerov_vec<double>(Func_Mat, yvec0_F1D_ch, yvec1_F1D_ch, x_F1D_x);
    Eigen::MatrixXd errvec_F2D_x_ch = yvecnum_F2D_ch_x;
    errvec_F2D_x_ch.row(0) = yvecnum_F2D_ch_x.row(0) - ycos_F1D_x.transpose();
    errvec_F2D_x_ch.row(1) = yvecnum_F2D_ch_x.row(1) - yexpp_F1D_x.transpose();

    std::cout << std::setw(15) << "x" 
    << std::setw(15) << "cos(2x)" << std::setw(15) << "err[cos(2x)]" 
    << std::setw(15) << "exp(2x)" << std::setw(15) << "err[exp(2x)]" << "\n";
    for (int i = 0; i < nxI; ++i) {
        std::cout << std::setw(15) << x_F1D_x(i) 
        << std::setw(15) << ycos_F1D_x(i) << std::setw(15) << errvec_F2D_x_ch(0, i) 
        << std::setw(15) << yexpp_F1D_x(i) << std::setw(15) << errvec_F2D_x_ch(1, i) << "\n";
    }


    // 测试Numerov矩阵方法
    std::cout << "\n Testing ivp_numerov_mat with F(x) = diag(-4.0, 4.0), y(x) = [[cos(2x), sin(2x)], [exp(-2x), exp(2x)]]:\n";
    Eigen::MatrixXd ymat0_F2D_ch_sol(2, 2), ymat1_F2D_ch_sol(2, 2);
    ymat0_F2D_ch_sol << ycos0, ysin0,
                        yexpn0, yexpp0;
    ymat1_F2D_ch_sol << ycos1, ysin1,
                        yexpn1, yexpp1;
    auto ymat_F3D_ch_sol_x = ivp_numerov_mat<double>(Func_Mat, ymat0_F2D_ch_sol, ymat1_F2D_ch_sol, x_F1D_x);
    auto errmat_F3D_ch_sol_x = ymat_F3D_ch_sol_x;
    for (int i = 0; i < nxI; ++i) {
        errmat_F3D_ch_sol_x(0, 0, i) = ymat_F3D_ch_sol_x(0, 0, i) - ycos_F1D_x(i);
        errmat_F3D_ch_sol_x(0, 1, i) = ymat_F3D_ch_sol_x(0, 1, i) - ysin_F1D_x(i);
        errmat_F3D_ch_sol_x(1, 0, i) = ymat_F3D_ch_sol_x(1, 0, i) - yexpn_F1D_x(i);
        errmat_F3D_ch_sol_x(1, 1, i) = ymat_F3D_ch_sol_x(1, 1, i) - yexpp_F1D_x(i);
    }

    std::cout << std::setw(15) << "x" 
    << std::setw(15) << "cos(2x)" << std::setw(15) << "err[cos(2x)]" 
    << std::setw(15) << "sin(2x)" << std::setw(15) << "err[sin(2x)]" 
    << std::setw(15) << "exp(-2x)" << std::setw(15) << "err[exp(-2x)]" 
    << std::setw(15) << "exp(+2x)" << std::setw(15) << "err[exp(+2x)]" << "\n";
    for (int i = 0; i < nxI; ++i) {
        std::cout << std::setw(15) << x_F1D_x(i) 
        << std::setw(15) << ymat_F3D_ch_sol_x(0, 0, i) << std::setw(15) << errmat_F3D_ch_sol_x(0, 0, i) 
        << std::setw(15) << ymat_F3D_ch_sol_x(0, 1, i) << std::setw(15) << errmat_F3D_ch_sol_x(0, 1, i)
        << std::setw(15) << ymat_F3D_ch_sol_x(1, 0, i) << std::setw(15) << errmat_F3D_ch_sol_x(1, 0, i) 
        << std::setw(15) << ymat_F3D_ch_sol_x(1, 1, i) << std::setw(15) << errmat_F3D_ch_sol_x(1, 1, i) << "\n";
    }
}

