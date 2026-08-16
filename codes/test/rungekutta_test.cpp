
#include <iostream>
#include <iomanip>
#include <Eigen/Core>
#include "ivp_rungekutta.hpp"


/**
 * y' = f = cos(x), y = sin(x)
 */
double sinfFunc(double x, double y) {
    return std::cos(x);
}

/**
 * [y1, y2]' = [y2, y1]
 * y1'' = y1, y1 = A1 exp(x) + B1 exp(-x)
 * y2'' = y2, y2 = A2 exp(x) + B2 exp(-x)
 */
Eigen::VectorXd expfFunc(double x, const Eigen::VectorXd& y) {
    Eigen::VectorXd dydx(2);
    dydx(0) = y(1); 
    dydx(1) = y(0); 
    return dydx;
}

int main() {
    int nxI = 100;
    double xmin_F = 0.0;
    double xmax_F = 10.0;
    Eigen::VectorXd x_F1D_x = Eigen::VectorXd::LinSpaced(nxI, xmin_F, xmax_F);
    double ysin0 = std::sin(x_F1D_x(0)); 
    double yexpp0 = std::exp(x_F1D_x(0)); 
    double yexpn0 = std::exp(-x_F1D_x(0)); 
    Eigen::VectorXd yexppvec0(2);
    yexppvec0(0) = yexpp0;
    yexppvec0(1) = yexpp0;
    Eigen::VectorXd yexpnvec0(2);
    yexpnvec0(0) = yexpn0;
    yexpnvec0(1) = - yexpn0;
    Eigen::VectorXd sin_F1D_x = x_F1D_x.array().sin();
    Eigen::VectorXd expp_F1D_x = x_F1D_x.array().exp();
    Eigen::VectorXd expn_F1D_x = (-x_F1D_x).array().exp();

    std::cout << std::scientific << std::setprecision(4) << std::right;

    // 测试标量求解器
    std::cout << "\n Testing ivp_rungekutta with f(x, y) = cos(x), y(x) = sin(x):\n";
    Eigen::VectorXd ysin_F1D_x = ivp_rk4<double>(sinfFunc, ysin0, x_F1D_x);
    Eigen::VectorXd err_sin_F1D_x = ysin_F1D_x - sin_F1D_x;
    std::cout << std::setw(15) << "x" 
    << std::setw(15) << "sin(x)"  << std::setw(15) << "err[sin(x)]" << "\n";
    for (int i = 0; i < nxI; ++i) {
        std::cout << std::setw(15) << x_F1D_x(i) 
        << std::setw(15) << sin_F1D_x(i) << std::setw(15) << err_sin_F1D_x(i) << "\n";
    }

    // 测试向量求解器
    std::cout << "\n Testing ivp_rk4_vec with [y1, y2]' = [y2, y1], y1(x) = A1 exp(x) + B1 exp(-x), y2(x) = A2 exp(x) + B2 exp(-x):\n";
    Eigen::MatrixXd yexpp_F2D_ch_x = ivp_rk4_vec<double>(expfFunc, yexppvec0, x_F1D_x);
    Eigen::MatrixXd yexpn_F2D_ch_x = ivp_rk4_vec<double>(expfFunc, yexpnvec0, x_F1D_x);
    Eigen::VectorXd err_expp_F1D_x = yexpp_F2D_ch_x.row(0) - expp_F1D_x.transpose();
    Eigen::VectorXd err_expn_F1D_x = yexpn_F2D_ch_x.row(0) - expn_F1D_x.transpose();
    std::cout << std::setw(15) << "x" 
    << std::setw(15) << "exp(x)"  << std::setw(15) << "err[exp(x)]"
    << std::setw(15) << "exp(-x)"  << std::setw(15) << "err[exp(-x)]" << "\n";
    for (int i = 0; i < nxI; ++i) {
        std::cout << std::setw(15) << x_F1D_x(i) 
        << std::setw(15) << expp_F1D_x(i) << std::setw(15) << err_expp_F1D_x(i)
        << std::setw(15) << expn_F1D_x(i) << std::setw(15) << err_expn_F1D_x(i) << "\n";
    }
    return 0;
}
