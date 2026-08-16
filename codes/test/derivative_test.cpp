#include <ios>
#include <cassert>
#include <cmath>
#include <complex>
#include <iomanip>
#include <ostream>
#include <iostream>
#include <Eigen/Dense>
#include "derivative.hpp"

double V0_F = 1.0;
double a0_F = 2.0;
double R0_F = 6.0;

/*
* Woods-Sakon potential and its derivatives
* math: WS(x) = V0 / (1 + exp((x - R0) / a0))
*/
double WS_func(double x_F) {
    double exp_F = std::exp((x_F - R0_F) / a0_F);
    return V0_F / (1 + exp_F);
}

/**
 * Derivative of Woods-Sakon potential
 * math: dWS(x) = V0 * (-1) / (1 + exp)^2 * exp * (1/a0)
 */
double dWS_func(double x_F) {
    double exp_F = std::exp((x_F - R0_F) / a0_F);
    return V0_F * (-1) * std::pow(1 + exp_F, -2) * exp_F * (1 / a0_F);
}


/**
 * Second derivative of Woods-Sakon potential
 * math: ddWS(x) = V0 * [ (-2) * (-1) / (1 + exp)^3 * exp^2 * (1/a0)^2 + (-1) / (1 + exp)^2 * exp * (1/a0)^2 ]
 */
double ddWS_func(double x_F) {
    double exp_F = std::exp((x_F - R0_F)/a0_F);
    double term1_F = V0_F * (-2) * (-1) * std::pow(1 + exp_F, -3) * std::pow(exp_F, 2) * std::pow(1 / a0_F, 2);
    double term2_F = V0_F * (-1) * std::pow(1 + exp_F, -2) * exp_F * std::pow(1 / a0_F, 2);
    return term1_F + term2_F;
}

int main() {
    int nx_I = 100;
    double xmin_F = 1e-10;
    double xmax_F = 15.0;
    Eigen::VectorXd x_F1D_x = Eigen::VectorXd::LinSpaced(nx_I, xmin_F, xmax_F);
    Eigen::VectorXd y_F1D_x(nx_I);
    for (int i = 0; i < nx_I; ++i) {
        y_F1D_x(i) = WS_func(x_F1D_x(i));
    }
    Eigen::VectorXd dydx_F1D_x = derivative1<double>(y_F1D_x, x_F1D_x);
    Eigen::VectorXd d2ydx2_F1D_x = derivative2<double>(y_F1D_x, x_F1D_x);
    std::cout << std::scientific << std::setprecision(2);
    std::cout << std::setw(15) << "x";
    std::cout << std::setw(15) << "WS(x)";
    std::cout << std::setw(15) << "dWX(x)";
    std::cout << std::setw(15) << "err(dWS)";
    std::cout << std::setw(15) << "err(dWSd)";
    std::cout << std::setw(15) << "ddWX(x)";
    std::cout << std::setw(15) << "err(ddWS)";
    std::cout << std::setw(15) << "err(ddWSd)";
    std::cout << std::endl;
    for (int i = 0; i < nx_I; ++i) {
        double x = x_F1D_x(i);
        std::cout << std::setw(15) << x;
        std::cout << std::setw(15) << WS_func(x);
        std::cout << std::setw(15) << dWS_func(x);
        std::cout << std::setw(15) << derivative1(WS_func, x) - dWS_func(x);
        std::cout << std::setw(15) << dydx_F1D_x(i) - dWS_func(x);
        std::cout << std::setw(15) << ddWS_func(x);
        std::cout << std::setw(15) << derivative2(WS_func, x) - ddWS_func(x);
        std::cout << std::setw(15) << d2ydx2_F1D_x(i) - ddWS_func(x);
        std::cout << std::endl;
    }

    int Nx_check_I = 7;
    Eigen::VectorXd x_check_F1D_x = Eigen::VectorXd::LinSpaced(Nx_check_I, -1.0, 1.0);
    Eigen::VectorXcd f_check_C1D_x(Nx_check_I);
    std::complex<double> coeff_C(1.0, 1.0);
    for (int x_I = 0; x_I < Nx_check_I; ++x_I) {
        f_check_C1D_x(x_I) = coeff_C * x_check_F1D_x(x_I) * x_check_F1D_x(x_I);
    }
    Eigen::VectorXcd d2f_check_C1D_x = derivative2(f_check_C1D_x, x_check_F1D_x);
    std::complex<double> d2f_expected_C = 2.0 * coeff_C;
    for (int x_I = 0; x_I < Nx_check_I; ++x_I) {
        assert(std::abs(d2f_check_C1D_x(x_I) - d2f_expected_C) < 1.0e-12);
        assert(std::abs(derivative2(f_check_C1D_x, x_check_F1D_x, x_I) - d2f_expected_C) < 1.0e-12);
    }
}
