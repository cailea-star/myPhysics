#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include "b_spline_basis.hpp"


double harmonic_oscillator_energy(int n_I) {
    return static_cast<double>(n_I) + 0.5;
}

double harmonic_oscillator_wavefunction(int n_I, double x_F) {
    double factorial_F = 1.0;
    for (int i_I = 2; i_I <= n_I; ++i_I) factorial_F *= static_cast<double>(i_I);

    double Hprev_F = 1.0;
    double Hcurr_F = 2.0 * x_F;
    double Hn_F = (n_I == 0) ? Hprev_F : Hcurr_F;
    for (int n_rec_I = 2; n_rec_I <= n_I; ++n_rec_I) {
        Hn_F = 2.0 * x_F * Hcurr_F - 2.0 * static_cast<double>(n_rec_I - 1) * Hprev_F;
        Hprev_F = Hcurr_F;
        Hcurr_F = Hn_F;
    }

    double norm_F = 1.0 / std::sqrt(std::pow(2.0, n_I) * factorial_F * std::sqrt(M_PI));
    return norm_F * Hn_F * std::exp(-0.5 * x_F * x_F);
}

int main() {
    int Nx_I = 121;
    int Nquad_I = 5;
    int Ncheck_I = 4;
    double xmin_F = -8.0;
    double xmax_F = 8.0;

    Eigen::VectorXd x_F1D_x = Eigen::VectorXd::LinSpaced(Nx_I, xmin_F, xmax_F);
    BSplineBasisFunction b_basis_func(x_F1D_x);
    const bool isZeroBound_B = true;
    BSplineBasis b_basis(b_basis_func, Nquad_I, isZeroBound_B);
    BSplineBasis b_basis_full(b_basis_func, Nquad_I, false);
    const Eigen::VectorXd& x_F1D_g = b_basis.x_F1D_grid;
    const Eigen::VectorXd& w_F1D_g = b_basis.w_F1D_grid;
    int Nb_I = b_basis_func.size();
    int Nb_inner_I = Nb_I - 2;
    if (b_basis.B_F2D_grid_b.cols() != Nb_inner_I || b_basis_full.B_F2D_grid_b.cols() != Nb_I) return 1;
    
    Eigen::VectorXd one_F1D_g = Eigen::VectorXd::Ones(x_F1D_g.size());
    Eigen::VectorXd half_F1D_g = Eigen::VectorXd::Constant(x_F1D_g.size(), 0.5);
    Eigen::VectorXd V_F1D_g = (0.5 * x_F1D_g.array().square()).matrix();

    Eigen::MatrixXd S_F2D_i_j = b_basis.calc_B_O_B(one_F1D_g);
    Eigen::MatrixXd V_F2D_i_j = b_basis.calc_B_O_B(V_F1D_g);
    Eigen::MatrixXd T_F2D_i_j = b_basis.calc_dB_O_dB(half_F1D_g);
    Eigen::MatrixXd H_F2D_i_j = T_F2D_i_j + V_F2D_i_j;
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> solver(H_F2D_i_j, S_F2D_i_j);
    if (solver.info() != Eigen::Success) return 1;

    std::cout << std::scientific << std::setprecision(6) << std::left;
    std::cout << "\n Testing BSplineBasis1D with V(x) = x^2 / 2:\n";
    std::cout << std::setw(10) << "n";
    std::cout << std::setw(18) << "E_BSpline";
    std::cout << std::setw(18) << "E_exact";
    std::cout << std::setw(18) << "error" << "\n";

    double maxErr_F = 0.0;
    for (int n_I = 0; n_I < Ncheck_I; ++n_I) {
        double E_F = solver.eigenvalues()(n_I);
        double Eexact_F = harmonic_oscillator_energy(n_I);
        double err_F = std::abs(E_F - Eexact_F);
        maxErr_F = std::max(maxErr_F, err_F);
        std::cout << std::setw(10) << n_I;
        std::cout << std::setw(18) << E_F;
        std::cout << std::setw(18) << Eexact_F;
        std::cout << std::setw(18) << err_F << "\n";
    }

    std::cout << "\n Wave Function Overlaps:\n";
    std::cout << std::setw(10) << "n";
    std::cout << std::setw(18) << "|overlap|";
    std::cout << std::setw(18) << "1-|overlap|" << "\n";

    double maxWaveErr_F = 0.0;
    for (int n_I = 0; n_I < Ncheck_I; ++n_I) {
        Eigen::VectorXd coeff_F1D_b = Eigen::VectorXd::Zero(Nb_I);
        coeff_F1D_b.segment(1, Nb_inner_I) = solver.eigenvectors().col(n_I);
        double overlap_F = 0.0;
        double normNum_F = 0.0;
        double normExact_F = 0.0;
        for (int g_I = 0; g_I < x_F1D_g.size(); ++g_I) {
            double psi_F = b_basis.b_funcs.WaveFunction(x_F1D_g(g_I), coeff_F1D_b)(0);
            double exact_F = harmonic_oscillator_wavefunction(n_I, x_F1D_g(g_I));
            overlap_F += w_F1D_g(g_I) * psi_F * exact_F;
            normNum_F += w_F1D_g(g_I) * psi_F * psi_F;
            normExact_F += w_F1D_g(g_I) * exact_F * exact_F;
        }

        double overlapAbs_F = std::abs(overlap_F) / std::sqrt(normNum_F * normExact_F);
        double waveErr_F = 1.0 - overlapAbs_F;
        maxWaveErr_F = std::max(maxWaveErr_F, waveErr_F);
        std::cout << std::setw(10) << n_I;
        std::cout << std::setw(18) << overlapAbs_F;
        std::cout << std::setw(18) << waveErr_F << "\n";
    }

    if (maxErr_F > 5.0e-4 || maxWaveErr_F > 5.0e-4) return 1;
    return 0;
}
