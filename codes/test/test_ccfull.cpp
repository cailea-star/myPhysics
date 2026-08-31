/**
 * @file    test_ccfull.cpp
 * @author  cailea
 * @date    2026-08-25
 * @brief   Coupled-channel fusion regression test.
 */

#include <Eigen/Core>

#include <array>
#include <cassert>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <tuple>

#include "basis_bspline.hpp"
#include "cc_fusion.hpp"
#include "cc_params.hpp"
#include "cc_potentials.hpp"

struct FusionReferencePoint {
    double E_F;
    double sigma_F;
    double Lavg_F;
    double P0_F;
};

/**
 * @brief  Sum grid-propagated fusion partial waves.
 * @math   (σ,⟨L⟩,P₀)=Σ_L Numerov(F_L)
 * @output Grid-propagated fusion observables.
 */
std::tuple<double, double, double> fusion_cross_section_grid(const Eigen::Ref<const Eigen::VectorXd>& r_F1D_r, const CCParams& params, const SphericalRadialBoundaryFunc& boundary_Func = fusion_inner_boundary, bool useNumerov_B = true) {
    double E_F = params.channel_1D_ch[0].Ech_F;
    return fusion_cross_section(params, [&](int L_I) {
        auto [rin_F, rwell_F, rmid_F, rbarrier_F, rout_F] = params.find_barrier(L_I, E_F);
        (void) rin_F;
        (void) rmid_F;
        (void) rbarrier_F;
        (void) rout_F;
        if (rwell_F < 0.0 || params.calc_Vall_scalar(rwell_F) + calc_Vcent(rwell_F, params.Varg_params, L_I) > E_F) {return 0.0;}
        Eigen::Tensor<std::complex<double>, 3, Eigen::ColMajor> u_C3D_ch_sol_r = useNumerov_B ? fusion_wavefunction_numerov(r_F1D_r, params, L_I, boundary_Func) : fusion_wavefunction_rk4(r_F1D_r, params, L_I, boundary_Func);
        return fusion_probability_u(r_F1D_r, params, L_I, u_C3D_ch_sol_r)(0);
    });
}

/**
 * @brief  Sum B-spline R-matrix fusion partial waves.
 * @math   (σ,⟨L⟩,P₀)=Σ_L R_L
 * @output R-matrix fusion observables.
 */
std::tuple<double, double, double> fusion_cross_section_rmatrix(RMatrix& r_matrix, const CCParams& params) {
    double E_F = params.channel_1D_ch[0].Ech_F;
    return fusion_cross_section(params, [&](int L_I) {
        auto [rin_F, rwell_F, rmid_F, rbarrier_F, rout_F] = params.find_barrier(L_I, E_F);
        (void) rin_F;
        (void) rmid_F;
        (void) rbarrier_F;
        (void) rout_F;
        if (rwell_F < 0.0 || params.calc_Vall_scalar(rwell_F) + calc_Vcent(rwell_F, params.Varg_params, L_I) > E_F) {return 0.0;}
        return fusion_probability_R(params, L_I, r_matrix)(0);
    });
}

/**
 * @brief  Validate grid and R-matrix fusion observables.
 * @math   |O-O_ref|<ε
 * @output Process status.
 */
int main() {
    constexpr double dr_F = 0.05;
    constexpr std::array<FusionReferencePoint, 5> reference_1D_point{{
        {58.0, 1.1329446309, 6.35460, 0.0239116754},
        {60.0, 15.8106006497, 6.99390, 0.2605522539},
        {62.0, 69.9070152939, 9.21422, 0.5962191065},
        {64.0, 156.6168982566, 11.31130, 0.9172747590},
        {66.0, 257.7196558111, 13.44984, 0.9886207192},
    }};
    auto make_params_Func = [](double E_F) {
        VARG Varg_params(144, 62, 16, 8);
        Varg_params.V0_F = 105.1;
        Varg_params.a0_F = 0.75;
        Varg_params.Rp0_F = 1.1 * std::cbrt(16.0);
        Varg_params.Rt0_F = 1.1 * std::cbrt(144.0);
        double Rt_reference_F = 1.06 * std::cbrt(144.0);
        Varg_params.Rp0_F += Varg_params.Rt0_F - Rt_reference_F;
        Varg_params.Rt0_F = Rt_reference_F;
        Varg_params.Rcoul_p_F = 0.0;
        Varg_params.Rcoul_t_F = Rt_reference_F;
        Varg_params.beta2t_F = 0.205;
        Varg_params.beta4t_F = 0.0;
        Varg_params.beta2tcoul_F = Varg_params.beta2t_F + 2.0 * std::sqrt(5.0 / std::numbers::pi) * Varg_params.beta2t_F * Varg_params.beta2t_F / 7.0;
        Varg_params.beta4tcoul_F = Varg_params.beta4t_F + 9.0 * Varg_params.beta2t_F * Varg_params.beta2t_F / (7.0 * std::sqrt(std::numbers::pi));
        return CCParams(144, 62, 16, 8, 4, Varg_params, {Channel(0, 0, 0, 0, E_F), Channel(0, 4, 0, 4, E_F - 1.81), Channel(0, 8, 0, 8, E_F - 1.81 * 20.0 / 6.0)});
    };
    CCParams params = make_params_Func(62.0);
    assert(params.is_valid());

    auto [rin_F, rwell_F, rmid_F, rbarrier_F, rout_F] = params.find_barrier(0, 0.0);
    (void) rin_F;
    (void) rmid_F;
    (void) rout_F;
    assert(std::abs(rbarrier_F - 10.8175) < 1.0e-4);
    int Nrstep_I = static_cast<int>(std::floor((30.0 - rwell_F) / dr_F));
    Eigen::VectorXd r_F1D_r = Eigen::VectorXd::LinSpaced(Nrstep_I + 2, rwell_F, rwell_F + dr_F * (Nrstep_I + 1));
    assert(r_F1D_r.size() > 3);

    params.print_info();
    std::cout << std::right << std::fixed << std::setprecision(10);
    std::cout << std::setw(12) << "E" << " | " << std::setw(16) << "sigma" << " | " << std::setw(16) << "err[sigma]" << " | " << std::setw(14) << "<L>" << " | " << std::setw(16) << "err[<L>]" << " | " << std::setw(14) << "P(L=0)" << " | " << std::setw(16) << "err[P(L=0)]" << '\n';
    for (const FusionReferencePoint& reference_point : reference_1D_point) {
        CCParams energy_params = make_params_Func(reference_point.E_F);
        auto [sigma_F, Lavg_F, P0_F] = fusion_cross_section_grid(r_F1D_r, energy_params);
        std::cout << std::setw(12) << reference_point.E_F << " | " << std::setw(16) << sigma_F << " | " << std::showpos << std::setw(16) << sigma_F - reference_point.sigma_F << std::noshowpos << " | " << std::setw(14) << Lavg_F << " | " << std::showpos << std::setw(16) << Lavg_F - reference_point.Lavg_F << std::noshowpos << " | " << std::setw(14) << P0_F << " | " << std::showpos << std::setw(16) << P0_F - reference_point.P0_F << std::noshowpos << '\n';
    }

    Eigen::VectorXd rR_F1D_r = Eigen::VectorXd::LinSpaced(81, r_F1D_r(0), r_F1D_r(r_F1D_r.size() - 1));
    BSplineBasisFunction b_basis_func(rR_F1D_r);
    BSplineBasis b_basis(b_basis_func, 5, false);
    RMatrix r_matrix(b_basis, static_cast<int>(params.channel_1D_ch.size()));
    for (const FusionReferencePoint& reference_point : reference_1D_point) {
        auto [sigma_F, Lavg_F, P0_F] = fusion_cross_section_rmatrix(r_matrix, make_params_Func(reference_point.E_F));
        std::cout << "R " << std::setw(9) << reference_point.E_F << " | " << std::setw(16) << sigma_F << " | " << std::showpos << std::setw(16) << sigma_F - reference_point.sigma_F << std::noshowpos << " | " << std::setw(14) << Lavg_F << " | " << std::showpos << std::setw(16) << Lavg_F - reference_point.Lavg_F << std::noshowpos << " | " << std::setw(14) << P0_F << " | " << std::showpos << std::setw(16) << P0_F - reference_point.P0_F << std::noshowpos << '\n';
        assert(std::abs(sigma_F - reference_point.sigma_F) < 1.0e-2);
        assert(std::abs(Lavg_F - reference_point.Lavg_F) < 1.0e-3);
        assert(std::abs(P0_F - reference_point.P0_F) < 2.0e-4);
    }
}
