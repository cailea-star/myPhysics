/**
 * @file    test_psm_eo.cpp
 * @author  cailea
 * @date    2026-09-28
 * @brief   Print EO PSM results.
 */

#include <utility>

#include "psm_spherical.hpp"
#include "psm_spherical_debug.hpp"

/**
 * @brief Solve Nilsson-BCS and projected configuration mixing.
 * @math (N,Z,ε₂,ε₄,Gₙ,Gₚ) → f,u,v,Ehf,Eqp,Δ.
 * @output HFBCS solutions, PSM energies, and collective amplitudes.
 * @note References: PSM/split/EO/EO_DATA and EO_OUT1.
 * @note Generic configuration cuts; standard-PSM comparison pending.
 */
int main() {
    const int TargetN_I = 96;
    const int TargetZ_I = 69;
    const double epsilon2_F = 0.250;
    const double epsilon4_F = 0.002;
    const double bn_F = 2.27512;
    const double bp_F = 2.40386;

    PSMSphericalSetting psm_setting;
    psm_setting.set_sp_neutron(bn_F, {4, 5, 6});
    psm_setting.set_sp_proton(bp_F, {3, 4, 5});
    psm_setting.set_hfbcs(TargetN_I, TargetZ_I, epsilon2_F, epsilon4_F);
    psm_setting.set_config_cut_neutron((Eigen::Vector2i() << 0, 2).finished(), (Eigen::Vector2d() << 0.0, 3.5).finished(), (Eigen::Vector2i() << 1, 20).finished(), (Eigen::Vector2i() << 0, 24).finished());
    psm_setting.set_config_cut_proton(Eigen::VectorXi::Constant(1, 1), Eigen::VectorXd::Constant(1, 2.5), Eigen::VectorXi::Constant(1, 12), Eigen::VectorXi::Constant(1, 11));

    PSMSpherical psm(std::move(psm_setting));
    psm.solve_hfbcs();
    print_hfbcs(psm);

    // Full signed-K configurations; Nphi = 1 disables PNP.
    // Nalpha = Ngamma = 65 resolves |Kcfg - K| through Imax.
    for (int TargetTwoI_I = 1; TargetTwoI_I <= 59; TargetTwoI_I += 2) {
        std::cout << "\nI = " << 0.5 * TargetTwoI_I << std::endl;
        psm.build_projection(TargetTwoI_I, 65, 40, 65, 1, 1);
        psm.solve_ci(0.16);
        assert(psm.Eci_F1D_eigenH.size() > 0 && psm.Eci_F1D_eigenH.allFinite());
        print_ci(psm);
        std::cout << "[Yrast] I = " << 0.5 * TargetTwoI_I << ", E [MeV] = " << psm.Eci_F1D_eigenH(0) << std::endl;
    }
    return 0;
}
