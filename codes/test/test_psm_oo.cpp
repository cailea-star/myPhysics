/**
 * @file    test_psm_oo.cpp
 * @author  cailea
 * @date    2026-09-28
 * @brief   Print OO PSM results.
 */

#include <utility>
#include <vector>

#include "psm_spherical.hpp"
#include "psm_spherical_debug.hpp"

/**
 * @brief Solve Nilsson-BCS and projected configuration mixing.
 * @math (N,Z,ε₂,ε₄,Gₙ,Gₚ) → f,u,v,Ehf,Eqp,Δ.
 * @output HFBCS solutions, PSM energies, and collective amplitudes.
 * @note References: PSM/split/OO/OO_DATA and OO_OUT1.
 * @note Generic configuration cuts; standard-PSM comparison pending.
 */
int main() {
    const int TargetN_I = 95;
    const int TargetZ_I = 69;
    const double epsilon2_F = 0.245;
    const double epsilon4_F = -0.002;
    const double bn_F = 2.27479;
    const double bp_F = 2.39932;

    const std::vector<int> Nshelln_I1D_Nshell{4, 5, 6};
    const std::vector<int> Nshellp_I1D_Nshell{3, 4, 5};
    const std::vector<int> Nshelln_I1D_Nactive{6};
    const std::vector<int> Nshellp_I1D_Nactive{5};

    const Eigen::VectorXi Ncqpn_I1D_Ncqp = Eigen::VectorXi::Constant(1, 1);
    const Eigen::VectorXd ECutn_F1D_Ncqp = Eigen::VectorXd::Constant(1, 2.5);
    const Eigen::VectorXi NCutn_I1D_Ncqp = Eigen::VectorXi::Constant(1, 8);
    const Eigen::VectorXi TwoKCutn_I1D_Ncqp = Eigen::VectorXi::Constant(1, 13);
    const Eigen::VectorXi Ncqpp_I1D_Ncqp = Eigen::VectorXi::Constant(1, 1);
    const Eigen::VectorXd ECutp_F1D_Ncqp = Eigen::VectorXd::Constant(1, 2.6);
    const Eigen::VectorXi NCutp_I1D_Ncqp = Eigen::VectorXi::Constant(1, 8);
    const Eigen::VectorXi TwoKCutp_I1D_Ncqp = Eigen::VectorXi::Constant(1, 11);

    PSMSphericalSetting psm_setting;
    psm_setting.set_sp_neutron(bn_F, Nshelln_I1D_Nshell, Nshelln_I1D_Nactive);
    psm_setting.set_sp_proton(bp_F, Nshellp_I1D_Nshell, Nshellp_I1D_Nactive);
    psm_setting.set_hfbcs(TargetN_I, TargetZ_I, epsilon2_F, epsilon4_F);
    psm_setting.set_config_cut_neutron(Ncqpn_I1D_Ncqp, ECutn_F1D_Ncqp, NCutn_I1D_Ncqp, TwoKCutn_I1D_Ncqp);
    psm_setting.set_config_cut_proton(Ncqpp_I1D_Ncqp, ECutp_F1D_Ncqp, NCutp_I1D_Ncqp, TwoKCutp_I1D_Ncqp);

    PSMSpherical psm(std::move(psm_setting));
    psm.solve_hfbcs();
    print_hfbcs(psm);

    // Full signed-K configurations; Nphi = 1 disables PNP.
    // Axial K selection: Nalpha = Ngamma = 1.
    const std::vector<int> TargetTwoI_I1D_twoI{0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50};
    psm.build_projection(40, 1, 1);
    psm.build_ci(0.16);
    for (const int TargetTwoI_I : TargetTwoI_I1D_twoI) {
        psm.solve_ci(TargetTwoI_I);
        print_ci(psm, TargetTwoI_I);
    }
    return 0;
}
