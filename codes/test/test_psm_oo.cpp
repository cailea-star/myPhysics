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
 * @note Standard PSM nuclear configuration selection; PNP disabled.
 */
int main() {
    const int TargetN_I = 95;
    const int TargetZ_I = 69;
    const double epsilon2_F = 0.245;
    const double epsilon4_F = -0.002;
    const double bn_F = 2.2747931109172965;
    const double bp_F = 2.3993171453537059;

    const std::vector<int> Nshelln_I1D_Nshell{4, 5, 6};
    const std::vector<int> Nshellp_I1D_Nshell{3, 4, 5};
    const std::vector<int> Nshelln_I1D_Nactive{6};
    const std::vector<int> Nshellp_I1D_Nactive{5};

    Eigen::MatrixXi Ncqp_I2D_Ncqp_np(3, 2);
    Ncqp_I2D_Ncqp_np.row(0) << 1, 0;
    Ncqp_I2D_Ncqp_np.row(1) << 0, 1;
    Ncqp_I2D_Ncqp_np.row(2) << 1, 1;
    const Eigen::VectorXd ECut_F1D_Ncqp = (Eigen::VectorXd(3) << 2.5, 2.6, 3.3).finished();
    const Eigen::VectorXi NCut_I1D_Ncqp = (Eigen::VectorXi(3) << 8, 8, 32).finished();
    const Eigen::VectorXi TwoKCut_I1D_Ncqp = (Eigen::VectorXi(3) << 13, 11, 24).finished();

    PSMSphericalSetting psm_setting;
    psm_setting.set_sp_neutron(bn_F, Nshelln_I1D_Nshell, Nshelln_I1D_Nactive);
    psm_setting.set_sp_proton(bp_F, Nshellp_I1D_Nshell, Nshellp_I1D_Nactive);
    psm_setting.set_hfbcs(TargetN_I, TargetZ_I, epsilon2_F, epsilon4_F);
    psm_setting.set_ci_cut(Ncqp_I2D_Ncqp_np, ECut_F1D_Ncqp, NCut_I1D_Ncqp, TwoKCut_I1D_Ncqp);

    PSMSpherical psm(std::move(psm_setting));
    psm.solve_hfbcs();
    print_hfbcs(psm);

    // Standard Kramers configurations; Nphi = 1 disables PNP.
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
