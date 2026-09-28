/**
 * @file    test_psm_oe.cpp
 * @author  cailea
 * @date    2026-09-28
 * @brief   Print OE Nilsson and BCS solutions.
 */

#include "psm_spherical.hpp"
#include "psm_spherical_debug.hpp"

/**
 * @brief Solve and print Nilsson-BCS results.
 * @math (N,Z,ε₂,ε₄,Gₙ,Gₚ) → f,u,v,Ehf,Eqp,Δ.
 * @output Neutron and proton HFBCS solutions.
 * @note References: PSM/split/OE/OE_DATA and OE_OUT1.
 * @note Unblocked BCS vacuum; no projected configuration mixing.
 */
int main() {
    const int TargetN_I = 95;
    const int TargetZ_I = 70;
    const double epsilon2_F = 0.242;
    const double epsilon4_F = 0.000;
    const double bn_F = 2.27963;
    const double bp_F = 2.39866;

    // G_q = [20.12 − 13.13(N_q−N_other)/A]/A.
    const double A_F = TargetN_I + TargetZ_I;
    const double Gn_F = (20.12 - 13.13 * (TargetN_I - TargetZ_I) / A_F) / A_F;
    const double Gp_F = (20.12 - 13.13 * (TargetZ_I - TargetN_I) / A_F) / A_F;

    PSMSpherical psm(bn_F, bp_F, {4, 5, 6}, {3, 4, 5});
    psm.build_hfbcs(TargetN_I, TargetZ_I, Gn_F, Gp_F, epsilon2_F, epsilon4_F);
    psm_debug::print_hfbcs(psm);
    return 0;
}
