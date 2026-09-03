/**
 * @file    hfb_axial_initialize.cpp
 * @author  cailea
 * @date    2026-05-05
 * @brief   Initialize axial HFB Woods-Saxon fields.
 */

#include "hfb_axial.hpp"

#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>

namespace {

struct WSGeometry {
    double distance_F = 0.0;
    double surface_scale_F = 0.0;
};

struct WSShapeGeometry {
    double b2_ws_F = 0.0;
    double b3_ws_F = 0.0;
    double b4_ws_F = 0.0;
    double surface_scale0_F = 0.0;
    double zcm_F = 0.0;

    /**
     * @brief  Calculate volume-preserving Woods-Saxon geometry.
     * @math   (β_2,β_3,β_4,b_0,A) → (b_2,b_3,b_4,f_0,z_{cm})
     * @output Woods-Saxon shape geometry.
     */
    static WSShapeGeometry from_beta(double beta2_F, double beta3_F, double beta4_F, double b0_basis_F, double hbzero_F, double Atarget_F) {
        const double pi_F = std::acos(-1.0);
        const double b2_F = beta2_F * std::sqrt(5.0 / (4.0 * pi_F));
        const double b3_F = beta3_F * std::sqrt(7.0 / (4.0 * pi_F));
        const double b4_F = beta4_F * std::sqrt(9.0 / (4.0 * pi_F));
        const double b2_squared_F = b2_F * b2_F;
        const double b2_cubed_F = b2_squared_F * b2_F;
        const double b3_squared_F = b3_F * b3_F;
        const double b3_cubed_F = b3_squared_F * b3_F;
        const double b4_squared_F = b4_F * b4_F;
        const double b4_cubed_F = b4_squared_F * b4_F;

        // (b_2,b_3,b_4) → f_V.
        const double volume_factor_F = 2.0 + (2.0 * b4_squared_F) / 3.0 + (40.0 * b2_F * b4_squared_F) / 231.0 + (36.0 * b4_cubed_F) / 1001.0 + (6.0 * b2_squared_F) / 5.0 + (4.0 * b2_cubed_F) / 35.0 + (12.0 * b2_squared_F * b4_F) / 35.0 + (6.0 * b3_squared_F) / 7.0 + (8.0 * b2_F * b3_squared_F) / 35.0 + (12.0 * b3_squared_F * b4_F) / 77.0;
        const double surface_scale0_F = std::cbrt(2.0 / volume_factor_F);

        // (b_2,b_3,b_4) → P_{cm}.
        const double zcm_polynomial_F = (72.0 * b2_F * b3_F) / 35.0 + (48.0 * b2_squared_F * b3_F) / 35.0 + (24.0 * b2_cubed_F * b3_F) / 55.0 + (96.0 * b3_cubed_F) / 385.0 + (1656.0 * b2_F * b3_cubed_F) / 5005.0 + (32.0 * b3_F * b4_F) / 21.0 + (2272.0 * b2_F * b3_F * b4_F) / 1155.0 + (1856.0 * b2_squared_F * b3_F * b4_F) / 2145.0 + (3904.0 * b3_cubed_F * b4_F) / 15015.0 + (544.0 * b3_F * b4_squared_F) / 1001.0 + (488.0 * b2_F * b3_F * b4_squared_F) / 715.0 + (512.0 * b3_F * b4_cubed_F) / 2431.0;
        const double r0_basis_F = (2.0 * hbzero_F / (41.0 * b0_basis_F * b0_basis_F)) * std::cbrt(Atarget_F);
        const double zcm_F = 3.0 * std::pow(1.0 / volume_factor_F, 4.0 / 3.0) * zcm_polynomial_F * r0_basis_F / (2.0 * std::cbrt(2.0));
        return {b2_F, b3_F, b4_F, surface_scale0_F, zcm_F};
    }

    /**
     * @brief  Evaluate Woods-Saxon geometry at one point.
     * @math   (r_⊥²,z) → (r,f)
     * @output Radius and deformed surface scale.
     */
    WSGeometry calc_geometry(double r2_F, double z_F) const {
        const double z_shift_F = z_F - zcm_F;
        const double z2_F = z_shift_F * z_shift_F;
        const double distance2_F = r2_F + z2_F;
        assert(distance2_F > 0.0);
        const double cos2_F = z2_F / distance2_F;
        const double abs_cos_F = std::sqrt(cos2_F);
        const double P2_F = 0.5 * (3.0 * cos2_F - 1.0);
        const double P3_F = 0.5 * (5.0 * abs_cos_F * abs_cos_F * abs_cos_F - 3.0 * abs_cos_F);
        const double P4_F = (35.0 * cos2_F * cos2_F - 30.0 * cos2_F + 3.0) / 8.0;
        const double surface_scale_F = surface_scale0_F * (1.0 + b2_ws_F * P2_F + b3_ws_F * P3_F + b4_ws_F * P4_F);
        return {std::sqrt(distance2_F), surface_scale_F};
    }
};

} // namespace

void AxialHFB::initialize_WS_field(int Ntarget_I, int Ztarget_I, double beta2_F, double beta3_F, double beta4_F) {
    assert(Ntarget_I >= 0);
    assert(Ztarget_I >= 0);
    assert(Ntarget_I + Ztarget_I > 0);

    // P_{WS} = (V_0,r_0,a_0,V_{LS},r_{LS},a_{LS},κ_V).
    const double V0WS_F = -71.28;
    const double r0WS_F = 1.2334;
    const double a0WS_F = 0.6150;
    const double V0LS_F = 11.1175;
    const double r0LS_F = 1.1443;
    const double a0LS_F = 0.6476;
    const double akv_F = 0.4616;

    // (ℏc,m_u) in MeV-fm units.
    const double hbarc_F = 197.32891;
    const double amu_F = 931.494013;

    // (N,Z) → (A,R_{WS},R_{LS}).
    const int Nz_I = axialconfig.Nz_I;
    const int Nr_I = axialconfig.Nr_I;
    const int Atarget_I = Ntarget_I + Ztarget_I;
    const double Atarget_F = static_cast<double>(Atarget_I);
    const EDFParamsSkyrme active_edf_ = hfbsettings.make_active_edf(edf_skyrme, Atarget_I);
    const double R0WS_F = r0WS_F * std::cbrt(Atarget_F);
    const double R0LS_F = r0LS_F * std::cbrt(Atarget_F);

    // (b_z,b_r,β_2,β_3,β_4) → G_{WS}.
    const double b0_basis_F = AxialConfig::bzbr_to_b0beta20(axialconfig.bz_F, axialconfig.br_F).first;
    const WSShapeGeometry ws_shape_ = WSShapeGeometry::from_beta(beta2_F, beta3_F, beta4_F, b0_basis_F, edf_skyrme.hbzero_F, Atarget_F);

    // I_n=(N-Z)/A; I_p=(Z-N)/A.
    const double asymmetry_n_F = static_cast<double>(Ntarget_I - Ztarget_I) / Atarget_F;
    const double asymmetry_p_F = static_cast<double>(Ztarget_I - Ntarget_I) / Atarget_F;
    const double V0WS_n_F = V0WS_F * (1.0 - akv_F * asymmetry_n_F);
    const double V0WS_p_F = V0WS_F * (1.0 - akv_F * asymmetry_p_F);
    const double VLS_prefactor_F = 0.5 * std::pow(hbarc_F / amu_F, 2.0);
    const double V0LS_n_F = VLS_prefactor_F * V0WS_n_F * V0LS_F;
    const double V0LS_p_F = VLS_prefactor_F * V0WS_p_F * V0LS_F;

    const auto& z_F1D_z = global_basis.z_F1D_z;
    const auto& r_F1D_r = global_basis.r_F1D_r;

    field_n.set_zero();
    field_p.set_zero();

    // (V_0,G_{WS}) → (v_{cent},v_{mass},v_{∇J},Δ).
    const auto initialize_field_Func = [&](AxialHFBField& field_, double V0WS_q_F, double hbzero_q_F, double V0LS_q_F) {
        for (int r_I = 0; r_I < Nr_I; ++r_I) {
            for (int z_I = 0; z_I < Nz_I; ++z_I) {
                const double r_F = r_F1D_r(r_I);
                const double r2_F = r_F * r_F;
                const double z_F = z_F1D_z(z_I);
                const WSGeometry geometry_ = ws_shape_.calc_geometry(r2_F, z_F);
                const double vcent_F = V0WS_q_F / (1.0 + std::exp((geometry_.distance_F - R0WS_F * geometry_.surface_scale_F) / a0WS_F));
                const double vdJ_F = -V0LS_q_F / (1.0 + std::exp((geometry_.distance_F - R0LS_F * geometry_.surface_scale_F) / a0LS_F));
                const double kappa_aux_F = 5.0e-3 * std::exp((geometry_.distance_F - R0WS_F * geometry_.surface_scale_F) / 2.0);
                field_.vmass_F2D_z_r(z_I, r_I) = hbzero_q_F;
                field_.vcent_F2D_z_r(z_I, r_I) = vcent_F;
                field_.vdJ_F2D_z_r(z_I, r_I) = vdJ_F;
                field_.vpair_F2D_z_r(z_I, r_I) = -100.0 * kappa_aux_F;
            }
        }
    };
    initialize_field_Func(field_n, V0WS_n_F, active_edf_.hbzeron_F, V0LS_n_F);
    initialize_field_Func(field_p, V0WS_p_F, active_edf_.hbzerop_F, V0LS_p_F);

    // P_{WS} → stdout.
    std::cout << "[AxialHFB::initialize_WS_field] Woods-Saxon initialization\n";
    std::cout << std::scientific << std::setprecision(3) << std::right;
    std::cout << std::setw(20) << "[WS] params:" << std::setw(10) << " V0ws =" << std::setw(10) << V0WS_F << std::setw(10) << " r0ws =" << std::setw(10) << r0WS_F << std::setw(10) << " a0ws =" << std::setw(10) << a0WS_F << "\n";
    std::cout << std::setw(20) << "[LS] params:" << std::setw(10) << " V0ls =" << std::setw(10) << V0LS_F << std::setw(10) << " r0ls =" << std::setw(10) << r0LS_F << std::setw(10) << " a0ls =" << std::setw(10) << a0LS_F << "\n";
    std::cout << std::setw(20) << "[Nucleus] params:" << std::setw(10) << " A =" << std::setw(10) << Atarget_I << std::setw(10) << " Z =" << std::setw(10) << Ztarget_I << std::setw(10) << " N =" << std::setw(10) << Ntarget_I << std::setw(10) << " asym_n =" << std::setw(10) << asymmetry_n_F << "\n";
    std::cout << std::setw(20) << "[Derived] params:" << std::setw(10) << " V0ws_n =" << std::setw(10) << V0WS_n_F << std::setw(10) << " V0ws_p =" << std::setw(10) << V0WS_p_F << std::setw(10) << " V0ls_n =" << std::setw(10) << V0LS_n_F << std::setw(10) << " V0ls_p =" << std::setw(10) << V0LS_p_F << "\n";
    std::cout << std::setw(20) << "[Shape] params:" << std::setw(10) << " beta2 =" << std::setw(10) << beta2_F << std::setw(10) << " beta3 =" << std::setw(10) << beta3_F << std::setw(10) << " beta4 =" << std::setw(10) << beta4_F << "\n";
    std::cout << std::setw(20) << "[Shape] derived:" << std::setw(10) << " zcm =" << std::setw(10) << ws_shape_.zcm_F << std::setw(10) << " fac =" << std::setw(10) << ws_shape_.surface_scale0_F << "\n";
}
