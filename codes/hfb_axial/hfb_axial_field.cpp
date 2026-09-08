/**
 * @file    hfb_axial_field.cpp
 * @author  cailea
 * @date    2026-05-05
 * @brief   Calculate axial Skyrme and pairing fields.
 */

#include <algorithm>
#include <cassert>
#include <cmath>

#include "hfb_axial.hpp"

namespace {

const double pi_F = std::acos(-1.0);

/**
 * @brief  Apply piecewise local pairing regularization.
 * @math   g_{reg}^{-1}=g^{-1}-I_{low}-I_{high}
 * @output Regularized local pairing coupling.
 */
double regularize_gr_at_point(double EspCut_F, double lambda_F, double vcent_F, double vmass_F, double gr_F) {
    if (gr_F == 0.0) {return gr_F;}

    const double Ewindow_F = EspCut_F - lambda_F;
    const double kFermi2_F = (lambda_F - vcent_F) / vmass_F;
    const double kHigh2_F = (lambda_F - vcent_F + Ewindow_F) / vmass_F;
    const double kLow2_F = (lambda_F - vcent_F - Ewindow_F) / vmass_F;

    // k_{low}²>0 → I_{low}+I_{high}.
    if (kLow2_F > 0.0) {
        const double kFermi_F = std::sqrt(kFermi2_F);
        const double kHigh_F = std::sqrt(kHigh2_F);
        const double kLow_F = std::sqrt(kLow2_F);
        const double low_term_F = 0.25 * (kLow_F / (pi_F * pi_F * vmass_F)) * (1.0 - 0.5 * (kFermi_F / kLow_F) * std::log((kFermi_F + kLow_F) / (kFermi_F - kLow_F)));
        const double high_term_F = 0.25 * (kHigh_F / (pi_F * pi_F * vmass_F)) * (1.0 - 0.5 * (kFermi_F / kHigh_F) * std::log((kHigh_F + kFermi_F) / (kHigh_F - kFermi_F)));
        gr_F = 1.0 / (1.0 / gr_F - low_term_F - high_term_F);
    } else if (kFermi2_F > 0.0) {
        // k_F²>0≥k_{low}² → I_{high}.
        const double kFermi_F = std::sqrt(kFermi2_F);
        const double kHigh_F = std::sqrt(kHigh2_F);
        const double high_term_F = 0.25 * (kHigh_F / (pi_F * pi_F * vmass_F)) * (1.0 - 0.5 * (kFermi_F / kHigh_F) * std::log((kHigh_F + kFermi_F) / (kHigh_F - kFermi_F)));
        gr_F = 1.0 / (1.0 / gr_F - high_term_F);
    } else if (kHigh2_F > 0.0) {
        // k_{high}²>0≥k_F² → I_{high}.
        const double kHigh_F = std::sqrt(kHigh2_F);
        const double kFermi_F = std::sqrt(std::abs(kFermi2_F));
        const double high_term_F = 0.25 * (kHigh_F / (pi_F * pi_F * vmass_F)) * (1.0 + (kFermi_F / kHigh_F) * std::atan(kFermi_F / kHigh_F));
        gr_F = 1.0 / (1.0 / gr_F - high_term_F);
    }
    return gr_F;
}

} // namespace

void AxialHFBField::set_zero() {
    vcent_F2D_z_r.setZero();
    vmass_F2D_z_r.setZero();
    vpair_F2D_z_r.setZero();
    vD2_F2D_z_r.setZero();
    vDr_F2D_z_r.setZero();
    vDz_F2D_z_r.setZero();
    vJzphi_F2D_z_r.setZero();
    vJphiz_F2D_z_r.setZero();
    vJphir_F2D_z_r.setZero();
    vJrphi_F2D_z_r.setZero();
    vdJ_F2D_z_r.setZero();
}

/**
 * @brief  Add local Skyrme functional derivatives.
 * @math   F_q(r) → F_q(r)+δE_{Skyrme}/δD_q(r)
 * @output Accumulated neutron and proton fields.
 */
void AxialHFBField::add_nuclei_fields(AxialHFBField& field_p_, AxialHFBField& field_n_, const AxialHFBDensity& density_p_, const AxialHFBDensity& density_n_, const EDFParamsSkyrme& edf_skyrme_) {
    const int Nz_I = static_cast<int>(field_n_.vcent_F2D_z_r.rows());
    const int Nr_I = static_cast<int>(field_n_.vcent_F2D_z_r.cols());

    // (D_n,D_p) → (D_0,D_1).
    const auto combine_rho01_Func = [](const Eigen::MatrixXd& rho_n_F2D_z_r, const Eigen::MatrixXd& rho_p_F2D_z_r, int z_I, int r_I) {
        const double rho_n_F = rho_n_F2D_z_r(z_I, r_I);
        const double rho_p_F = rho_p_F2D_z_r(z_I, r_I);
        const double rho0_F = rho_n_F + rho_p_F;
        const double rho1_F = rho_n_F - rho_p_F;
        return std::make_pair(rho0_F, rho1_F);
    };

    // D_t → δ𝓔/δD_q.
    #pragma omp parallel for collapse(2) schedule(static)
    for (int r_I = 0; r_I < Nr_I; ++r_I) {
        for (int z_I = 0; z_I < Nz_I; ++z_I) {
            const auto [dJ0_F, dJ1_F] = combine_rho01_Func(density_n_.dJ_F2D_z_r, density_p_.dJ_F2D_z_r, z_I, r_I);
            const auto [rho0_F, rho1_F] = combine_rho01_Func(density_n_.rho_F2D_z_r, density_p_.rho_F2D_z_r, z_I, r_I);
            const auto [tau0_F, tau1_F] = combine_rho01_Func(density_n_.tau_F2D_z_r, density_p_.tau_F2D_z_r, z_I, r_I);
            const auto [rhoD20_F, rhoD21_F] = combine_rho01_Func(density_n_.rhoD2_F2D_z_r, density_p_.rhoD2_F2D_z_r, z_I, r_I);
            const auto [rhoDr0_F, rhoDr1_F] = combine_rho01_Func(density_n_.rhoDr_F2D_z_r, density_p_.rhoDr_F2D_z_r, z_I, r_I);
            const auto [rhoDz0_F, rhoDz1_F] = combine_rho01_Func(density_n_.rhoDz_F2D_z_r, density_p_.rhoDz_F2D_z_r, z_I, r_I);
            const auto [Jphiz0_F, Jphiz1_F] = combine_rho01_Func(density_n_.Jphiz_F2D_z_r, density_p_.Jphiz_F2D_z_r, z_I, r_I);
            const auto [Jzphi0_F, Jzphi1_F] = combine_rho01_Func(density_n_.Jzphi_F2D_z_r, density_p_.Jzphi_F2D_z_r, z_I, r_I);
            const auto [Jphir0_F, Jphir1_F] = combine_rho01_Func(density_n_.Jphir_F2D_z_r, density_p_.Jphir_F2D_z_r, z_I, r_I);
            const auto [Jrphi0_F, Jrphi1_F] = combine_rho01_Func(density_n_.Jrphi_F2D_z_r, density_p_.Jrphi_F2D_z_r, z_I, r_I);

            // δ𝓔/δD_t ← 0.
            double dH_drho0_F = 0.0;
            double dH_drho1_F = 0.0;
            double dH_dtau0_F = 0.0;
            double dH_dtau1_F = 0.0;
            double dH_drhoD20_F = 0.0;
            double dH_drhoD21_F = 0.0;
            double dH_drhoDr0_F = 0.0;
            double dH_drhoDr1_F = 0.0;
            double dH_drhoDz0_F = 0.0;
            double dH_drhoDz1_F = 0.0;
            double dH_ddJ0_F = 0.0;
            double dH_ddJ1_F = 0.0;
            double dH_dJphiz0_F = 0.0;
            double dH_dJphiz1_F = 0.0;
            double dH_dJzphi0_F = 0.0;
            double dH_dJzphi1_F = 0.0;
            double dH_dJphir0_F = 0.0;
            double dH_dJphir1_F = 0.0;
            double dH_dJrphi0_F = 0.0;
            double dH_dJrphi1_F = 0.0;

            // 𝓔_ρ=(C_t^ρ+C_t^{ρ,α}ρ_0^α)ρ_t².
            const double rho0Sigma_F = std::pow(rho0_F, edf_skyrme_.sigma_F);
            const double drho0Sigma_drho0_F = edf_skyrme_.sigma_F * std::pow(rho0_F, edf_skyrme_.sigma_F - 1.0);
            dH_drho0_F += 2.0 * (edf_skyrme_.Crho_0_F + edf_skyrme_.Cdrho_0_F * rho0Sigma_F) * rho0_F + drho0Sigma_drho0_F * (edf_skyrme_.Cdrho_0_F * rho0_F * rho0_F + edf_skyrme_.Cdrho_1_F * rho1_F * rho1_F);
            dH_drho1_F += 2.0 * (edf_skyrme_.Crho_1_F + edf_skyrme_.Cdrho_1_F * rho0Sigma_F) * rho1_F;

            // 𝓔_τ=Σ_t C_t^τρ_tτ_t.
            dH_drho0_F += edf_skyrme_.Ctau_0_F * tau0_F;
            dH_drho1_F += edf_skyrme_.Ctau_1_F * tau1_F;
            dH_dtau0_F += edf_skyrme_.Ctau_0_F * rho0_F;
            dH_dtau1_F += edf_skyrme_.Ctau_1_F * rho1_F;

            // 𝓔_{Δρ}=Σ_t C_t^{Δρ}ρ_tΔρ_t.
            dH_drho0_F += edf_skyrme_.CrDr_0_F * rhoD20_F;
            dH_drho1_F += edf_skyrme_.CrDr_1_F * rhoD21_F;
            dH_drhoD20_F += edf_skyrme_.CrDr_0_F * rho0_F;
            dH_drhoD21_F += edf_skyrme_.CrDr_1_F * rho1_F;

            // 𝓔_{∇ρ}=Σ_t C_t^{∇ρ}|∇ρ_t|².
            dH_drhoDr0_F += 2.0 * edf_skyrme_.Cnrho_0_F * rhoDr0_F;
            dH_drhoDr1_F += 2.0 * edf_skyrme_.Cnrho_1_F * rhoDr1_F;
            dH_drhoDz0_F += 2.0 * edf_skyrme_.Cnrho_0_F * rhoDz0_F;
            dH_drhoDz1_F += 2.0 * edf_skyrme_.Cnrho_1_F * rhoDz1_F;

            // 𝓔_{∇J}=Σ_t C_t^{∇J}ρ_t∇·J_t.
            dH_drho0_F += edf_skyrme_.CrdJ_0_F * dJ0_F;
            dH_drho1_F += edf_skyrme_.CrdJ_1_F * dJ1_F;
            dH_ddJ0_F += edf_skyrme_.CrdJ_0_F * rho0_F;
            dH_ddJ1_F += edf_skyrme_.CrdJ_1_F * rho1_F;

            // J_r=J_{φz}-J_{zφ}; J_z=J_{rφ}-J_{φr}.
            const double Jr0_F = Jphiz0_F - Jzphi0_F;
            const double Jr1_F = Jphiz1_F - Jzphi1_F;
            const double Jz0_F = Jrphi0_F - Jphir0_F;
            const double Jz1_F = Jrphi1_F - Jphir1_F;
            dH_drhoDr0_F += edf_skyrme_.CJdr_0_F * Jr0_F;
            dH_drhoDr1_F += edf_skyrme_.CJdr_1_F * Jr1_F;
            dH_drhoDz0_F += edf_skyrme_.CJdr_0_F * Jz0_F;
            dH_drhoDz1_F += edf_skyrme_.CJdr_1_F * Jz1_F;
            const double dHJdr_dJr0_F = edf_skyrme_.CJdr_0_F * rhoDr0_F;
            const double dHJdr_dJr1_F = edf_skyrme_.CJdr_1_F * rhoDr1_F;
            const double dHJdr_dJz0_F = edf_skyrme_.CJdr_0_F * rhoDz0_F;
            const double dHJdr_dJz1_F = edf_skyrme_.CJdr_1_F * rhoDz1_F;
            dH_dJphiz0_F += dHJdr_dJr0_F;
            dH_dJphiz1_F += dHJdr_dJr1_F;
            dH_dJzphi0_F -= dHJdr_dJr0_F;
            dH_dJzphi1_F -= dHJdr_dJr1_F;
            dH_dJphir0_F -= dHJdr_dJz0_F;
            dH_dJphir1_F -= dHJdr_dJz1_F;
            dH_dJrphi0_F += dHJdr_dJz0_F;
            dH_dJrphi1_F += dHJdr_dJz1_F;

            // 𝓔_J=Σ_t C_t^JΣ_{ab}J_{t,ab}².
            dH_dJphiz0_F += 2.0 * edf_skyrme_.CJ_0_F * Jphiz0_F;
            dH_dJphiz1_F += 2.0 * edf_skyrme_.CJ_1_F * Jphiz1_F;
            dH_dJzphi0_F += 2.0 * edf_skyrme_.CJ_0_F * Jzphi0_F;
            dH_dJzphi1_F += 2.0 * edf_skyrme_.CJ_1_F * Jzphi1_F;
            dH_dJphir0_F += 2.0 * edf_skyrme_.CJ_0_F * Jphir0_F;
            dH_dJphir1_F += 2.0 * edf_skyrme_.CJ_1_F * Jphir1_F;
            dH_dJrphi0_F += 2.0 * edf_skyrme_.CJ_0_F * Jrphi0_F;
            dH_dJrphi1_F += 2.0 * edf_skyrme_.CJ_1_F * Jrphi1_F;

            // 𝓔_{J̄}=Σ_t C_t^{J̄}Σ_{ab}J_{t,ab}J_{t,ba}.
            dH_dJphiz0_F += 2.0 * edf_skyrme_.CJbar_0_F * Jzphi0_F;
            dH_dJphiz1_F += 2.0 * edf_skyrme_.CJbar_1_F * Jzphi1_F;
            dH_dJzphi0_F += 2.0 * edf_skyrme_.CJbar_0_F * Jphiz0_F;
            dH_dJzphi1_F += 2.0 * edf_skyrme_.CJbar_1_F * Jphiz1_F;
            dH_dJphir0_F += 2.0 * edf_skyrme_.CJbar_0_F * Jrphi0_F;
            dH_dJphir1_F += 2.0 * edf_skyrme_.CJbar_1_F * Jrphi1_F;
            dH_dJrphi0_F += 2.0 * edf_skyrme_.CJbar_0_F * Jphir0_F;
            dH_dJrphi1_F += 2.0 * edf_skyrme_.CJbar_1_F * Jphir1_F;

            // δ/δD_n=δ/δD_0+δ/δD_1.
            field_n_.vdJ_F2D_z_r(z_I, r_I) += dH_ddJ0_F + dH_ddJ1_F;
            field_n_.vcent_F2D_z_r(z_I, r_I) += dH_drho0_F + dH_drho1_F;
            field_n_.vmass_F2D_z_r(z_I, r_I) += dH_dtau0_F + dH_dtau1_F + edf_skyrme_.hbzeron_F;
            field_n_.vD2_F2D_z_r(z_I, r_I) += dH_drhoD20_F + dH_drhoD21_F;
            field_n_.vDr_F2D_z_r(z_I, r_I) += dH_drhoDr0_F + dH_drhoDr1_F;
            field_n_.vDz_F2D_z_r(z_I, r_I) += dH_drhoDz0_F + dH_drhoDz1_F;
            field_n_.vJphiz_F2D_z_r(z_I, r_I) += dH_dJphiz0_F + dH_dJphiz1_F;
            field_n_.vJzphi_F2D_z_r(z_I, r_I) += dH_dJzphi0_F + dH_dJzphi1_F;
            field_n_.vJphir_F2D_z_r(z_I, r_I) += dH_dJphir0_F + dH_dJphir1_F;
            field_n_.vJrphi_F2D_z_r(z_I, r_I) += dH_dJrphi0_F + dH_dJrphi1_F;

            // δ/δD_p=δ/δD_0-δ/δD_1.
            field_p_.vdJ_F2D_z_r(z_I, r_I) += dH_ddJ0_F - dH_ddJ1_F;
            field_p_.vcent_F2D_z_r(z_I, r_I) += dH_drho0_F - dH_drho1_F;
            field_p_.vmass_F2D_z_r(z_I, r_I) += dH_dtau0_F - dH_dtau1_F + edf_skyrme_.hbzerop_F;
            field_p_.vD2_F2D_z_r(z_I, r_I) += dH_drhoD20_F - dH_drhoD21_F;
            field_p_.vDr_F2D_z_r(z_I, r_I) += dH_drhoDr0_F - dH_drhoDr1_F;
            field_p_.vDz_F2D_z_r(z_I, r_I) += dH_drhoDz0_F - dH_drhoDz1_F;
            field_p_.vJphiz_F2D_z_r(z_I, r_I) += dH_dJphiz0_F - dH_dJphiz1_F;
            field_p_.vJzphi_F2D_z_r(z_I, r_I) += dH_dJzphi0_F - dH_dJzphi1_F;
            field_p_.vJphir_F2D_z_r(z_I, r_I) += dH_dJphir0_F - dH_dJphir1_F;
            field_p_.vJrphi_F2D_z_r(z_I, r_I) += dH_dJrphi0_F - dH_dJrphi1_F;
        }
    }
}

/**
 * @brief  Add direct and Slater-exchange Coulomb fields.
 * @math   v_C=K_Cρ_p-e²C_{ex}(3ρ_p/π)^{1/3}
 * @output Updated proton central field.
 */
void AxialHFBField::add_coulomb_field(AxialHFBField& field_p_, const AxialHFBDensity& density_p_, const AxialCoulombField& coulomb_field_, const EDFParamsSkyrme& edf_skyrme_) {
    assert(coulomb_field_.isBuilt_B);

    const int Nz_I = static_cast<int>(field_p_.vcent_F2D_z_r.rows());
    const int Nr_I = static_cast<int>(field_p_.vcent_F2D_z_r.cols());

    // ρ_p → v_C^{dir}.
    const Eigen::MatrixXd vCoulombDirect_F2D_z_r = coulomb_field_.calc_direct_field(density_p_.rho_F2D_z_r);
    const double coex_F = -edf_skyrme_.e2charg_F * std::cbrt(3.0 / pi_F);

    // v_C=v_C^{dir}+v_C^{Slater}.
    #pragma omp parallel for collapse(2) schedule(static)
    for (int r_I = 0; r_I < Nr_I; ++r_I) {
        for (int z_I = 0; z_I < Nz_I; ++z_I) {
            const double rhoPositive_p_F = std::max(density_p_.rho_F2D_z_r(z_I, r_I), 0.0);
            const double vCoulomb_F = vCoulombDirect_F2D_z_r(z_I, r_I) + edf_skyrme_.CExPar_F * coex_F * std::cbrt(rhoPositive_p_F);
            field_p_.vcent_F2D_z_r(z_I, r_I) += vCoulomb_F;
        }
    }
}

/**
 * @brief  Calculate local pairing and rearrangement fields.
 * @math   Δ_q=g_qκ_q; g_q=C_q^{V0}(1-C_q^{V1}ρ_0/ρ_c)
 * @output Updated pairing and central fields.
 */
void AxialHFBField::add_pairing_fields(AxialHFBField& field_p_, AxialHFBField& field_n_, const AxialHFBDensity& density_p_, const AxialHFBDensity& density_n_, const EDFParamsSkyrme& edf_skyrme_, const HFBSettings& hfbsettings_, double lambda_n_F, double lambda_p_F) {
    assert(hfbsettings_.termSwitches.addLocalPair_B);

    const int Nz_I = static_cast<int>(field_n_.vcent_F2D_z_r.rows());
    const int Nr_I = static_cast<int>(field_n_.vcent_F2D_z_r.cols());
    const double rhoc_F = 0.16;

    // (ρ_0,g_q,κ_q) → (g_{reg,q},Δ_q).
    if (hfbsettings_.termSwitches.useLocalPairRegularization_B) {
        #pragma omp parallel for collapse(2) schedule(static)
        for (int r_I = 0; r_I < Nr_I; ++r_I) {
            for (int z_I = 0; z_I < Nz_I; ++z_I) {
                const double rho0_F = density_n_.rho_F2D_z_r(z_I, r_I) + density_p_.rho_F2D_z_r(z_I, r_I);
                const double kappa_n_F = density_n_.kappa_F2D_z_r(z_I, r_I);
                const double kappa_p_F = density_p_.kappa_F2D_z_r(z_I, r_I);
                double gr_n_F = edf_skyrme_.CpV0_0_F * (1.0 - (rho0_F / rhoc_F) * edf_skyrme_.CpV1_0_F);
                double gr_p_F = edf_skyrme_.CpV0_1_F * (1.0 - (rho0_F / rhoc_F) * edf_skyrme_.CpV1_1_F);
                gr_n_F = regularize_gr_at_point(hfbsettings_.EspCut_F, lambda_n_F, field_n_.vcent_F2D_z_r(z_I, r_I), field_n_.vmass_F2D_z_r(z_I, r_I), gr_n_F);
                gr_p_F = regularize_gr_at_point(hfbsettings_.EspCut_F, lambda_p_F, field_p_.vcent_F2D_z_r(z_I, r_I), field_p_.vmass_F2D_z_r(z_I, r_I), gr_p_F);
                field_n_.vpair_F2D_z_r(z_I, r_I) += kappa_n_F * gr_n_F;
                field_p_.vpair_F2D_z_r(z_I, r_I) += kappa_p_F * gr_p_F;
            }
        }
        return;
    }

    // (g_q,κ_q) → (Δ_q,δ𝓔_pair/δρ_0).
    #pragma omp parallel for collapse(2) schedule(static)
    for (int r_I = 0; r_I < Nr_I; ++r_I) {
        for (int z_I = 0; z_I < Nz_I; ++z_I) {
            const double rho0_F = density_n_.rho_F2D_z_r(z_I, r_I) + density_p_.rho_F2D_z_r(z_I, r_I);
            const double kappa_n_F = density_n_.kappa_F2D_z_r(z_I, r_I);
            const double kappa_p_F = density_p_.kappa_F2D_z_r(z_I, r_I);
            const double gr_n_F = edf_skyrme_.CpV0_0_F * (1.0 - (rho0_F / rhoc_F) * edf_skyrme_.CpV1_0_F);
            const double gr_p_F = edf_skyrme_.CpV0_1_F * (1.0 - (rho0_F / rhoc_F) * edf_skyrme_.CpV1_1_F);
            field_n_.vpair_F2D_z_r(z_I, r_I) += kappa_n_F * gr_n_F;
            field_p_.vpair_F2D_z_r(z_I, r_I) += kappa_p_F * gr_p_F;
            const double dHpair_drho0_F = -(edf_skyrme_.CpV0_0_F * edf_skyrme_.CpV1_0_F / rhoc_F) * kappa_n_F * kappa_n_F - (edf_skyrme_.CpV0_1_F * edf_skyrme_.CpV1_1_F / rhoc_F) * kappa_p_F * kappa_p_F;
            field_n_.vcent_F2D_z_r(z_I, r_I) += dHpair_drho0_F;
            field_p_.vcent_F2D_z_r(z_I, r_I) += dHpair_drho0_F;
        }
    }
}
