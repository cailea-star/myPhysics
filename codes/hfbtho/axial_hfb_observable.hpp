/**
 * @file    axial_hfb_observable.hpp
 * @author  cailea
 * @date    2026-05-02
 * @brief   Calculate axial HFB observables.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

#include "axial_hfb.hpp"

class AxialHFB;
class AxialHFBBlocking;

/**
 * @brief Store axial HFB observables.
 */
class AxialHFBObservable {
public:
    double A_F = 0.0; // A = N + Z.
    double Etot_F = 0.0; // E_{tot} [MeV].
    double beta2_F = 0.0; // β_2.
    double rms_n_F = 0.0; // √⟨r_n²⟩ [fm].
    double rms_p_F = 0.0; // √⟨r_p²⟩ [fm].
    double Epair_n_F = 0.0; // E_{pair,n} [MeV].
    double Epair_p_F = 0.0; // E_{pair,p} [MeV].
    double Delta_n_F = 0.0; // ⟨Δ_n⟩ [MeV].
    double Delta_p_F = 0.0; // ⟨Δ_p⟩ [MeV].
    double lambda_n_F = 0.0; // λ_n [MeV].
    double lambda_p_F = 0.0; // λ_p [MeV].
    double lambda2_n_F = 0.0; // λ_{2,n} [MeV].
    double lambda2_p_F = 0.0; // λ_{2,p} [MeV].

public:
    /**
     * @brief Calculate unblocked HFB observables.
     * @math HFB → O
     * @output Updated observables.
     */
    void update_observable(const AxialHFB& hfb_);

    /**
     * @brief Calculate singly blocked HFB observables.
     * @math (HFB,B_μ) → O
     * @output Updated observables.
     */
    void update_observable(const AxialHFB& hfb_, const AxialHFBBlocking& activeBlocking_);

    /**
     * @brief Calculate multiply blocked HFB observables.
     * @math (HFB,{B_μ}) → O
     * @output Updated observables.
     */
    void update_observable(const AxialHFB& hfb_, const std::vector<AxialHFBBlocking>& activeBlockings_);
};

namespace {

struct MatrixTrace {
    double Eph_F = 0.0;
    double Epair_F = 0.0;
    double DeltaRho_F = 0.0;
};

struct LocalEnergyTrace {
    double Ephysical_F = 0.0;
    double EphTrace_F = 0.0;
};

/**
 * @brief Calculate field-density matrix traces.
 * @math (Γ,ρ,Δ,κ) → (TrΓρ,TrΔκ,TrΔρ)
 * @output Particle-hole and pairing traces.
 */
MatrixTrace calc_matrix_trace(const AxialHFBBlockList& blocklist_) {
    MatrixTrace trace_;

    // Σ_b Tr(Γ_bρ_b,Δ_bκ_b,Δ_bρ_b).
    for (const auto& block_ : blocklist_.blocks_S1D_block) {
        trace_.Eph_F += block_.Gamma_F2D_bsp_bsp.cwiseProduct(block_.rho_F2D_bsp_bsp).sum();
        trace_.Epair_F += block_.Delta_F2D_bsp_bsp.cwiseProduct(block_.kappa_F2D_bsp_bsp).sum();
        trace_.DeltaRho_F += block_.Delta_F2D_bsp_bsp.cwiseProduct(block_.rho_F2D_bsp_bsp).sum();
    }

    return trace_;
}

/**
 * @brief Calculate local EDF energy traces.
 * @math (ρ,τ,J,h) → (E_{local},Tr[hρ])
 * @output Physical energy and field trace.
 */
LocalEnergyTrace calc_local_energy_trace(const AxialHFB& hfb_, const EDFParamsSkyrme& edf_skyrme_) {
    LocalEnergyTrace trace_;
    const auto& density_n_ = hfb_.densities.density_n;
    const auto& density_p_ = hfb_.densities.density_p;
    const auto& field_n_ = hfb_.fields.field_n;
    const auto& field_p_ = hfb_.fields.field_p;
    const auto& w_F2D_z_r = hfb_.global_basis.w_F2D_z_r;
    const int Nr_I = hfb_.axialconfig.Nr_I;
    const int Nz_I = hfb_.axialconfig.Nz_I;

    // E_C^{Slater} = C_Cρ_p^{4/3}.
    constexpr double pi_F = 3.141592653589793238462643383279502884;
    const double coulombExchangeCoeff_F = -0.75 * edf_skyrme_.CExPar_F * edf_skyrme_.e2charg_F * std::cbrt(3.0 / pi_F);

    // Tr[h_qρ_q](z,r).
    const auto calc_species_ph_trace_Func = [](const AxialHFBDensity& density_, const AxialHFBField& field_, int z_I, int r_I) {
        return field_.vcent_F2D_z_r(z_I, r_I) * density_.rho_F2D_z_r(z_I, r_I)
             + field_.vmass_F2D_z_r(z_I, r_I) * density_.tau_F2D_z_r(z_I, r_I)
             + field_.vD2_F2D_z_r(z_I, r_I) * density_.rhoD2_F2D_z_r(z_I, r_I)
             + field_.vDr_F2D_z_r(z_I, r_I) * density_.rhoDr_F2D_z_r(z_I, r_I)
             + field_.vDz_F2D_z_r(z_I, r_I) * density_.rhoDz_F2D_z_r(z_I, r_I)
             + field_.vdJ_F2D_z_r(z_I, r_I) * density_.dJ_F2D_z_r(z_I, r_I)
             + field_.vJzphi_F2D_z_r(z_I, r_I) * density_.Jzphi_F2D_z_r(z_I, r_I)
             + field_.vJphiz_F2D_z_r(z_I, r_I) * density_.Jphiz_F2D_z_r(z_I, r_I)
             + field_.vJphir_F2D_z_r(z_I, r_I) * density_.Jphir_F2D_z_r(z_I, r_I)
             + field_.vJrphi_F2D_z_r(z_I, r_I) * density_.Jrphi_F2D_z_r(z_I, r_I);
    };

    // E_{local} = ∫d³r ℰ_{Skyrme}.
    for (int r_I = 0; r_I < Nr_I; ++r_I) {
        for (int z_I = 0; z_I < Nz_I; ++z_I) {
            const double w_F = w_F2D_z_r(z_I, r_I);
            const double rho_n_F = density_n_.rho_F2D_z_r(z_I, r_I);
            const double rho_p_F = density_p_.rho_F2D_z_r(z_I, r_I);
            const double tau_n_F = density_n_.tau_F2D_z_r(z_I, r_I);
            const double tau_p_F = density_p_.tau_F2D_z_r(z_I, r_I);
            const double rho0_F = rho_n_F + rho_p_F;
            const double rho1_F = rho_n_F - rho_p_F;
            const double tau0_F = tau_n_F + tau_p_F;
            const double tau1_F = tau_n_F - tau_p_F;
            const double rhoD20_F = density_n_.rhoD2_F2D_z_r(z_I, r_I) + density_p_.rhoD2_F2D_z_r(z_I, r_I);
            const double rhoD21_F = density_n_.rhoD2_F2D_z_r(z_I, r_I) - density_p_.rhoD2_F2D_z_r(z_I, r_I);
            const double rhoDr0_F = density_n_.rhoDr_F2D_z_r(z_I, r_I) + density_p_.rhoDr_F2D_z_r(z_I, r_I);
            const double rhoDr1_F = density_n_.rhoDr_F2D_z_r(z_I, r_I) - density_p_.rhoDr_F2D_z_r(z_I, r_I);
            const double rhoDz0_F = density_n_.rhoDz_F2D_z_r(z_I, r_I) + density_p_.rhoDz_F2D_z_r(z_I, r_I);
            const double rhoDz1_F = density_n_.rhoDz_F2D_z_r(z_I, r_I) - density_p_.rhoDz_F2D_z_r(z_I, r_I);
            const double dJ0_F = density_n_.dJ_F2D_z_r(z_I, r_I) + density_p_.dJ_F2D_z_r(z_I, r_I);
            const double dJ1_F = density_n_.dJ_F2D_z_r(z_I, r_I) - density_p_.dJ_F2D_z_r(z_I, r_I);

            // J_t = (J_{φz}-J_{zφ},J_{rφ}-J_{φr}).
            const double Jzphi0_F = density_n_.Jzphi_F2D_z_r(z_I, r_I) + density_p_.Jzphi_F2D_z_r(z_I, r_I);
            const double Jzphi1_F = density_n_.Jzphi_F2D_z_r(z_I, r_I) - density_p_.Jzphi_F2D_z_r(z_I, r_I);
            const double Jphiz0_F = density_n_.Jphiz_F2D_z_r(z_I, r_I) + density_p_.Jphiz_F2D_z_r(z_I, r_I);
            const double Jphiz1_F = density_n_.Jphiz_F2D_z_r(z_I, r_I) - density_p_.Jphiz_F2D_z_r(z_I, r_I);
            const double Jphir0_F = density_n_.Jphir_F2D_z_r(z_I, r_I) + density_p_.Jphir_F2D_z_r(z_I, r_I);
            const double Jphir1_F = density_n_.Jphir_F2D_z_r(z_I, r_I) - density_p_.Jphir_F2D_z_r(z_I, r_I);
            const double Jrphi0_F = density_n_.Jrphi_F2D_z_r(z_I, r_I) + density_p_.Jrphi_F2D_z_r(z_I, r_I);
            const double Jrphi1_F = density_n_.Jrphi_F2D_z_r(z_I, r_I) - density_p_.Jrphi_F2D_z_r(z_I, r_I);
            const double Jr0_F = Jphiz0_F - Jzphi0_F;
            const double Jr1_F = Jphiz1_F - Jzphi1_F;
            const double Jz0_F = Jrphi0_F - Jphir0_F;
            const double Jz1_F = Jrphi1_F - Jphir1_F;

            // ρ_0^σ = max(ρ_0,0)^σ.
            const double rho0Positive_F = std::max(rho0_F, 0.0);
            const double rho0Sigma_F = rho0Positive_F > 0.0 ? std::pow(rho0Positive_F, edf_skyrme_.sigma_F) : 0.0;

            // ℰ_{Skyrme} = ℰ_{even}+ℰ_C^{Slater}.
            double energyDensity_F = 0.0;
            energyDensity_F += edf_skyrme_.hbzeron_F * tau_n_F + edf_skyrme_.hbzerop_F * tau_p_F;
            energyDensity_F += edf_skyrme_.Ctau_0_F * rho0_F * tau0_F + edf_skyrme_.Ctau_1_F * rho1_F * tau1_F;
            energyDensity_F += edf_skyrme_.Crho_0_F * rho0_F * rho0_F + edf_skyrme_.Crho_1_F * rho1_F * rho1_F;
            energyDensity_F += rho0Sigma_F * (edf_skyrme_.Cdrho_0_F * rho0_F * rho0_F + edf_skyrme_.Cdrho_1_F * rho1_F * rho1_F);
            energyDensity_F += edf_skyrme_.CrDr_0_F * rho0_F * rhoD20_F + edf_skyrme_.CrDr_1_F * rho1_F * rhoD21_F;
            energyDensity_F += edf_skyrme_.Cnrho_0_F * (rhoDr0_F * rhoDr0_F + rhoDz0_F * rhoDz0_F);
            energyDensity_F += edf_skyrme_.Cnrho_1_F * (rhoDr1_F * rhoDr1_F + rhoDz1_F * rhoDz1_F);
            energyDensity_F += edf_skyrme_.CrdJ_0_F * rho0_F * dJ0_F + edf_skyrme_.CrdJ_1_F * rho1_F * dJ1_F;
            energyDensity_F += edf_skyrme_.CJdr_0_F * (rhoDr0_F * Jr0_F + rhoDz0_F * Jz0_F);
            energyDensity_F += edf_skyrme_.CJdr_1_F * (rhoDr1_F * Jr1_F + rhoDz1_F * Jz1_F);
            energyDensity_F += edf_skyrme_.CJ_0_F * (Jphiz0_F * Jphiz0_F + Jzphi0_F * Jzphi0_F + Jphir0_F * Jphir0_F + Jrphi0_F * Jrphi0_F);
            energyDensity_F += edf_skyrme_.CJ_1_F * (Jphiz1_F * Jphiz1_F + Jzphi1_F * Jzphi1_F + Jphir1_F * Jphir1_F + Jrphi1_F * Jrphi1_F);
            energyDensity_F += 2.0 * edf_skyrme_.CJbar_0_F * (Jphiz0_F * Jzphi0_F + Jphir0_F * Jrphi0_F);
            energyDensity_F += 2.0 * edf_skyrme_.CJbar_1_F * (Jphiz1_F * Jzphi1_F + Jphir1_F * Jrphi1_F);

            if (hfb_.hfbsettings.termSwitches.addLocalCoulomb_B) {
                const double rho_pPositive_F = std::max(rho_p_F, 0.0);
                energyDensity_F += coulombExchangeCoeff_F * std::pow(rho_pPositive_F, 4.0 / 3.0);
            }

            trace_.Ephysical_F += w_F * energyDensity_F;
            trace_.EphTrace_F += w_F * (calc_species_ph_trace_Func(density_n_, field_n_, z_I, r_I) + calc_species_ph_trace_Func(density_p_, field_p_, z_I, r_I));
        }
    }

    // E_C^{dir} = 1/2 ∫ρ_pV_C^{dir}d³r.
    if (hfb_.hfbsettings.termSwitches.addLocalCoulomb_B && hfb_.fields.isCoulombKernelBuilt_B) {
        const int Ngrid_I = Nz_I * Nr_I;
        Eigen::Map<const Eigen::VectorXd> rho_F1D_grid(density_p_.rho_F2D_z_r.data(), Ngrid_I);
        Eigen::Map<const Eigen::MatrixXd> coulombKernel_F2D_source_target(hfb_.fields.coulomb_F4D_zs_rs_zt_rt.data(), Ngrid_I, Ngrid_I);
        Eigen::VectorXd Vcoulomb_F1D_grid = coulombKernel_F2D_source_target.transpose() * rho_F1D_grid;
        Eigen::Map<const Eigen::MatrixXd> Vcoulomb_F2D_z_r(Vcoulomb_F1D_grid.data(), Nz_I, Nr_I);

        for (int r_I = 0; r_I < Nr_I; ++r_I) {
            for (int z_I = 0; z_I < Nz_I; ++z_I) {
                trace_.Ephysical_F += 0.5 * w_F2D_z_r(z_I, r_I) * density_p_.rho_F2D_z_r(z_I, r_I) * Vcoulomb_F2D_z_r(z_I, r_I);
            }
        }
    }

    return trace_;
}

/**
 * @brief Estimate the last occupied single-particle energy.
 * @math λ_{last}=E_{sp,N/2}
 * @output Estimated last occupied energy.
 */
double estimate_last_lambda(const AxialHFBBlockList& blocklist_, const HFBSettings& hfbsettings_, const std::vector<AxialHFBBlocking>& activeBlockings_, bool isNeutron_B_, double Nparticle_F_, double pairTrace_F_) {
    double lambdaLast_F = blocklist_.lambda_F;

    // Pairing or T>0 → λ_{last}=λ.
    if (std::abs(pairTrace_F_) >= 1.0e-4 || hfbsettings_.temperature_F > 0.0) {
        return lambdaLast_F;
    }

    constexpr double EspCutTolerance_F = 1.0e-6;
    static const double EspCutTail_F = std::log(1.0 / EspCutTolerance_F - 1.0) / 100.0;
    std::vector<double> Esp_F1D_bqp;

    // E_{sp}=λ+E_{bqp}(1-2V²).
    for (int block_I = 0; block_I < static_cast<int>(blocklist_.blocks_S1D_block.size()); ++block_I) {
        const AxialHFBBlock& block_ = blocklist_.blocks_S1D_block[block_I];

        for (int bqp_I = 0; bqp_I < static_cast<int>(block_.Eqp_F1D_bqp.size()); ++bqp_I) {
            const double Eqp_F = block_.Eqp_F1D_bqp(bqp_I);
            if (!(Eqp_F > 0.0)) {
                continue;
            }

            double V2_F = block_.V_F2D_bsp_bqp.col(bqp_I).squaredNorm();

            // B_μ → V_μ²=1/2.
            for (const AxialHFBBlocking& activeBlocking_ : activeBlockings_) {
                if (activeBlocking_.isNeutron_B == isNeutron_B_ && activeBlocking_.block_I == block_I && activeBlocking_.bqp_I == bqp_I) {
                    V2_F = 0.5;
                }
            }

            const double Esp_F = blocklist_.lambda_F + Eqp_F * (1.0 - 2.0 * V2_F);
            if (hfbsettings_.useEspCut_B && Esp_F > hfbsettings_.EspCut_F + EspCutTail_F) {
                continue;
            }
            Esp_F1D_bqp.push_back(Esp_F);
        }
    }

    // λ_{last}=E_{sp,N/2}.
    const int Npair_I = static_cast<int>(std::floor(Nparticle_F_ + 0.1)) / 2;
    if (Npair_I <= 0 || Npair_I > static_cast<int>(Esp_F1D_bqp.size())) {
        return lambdaLast_F;
    }
    std::sort(Esp_F1D_bqp.begin(), Esp_F1D_bqp.end());
    lambdaLast_F = Esp_F1D_bqp[Npair_I - 1];

    return lambdaLast_F;
}

} // namespace

inline void AxialHFBObservable::update_observable(const AxialHFB& hfb_) {
    const std::vector<AxialHFBBlocking> noBlockings_;
    update_observable(hfb_, noBlockings_);
}

inline void AxialHFBObservable::update_observable(const AxialHFB& hfb_, const AxialHFBBlocking& activeBlocking_) {
    const std::vector<AxialHFBBlocking> activeBlockings_{activeBlocking_};
    update_observable(hfb_, activeBlockings_);
}

inline void AxialHFBObservable::update_observable(const AxialHFB& hfb_, const std::vector<AxialHFBBlocking>& activeBlockings_) {
    const auto& density_n_ = hfb_.densities.density_n;
    const auto& density_p_ = hfb_.densities.density_p;
    const int Nr_I = hfb_.axialconfig.Nr_I;
    const int Nz_I = hfb_.axialconfig.Nz_I;
    const auto& blocklist_n_ = hfb_.blocks.blocklist_n;
    const auto& blocklist_p_ = hfb_.blocks.blocklist_p;
    const auto& hfbsettings_ = hfb_.hfbsettings;
    const auto& w_F2D_z_r = hfb_.global_basis.w_F2D_z_r;
    const auto& r_F1D_r = hfb_.global_basis.r_F1D_r;
    const auto& z_F1D_z = hfb_.global_basis.z_F1D_z;

    // ∫d³r = Σ_{z,r}w_{zr}.
    const auto weight3D_Func = [&](int r_I, int z_I) {
        return w_F2D_z_r(z_I, r_I);
    };

    double r2Sum_n_F = 0.0;
    double r2Sum_p_F = 0.0;
    double r2P2Sum_n_F = 0.0;
    double r2P2Sum_p_F = 0.0;
    double NSum_n_F = 0.0;
    double NSum_p_F = 0.0;

    // (N_q,⟨r²⟩_q,Q_{2q}) = ∫ρ_q(1,r²,2r²P₂)d³r.
    for (int r_I = 0; r_I < Nr_I; ++r_I) {
        for (int z_I = 0; z_I < Nz_I; ++z_I) {
            const double r_F = r_F1D_r(r_I);
            const double r2_F = r_F * r_F;
            const double z_F = z_F1D_z(z_I);
            const double z2_F = z_F * z_F;
            const double distance2_F = z2_F + r2_F;
            const double w_F = weight3D_Func(r_I, z_I);
            const double rho_n_F = density_n_.rho_F2D_z_r(z_I, r_I);
            const double rho_p_F = density_p_.rho_F2D_z_r(z_I, r_I);
            const double r2P2_F = 1.5 * z2_F - 0.5 * distance2_F;

            NSum_n_F += w_F * rho_n_F;
            NSum_p_F += w_F * rho_p_F;
            r2Sum_n_F += w_F * distance2_F * rho_n_F;
            r2Sum_p_F += w_F * distance2_F * rho_p_F;
            r2P2Sum_n_F += w_F * r2P2_F * rho_n_F;
            r2P2Sum_p_F += w_F * r2P2_F * rho_p_F;
        }
    }

    // E = E_{local}+E_{nonlocal}+E_{pair}.
    const MatrixTrace matrixTrace_n_ = calc_matrix_trace(blocklist_n_);
    const MatrixTrace matrixTrace_p_ = calc_matrix_trace(blocklist_p_);
    const double Epair_n_F = matrixTrace_n_.Epair_F;
    const double Epair_p_F = matrixTrace_p_.Epair_F;
    const int Atarget_I = static_cast<int>(std::floor(NSum_n_F + NSum_p_F + 0.5));
    const EDFParamsSkyrme activeEDF_ = hfbsettings_.make_active_edf(hfb_.edf_skyrme, Atarget_I);
    const LocalEnergyTrace localEnergyTrace_ = calc_local_energy_trace(hfb_, activeEDF_);
    const double EphMatrixTrace_F = matrixTrace_n_.Eph_F + matrixTrace_p_.Eph_F;
    double EnonlocalPh_F = 0.0;

    if (hfbsettings_.termSwitches.addFiniteRange_B || hfbsettings_.termSwitches.addFiniteRangeCoulomb_B) {
        EnonlocalPh_F = EphMatrixTrace_F - 0.5 * localEnergyTrace_.EphTrace_F;
    }

    const double DeltaSum_n_F = -2.0 * matrixTrace_n_.DeltaRho_F;
    const double DeltaSum_p_F = -2.0 * matrixTrace_p_.DeltaRho_F;
    double Etot_F_ = localEnergyTrace_.Ephysical_F + EnonlocalPh_F + Epair_n_F + Epair_p_F;

    // O = (A,r_{rms},Q_2,Δ,λ,λ_2,E).
    const double N_n_F = NSum_n_F;
    const double N_p_F = NSum_p_F;
    const double A_F_ = N_n_F + N_p_F;
    const double rms_n_F_ = N_n_F > 0.0 ? std::sqrt(r2Sum_n_F / N_n_F) : 0.0;
    const double rms_p_F_ = N_p_F > 0.0 ? std::sqrt(r2Sum_p_F / N_p_F) : 0.0;
    const double rmsTotal_F = A_F_ > 0.0 ? std::sqrt((N_n_F * rms_n_F_ * rms_n_F_ + N_p_F * rms_p_F_ * rms_p_F_) / A_F_) : 0.0;
    const double Q2_n_F = 2.0 * r2P2Sum_n_F;
    const double Q2_p_F = 2.0 * r2P2Sum_p_F;
    constexpr double pi_F = 3.141592653589793238462643383279502884;
    const double beta2_F_ = rmsTotal_F > 0.0 ? std::sqrt(pi_F / 5.0) * (Q2_n_F + Q2_p_F) / (rmsTotal_F * rmsTotal_F * A_F_) : 0.0;
    const double Delta_n_F_ = N_n_F > 0.0 ? DeltaSum_n_F / N_n_F : 0.0;
    const double Delta_p_F_ = N_p_F > 0.0 ? DeltaSum_p_F / N_p_F : 0.0;
    double lambda2_n_F_ = 0.0;
    double lambda2_p_F_ = 0.0;

    if (hfbsettings_.useLipkinNogami_B) {
        lambda2_n_F_ = blocklist_n_.lambda2_F;
        lambda2_p_F_ = blocklist_p_.lambda2_F;
        Etot_F_ += blocklist_n_.ELipkinNogami_F + blocklist_p_.ELipkinNogami_F;
    }

    const double lambda_n_F_ = estimate_last_lambda(blocklist_n_, hfbsettings_, activeBlockings_, true, N_n_F, Epair_n_F);
    const double lambda_p_F_ = estimate_last_lambda(blocklist_p_, hfbsettings_, activeBlockings_, false, N_p_F, Epair_p_F);

    // O_{new} → this.
    A_F = A_F_;
    Etot_F = Etot_F_;
    beta2_F = beta2_F_;
    rms_n_F = rms_n_F_;
    rms_p_F = rms_p_F_;
    this->Epair_n_F = Epair_n_F;
    this->Epair_p_F = Epair_p_F;
    Delta_n_F = std::abs(Delta_n_F_);
    Delta_p_F = std::abs(Delta_p_F_);
    lambda_n_F = lambda_n_F_;
    lambda_p_F = lambda_p_F_;
    lambda2_n_F = lambda2_n_F_;
    lambda2_p_F = lambda2_p_F_;
}
