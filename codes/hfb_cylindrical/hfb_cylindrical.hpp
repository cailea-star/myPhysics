/**
 * @file    hfb_cylindrical.hpp
 * @author  cailea
 * @date    2026-09-12
 * @brief   Cylindrical Kramers HFB model interface.
 */

#pragma once

#include <cassert>
#include <limits>

#include "hfb_kramers_nucleus.hpp"
#include "hfb_edf_setting.hpp"
#include "hfb_edf_skyrme.hpp"
#include "hfb_edf_gogny.hpp"
#include "cylindrical_config.hpp"
#include "cylindrical_basis.hpp"
#include "cylindrical_density.hpp"
#include "cylindrical_coulomb_field.hpp"
#include "cylindrical_gaussian_gogny.hpp"
#include "cylindrical_gaussian_coulomb.hpp"

class HFBKramersBlocking;

/**
 * @brief Store coordinate-space fields for one species.
 */
class HFBCylindricalField {
public:
    Eigen::MatrixXd vcent_F2D_z_r{}; // δE/δρ(z,r_⊥) [MeV].
    Eigen::MatrixXd vmass_F2D_z_r{}; // δE/δτ(z,r_⊥) [MeV fm²].
    Eigen::MatrixXd vpair_F2D_z_r{}; // δE/δκ(z,r_⊥) [MeV].
    Eigen::MatrixXd vD2_F2D_z_r{}; // δE/δΔρ(z,r_⊥) [MeV fm²].
    Eigen::MatrixXd vDr_F2D_z_r{}; // δE/δ∂_{r⊥}ρ(z,r_⊥) [MeV fm].
    Eigen::MatrixXd vDz_F2D_z_r{}; // δE/δ∂_zρ(z,r_⊥) [MeV fm].
    Eigen::MatrixXd vJzphi_F2D_z_r{}; // δE/δJ_{zφ}(z,r_⊥) [MeV fm].
    Eigen::MatrixXd vJphiz_F2D_z_r{}; // δE/δJ_{φz}(z,r_⊥) [MeV fm].
    Eigen::MatrixXd vJphir_F2D_z_r{}; // δE/δJ_{φr}(z,r_⊥) [MeV fm].
    Eigen::MatrixXd vJrphi_F2D_z_r{}; // δE/δJ_{rφ}(z,r_⊥) [MeV fm].
    Eigen::MatrixXd vdJ_F2D_z_r{}; // δE/δ∇·J(z,r_⊥) [MeV fm²].

public:
    /**
     * @brief Construct zero coordinate-space fields.
     * @math (N_z,N_r) → 0_{N_z×N_r}
     * @output Zero-initialized field grids.
     */
    explicit HFBCylindricalField(const CylindricalSetting& cylindricalsetting_) {
        // (N_z,N_r) → {F(z,r)}=0.
        vcent_F2D_z_r.resize(cylindricalsetting_.Nz_I, cylindricalsetting_.Nr_I);
        vmass_F2D_z_r.resize(cylindricalsetting_.Nz_I, cylindricalsetting_.Nr_I);
        vpair_F2D_z_r.resize(cylindricalsetting_.Nz_I, cylindricalsetting_.Nr_I);
        vD2_F2D_z_r.resize(cylindricalsetting_.Nz_I, cylindricalsetting_.Nr_I);
        vDr_F2D_z_r.resize(cylindricalsetting_.Nz_I, cylindricalsetting_.Nr_I);
        vDz_F2D_z_r.resize(cylindricalsetting_.Nz_I, cylindricalsetting_.Nr_I);
        vJzphi_F2D_z_r.resize(cylindricalsetting_.Nz_I, cylindricalsetting_.Nr_I);
        vJphiz_F2D_z_r.resize(cylindricalsetting_.Nz_I, cylindricalsetting_.Nr_I);
        vJphir_F2D_z_r.resize(cylindricalsetting_.Nz_I, cylindricalsetting_.Nr_I);
        vJrphi_F2D_z_r.resize(cylindricalsetting_.Nz_I, cylindricalsetting_.Nr_I);
        vdJ_F2D_z_r.resize(cylindricalsetting_.Nz_I, cylindricalsetting_.Nr_I);

        set_zero();
    }

    /**
     * @brief Reset all field grids.
     * @math {v_{cent},v_{mass},v_{pair},v_{D²},v_{∇},v_J} → 0
     * @output Zeroed field grids.
     */
    void set_zero();

    /**
     * @brief Add direct and exchange Coulomb fields.
     * @math ρ_p → v_C^{dir}+v_C^{Slater}
     * @output Updated proton central field.
     */
    static void add_coulomb_field(HFBCylindricalField& field_p_, const CylindricalDensity& density_p_, const CylindricalCoulombField& coulomb_field_, const EDFParamsSkyrme& edf_skyrme_);

    /**
     * @brief Add local pairing fields.
     * @math (ρ_0,κ_q,λ_q) → Δ_q
     * @output Updated neutron and proton pairing fields.
     */
    static void add_pairing_fields(HFBCylindricalField& field_p_, HFBCylindricalField& field_n_, const CylindricalDensity& density_p_, const CylindricalDensity& density_n_, const EDFParamsSkyrme& edf_skyrme_, const HFBEDFSetting& hfbedfsetting_, double lambda_n_F, double lambda_p_F);

    /**
     * @brief Add local neutron and proton fields.
     * @math (F_p,F_n) → (F_p,F_n)+(δE/δD_p,δE/δD_n)
     * @output Accumulated coordinate-space fields.
     */
    static void add_nuclei_fields(HFBCylindricalField& field_p_, HFBCylindricalField& field_n_, const CylindricalDensity& density_p_, const CylindricalDensity& density_n_, const EDFParamsSkyrme& edf_skyrme_);
};

class HFBKramersNucleusCylindrical : public HFBKramersNucleus {
public:
    double beta2Initial_F = 0.0; // Initial quadrupole deformation.
    double beta3Initial_F = 0.0; // Initial octupole deformation.
    double beta4Initial_F = 0.0; // Initial hexadecapole deformation.
    const HFBEDFSetting hfbedfsetting;
    const CylindricalSetting cylindricalsetting;
    CylindricalBasis2D cylindricalbasis;
    EDFParamsSkyrme edf_skyrme; // Base local EDF.
    CylindricalGaussianGogny gaussian_gogny; // Finite-range Gogny interaction.
    CylindricalGaussianCoulomb gaussian_coulomb; // Finite-range Coulomb interaction.
    CylindricalCoulombField coulomb_field; // ρ_p → v_C^{dir}.
    CylindricalDensity density_neutron;
    CylindricalDensity density_proton;
    HFBCylindricalField field_neutron;
    HFBCylindricalField field_proton;

public:
    /**
     * @brief Construct cylindrical Kramers HFB workspaces and interactions.
     * @math (C_cyl,P_HFB,P_Skyrme,P_Gogny) → HFB; η = 2Σ.
     * @output Allocated species, basis, densities, fields, and interaction objects.
     */
    HFBKramersNucleusCylindrical(const CylindricalSetting& cylindricalsetting_, const HFBEDFSetting& hfbedfsetting_, const EDFParamsSkyrme& edf_skyrme_, const EDFParamsGogny& edf_gogny_ = EDFParamsGogny::D1S())
    : HFBKramersNucleus(cylindricalsetting_.eta_F2D_block_bsp), hfbedfsetting(hfbedfsetting_), cylindricalsetting(cylindricalsetting_), cylindricalbasis(cylindricalsetting), edf_skyrme(edf_skyrme_), gaussian_gogny(cylindricalsetting, edf_gogny_.forceName_Str, edf_gogny_.mu_F1D_g, edf_gogny_.W_F1D_g, edf_gogny_.B_F1D_g, edf_gogny_.H_F1D_g, edf_gogny_.M_F1D_g), gaussian_coulomb(cylindricalsetting, edf_skyrme.e2charg_F), coulomb_field(cylindricalbasis), density_neutron(cylindricalsetting, cylindricalbasis), density_proton(cylindricalsetting, cylindricalbasis), field_neutron(cylindricalsetting), field_proton(cylindricalsetting) {
        assert(cylindricalsetting.useTimeReversal_B);

        // P_HFB → common tolerance and species temperature/cutoff.
        accuracy_F = hfbedfsetting.accuracy_F;
        hfb_neutron.temperature_F = hfbedfsetting.temperature_F;
        hfb_proton.temperature_F = hfbedfsetting.temperature_F;
        hfb_neutron.EspCut_F = hfbedfsetting.useEspCut_B ? hfbedfsetting.EspCut_F : std::numeric_limits<double>::infinity();
        hfb_proton.EspCut_F = hfb_neutron.EspCut_F;
    }

    /**
     * @brief Initialize neutron and proton Woods-Saxon fields.
     * @math (N,Z,β₂,β₃,β₄) → (v_cent,v_mass,v_dJ,v_pair)_{n,p}.
     * @output Initialized neutron and proton coordinate-space fields.
     */
    void initialize_WS_field();

    /**
     * @brief Accumulate coordinate-space fields into HFB block matrices.
     * @math (Γ,Δ) → (Γ,Δ)+(Γ[field],Δ[field]).
     * @output Accumulated neutron and proton Gamma and Delta.
     */
    void add_Gamma_Delta_from_field();




    /**
     * @brief Initialize neutron and proton one-body fields.
     * @math (N,Z) → (h₀,n,h₀,p).
     * @output Initialized neutron and proton one-body fields.
     */
    void initialize_h0() override;

    /**
     * @brief Initialize neutron and proton HFB fields.
     * @math (N,Z) → (Γ_n,Δ_n,Γ_p,Δ_p)_initial.
     * @output Initialized neutron and proton block fields.
     */
    void initialize_GammaDelta() override;

    /**
     * @brief Update neutron and proton HFB fields.
     * @math (ρ_n,κ_n,ρ_p,κ_p) → (Γ_n,Δ_n,Γ_p,Δ_p).
     * @output Updated neutron and proton block fields.
     */
    void update_Gamma_Delta() override;

    /**
     * @brief Print the current iteration summary.
     * @math (i,ε,α) → stdout.
     * @output Iteration diagnostics.
     */
    void print_abstract(int iteration_I, double error_F, double mixing_F) const override;

    /**
     * @brief Print iteration observables and supplied blocking trackers.
     * @math (i,ε,α,O,B_n,B_p) → stdout.
     * @output Formatted iteration and blocking rows.
     */
    void print_abstract(int iteration_I, double error_F, double mixing_F, const std::vector<HFBKramersBlocking>& neutronBlockings_, const std::vector<HFBKramersBlocking>& protonBlockings_) const;
};
