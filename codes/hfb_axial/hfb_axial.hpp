/**
 * @file    hfb_axial.hpp
 * @author  cailea
 * @date    2026-05-02
 * @brief   Define axial HFB blocks, densities, fields, and solver state.
 */

#pragma once

#include <cassert>
#include <utility>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Eigenvalues>

#include "cylindrical_basis.hpp"
#include "cylindrical_config.hpp"
#include "cylindrical_coulomb_field.hpp"
#include "cylindrical_gaussian_coulomb.hpp"
#include "cylindrical_gaussian_gogny.hpp"
#include "hfb_edf_gogny.hpp"
#include "hfb_edf_skyrme.hpp"
#include "hfb_settings.hpp"

class AxialHFBBlock;
class AxialHFBBlockList;
class AxialHFBBlocking;
class AxialHFBDensity;
class AxialHFBField;

/**
 * @brief Store one time-reversal-invariant axial block.
 */
class AxialHFBBlock {
public:
    int twoOmega_I = 0; // 2Ω = 2Λ + 2Σ.
    std::vector<CylindricalSPLabel> labels_S1D_bsp{}; // {α_bsp} within one block.
    std::vector<int> indices_I1D_bsp{}; // sp(block,bsp): global indices.
    std::vector<int> indices_I1D_bup{}; // {bsp | Σ=+1/2}.
    std::vector<int> indices_I1D_bdn{}; // {bsp | Σ=-1/2}.
    Eigen::VectorXd twoSigma_F1D_bsp{}; // 2Σ_bsp=±1.
    Eigen::MatrixXd UPos_F2D_bsp_bqp{}; // U⁺
    Eigen::MatrixXd VNeg_F2D_bsp_bqp{}; // V⁻.
    Eigen::MatrixXd rhoNegNeg_F2D_bsp_bsp{}; // Normal density ρ⁻⁻_{bsp,bsp'}.
    Eigen::MatrixXd rhoPosPos_F2D_bsp_bsp{}; // Derived ρ⁺⁺=Dρ⁻⁻D; D=diag(2Σ).
    Eigen::MatrixXd kappaPosNeg_F2D_bsp_bsp{}; // Pair density κ⁺⁻_{bsp,bsp'}.
    Eigen::MatrixXd kappaNegPos_F2D_bsp_bsp{}; // Derived κ⁻⁺=-(κ⁺⁻)ᵀ.
    Eigen::MatrixXd DeltaPosNeg_F2D_bsp_bsp{}; // Pairing field Δ⁺⁻_{bsp,bsp'} [MeV].
    Eigen::MatrixXd DeltaNegPos_F2D_bsp_bsp{}; // Pairing field Δ⁻⁺=-(Δ⁺⁻)ᵀ [MeV].
    Eigen::MatrixXd GammaPosPos_F2D_bsp_bsp{}; // h⁺⁺ including kinetic energy [MeV].
    Eigen::MatrixXd GammaNegNeg_F2D_bsp_bsp{}; // h⁻⁻ including kinetic energy; h⁻⁻=Dh⁺⁺D [MeV]; D=diag(2Σ).
    Eigen::VectorXd Eqp_F1D_bqp{}; // E⁺=E⁻; columns of U⁺,V⁻ [MeV].
    Eigen::MatrixXd HPos_F2D_2bsp_2bsp{}; // ℋ⁺ workspace [MeV].
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> HPos_eigensolver{}; // ℋ⁺ diagonalization workspace.

public:
    /**
     * @brief Construct one axial HFB block.
     * @math {α_{bsp},sp_{global}} → B_{2Ω}
     * @output Zero-initialized block matrices.
    */
    AxialHFBBlock(const std::vector<CylindricalSPLabel>& labels_S1D_bsp_, const std::vector<int>& indices_I1D_bsp_) {
        // ({α_bsp},{sp_global}) → block metadata.
        labels_S1D_bsp = labels_S1D_bsp_;
        indices_I1D_bsp = indices_I1D_bsp_;
        assert(!labels_S1D_bsp.empty());
        assert(labels_S1D_bsp.size() == indices_I1D_bsp.size());
        twoOmega_I = labels_S1D_bsp.front().twoOmega_I;

        // N_bsp → {0_{bsp×bsp},0_bsp}.
        UPos_F2D_bsp_bqp.resize(labels_S1D_bsp.size(), labels_S1D_bsp.size());
        VNeg_F2D_bsp_bqp.resize(labels_S1D_bsp.size(), labels_S1D_bsp.size());
        rhoNegNeg_F2D_bsp_bsp.resize(labels_S1D_bsp.size(), labels_S1D_bsp.size());
        rhoPosPos_F2D_bsp_bsp.resize(labels_S1D_bsp.size(), labels_S1D_bsp.size());
        kappaPosNeg_F2D_bsp_bsp.resize(labels_S1D_bsp.size(), labels_S1D_bsp.size());
        kappaNegPos_F2D_bsp_bsp.resize(labels_S1D_bsp.size(), labels_S1D_bsp.size());
        DeltaPosNeg_F2D_bsp_bsp.resize(labels_S1D_bsp.size(), labels_S1D_bsp.size());
        DeltaNegPos_F2D_bsp_bsp.resize(labels_S1D_bsp.size(), labels_S1D_bsp.size());
        GammaPosPos_F2D_bsp_bsp.resize(labels_S1D_bsp.size(), labels_S1D_bsp.size());
        GammaNegNeg_F2D_bsp_bsp.resize(labels_S1D_bsp.size(), labels_S1D_bsp.size());
        Eqp_F1D_bqp.resize(labels_S1D_bsp.size());
        HPos_F2D_2bsp_2bsp.resize(2 * labels_S1D_bsp.size(), 2 * labels_S1D_bsp.size());
        HPos_eigensolver = Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>(HPos_F2D_2bsp_2bsp.rows());

        UPos_F2D_bsp_bqp.setZero();
        VNeg_F2D_bsp_bqp.setZero();
        rhoNegNeg_F2D_bsp_bsp.setZero();
        rhoPosPos_F2D_bsp_bsp.setZero();
        kappaPosNeg_F2D_bsp_bsp.setZero();
        kappaNegPos_F2D_bsp_bsp.setZero();
        DeltaPosNeg_F2D_bsp_bsp.setZero();
        DeltaNegPos_F2D_bsp_bsp.setZero();
        GammaPosPos_F2D_bsp_bsp.setZero();
        GammaNegNeg_F2D_bsp_bsp.setZero();
        Eqp_F1D_bqp.setZero();
        HPos_F2D_2bsp_2bsp.setZero();

        // Σ_bsp → ({bsp_↑},{bsp_↓}).
        twoSigma_F1D_bsp.resize(labels_S1D_bsp.size());
        indices_I1D_bup.reserve(labels_S1D_bsp.size());
        indices_I1D_bdn.reserve(labels_S1D_bsp.size());
        for (int bsp_I = 0; bsp_I < static_cast<int>(labels_S1D_bsp.size()); ++bsp_I) {
            twoSigma_F1D_bsp(bsp_I) = labels_S1D_bsp[bsp_I].twoSigma_I;
            if (labels_S1D_bsp[bsp_I].twoSigma_I == 1) {indices_I1D_bup.push_back(bsp_I);}
            else {indices_I1D_bdn.push_back(bsp_I);}
        }
    }

    /**
     * @brief Reset particle-hole and pairing fields.
     * @math (Γ⁺⁺,Γ⁻⁻,Δ⁺⁻,Δ⁻⁺) → 0
     * @output Zeroed block fields.
     */
    void set_zero_Gamma_Delta();

    /**
     * @brief Add coordinate-space fields to this block.
     * @math (Γ_q,Δ_q) → (Γ_q,Δ_q)+(Γ_q^{loc},Δ_q^{loc})
     * @output Accumulated block fields.
     */
    void add_Gamma_Delta_from_field(const AxialHFBField& field_, const CylindricalBasis2D& global_basis_);

    /**
     * @brief Add the Lipkin-Nogami field correction.
     * @math h' → h'-2λ_2(1-2ρ)
     * @output Updated particle-hole field.
     */
    void add_lipkin_nogami(double lambda2_F);

    /**
     * @brief Solve the block HFB eigenproblem.
     * @math ℋ(U,V)^T=E_{bqp}(U,V)^T
     * @output Updated U, V, E, ρ, and κ.
     */
    void update_UV_E_rho_kappa(const HFBSettings& hfbsettings_, double lambda_F);
};

/**
 * @brief Store all axial blocks for one nucleon species.
 */
class AxialHFBBlockList {
public:
    bool isNeutron_B = false; // true: neutron; false: proton.
    double lambda_F = -7.0; // Chemical potential λ [MeV].
    double lambda2_F = 0.0; // Lipkin-Nogami λ₂ [MeV].
    double ELipkinNogami_F = 0.0; // Lipkin-Nogami energy [MeV].
    std::vector<AxialHFBBlock> blocks_X1D_block{}; // {B_{2Ω,π}}.

public:
    /**
     * @brief Construct all blocks for one species.
     * @math C_{axial} → {B_{2Ω,π}}
     * @output Zero-initialized block list.
     */
    AxialHFBBlockList(const CylindricalSetting& cylindricalsetting_, const HFBSettings&, bool isNeutron_B_) {
        isNeutron_B = isNeutron_B_;
        // C_axial → {B_{2Ω,π}}.
        assert(cylindricalsetting_.labels_S2D_block_bsp.size() == cylindricalsetting_.indices_I2D_block_bsp.size());
        const int Nblock_I = static_cast<int>(cylindricalsetting_.labels_S2D_block_bsp.size());
        blocks_X1D_block.reserve(Nblock_I);
        for (int block_I = 0; block_I < Nblock_I; ++block_I) {
            blocks_X1D_block.emplace_back(cylindricalsetting_.labels_S2D_block_bsp[block_I], cylindricalsetting_.indices_I2D_block_bsp[block_I]);
        }
    }

    /**
     * @brief Reset all block fields.
     * @math {Γ⁺⁺,Γ⁻⁻,Δ⁺⁻,Δ⁻⁺}_{block} → 0
     * @output Zeroed block fields.
     */
    void set_zero_Gamma_Delta();

    /**
     * @brief Add coordinate-space fields to all blocks.
     * @math {Γ_q,Δ_q}_{block} → {Γ_q+Γ_q^{loc},Δ_q+Δ_q^{loc}}_{block}
     * @output Accumulated block fields.
     */
    void add_Gamma_Delta_from_field(const AxialHFBField& field_, const CylindricalBasis2D& global_basis_);

    /**
     * @brief Apply the Lipkin-Nogami correction.
     * @math λ_2 → ({Γ_q}_{block},E_{LN})
     * @output Updated blocks and correction energy.
     * @note Supplied average gap overrides canonical reconstruction.
     */
    void add_lipkin_nogami();

    /**
     * @brief Solve every HFB block.
     * @math {ℋ_{block}} → {U,V,E,ρ,κ}_{block}
     * @output Updated block solutions.
     */
    void update_UV_E_rho_kappa(const HFBSettings& hfbsettings_);

    /**
     * @brief Add finite-range Gogny particle-hole fields.
     * @math \bar v_Gρ → (Γ_n,Γ_p)
     * @output Updated particle-hole matrices.
     */
    static void add_Gamma_from_Gogny(AxialHFBBlockList& blocklist_p_, AxialHFBBlockList& blocklist_n_, const CylindricalGaussianGogny& gogny_);

    /**
     * @brief Add finite-range Gogny pairing fields.
     * @math \bar v_Gκ → (Δ_n,Δ_p)
     * @output Updated pairing matrices.
     */
    static void add_Delta_from_Gogny(AxialHFBBlockList& blocklist_p_, AxialHFBBlockList& blocklist_n_, const CylindricalGaussianGogny& gogny_);

    /**
     * @brief Add finite-range Coulomb matrix elements.
     * @math \bar v_Cρ_p → Γ_p
     * @output Updated proton particle-hole matrices.
     */
    static void add_coulomb_from_Gaussian(AxialHFBBlockList& blocklist_p_, const CylindricalGaussianCoulomb& coulomb_);
};

/**
 * @brief Track one blocked quasiparticle state.
 */
class AxialHFBBlocking {
public:
    int block_I = -1; // Symmetry-block index.
    int bqp_I = -1; // Tracked block quasiparticle; -1 when absent.
    bool isNeutron_B = false; // Neutron when true; proton otherwise.
    double overlap_F = 0.0; // Blocking-state overlap.
    Eigen::VectorXd blockedV_F1D_bsp{}; // Previous blocked V column.
    Eigen::VectorXd blockedU_F1D_bsp{}; // Previous blocked U column.

public:
    /**
     * @brief Construct a blocking-state tracker.
     * @math (block,bqp,U,V) → B_{blocked}
     * @output Initialized blocked-state reference.
     */
    AxialHFBBlocking(const AxialHFBBlockList& blocklist_, bool isNeutron_B_, int block_I_, int bqp_I_) {
        // (block,bqp,q) → tracker metadata.
        block_I = block_I_;
        bqp_I = bqp_I_;
        isNeutron_B = isNeutron_B_;

        // (block,bqp) → (V_μ,U_μ).
        assert(block_I >= 0 && block_I < static_cast<int>(blocklist_.blocks_X1D_block.size()));
        const AxialHFBBlock& block_ = blocklist_.blocks_X1D_block[block_I];
        assert(bqp_I >= 0 && bqp_I < block_.VNeg_F2D_bsp_bqp.cols());
        assert(bqp_I < block_.UPos_F2D_bsp_bqp.cols());
        blockedV_F1D_bsp = block_.VNeg_F2D_bsp_bqp.col(bqp_I);
        blockedU_F1D_bsp = block_.UPos_F2D_bsp_bqp.col(bqp_I);
    }

    /**
     * @brief Apply equal-filling blocking to one block list.
     * @math (ρ,κ,U_μ,V_μ) → (ρ',κ')
     * @output Updated blocked densities and overlap.
     */
    void apply_blocking(AxialHFBBlockList& blocklist_);

    /**
     * @brief Enumerate low-energy blocking candidates.
     * @math E_{bqp}>0, |E_{bqp}-E_{min}|≤E_{cut} → {B_{blocked}}
     * @output Ranked blocking candidates.
     */
    static std::vector<AxialHFBBlocking> list_candidates(const HFBSettings& hfbsettings_, const AxialHFBBlockList& blocklist_, bool isNeutron_B_);
};

/**
 * @brief Store coordinate-space densities for one species.
 */
class AxialHFBDensity {
public:
    Eigen::MatrixXd rho_F2D_z_r{}; // ρ(z,r_⊥) [fm⁻³].
    Eigen::MatrixXd tau_F2D_z_r{}; // τ(z,r_⊥) [fm⁻⁵].
    Eigen::MatrixXd kappa_F2D_z_r{}; // κ(z,r_⊥) [fm⁻³].
    Eigen::MatrixXd rhoD2_F2D_z_r{}; // Δρ(z,r_⊥) [fm⁻⁵].
    Eigen::MatrixXd rhoDr_F2D_z_r{}; // ∂_{r⊥}ρ(z,r_⊥) [fm⁻⁴].
    Eigen::MatrixXd rhoDz_F2D_z_r{}; // ∂_zρ(z,r_⊥) [fm⁻⁴].
    Eigen::MatrixXd Jzphi_F2D_z_r{}; // J_{zφ}(z,r_⊥) [fm⁻⁴].
    Eigen::MatrixXd Jphiz_F2D_z_r{}; // J_{φz}(z,r_⊥) [fm⁻⁴].
    Eigen::MatrixXd Jphir_F2D_z_r{}; // J_{φr}(z,r_⊥) [fm⁻⁴].
    Eigen::MatrixXd Jrphi_F2D_z_r{}; // J_{rφ}(z,r_⊥) [fm⁻⁴].
    Eigen::MatrixXd dJ_F2D_z_r{}; // ∇·J(z,r_⊥) [fm⁻⁵].

public:
    /**
     * @brief Construct zero coordinate-space densities.
     * @math (N_z,N_r) → 0_{N_z×N_r}
     * @output Zero-initialized density grids.
     */
    explicit AxialHFBDensity(const CylindricalSetting& cylindricalsetting_) {
        // (N_z,N_r) → {D(z,r)}=0.
        rho_F2D_z_r.resize(cylindricalsetting_.Nz_I, cylindricalsetting_.Nr_I);
        tau_F2D_z_r.resize(cylindricalsetting_.Nz_I, cylindricalsetting_.Nr_I);
        kappa_F2D_z_r.resize(cylindricalsetting_.Nz_I, cylindricalsetting_.Nr_I);
        rhoD2_F2D_z_r.resize(cylindricalsetting_.Nz_I, cylindricalsetting_.Nr_I);
        rhoDr_F2D_z_r.resize(cylindricalsetting_.Nz_I, cylindricalsetting_.Nr_I);
        rhoDz_F2D_z_r.resize(cylindricalsetting_.Nz_I, cylindricalsetting_.Nr_I);
        Jzphi_F2D_z_r.resize(cylindricalsetting_.Nz_I, cylindricalsetting_.Nr_I);
        Jphiz_F2D_z_r.resize(cylindricalsetting_.Nz_I, cylindricalsetting_.Nr_I);
        Jphir_F2D_z_r.resize(cylindricalsetting_.Nz_I, cylindricalsetting_.Nr_I);
        Jrphi_F2D_z_r.resize(cylindricalsetting_.Nz_I, cylindricalsetting_.Nr_I);
        dJ_F2D_z_r.resize(cylindricalsetting_.Nz_I, cylindricalsetting_.Nr_I);

        set_zero();
    }

    /**
     * @brief Reset all density grids.
     * @math {ρ,τ,κ,Δρ,∇ρ,J} → 0
     * @output Zeroed density grids.
     */
    void set_zero();

    /**
     * @brief Rebuild densities from all species blocks.
     * @math {B_{block}} → (ρ,τ,κ,Δρ,∇ρ,J)
     * @output Updated density grids.
     */
    void update_density(const CylindricalBasis2D& global_basis_, const AxialHFBBlockList& blocklist_);
};

/**
 * @brief Store coordinate-space fields for one species.
 */
class AxialHFBField {
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
    explicit AxialHFBField(const CylindricalSetting& cylindricalsetting_) {
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
    static void add_coulomb_field(AxialHFBField& field_p_, const AxialHFBDensity& density_p_, const CylindricalCoulombField& coulomb_field_, const EDFParamsSkyrme& edf_skyrme_);

    /**
     * @brief Add local pairing fields.
     * @math (ρ_0,κ_q,λ_q) → Δ_q
     * @output Updated neutron and proton pairing fields.
     */
    static void add_pairing_fields(AxialHFBField& field_p_, AxialHFBField& field_n_, const AxialHFBDensity& density_p_, const AxialHFBDensity& density_n_, const EDFParamsSkyrme& edf_skyrme_, const HFBSettings& hfbsettings_, double lambda_n_F, double lambda_p_F);

    /**
     * @brief Add local neutron and proton fields.
     * @math (F_p,F_n) → (F_p,F_n)+(δE/δD_p,δE/δD_n)
     * @output Accumulated coordinate-space fields.
     */
    static void add_nuclei_fields(AxialHFBField& field_p_, AxialHFBField& field_n_, const AxialHFBDensity& density_p_, const AxialHFBDensity& density_n_, const EDFParamsSkyrme& edf_skyrme_);
};

/**
 * @brief Own one complete axial HFB calculation.
 */
class AxialHFB {
public:
    CylindricalSetting cylindricalsetting; // Axial basis configuration.
    HFBSettings hfbsettings; // HFB iteration controls.
    EDFParamsSkyrme edf_skyrme; // Base local EDF.
    CylindricalGaussianGogny gogny; // Finite-range Gogny interaction.
    CylindricalGaussianCoulomb coulomb; // Finite-range Coulomb interaction.
    CylindricalBasis2D global_basis; // Global axial basis.
    AxialHFBBlockList blocklist_n; // Neutron HFB blocks.
    AxialHFBBlockList blocklist_p; // Proton HFB blocks.
    AxialHFBDensity density_p; // Proton local densities.
    AxialHFBDensity density_n; // Neutron local densities.
    AxialHFBField field_p; // Proton local fields.
    AxialHFBField field_n; // Neutron local fields.
    CylindricalCoulombField coulombField; // ρ_p → v_C^{dir}.

public:
    /**
     * @brief Construct an axial HFB solver with D1S.
     * @math (C_{axial},P_{HFB},C_{Skyrme}) → HFB[D1S]
     * @output Initialized HFB solver.
     */
    AxialHFB(const CylindricalSetting& cylindricalsetting_, const HFBSettings& hfbsettings_, const EDFParamsSkyrme& edf_skyrme_)
    : AxialHFB(cylindricalsetting_, hfbsettings_, edf_skyrme_, EDFParamsGogny::D1S()) {}

    /**
     * @brief Construct an axial HFB solver from Gogny parameters.
     * @math (C_{axial},P_{HFB},C_{Skyrme},P_G) → HFB
     * @output Initialized HFB solver.
     */
    AxialHFB(const CylindricalSetting& cylindricalsetting_, const HFBSettings& hfbsettings_, const EDFParamsSkyrme& edf_skyrme_, const EDFParamsGogny& edf_gogny_)
    : AxialHFB(cylindricalsetting_, hfbsettings_, edf_skyrme_, CylindricalGaussianGogny(cylindricalsetting_, edf_gogny_.forceName_Str, edf_gogny_.mu_F1D_g, edf_gogny_.W_F1D_g, edf_gogny_.B_F1D_g, edf_gogny_.H_F1D_g, edf_gogny_.M_F1D_g)) {}

    /**
     * @brief Construct an axial HFB solver from interactions.
     * @math (C_{axial},P_{HFB},C_{Skyrme},V_G) → HFB
     * @output Initialized HFB solver.
     */
    AxialHFB(const CylindricalSetting& cylindricalsetting_, const HFBSettings& hfbsettings_, const EDFParamsSkyrme& edf_skyrme_, CylindricalGaussianGogny gogny_)
    : cylindricalsetting(cylindricalsetting_), hfbsettings(hfbsettings_), edf_skyrme(edf_skyrme_), gogny(std::move(gogny_)), coulomb(cylindricalsetting_, edf_skyrme_.e2charg_F), global_basis(cylindricalsetting_), blocklist_n(cylindricalsetting_, hfbsettings_, true), blocklist_p(cylindricalsetting_, hfbsettings_, false), density_p(cylindricalsetting_), density_n(cylindricalsetting_), field_p(cylindricalsetting_), field_n(cylindricalsetting_), coulombField(global_basis) {}

    /**
     * @brief Initialize deformed Woods-Saxon fields.
     * @math (N,Z,β_2,β_3,β_4) → (h_n^{WS},h_p^{WS})
     * @output Initialized coordinate-space fields.
     */
    void initialize_WS_field(int Ntarget_I, int Ztarget_I, double beta2_F = 0.0, double beta3_F = 0.0, double beta4_F = 0.0);

    /**
     * @brief Iterate the unblocked axial HFB equations.
     * @math (N,Z) → HFB_{converged}
     * @output Updated converged solver state.
     */
    void iterate(int Ntarget_I, int Ztarget_I);

    /**
     * @brief Iterate one blocked axial HFB state.
     * @math (N,Z,B_μ) → HFB_{blocked}
     * @output Updated blocked solver state.
     */
    void iterate(int Ntarget_I, int Ztarget_I, AxialHFBBlocking& activeBlocking_);

    /**
     * @brief Iterate multiple blocked axial HFB states.
     * @math (N,Z,{B_μ}) → {HFB_{blocked}}
     * @output Updated blocked solver states.
     */
    void iterate(int Ntarget_I, int Ztarget_I, std::vector<AxialHFBBlocking>& activeBlockings_);

private:
    /**
     * @brief Update chemical potential using Brent root search.
     * @math N_q(λ_q)-N_q^{target}=0
     * @output Updated λ, U, V, E, ρ, and κ.
     * @note Runtime branches implement root bracketing.
     */
    void update_blocklist_lambda(AxialHFBBlockList& blocklist_, int Ntarget_I, std::vector<AxialHFBBlocking>& activeBlockings_, double lambdaTolerance_F);
};
