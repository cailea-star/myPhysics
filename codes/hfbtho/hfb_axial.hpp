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
#include <unsupported/Eigen/CXX11/Tensor>

#include "axial_basis.hpp"
#include "axial_config.hpp"
#include "axial_gaussian_coulomb.hpp"
#include "axial_gaussian_gogny.hpp"
#include "hfb_edf_gogny.hpp"
#include "hfb_edf_skyrme.hpp"

/**
 * @brief Store HFB iteration controls and term switches.
 */
class HFBSettings {
public:
    /**
     * @brief Store enabled HFB energy terms.
     */
    class TermSwitches {
    public:
        bool addKinetic_B; // E_kin = (ℏ²/2m)τ.
        bool addLocalRhoRho_B; // E_ρρ = Σ_t C_t^ρρ_t².
        bool addLocalRhoAlpha_B; // E_ρα = Σ_t C_t^{ρ,α}ρ_0^αρ_t².
        bool addLocalRhoTau_B; // E_ρτ = Σ_t C_t^τρ_tτ_t.
        bool addLocalSurface_B; // E_surf = Σ_t(C_t^{Δρ}ρ_tΔρ_t+C_t^{∇ρ}|∇ρ_t|²).
        bool addLocalSpinOrbit_B; // E_SO = Σ_t(C_t^{∇J}ρ_t∇·J_t+C_t^{J∇ρ}J_t·∇ρ_t).
        bool addLocalTensor_B; // E_tensor = Σ_t(C_t^J J_t²+C_t^{J̄}J_{t,ab}J_{t,ba}).
        bool addLocalCoulomb_B; // E_C = E_C^dir+E_C^Slater.
        bool addLocalPair_B; // E_pair = Σ_q∫κ_q*Δ_q d³r.
        bool useLocalPairRegularization_B; // C_pair → C_pair^reg.
        bool addFiniteRange_B; // V_G = Σ_g(W_g+B_gP_σ-H_gP_τ-M_gP_σP_τ)e^{-r₁₂²/μ_g²}.
        bool addFiniteRangeCoulomb_B; // V_C = e²/r₁₂.

    public:
        /**
         * @brief Construct default local-EDF switches.
         * @math ∅ → S_{Skyrme}
         * @output Enabled standard Skyrme terms.
         */
        TermSwitches()
        : addKinetic_B(true), addLocalRhoRho_B(true), addLocalRhoAlpha_B(true), addLocalRhoTau_B(true), addLocalSurface_B(true), addLocalSpinOrbit_B(true), addLocalTensor_B(false), addLocalCoulomb_B(true), addLocalPair_B(true), useLocalPairRegularization_B(false), addFiniteRange_B(false), addFiniteRangeCoulomb_B(false) {}

        /**
         * @brief Build local Skyrme term switches.
         * @math ∅ → S_{Skyrme}
         * @output Standard local-EDF switches.
         */
        static TermSwitches skyrme() {
            return TermSwitches();
        }

        /**
         * @brief Build finite-range Gogny term switches.
         * @math S_{Skyrme} → S_{Gogny}
         * @output Gogny interaction switches.
         */
        static TermSwitches gogny() {
            TermSwitches termSwitches_;
            termSwitches_.addLocalRhoRho_B = false;
            termSwitches_.addLocalRhoTau_B = false;
            termSwitches_.addLocalSurface_B = false;
            termSwitches_.addLocalTensor_B = false;
            termSwitches_.addLocalPair_B = false;
            termSwitches_.addFiniteRange_B = true;
            termSwitches_.addLocalCoulomb_B = false;
            termSwitches_.addFiniteRangeCoulomb_B = true;
            return termSwitches_;
        }
    };

    int Nblocking_I; // Number of blocking candidates.
    double EblockingCut_F; // Blocking window [MeV].
    double EspCut_F; // Active quasiparticle cutoff [MeV].
    double temperature_F; // Temperature [MeV].
    double accuracy_F; // SCF convergence tolerance.
    bool useEspCut_B; // E_qp ≤ E_sp^cut.
    bool useLipkinNogami_B; // Enable Lipkin-Nogami correction.
    bool useCmCorrection_B; // ℏ²/(2m) → [1-1/A]ℏ²/(2m).
    TermSwitches termSwitches; // Enabled energy terms.

public:
    /**
     * @brief Construct default HFB controls.
     * @math ∅ → P_{HFB}
     * @output Initialized solver settings.
     */
    HFBSettings()
    : Nblocking_I(6), EblockingCut_F(1.0), EspCut_F(60.0), temperature_F(0.0), accuracy_F(1.0e-5), useEspCut_B(true), useLipkinNogami_B(true), useCmCorrection_B(false), termSwitches(TermSwitches::skyrme()) {}

    /**
     * @brief Calculate the spherical oscillator length.
     * @math A → b_0
     * @output Oscillator length b_0 [fm].
     */
    static double calc_b0(int Atarget_I);

    /**
     * @brief Apply nucleus-dependent EDF corrections.
     * @math (C,A,P_{HFB}) → C_{active}
     * @output Active Skyrme EDF parameters.
     */
    EDFParamsSkyrme make_active_edf(const EDFParamsSkyrme& base_edf_, int Atarget_I) const;
};

class AxialHFBBlock;
class AxialHFBBlockList;
class AxialHFBBlockListSystem;
class AxialHFBBlocking;
class AxialHFBDensity;
class AxialHFBDensitySystem;
class AxialHFBField;
class AxialHFBFieldSystem;

/**
 * @brief Store one axial symmetry block.
 */
class AxialHFBBlock {
public:
    int twoOmega_I; // 2Ω = 2Λ + 2Σ.
    std::vector<AxialSPLabel> labels_S1D_bsp; // {α_bsp} within one block.
    std::vector<int> indices_I1D_bsp; // sp(block,bsp): global indices.
    std::vector<int> indices_I1D_bup; // {bsp | Σ=+1/2}.
    std::vector<int> indices_I1D_bdn; // {bsp | Σ=-1/2}.
    Eigen::MatrixXd U_F2D_bsp_bqp; // Bogoliubov U_{bsp,bqp}.
    Eigen::MatrixXd V_F2D_bsp_bqp; // Bogoliubov V_{bsp,bqp}.
    Eigen::MatrixXd rho_F2D_bsp_bsp; // Normal density ρ_{bsp,bsp'}.
    Eigen::MatrixXd kappa_F2D_bsp_bsp; // Pair density κ_{bsp,bsp'}.
    Eigen::MatrixXd Delta_F2D_bsp_bsp; // Pairing field Δ_{bsp,bsp'} [MeV].
    Eigen::MatrixXd Gamma_F2D_bsp_bsp; // Particle-hole field Γ_{bsp,bsp'} [MeV].
    Eigen::VectorXd Eqp_F1D_bqp; // Quasiparticle energies E_bqp [MeV].

public:
    /**
     * @brief Construct one axial HFB block.
     * @math {α_{bsp},sp_{global}} → B_{2Ω}
     * @output Zero-initialized block matrices.
    */
    AxialHFBBlock(const std::vector<AxialSPLabel>& labels_S1D_bsp_, const std::vector<int>& indices_I1D_bsp_)
    : twoOmega_I(0), labels_S1D_bsp(labels_S1D_bsp_), indices_I1D_bsp(indices_I1D_bsp_), indices_I1D_bup(), indices_I1D_bdn(), U_F2D_bsp_bqp(Eigen::MatrixXd::Zero(labels_S1D_bsp_.size(), labels_S1D_bsp_.size())), V_F2D_bsp_bqp(Eigen::MatrixXd::Zero(labels_S1D_bsp_.size(), labels_S1D_bsp_.size())), rho_F2D_bsp_bsp(Eigen::MatrixXd::Zero(labels_S1D_bsp_.size(), labels_S1D_bsp_.size())), kappa_F2D_bsp_bsp(Eigen::MatrixXd::Zero(labels_S1D_bsp_.size(), labels_S1D_bsp_.size())), Delta_F2D_bsp_bsp(Eigen::MatrixXd::Zero(labels_S1D_bsp_.size(), labels_S1D_bsp_.size())), Gamma_F2D_bsp_bsp(Eigen::MatrixXd::Zero(labels_S1D_bsp_.size(), labels_S1D_bsp_.size())), Eqp_F1D_bqp(Eigen::VectorXd::Zero(labels_S1D_bsp_.size())) {
        assert(!labels_S1D_bsp.empty());
        assert(labels_S1D_bsp.size() == indices_I1D_bsp.size());
        twoOmega_I = labels_S1D_bsp.front().twoOmega_I;
        indices_I1D_bup.reserve(labels_S1D_bsp.size());
        indices_I1D_bdn.reserve(labels_S1D_bsp.size());
        for (int bsp_I = 0; bsp_I < static_cast<int>(labels_S1D_bsp.size()); ++bsp_I) {
            if (labels_S1D_bsp[bsp_I].twoSigma_I == 1) {indices_I1D_bup.push_back(bsp_I);}
            else {indices_I1D_bdn.push_back(bsp_I);}
        }
    }

    /**
     * @brief Project coordinate-space fields into this block.
     * @math (h_q(r),Δ_q(r),Φ) → (Γ_q,Δ_q)
     * @output Updated block fields.
     */
    void update_Gamma_Delta_from_field(const AxialHFBField& field_, const AxialBasis& global_basis_);

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
    double lambda_F; // Chemical potential λ [MeV].
    double lambda2_F; // Lipkin-Nogami λ₂ [MeV].
    double ELipkinNogami_F; // Lipkin-Nogami energy [MeV].
    std::vector<AxialHFBBlock> blocks_S1D_block; // {B_{2Ω,π}}.

public:
    /**
     * @brief Construct all blocks for one species.
     * @math C_{axial} → {B_{2Ω,π}}
     * @output Zero-initialized block list.
     */
    AxialHFBBlockList(const AxialConfig& axialconfig_, const HFBSettings&)
    : lambda_F(0.0), lambda2_F(0.0), ELipkinNogami_F(0.0), blocks_S1D_block() {
        assert(axialconfig_.labels_S2D_block_bsp.size() == axialconfig_.indices_I2D_block_bsp.size());
        const int Nblock_I = static_cast<int>(axialconfig_.labels_S2D_block_bsp.size());
        blocks_S1D_block.reserve(Nblock_I);
        for (int block_I = 0; block_I < Nblock_I; ++block_I) {
            blocks_S1D_block.emplace_back(axialconfig_.labels_S2D_block_bsp[block_I], axialconfig_.indices_I2D_block_bsp[block_I]);
        }
    }

    /**
     * @brief Project coordinate-space fields into all blocks.
     * @math (h_q(r),Δ_q(r),Φ) → {Γ_q,Δ_q}_{block}
     * @output Updated block fields.
     */
    void update_Gamma_Delta_from_field(const AxialHFBField& field_, const AxialBasis& global_basis_);

    /**
     * @brief Apply the Lipkin-Nogami correction.
     * @math λ_2 → ({Γ_q}_{block},E_{LN})
     * @output Updated blocks and correction energy.
     */
    void add_lipkin_nogami();

    /**
     * @brief Solve every HFB block.
     * @math {ℋ_{block}} → {U,V,E,ρ,κ}_{block}
     * @output Updated block solutions.
     */
    void update_UV_E_rho_kappa(const HFBSettings& hfbsettings_);
};

/**
 * @brief Manage coupled neutron and proton HFB blocks.
 */
class AxialHFBBlockListSystem {
public:
    AxialHFBBlockList blocklist_n; // Neutron blocks.
    AxialHFBBlockList blocklist_p; // Proton blocks.

public:
    /**
     * @brief Construct neutron and proton block lists.
     * @math C_{axial} → ({B_n},{B_p})
     * @output Zero-initialized species blocks.
     */
    AxialHFBBlockListSystem(const AxialConfig& axialconfig_, const HFBSettings& hfbsettings_)
    : blocklist_n(axialconfig_, hfbsettings_), blocklist_p(axialconfig_, hfbsettings_) {}

    /**
     * @brief Project neutron and proton fields into blocks.
     * @math (h_n,h_p,Δ_n,Δ_p,Φ) → ({Γ_n,Δ_n},{Γ_p,Δ_p})
     * @output Updated neutron and proton blocks.
     */
    void update_Gamma_Delta_from_fields(const AxialHFBFieldSystem& fields_, const AxialBasis& global_basis_);

    /**
     * @brief Add finite-range Gogny particle-hole fields.
     * @math \bar v_Gρ → (Γ_n,Γ_p)
     * @output Updated particle-hole matrices.
     */
    void add_Gamma_from_Gogny(const AxialGaussianGogny& gogny_);

    /**
     * @brief Add finite-range Gogny pairing fields.
     * @math \bar v_Gκ → (Δ_n,Δ_p)
     * @output Updated pairing matrices.
     */
    void add_Delta_from_Gogny(const AxialGaussianGogny& gogny_);

    /**
     * @brief Add finite-range Coulomb matrix elements.
     * @math \bar v_Cρ_p → Γ_p
     * @output Updated proton particle-hole matrices.
     */
    void add_coulomb_from_Gaussian(const AxialGaussianCoulomb& coulomb_);

    /**
     * @brief Apply species-wise Lipkin-Nogami corrections.
     * @math (λ_{2,n},λ_{2,p}) → ({Γ_n},{Γ_p})
     * @output Updated neutron and proton blocks.
     */
    void add_lipkin_nogami();
};

/**
 * @brief Track one blocked quasiparticle state.
 */
class AxialHFBBlocking {
public:
    int block_I; // Symmetry-block index.
    int bqp_I; // Tracked block quasiparticle; -1 when absent.
    bool isNeutron_B; // Neutron when true; proton otherwise.
    double overlap_F; // Blocking-state overlap.
    Eigen::VectorXd blockedV_F1D_bsp; // Previous blocked V column.
    Eigen::VectorXd blockedU_F1D_bsp; // Previous blocked U column.

public:
    /**
     * @brief Construct a blocking-state tracker.
     * @math (block,bqp,U,V) → B_{blocked}
     * @output Initialized blocked-state reference.
     */
    AxialHFBBlocking(const AxialHFBBlockList& blocklist_, bool isNeutron_B_, int block_I_, int bqp_I_)
    : block_I(block_I_), bqp_I(bqp_I_), isNeutron_B(isNeutron_B_), overlap_F(0.0), blockedV_F1D_bsp(), blockedU_F1D_bsp() {
        assert(block_I >= 0 && block_I < static_cast<int>(blocklist_.blocks_S1D_block.size()));
        const AxialHFBBlock& block_ = blocklist_.blocks_S1D_block[block_I];
        assert(bqp_I >= 0 && bqp_I < block_.V_F2D_bsp_bqp.cols());
        assert(bqp_I < block_.U_F2D_bsp_bqp.cols());
        blockedV_F1D_bsp = block_.V_F2D_bsp_bqp.col(bqp_I);
        blockedU_F1D_bsp = block_.U_F2D_bsp_bqp.col(bqp_I);
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
    Eigen::MatrixXd rho_F2D_z_r; // ρ(z,r_⊥) [fm⁻³].
    Eigen::MatrixXd tau_F2D_z_r; // τ(z,r_⊥) [fm⁻⁵].
    Eigen::MatrixXd kappa_F2D_z_r; // κ(z,r_⊥) [fm⁻³].
    Eigen::MatrixXd rhoD2_F2D_z_r; // Δρ(z,r_⊥) [fm⁻⁵].
    Eigen::MatrixXd rhoDr_F2D_z_r; // ∂_{r⊥}ρ(z,r_⊥) [fm⁻⁴].
    Eigen::MatrixXd rhoDz_F2D_z_r; // ∂_zρ(z,r_⊥) [fm⁻⁴].
    Eigen::MatrixXd Jzphi_F2D_z_r; // J_{zφ}(z,r_⊥) [fm⁻⁴].
    Eigen::MatrixXd Jphiz_F2D_z_r; // J_{φz}(z,r_⊥) [fm⁻⁴].
    Eigen::MatrixXd Jphir_F2D_z_r; // J_{φr}(z,r_⊥) [fm⁻⁴].
    Eigen::MatrixXd Jrphi_F2D_z_r; // J_{rφ}(z,r_⊥) [fm⁻⁴].
    Eigen::MatrixXd dJ_F2D_z_r; // ∇·J(z,r_⊥) [fm⁻⁵].

public:
    /**
     * @brief Construct zero coordinate-space densities.
     * @math (N_z,N_r) → 0_{N_z×N_r}
     * @output Zero-initialized density grids.
     */
    explicit AxialHFBDensity(const AxialConfig& axialconfig_)
    : rho_F2D_z_r(Eigen::MatrixXd::Zero(axialconfig_.Nz_I, axialconfig_.Nr_I)), tau_F2D_z_r(Eigen::MatrixXd::Zero(axialconfig_.Nz_I, axialconfig_.Nr_I)), kappa_F2D_z_r(Eigen::MatrixXd::Zero(axialconfig_.Nz_I, axialconfig_.Nr_I)), rhoD2_F2D_z_r(Eigen::MatrixXd::Zero(axialconfig_.Nz_I, axialconfig_.Nr_I)), rhoDr_F2D_z_r(Eigen::MatrixXd::Zero(axialconfig_.Nz_I, axialconfig_.Nr_I)), rhoDz_F2D_z_r(Eigen::MatrixXd::Zero(axialconfig_.Nz_I, axialconfig_.Nr_I)), Jzphi_F2D_z_r(Eigen::MatrixXd::Zero(axialconfig_.Nz_I, axialconfig_.Nr_I)), Jphiz_F2D_z_r(Eigen::MatrixXd::Zero(axialconfig_.Nz_I, axialconfig_.Nr_I)), Jphir_F2D_z_r(Eigen::MatrixXd::Zero(axialconfig_.Nz_I, axialconfig_.Nr_I)), Jrphi_F2D_z_r(Eigen::MatrixXd::Zero(axialconfig_.Nz_I, axialconfig_.Nr_I)), dJ_F2D_z_r(Eigen::MatrixXd::Zero(axialconfig_.Nz_I, axialconfig_.Nr_I)) {}

    /**
     * @brief Accumulate one block into coordinate space.
     * @math (ρ_{ij},κ_{ij},Φ_i,Φ_j) → (ρ,τ,κ,Δρ,∇ρ,J)
     * @output Accumulated density grids.
     */
    void add_density(const AxialBasis& global_basis_, const AxialHFBBlock& block_);

    /**
     * @brief Rebuild densities from all species blocks.
     * @math {B_{block}} → (ρ,τ,κ,Δρ,∇ρ,J)
     * @output Updated density grids.
     */
    void update_density(const AxialBasis& global_basis_, const AxialHFBBlockList& blocklist_);
};

/**
 * @brief Store coupled neutron and proton densities.
 */
class AxialHFBDensitySystem {
public:
    AxialHFBDensity density_p; // Proton densities.
    AxialHFBDensity density_n; // Neutron densities.

public:
    /**
     * @brief Construct proton and neutron density grids.
     * @math C_{axial} → (D_p,D_n)
     * @output Zero-initialized species densities.
     */
    explicit AxialHFBDensitySystem(const AxialConfig& axialconfig_)
    : density_p(axialconfig_), density_n(axialconfig_) {}
};

/**
 * @brief Store coordinate-space fields for one species.
 */
class AxialHFBField {
public:
    Eigen::MatrixXd vcent_F2D_z_r; // δE/δρ(z,r_⊥) [MeV].
    Eigen::MatrixXd vmass_F2D_z_r; // δE/δτ(z,r_⊥) [MeV fm²].
    Eigen::MatrixXd vpair_F2D_z_r; // δE/δκ(z,r_⊥) [MeV].
    Eigen::MatrixXd vD2_F2D_z_r; // δE/δΔρ(z,r_⊥) [MeV fm²].
    Eigen::MatrixXd vDr_F2D_z_r; // δE/δ∂_{r⊥}ρ(z,r_⊥) [MeV fm].
    Eigen::MatrixXd vDz_F2D_z_r; // δE/δ∂_zρ(z,r_⊥) [MeV fm].
    Eigen::MatrixXd vJzphi_F2D_z_r; // δE/δJ_{zφ}(z,r_⊥) [MeV fm].
    Eigen::MatrixXd vJphiz_F2D_z_r; // δE/δJ_{φz}(z,r_⊥) [MeV fm].
    Eigen::MatrixXd vJphir_F2D_z_r; // δE/δJ_{φr}(z,r_⊥) [MeV fm].
    Eigen::MatrixXd vJrphi_F2D_z_r; // δE/δJ_{rφ}(z,r_⊥) [MeV fm].
    Eigen::MatrixXd vdJ_F2D_z_r; // δE/δ∇·J(z,r_⊥) [MeV fm²].

public:
    /**
     * @brief Construct zero coordinate-space fields.
     * @math (N_z,N_r) → 0_{N_z×N_r}
     * @output Zero-initialized field grids.
     */
    explicit AxialHFBField(const AxialConfig& axialconfig_)
    : vcent_F2D_z_r(Eigen::MatrixXd::Zero(axialconfig_.Nz_I, axialconfig_.Nr_I)), vmass_F2D_z_r(Eigen::MatrixXd::Zero(axialconfig_.Nz_I, axialconfig_.Nr_I)), vpair_F2D_z_r(Eigen::MatrixXd::Zero(axialconfig_.Nz_I, axialconfig_.Nr_I)), vD2_F2D_z_r(Eigen::MatrixXd::Zero(axialconfig_.Nz_I, axialconfig_.Nr_I)), vDr_F2D_z_r(Eigen::MatrixXd::Zero(axialconfig_.Nz_I, axialconfig_.Nr_I)), vDz_F2D_z_r(Eigen::MatrixXd::Zero(axialconfig_.Nz_I, axialconfig_.Nr_I)), vJzphi_F2D_z_r(Eigen::MatrixXd::Zero(axialconfig_.Nz_I, axialconfig_.Nr_I)), vJphiz_F2D_z_r(Eigen::MatrixXd::Zero(axialconfig_.Nz_I, axialconfig_.Nr_I)), vJphir_F2D_z_r(Eigen::MatrixXd::Zero(axialconfig_.Nz_I, axialconfig_.Nr_I)), vJrphi_F2D_z_r(Eigen::MatrixXd::Zero(axialconfig_.Nz_I, axialconfig_.Nr_I)), vdJ_F2D_z_r(Eigen::MatrixXd::Zero(axialconfig_.Nz_I, axialconfig_.Nr_I)) {}

    /**
     * @brief Calculate the regularized pairing coupling.
     * @math (E_{cut},λ,g) → g_{reg}
     * @output Regularized coupling grid.
     */
    void regularize_gr(double EspCut_F, double lambda_F, Eigen::MatrixXd& gr_F2D_z_r) const;
};

/**
 * @brief Manage coupled neutron and proton fields.
 */
class AxialHFBFieldSystem {
public:
    AxialHFBField field_p; // Proton fields.
    AxialHFBField field_n; // Neutron fields.
    bool isCoulombKernelBuilt_B; // Coulomb kernel readiness.
    Eigen::Tensor<double, 4, Eigen::ColMajor> coulomb_F4D_zs_rs_zt_rt; // K_C(z_s,r_s,z_t,r_t).

public:
    /**
     * @brief Construct species fields and Coulomb storage.
     * @math C_{axial} → (F_p,F_n,K_C)
     * @output Zero-initialized fields and kernel.
     */
    AxialHFBFieldSystem(const AxialConfig& axialconfig_, const HFBSettings&, const EDFParamsSkyrme&)
    : field_p(axialconfig_), field_n(axialconfig_), isCoulombKernelBuilt_B(false), coulomb_F4D_zs_rs_zt_rt(axialconfig_.Nz_I, axialconfig_.Nr_I, axialconfig_.Nz_I, axialconfig_.Nr_I) {}

    /**
     * @brief Build the axial direct-Coulomb kernel.
     * @math (Φ,e²,reflection) → K_C
     * @output Initialized Coulomb kernel.
     */
    void build_coulomb_kernel(const AxialBasis& basis_, const EDFParamsSkyrme& edf_skyrme_, bool useReflection_B);

    /**
     * @brief Add direct and exchange Coulomb fields.
     * @math ρ_p → v_C^{dir}+v_C^{Slater}
     * @output Updated proton central field.
     */
    void add_coulomb_field(const AxialHFBDensitySystem& densities_, const EDFParamsSkyrme& edf_skyrme_, const HFBSettings& hfbsettings_);

    /**
     * @brief Add local pairing fields.
     * @math (ρ_0,κ_q,λ_q) → Δ_q
     * @output Updated neutron and proton pairing fields.
     */
    void add_pairing_fields(const AxialHFBDensitySystem& densities_, const EDFParamsSkyrme& edf_skyrme_, const HFBSettings& hfbsettings_, double lambda_n_F, double lambda_p_F);

    /**
     * @brief Rebuild local neutron and proton fields.
     * @math (D_n,D_p,C,S) → (F_n,F_p)
     * @output Updated coordinate-space fields.
     */
    void update_nuclei_fields(const AxialHFBDensitySystem& densities_, const EDFParamsSkyrme& edf_skyrme_, const HFBSettings& hfbsettings_);
};

/**
 * @brief Own one complete axial HFB calculation.
 */
class AxialHFB {
public:
    AxialConfig axialconfig; // Axial basis configuration.
    HFBSettings hfbsettings; // HFB iteration controls.
    EDFParamsSkyrme edf_skyrme; // Base local EDF.
    AxialGaussianGogny gogny; // Finite-range Gogny interaction.
    AxialGaussianCoulomb coulomb; // Finite-range Coulomb interaction.
    AxialBasis global_basis; // Global axial basis.
    AxialHFBBlockListSystem blocks; // q ∈ {n,p}: HFB blocks.
    AxialHFBDensitySystem densities; // q ∈ {n,p}: local densities.
    AxialHFBFieldSystem fields; // q ∈ {n,p}: local fields.

public:
    /**
     * @brief Construct an axial HFB solver with D1S.
     * @math (C_{axial},P_{HFB},C_{Skyrme}) → HFB[D1S]
     * @output Initialized HFB solver.
     */
    AxialHFB(const AxialConfig& axialconfig_, const HFBSettings& hfbsettings_, const EDFParamsSkyrme& edf_skyrme_)
    : AxialHFB(axialconfig_, hfbsettings_, edf_skyrme_, EDFParamsGogny::D1S()) {}

    /**
     * @brief Construct an axial HFB solver from Gogny parameters.
     * @math (C_{axial},P_{HFB},C_{Skyrme},P_G) → HFB
     * @output Initialized HFB solver.
     */
    AxialHFB(const AxialConfig& axialconfig_, const HFBSettings& hfbsettings_, const EDFParamsSkyrme& edf_skyrme_, const EDFParamsGogny& edf_gogny_)
    : AxialHFB(axialconfig_, hfbsettings_, edf_skyrme_, AxialGaussianGogny(axialconfig_, edf_gogny_.forceName_Str, edf_gogny_.mu_F1D_g, edf_gogny_.W_F1D_g, edf_gogny_.B_F1D_g, edf_gogny_.H_F1D_g, edf_gogny_.M_F1D_g)) {}

    /**
     * @brief Construct an axial HFB solver from interactions.
     * @math (C_{axial},P_{HFB},C_{Skyrme},V_G) → HFB
     * @output Initialized HFB solver.
     */
    AxialHFB(const AxialConfig& axialconfig_, const HFBSettings& hfbsettings_, const EDFParamsSkyrme& edf_skyrme_, AxialGaussianGogny gogny_)
    : axialconfig(axialconfig_), hfbsettings(hfbsettings_), edf_skyrme(edf_skyrme_), gogny(std::move(gogny_)), coulomb(axialconfig_, edf_skyrme_.e2charg_F), global_basis(axialconfig_, axialconfig_.labels_S1D_sp), blocks(axialconfig_, hfbsettings_), densities(axialconfig_), fields(axialconfig_, hfbsettings_, edf_skyrme_) {}

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
};
