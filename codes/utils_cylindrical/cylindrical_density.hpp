/**
 * @file    cylindrical_density.hpp
 * @author  cailea
 * @date    2026-09-11
 * @brief   Cylindrical coordinate-space density contractions.
 */

#pragma once

#include <numbers>
#include <vector>

#include <Eigen/Core>

#include "cylindrical_basis.hpp"
#include "cylindrical_config.hpp"

class CylindricalDensity {
public:
    const CylindricalBasis2D& basis;

    std::vector<std::vector<CylindricalSPLabel>> labels_S2D_block_bsp{};
    std::vector<std::vector<int>> indices_I2D_block_bsp{};
    std::vector<std::vector<int>> bup_I2D_block_bup{}; // bsp(block,bup): 2Σ = +1.
    std::vector<std::vector<int>> bdn_I2D_block_bdn{}; // bsp(block,bdn): 2Σ = -1.

    // Im: real density from Im(ρ).
    Eigen::MatrixXd rho_F2D_z_r{};
    Eigen::MatrixXd tau_F2D_z_r{};
    Eigen::MatrixXd rhoLaplacian_F2D_z_r{};
    Eigen::MatrixXd rhoDr_F2D_z_r{};
    Eigen::MatrixXd rhoDz_F2D_z_r{};
    Eigen::MatrixXd jz_F2D_z_r{};
    Eigen::MatrixXd jr_F2D_z_r{};
    Eigen::MatrixXd jphi_F2D_z_r{};
    Eigen::MatrixXd sz_F2D_z_r{};
    Eigen::MatrixXd sr_F2D_z_r{};
    Eigen::MatrixXd sphi_F2D_z_r{};
    Eigen::MatrixXd Jphiz_F2D_z_r{};
    Eigen::MatrixXd Jzphi_F2D_z_r{};
    Eigen::MatrixXd Jphir_F2D_z_r{};
    Eigen::MatrixXd Jrphi_F2D_z_r{};
    Eigen::MatrixXd Jrz_F2D_z_r{};
    Eigen::MatrixXd Jzr_F2D_z_r{};
    Eigen::MatrixXd Jzz_F2D_z_r{};
    Eigen::MatrixXd Jrr_F2D_z_r{};
    Eigen::MatrixXd Jphiphi_F2D_z_r{};
    Eigen::MatrixXd dJ_F2D_z_r{};
    Eigen::MatrixXd kappa_F2D_z_r{};

public:
    /**
     * @brief  Copy configured blocks and allocate density buffers.
     * @math   α_(b,i) → sp(b,i); D(z,r) = 0.
     * @output Labels, basis indices, spin groups, and zeroed buffers.
     * @note   The referenced basis must outlive this object.
     */
    CylindricalDensity(const CylindricalSetting& cylindricalsetting_, const CylindricalBasis2D& basis_)
    : basis(basis_) {
        // α_(b,i), sp(b,i) ← setting.
        labels_S2D_block_bsp = cylindricalsetting_.labels_S2D_block_bsp;
        indices_I2D_block_bsp = cylindricalsetting_.indices_I2D_block_bsp;

        // α_(block,bsp) → {bsp_↑,bsp_↓}.
        bup_I2D_block_bup.resize(labels_S2D_block_bsp.size());
        bdn_I2D_block_bdn.resize(labels_S2D_block_bsp.size());
        for (int block_I = 0; block_I < static_cast<int>(labels_S2D_block_bsp.size()); ++block_I) {
            for (int bsp_I = 0; bsp_I < static_cast<int>(labels_S2D_block_bsp[block_I].size()); ++bsp_I) {
                if (labels_S2D_block_bsp[block_I][bsp_I].twoSigma_I == 1) {bup_I2D_block_bup[block_I].push_back(bsp_I);}
                if (labels_S2D_block_bsp[block_I][bsp_I].twoSigma_I == -1) {bdn_I2D_block_bdn[block_I].push_back(bsp_I);}
            }
        }

        // D(z,r) ∈ ℝ^{Nz×Nr}.
        rho_F2D_z_r.resize(cylindricalsetting_.Nz_I, cylindricalsetting_.Nr_I);
        tau_F2D_z_r.resize(cylindricalsetting_.Nz_I, cylindricalsetting_.Nr_I);
        rhoLaplacian_F2D_z_r.resize(cylindricalsetting_.Nz_I, cylindricalsetting_.Nr_I);
        rhoDr_F2D_z_r.resize(cylindricalsetting_.Nz_I, cylindricalsetting_.Nr_I);
        rhoDz_F2D_z_r.resize(cylindricalsetting_.Nz_I, cylindricalsetting_.Nr_I);
        jz_F2D_z_r.resize(cylindricalsetting_.Nz_I, cylindricalsetting_.Nr_I);
        jr_F2D_z_r.resize(cylindricalsetting_.Nz_I, cylindricalsetting_.Nr_I);
        jphi_F2D_z_r.resize(cylindricalsetting_.Nz_I, cylindricalsetting_.Nr_I);
        sz_F2D_z_r.resize(cylindricalsetting_.Nz_I, cylindricalsetting_.Nr_I);
        sr_F2D_z_r.resize(cylindricalsetting_.Nz_I, cylindricalsetting_.Nr_I);
        sphi_F2D_z_r.resize(cylindricalsetting_.Nz_I, cylindricalsetting_.Nr_I);
        Jphiz_F2D_z_r.resize(cylindricalsetting_.Nz_I, cylindricalsetting_.Nr_I);
        Jzphi_F2D_z_r.resize(cylindricalsetting_.Nz_I, cylindricalsetting_.Nr_I);
        Jphir_F2D_z_r.resize(cylindricalsetting_.Nz_I, cylindricalsetting_.Nr_I);
        Jrphi_F2D_z_r.resize(cylindricalsetting_.Nz_I, cylindricalsetting_.Nr_I);
        Jrz_F2D_z_r.resize(cylindricalsetting_.Nz_I, cylindricalsetting_.Nr_I);
        Jzr_F2D_z_r.resize(cylindricalsetting_.Nz_I, cylindricalsetting_.Nr_I);
        Jzz_F2D_z_r.resize(cylindricalsetting_.Nz_I, cylindricalsetting_.Nr_I);
        Jrr_F2D_z_r.resize(cylindricalsetting_.Nz_I, cylindricalsetting_.Nr_I);
        Jphiphi_F2D_z_r.resize(cylindricalsetting_.Nz_I, cylindricalsetting_.Nr_I);
        dJ_F2D_z_r.resize(cylindricalsetting_.Nz_I, cylindricalsetting_.Nr_I);
        kappa_F2D_z_r.resize(cylindricalsetting_.Nz_I, cylindricalsetting_.Nr_I);

        set_zero();
    }

    /**
     * @brief  Reset coordinate-space density grids.
     * @math   D(z,r) → 0.
     * @output Zeroed coordinate-space density grids.
     */
    void set_zero();

    /**
     * @brief  Contract block matrices with cylindrical basis functions.
     * @math   {ρ⁺⁺,ρ⁻⁻,κ⁺⁻} → {ρ,τ,∇ρ,∇²ρ,j,s,J,∇·J,κ}.
     * @math   κ⁻⁺ = −(κ⁺⁻)ᵀ.
     * @output Rebuilt grids summed over all blocks.
     * @note   Real symmetric ρ inputs; imaginary-part contributions vanish.
     */
    void update_density(const std::vector<Eigen::MatrixXd>& rhoPosPos_F3D_block_bsp_bsp, const std::vector<Eigen::MatrixXd>& rhoNegNeg_F3D_block_bsp_bsp, const std::vector<Eigen::MatrixXd>& kappaPosNeg_F3D_block_bsp_bsp);
};

/**
 * @brief  Reset coordinate-space density grids.
 * @math   D(z,r) → 0.
 * @output Zeroed coordinate-space density grids.
 */
inline void CylindricalDensity::set_zero() {
    // {ρ,τ,∇ρ,∇²ρ,j,s,J,∇·J,κ} = 0.
    rho_F2D_z_r.setZero();
    tau_F2D_z_r.setZero();
    rhoLaplacian_F2D_z_r.setZero();
    rhoDr_F2D_z_r.setZero();
    rhoDz_F2D_z_r.setZero();
    jz_F2D_z_r.setZero();
    jr_F2D_z_r.setZero();
    jphi_F2D_z_r.setZero();
    sz_F2D_z_r.setZero();
    sr_F2D_z_r.setZero();
    sphi_F2D_z_r.setZero();
    Jphiz_F2D_z_r.setZero();
    Jzphi_F2D_z_r.setZero();
    Jphir_F2D_z_r.setZero();
    Jrphi_F2D_z_r.setZero();
    Jrz_F2D_z_r.setZero();
    Jzr_F2D_z_r.setZero();
    Jzz_F2D_z_r.setZero();
    Jrr_F2D_z_r.setZero();
    Jphiphi_F2D_z_r.setZero();
    dJ_F2D_z_r.setZero();
    kappa_F2D_z_r.setZero();
}

/**
 * @brief  Contract real block matrices over spin groups.
 * @math   (r,z,block,Σ_α,Σ_β) → (bsp_α,bsp_β,sp_α,sp_β).
 * @output Rebuilt coordinate-space densities from all four spin groups.
 */
inline void CylindricalDensity::update_density(const std::vector<Eigen::MatrixXd>& rhoPosPos_F3D_block_bsp_bsp, const std::vector<Eigen::MatrixXd>& rhoNegNeg_F3D_block_bsp_bsp, const std::vector<Eigen::MatrixXd>& kappaPosNeg_F3D_block_bsp_bsp) {
    set_zero();

    // Σ_α = Σ_β = +1/2.
    const auto add_upup_Func = [&](int r_I, int z_I, int block_I) {
        const auto& rhoPosPos_F2D_bsp_bsp = rhoPosPos_F3D_block_bsp_bsp[block_I];
        const auto& rhoNegNeg_F2D_bsp_bsp = rhoNegNeg_F3D_block_bsp_bsp[block_I];
        const auto& kappaPosNeg_F2D_bsp_bsp = kappaPosNeg_F3D_block_bsp_bsp[block_I];

        const double rInv_F = 1.0 / basis.r_F1D_r(r_I);
        const double rInv2_F = rInv_F * rInv_F;
        const double norm_F = 1.0 / (2.0 * std::numbers::pi);

        for (int bUp1_I : bup_I2D_block_bup[block_I]) {
            const int spUp1_I = indices_I2D_block_bsp[block_I][bUp1_I];
            const double Lambda_F = labels_S2D_block_bsp[block_I][bUp1_I].Lambda_I;
            const double phi1_F = basis.phi_F3D_sp_z_r(spUp1_I, z_I, r_I);
            const double dphi1dr_F = basis.dphidr_F3D_sp_z_r(spUp1_I, z_I, r_I);
            const double dphi1dz_F = basis.dphidz_F3D_sp_z_r(spUp1_I, z_I, r_I);

            // L_α = (∂²_z + ∂²_r + r⁻¹∂_r − Λ²r⁻²)φ_α.
            const double L1_F = basis.ddphidz_F3D_sp_z_r(spUp1_I, z_I, r_I) + basis.ddphidr_F3D_sp_z_r(spUp1_I, z_I, r_I) + rInv_F * dphi1dr_F - Lambda_F * Lambda_F * rInv2_F * phi1_F;

            for (int bUp2_I : bup_I2D_block_bup[block_I]) {
                const int spUp2_I = indices_I2D_block_bsp[block_I][bUp2_I];
                const double phi2_F = basis.phi_F3D_sp_z_r(spUp2_I, z_I, r_I);
                const double dphi2dr_F = basis.dphidr_F3D_sp_z_r(spUp2_I, z_I, r_I);
                const double dphi2dz_F = basis.dphidz_F3D_sp_z_r(spUp2_I, z_I, r_I);

                // Λ_β = Λ_α within the up-spin group.
                const double L2_F = basis.ddphidz_F3D_sp_z_r(spUp2_I, z_I, r_I) + basis.ddphidr_F3D_sp_z_r(spUp2_I, z_I, r_I) + rInv_F * dphi2dr_F - Lambda_F * Lambda_F * rInv2_F * phi2_F;

                // R± = ρ⁺⁺ ± ρ⁻⁻.
                const double RPlus_F = rhoPosPos_F2D_bsp_bsp(bUp1_I, bUp2_I) + rhoNegNeg_F2D_bsp_bsp(bUp1_I, bUp2_I);
                const double RMinus_F = rhoPosPos_F2D_bsp_bsp(bUp1_I, bUp2_I) - rhoNegNeg_F2D_bsp_bsp(bUp1_I, bUp2_I);
                const double rhoKernel_F = phi1_F * phi2_F;

                // τ_αβ = ∂_zφ_α∂_zφ_β + ∂_rφ_α∂_rφ_β + Λ²r⁻²φ_αφ_β.
                const double tauKernel_F = dphi1dz_F * dphi2dz_F + dphi1dr_F * dphi2dr_F + Lambda_F * Lambda_F * rInv2_F * rhoKernel_F;
                const double rhoLaplacianKernel_F = L1_F * phi2_F + phi1_F * L2_F + 2.0 * tauKernel_F;

                // ρ, τ, ∇²ρ ← R⁺.
                rho_F2D_z_r(z_I, r_I) += norm_F * RPlus_F * rhoKernel_F;
                tau_F2D_z_r(z_I, r_I) += norm_F * RPlus_F * tauKernel_F;
                rhoLaplacian_F2D_z_r(z_I, r_I) += norm_F * RPlus_F * rhoLaplacianKernel_F;

                // ∇ρ: ∂(φ_αφ_β) = (∂φ_α)φ_β + φ_α(∂φ_β).
                rhoDr_F2D_z_r(z_I, r_I) += norm_F * RPlus_F * (dphi1dr_F * phi2_F + phi1_F * dphi2dr_F);
                rhoDz_F2D_z_r(z_I, r_I) += norm_F * RPlus_F * (dphi1dz_F * phi2_F + phi1_F * dphi2dz_F);

                // j_φ, s_z ← R⁻; J_φz ← R⁺.
                jphi_F2D_z_r(z_I, r_I) += norm_F * RMinus_F * Lambda_F * rInv_F * rhoKernel_F;
                sz_F2D_z_r(z_I, r_I) += norm_F * RMinus_F * rhoKernel_F;
                Jphiz_F2D_z_r(z_I, r_I) += norm_F * RPlus_F * Lambda_F * rInv_F * rhoKernel_F;

                // ∇·J: 2 Re[(∇_rφ_β†)σ_z(∇_φφ_α)/i].
                dJ_F2D_z_r(z_I, r_I) += 2.0 * norm_F * RPlus_F * Lambda_F * rInv_F * dphi2dr_F * phi1_F;

                // κ: representative ↑↑ contribution.
                kappa_F2D_z_r(z_I, r_I) += norm_F * kappaPosNeg_F2D_bsp_bsp(bUp1_I, bUp2_I) * rhoKernel_F;
            }
        }
    };

    // Σ_α = Σ_β = -1/2.
    const auto add_dndn_Func = [&](int r_I, int z_I, int block_I) {
        const auto& rhoPosPos_F2D_bsp_bsp = rhoPosPos_F3D_block_bsp_bsp[block_I];
        const auto& rhoNegNeg_F2D_bsp_bsp = rhoNegNeg_F3D_block_bsp_bsp[block_I];
        const auto& kappaPosNeg_F2D_bsp_bsp = kappaPosNeg_F3D_block_bsp_bsp[block_I];

        const double rInv_F = 1.0 / basis.r_F1D_r(r_I);
        const double rInv2_F = rInv_F * rInv_F;
        const double norm_F = 1.0 / (2.0 * std::numbers::pi);

        for (int bDown1_I : bdn_I2D_block_bdn[block_I]) {
            const int spDown1_I = indices_I2D_block_bsp[block_I][bDown1_I];
            const double Lambda_F = labels_S2D_block_bsp[block_I][bDown1_I].Lambda_I;
            const double phi1_F = basis.phi_F3D_sp_z_r(spDown1_I, z_I, r_I);
            const double dphi1dr_F = basis.dphidr_F3D_sp_z_r(spDown1_I, z_I, r_I);
            const double dphi1dz_F = basis.dphidz_F3D_sp_z_r(spDown1_I, z_I, r_I);

            // L_α = (∂²_z + ∂²_r + r⁻¹∂_r − Λ²r⁻²)φ_α.
            const double L1_F = basis.ddphidz_F3D_sp_z_r(spDown1_I, z_I, r_I) + basis.ddphidr_F3D_sp_z_r(spDown1_I, z_I, r_I) + rInv_F * dphi1dr_F - Lambda_F * Lambda_F * rInv2_F * phi1_F;

            for (int bDown2_I : bdn_I2D_block_bdn[block_I]) {
                const int spDown2_I = indices_I2D_block_bsp[block_I][bDown2_I];
                const double phi2_F = basis.phi_F3D_sp_z_r(spDown2_I, z_I, r_I);
                const double dphi2dr_F = basis.dphidr_F3D_sp_z_r(spDown2_I, z_I, r_I);
                const double dphi2dz_F = basis.dphidz_F3D_sp_z_r(spDown2_I, z_I, r_I);

                // Λ_β = Λ_α within the down-spin group.
                const double L2_F = basis.ddphidz_F3D_sp_z_r(spDown2_I, z_I, r_I) + basis.ddphidr_F3D_sp_z_r(spDown2_I, z_I, r_I) + rInv_F * dphi2dr_F - Lambda_F * Lambda_F * rInv2_F * phi2_F;

                // R± = ρ⁺⁺ ± ρ⁻⁻.
                const double RPlus_F = rhoPosPos_F2D_bsp_bsp(bDown1_I, bDown2_I) + rhoNegNeg_F2D_bsp_bsp(bDown1_I, bDown2_I);
                const double RMinus_F = rhoPosPos_F2D_bsp_bsp(bDown1_I, bDown2_I) - rhoNegNeg_F2D_bsp_bsp(bDown1_I, bDown2_I);
                const double rhoKernel_F = phi1_F * phi2_F;

                // τ_αβ = ∂_zφ_α∂_zφ_β + ∂_rφ_α∂_rφ_β + Λ²r⁻²φ_αφ_β.
                const double tauKernel_F = dphi1dz_F * dphi2dz_F + dphi1dr_F * dphi2dr_F + Lambda_F * Lambda_F * rInv2_F * rhoKernel_F;
                const double rhoLaplacianKernel_F = L1_F * phi2_F + phi1_F * L2_F + 2.0 * tauKernel_F;

                // ρ, τ, ∇²ρ ← R⁺.
                rho_F2D_z_r(z_I, r_I) += norm_F * RPlus_F * rhoKernel_F;
                tau_F2D_z_r(z_I, r_I) += norm_F * RPlus_F * tauKernel_F;
                rhoLaplacian_F2D_z_r(z_I, r_I) += norm_F * RPlus_F * rhoLaplacianKernel_F;

                // ∇ρ: ∂(φ_αφ_β) = (∂φ_α)φ_β + φ_α(∂φ_β).
                rhoDr_F2D_z_r(z_I, r_I) += norm_F * RPlus_F * (dphi1dr_F * phi2_F + phi1_F * dphi2dr_F);
                rhoDz_F2D_z_r(z_I, r_I) += norm_F * RPlus_F * (dphi1dz_F * phi2_F + phi1_F * dphi2dz_F);

                // j_φ ← R⁻; s_z ← −R⁻; J_φz ← −R⁺.
                jphi_F2D_z_r(z_I, r_I) += norm_F * RMinus_F * Lambda_F * rInv_F * rhoKernel_F;
                sz_F2D_z_r(z_I, r_I) -= norm_F * RMinus_F * rhoKernel_F;
                Jphiz_F2D_z_r(z_I, r_I) -= norm_F * RPlus_F * Lambda_F * rInv_F * rhoKernel_F;

                // ∇·J: 2 Re[(∇_rφ_β†)σ_z(∇_φφ_α)/i].
                dJ_F2D_z_r(z_I, r_I) -= 2.0 * norm_F * RPlus_F * Lambda_F * rInv_F * dphi2dr_F * phi1_F;

                // κ: representative ↓↓ contribution.
                kappa_F2D_z_r(z_I, r_I) -= norm_F * kappaPosNeg_F2D_bsp_bsp(bDown1_I, bDown2_I) * rhoKernel_F;
            }
        }
    };

    // Σ_α = +1/2, Σ_β = -1/2.
    const auto add_updn_Func = [&](int r_I, int z_I, int block_I) {
        const auto& rhoPosPos_F2D_bsp_bsp = rhoPosPos_F3D_block_bsp_bsp[block_I];
        const auto& rhoNegNeg_F2D_bsp_bsp = rhoNegNeg_F3D_block_bsp_bsp[block_I];

        const double rInv_F = 1.0 / basis.r_F1D_r(r_I);
        const double norm_F = 1.0 / (2.0 * std::numbers::pi);

        for (int bUp_I : bup_I2D_block_bup[block_I]) {
            const int spUp_I = indices_I2D_block_bsp[block_I][bUp_I];
            const double LambdaUp_F = labels_S2D_block_bsp[block_I][bUp_I].Lambda_I;
            const double phiUp_F = basis.phi_F3D_sp_z_r(spUp_I, z_I, r_I);
            const double dphiUpdr_F = basis.dphidr_F3D_sp_z_r(spUp_I, z_I, r_I);
            const double dphiUpdz_F = basis.dphidz_F3D_sp_z_r(spUp_I, z_I, r_I);

            for (int bDown_I : bdn_I2D_block_bdn[block_I]) {
                const int spDown_I = indices_I2D_block_bsp[block_I][bDown_I];
                const double LambdaDown_F = labels_S2D_block_bsp[block_I][bDown_I].Lambda_I;
                const double phiDown_F = basis.phi_F3D_sp_z_r(spDown_I, z_I, r_I);
                const double dphiDowndr_F = basis.dphidr_F3D_sp_z_r(spDown_I, z_I, r_I);
                const double dphiDowndz_F = basis.dphidz_F3D_sp_z_r(spDown_I, z_I, r_I);

                // R±_↑↓ = ρ⁺⁺_↑↓ ± ρ⁻⁻_↑↓.
                const double RPlus_F = rhoPosPos_F2D_bsp_bsp(bUp_I, bDown_I) + rhoNegNeg_F2D_bsp_bsp(bUp_I, bDown_I);
                const double RMinus_F = rhoPosPos_F2D_bsp_bsp(bUp_I, bDown_I) - rhoNegNeg_F2D_bsp_bsp(bUp_I, bDown_I);

                const double rhoKernel_F = phiUp_F * phiDown_F;
                const double JrphiKernel_F = phiDown_F * dphiUpdr_F - phiUp_F * dphiDowndr_F;
                const double JzphiKernel_F = phiDown_F * dphiUpdz_F - phiUp_F * dphiDowndz_F;

                // s_r ← R⁺; J_{rφ}, J_{zφ}, J_{φr} ← R⁻.
                sr_F2D_z_r(z_I, r_I) += norm_F * RPlus_F * rhoKernel_F;
                Jrphi_F2D_z_r(z_I, r_I) += 0.5 * norm_F * RMinus_F * JrphiKernel_F;
                Jzphi_F2D_z_r(z_I, r_I) += 0.5 * norm_F * RMinus_F * JzphiKernel_F;
                Jphir_F2D_z_r(z_I, r_I) += 0.5 * norm_F * RMinus_F * (LambdaUp_F + LambdaDown_F) * rInv_F * rhoKernel_F;

                // ∇·J: angular-bra and σ_φ kernels, including 2 Re.
                dJ_F2D_z_r(z_I, r_I) += 2.0 * norm_F * RMinus_F * (-LambdaDown_F * rInv_F * phiDown_F * dphiUpdz_F + dphiDowndz_F * dphiUpdr_F);
            }
        }
    };

    // Σ_α = -1/2, Σ_β = +1/2.
    const auto add_dnup_Func = [&](int r_I, int z_I, int block_I) {
        const auto& rhoPosPos_F2D_bsp_bsp = rhoPosPos_F3D_block_bsp_bsp[block_I];
        const auto& rhoNegNeg_F2D_bsp_bsp = rhoNegNeg_F3D_block_bsp_bsp[block_I];

        const double rInv_F = 1.0 / basis.r_F1D_r(r_I);
        const double norm_F = 1.0 / (2.0 * std::numbers::pi);

        for (int bDown_I : bdn_I2D_block_bdn[block_I]) {
            const int spDown_I = indices_I2D_block_bsp[block_I][bDown_I];
            const double LambdaDown_F = labels_S2D_block_bsp[block_I][bDown_I].Lambda_I;
            const double phiDown_F = basis.phi_F3D_sp_z_r(spDown_I, z_I, r_I);
            const double dphiDowndr_F = basis.dphidr_F3D_sp_z_r(spDown_I, z_I, r_I);
            const double dphiDowndz_F = basis.dphidz_F3D_sp_z_r(spDown_I, z_I, r_I);

            for (int bUp_I : bup_I2D_block_bup[block_I]) {
                const int spUp_I = indices_I2D_block_bsp[block_I][bUp_I];
                const double LambdaUp_F = labels_S2D_block_bsp[block_I][bUp_I].Lambda_I;
                const double phiUp_F = basis.phi_F3D_sp_z_r(spUp_I, z_I, r_I);
                const double dphiUpdr_F = basis.dphidr_F3D_sp_z_r(spUp_I, z_I, r_I);
                const double dphiUpdz_F = basis.dphidz_F3D_sp_z_r(spUp_I, z_I, r_I);

                // R±_↓↑ = ρ⁺⁺_↓↑ ± ρ⁻⁻_↓↑.
                const double RPlus_F = rhoPosPos_F2D_bsp_bsp(bDown_I, bUp_I) + rhoNegNeg_F2D_bsp_bsp(bDown_I, bUp_I);
                const double RMinus_F = rhoPosPos_F2D_bsp_bsp(bDown_I, bUp_I) - rhoNegNeg_F2D_bsp_bsp(bDown_I, bUp_I);

                const double rhoKernel_F = phiDown_F * phiUp_F;
                const double JrphiKernel_F = phiUp_F * dphiDowndr_F - phiDown_F * dphiUpdr_F;
                const double JzphiKernel_F = phiUp_F * dphiDowndz_F - phiDown_F * dphiUpdz_F;

                // σ_φ: ↓↑ contributes with the opposite sign.
                sr_F2D_z_r(z_I, r_I) += norm_F * RPlus_F * rhoKernel_F;
                Jrphi_F2D_z_r(z_I, r_I) -= 0.5 * norm_F * RMinus_F * JrphiKernel_F;
                Jzphi_F2D_z_r(z_I, r_I) -= 0.5 * norm_F * RMinus_F * JzphiKernel_F;
                Jphir_F2D_z_r(z_I, r_I) += 0.5 * norm_F * RMinus_F * (LambdaDown_F + LambdaUp_F) * rInv_F * rhoKernel_F;

                // ∇·J: angular-bra and σ_φ kernels, including 2 Re.
                dJ_F2D_z_r(z_I, r_I) += 2.0 * norm_F * RMinus_F * (-LambdaUp_F * rInv_F * phiUp_F * dphiDowndz_F - dphiUpdz_F * dphiDowndr_F);
            }
        }
    };

    for (int r_I = 0; r_I < basis.r_F1D_r.size(); ++r_I) {
        for (int z_I = 0; z_I < basis.z_F1D_z.size(); ++z_I) {
            for (int block_I = 0; block_I < labels_S2D_block_bsp.size(); ++block_I) {
                add_upup_Func(r_I, z_I, block_I);
                add_dndn_Func(r_I, z_I, block_I);
                add_updn_Func(r_I, z_I, block_I);
                add_dnup_Func(r_I, z_I, block_I);
            }
        }
    }
}
