/**
 * @file    hfb_axial_density.cpp
 * @author  cailea
 * @date    2026-05-05
 * @brief   Calculate axial HFB coordinate-space densities.
 */

#include "hfb_axial.hpp"

#include <cmath>
#include <vector>

void AxialHFBDensity::set_zero() {
    rho_F2D_z_r.setZero();
    tau_F2D_z_r.setZero();
    kappa_F2D_z_r.setZero();
    rhoD2_F2D_z_r.setZero();
    rhoDr_F2D_z_r.setZero();
    rhoDz_F2D_z_r.setZero();
    dJ_F2D_z_r.setZero();
    Jzphi_F2D_z_r.setZero();
    Jphiz_F2D_z_r.setZero();
    Jphir_F2D_z_r.setZero();
    Jrphi_F2D_z_r.setZero();
}

void AxialHFBDensity::update_density(const AxialBasis& global_basis_, const AxialHFBBlockList& blocklist_) {
    set_zero();
    const int Nz_I = static_cast<int>(rho_F2D_z_r.rows());
    const int Nr_I = static_cast<int>(rho_F2D_z_r.cols());
    const int Nblock_I = static_cast<int>(blocklist_.blocks_X1D_block.size());
    constexpr double kramersDegeneracy_F = 2.0;

    struct DensityValues {
        double rho_F = 0.0;
        double tau_F = 0.0;
        double kappa_F = 0.0;
        double rhoD2_F = 0.0;
        double rhoDr_F = 0.0;
        double rhoDz_F = 0.0;
        double dJ_F = 0.0;
        double Jzphi_F = 0.0;
        double Jphiz_F = 0.0;
        double Jphir_F = 0.0;
        double Jrphi_F = 0.0;
    };

    // Σ_1=Σ_2 → {ρ,τ,κ,Δρ,∇ρ,J_{φz},∇·J}.
    const auto add_same_spin_at_oneBlockPoint_Func = [&](DensityValues& densityValues_, const AxialHFBBlock& block_, int z_I, int r_I, const std::vector<int>& indices_I1D_bsp_, double Lambda_F, double spinSign_F) {
        const int Nbsp_I = static_cast<int>(block_.labels_S1D_bsp.size());
        const int spStart_I = block_.indices_I1D_bsp.front();
        const double rInv_F = 1.0 / global_basis_.r_F1D_r(r_I);
        const double rInv2_F = rInv_F * rInv_F;
        const double Lambda2_F = Lambda_F * Lambda_F;
        const int Nspin_I = static_cast<int>(indices_I1D_bsp_.size());

        // Φ(sp,z,r) → local Eigen views.
        Eigen::Map<const Eigen::VectorXd> phi_F1D_bsp(&global_basis_.phi_F3D_sp_z_r(spStart_I, z_I, r_I), Nbsp_I);
        Eigen::Map<const Eigen::VectorXd> dphidr_F1D_bsp(&global_basis_.dphidr_F3D_sp_z_r(spStart_I, z_I, r_I), Nbsp_I);
        Eigen::Map<const Eigen::VectorXd> dphidz_F1D_bsp(&global_basis_.dphidz_F3D_sp_z_r(spStart_I, z_I, r_I), Nbsp_I);
        Eigen::Map<const Eigen::VectorXd> ddphidr_F1D_bsp(&global_basis_.ddphidr_F3D_sp_z_r(spStart_I, z_I, r_I), Nbsp_I);
        Eigen::Map<const Eigen::VectorXd> ddphidz_F1D_bsp(&global_basis_.ddphidz_F3D_sp_z_r(spStart_I, z_I, r_I), Nbsp_I);

        for (int index2_I = 0; index2_I < Nspin_I; ++index2_I) {
            for (int index1_I = index2_I; index1_I < Nspin_I; ++index1_I) {
                const int bsp1_I = indices_I1D_bsp_[index1_I];
                const int bsp2_I = indices_I1D_bsp_[index2_I];
                const double rho12_F = block_.rho_F2D_bsp_bsp(bsp1_I, bsp2_I);
                const double kappa12_F = block_.kappa_F2D_bsp_bsp(bsp1_I, bsp2_I);
                if (std::abs(rho12_F) < 1.0e-14 && std::abs(kappa12_F) < 1.0e-14) {continue;}

                const double symmetryFactor_F = (bsp1_I == bsp2_I) ? 1.0 : 2.0;
                const double phi1_F = phi_F1D_bsp(bsp1_I);
                const double phi2_F = phi_F1D_bsp(bsp2_I);
                const double dphi1dr_F = dphidr_F1D_bsp(bsp1_I);
                const double dphi2dr_F = dphidr_F1D_bsp(bsp2_I);
                const double dphi1dz_F = dphidz_F1D_bsp(bsp1_I);
                const double dphi2dz_F = dphidz_F1D_bsp(bsp2_I);

                // ∇²φ = ∂²_rφ+r⁻¹∂_rφ+∂²_zφ-Λ²r⁻²φ.
                const double laplacianPhi1_F = ddphidr_F1D_bsp(bsp1_I) + rInv_F * dphi1dr_F + ddphidz_F1D_bsp(bsp1_I) - rInv2_F * Lambda2_F * phi1_F;
                const double laplacianPhi2_F = ddphidr_F1D_bsp(bsp2_I) + rInv_F * dphi2dr_F + ddphidz_F1D_bsp(bsp2_I) - rInv2_F * Lambda2_F * phi2_F;

                // (φ_1,φ_2,∇φ_1,∇φ_2) → density kernels.
                const double rhoKernel_F = phi1_F * phi2_F;
                const double tauKernel_F = dphi1dr_F * dphi2dr_F + dphi1dz_F * dphi2dz_F + rInv2_F * Lambda2_F * rhoKernel_F;
                const double rhoD2Kernel_F = laplacianPhi1_F * phi2_F + phi1_F * laplacianPhi2_F + 2.0 * tauKernel_F;
                const double rhoDrKernel_F = dphi1dr_F * phi2_F + phi1_F * dphi2dr_F;
                const double rhoDzKernel_F = dphi1dz_F * phi2_F + phi1_F * dphi2dz_F;
                const double JphizKernel_F = spinSign_F * Lambda_F * rInv_F * rhoKernel_F;
                const double dJKernel_F = spinSign_F * Lambda_F * rInv_F * rhoDrKernel_F;

                // ρ_12K_12 → local normal densities.
                const double rhoWeight_F = symmetryFactor_F * rho12_F;
                densityValues_.rho_F += rhoWeight_F * rhoKernel_F;
                densityValues_.tau_F += rhoWeight_F * tauKernel_F;
                densityValues_.rhoD2_F += rhoWeight_F * rhoD2Kernel_F;
                densityValues_.rhoDr_F += rhoWeight_F * rhoDrKernel_F;
                densityValues_.rhoDz_F += rhoWeight_F * rhoDzKernel_F;
                densityValues_.Jphiz_F += rhoWeight_F * JphizKernel_F;
                densityValues_.dJ_F += rhoWeight_F * dJKernel_F;
                densityValues_.kappa_F += 0.5 * symmetryFactor_F * kappa12_F * rhoKernel_F;
            }
        }
    };

    // Σ_1≠Σ_2 → {J_{rφ},J_{zφ},J_{φr},∇·J}.
    const auto add_opposite_spin_at_oneBlockPoint_Func = [&](DensityValues& densityValues_, const AxialHFBBlock& block_, int z_I, int r_I) {
        const int Nbsp_I = static_cast<int>(block_.labels_S1D_bsp.size());
        const int twoOmega_I = block_.twoOmega_I;
        const double Omega_F = 0.5 * twoOmega_I;
        const double LambdaUp_F = 0.5 * twoOmega_I - 0.5;
        const double LambdaDown_F = 0.5 * twoOmega_I + 0.5;
        const int spStart_I = block_.indices_I1D_bsp.front();
        const double rInv_F = 1.0 / global_basis_.r_F1D_r(r_I);

        // Φ(sp,z,r) → local Eigen views.
        Eigen::Map<const Eigen::VectorXd> phi_F1D_bsp(&global_basis_.phi_F3D_sp_z_r(spStart_I, z_I, r_I), Nbsp_I);
        Eigen::Map<const Eigen::VectorXd> dphidr_F1D_bsp(&global_basis_.dphidr_F3D_sp_z_r(spStart_I, z_I, r_I), Nbsp_I);
        Eigen::Map<const Eigen::VectorXd> dphidz_F1D_bsp(&global_basis_.dphidz_F3D_sp_z_r(spStart_I, z_I, r_I), Nbsp_I);

        for (int bspUp_I : block_.indices_I1D_bup) {
            for (int bspDown_I : block_.indices_I1D_bdn) {
                const double rhoUpDown_F = block_.rho_F2D_bsp_bsp(bspUp_I, bspDown_I);
                if (std::abs(rhoUpDown_F) < 1.0e-14) {continue;}

                const double phiUp_F = phi_F1D_bsp(bspUp_I);
                const double phiDown_F = phi_F1D_bsp(bspDown_I);
                const double dphiUpdr_F = dphidr_F1D_bsp(bspUp_I);
                const double dphiDowndr_F = dphidr_F1D_bsp(bspDown_I);
                const double dphiUpdz_F = dphidz_F1D_bsp(bspUp_I);
                const double dphiDowndz_F = dphidz_F1D_bsp(bspDown_I);

                densityValues_.Jrphi_F += rhoUpDown_F * (phiDown_F * dphiUpdr_F - phiUp_F * dphiDowndr_F);
                densityValues_.Jzphi_F += rhoUpDown_F * (phiDown_F * dphiUpdz_F - phiUp_F * dphiDowndz_F);
                densityValues_.Jphir_F += rhoUpDown_F * 2.0 * Omega_F * rInv_F * phiUp_F * phiDown_F;
                densityValues_.dJ_F += rhoUpDown_F * (2.0 * (-LambdaUp_F * rInv_F * phiUp_F * dphiDowndz_F - LambdaDown_F * rInv_F * phiDown_F * dphiUpdz_F) + 2.0 * (dphiDowndz_F * dphiUpdr_F - dphiUpdz_F * dphiDowndr_F));
            }
        }
    };

    // {B_block} → D(z,r).
    #pragma omp parallel for collapse(2) schedule(static)
    for (int r_I = 0; r_I < Nr_I; ++r_I) {
        for (int z_I = 0; z_I < Nz_I; ++z_I) {
            DensityValues densityValues_{};
            for (int block_I = 0; block_I < Nblock_I; ++block_I) {
                const AxialHFBBlock& block_ = blocklist_.blocks_X1D_block[block_I];
                const int twoOmega_I = block_.twoOmega_I;
                const double LambdaUp_F = 0.5 * twoOmega_I - 0.5;
                const double LambdaDown_F = 0.5 * twoOmega_I + 0.5;
                add_same_spin_at_oneBlockPoint_Func(densityValues_, block_, z_I, r_I, block_.indices_I1D_bup, LambdaUp_F, 1.0);
                add_same_spin_at_oneBlockPoint_Func(densityValues_, block_, z_I, r_I, block_.indices_I1D_bdn, LambdaDown_F, -1.0);
                add_opposite_spin_at_oneBlockPoint_Func(densityValues_, block_, z_I, r_I);
            }

            rho_F2D_z_r(z_I, r_I) = kramersDegeneracy_F * densityValues_.rho_F;
            tau_F2D_z_r(z_I, r_I) = kramersDegeneracy_F * densityValues_.tau_F;
            kappa_F2D_z_r(z_I, r_I) = kramersDegeneracy_F * densityValues_.kappa_F;
            rhoD2_F2D_z_r(z_I, r_I) = kramersDegeneracy_F * densityValues_.rhoD2_F;
            rhoDr_F2D_z_r(z_I, r_I) = kramersDegeneracy_F * densityValues_.rhoDr_F;
            rhoDz_F2D_z_r(z_I, r_I) = kramersDegeneracy_F * densityValues_.rhoDz_F;
            dJ_F2D_z_r(z_I, r_I) = kramersDegeneracy_F * densityValues_.dJ_F;
            Jzphi_F2D_z_r(z_I, r_I) = kramersDegeneracy_F * densityValues_.Jzphi_F;
            Jphiz_F2D_z_r(z_I, r_I) = kramersDegeneracy_F * densityValues_.Jphiz_F;
            Jphir_F2D_z_r(z_I, r_I) = kramersDegeneracy_F * densityValues_.Jphir_F;
            Jrphi_F2D_z_r(z_I, r_I) = kramersDegeneracy_F * densityValues_.Jrphi_F;
        }
    }
}
