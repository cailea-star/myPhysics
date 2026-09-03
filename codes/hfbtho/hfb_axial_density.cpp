/**
 * @file    hfb_axial_density.cpp
 * @author  cailea
 * @date    2026-05-05
 * @brief   Calculate axial HFB coordinate-space densities.
 */

#include "hfb_axial.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

#include <omp.h>

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

void AxialHFBDensity::add_density(const AxialBasis& global_basis_, const AxialHFBBlock& block_) {
    const int Nz_I = static_cast<int>(rho_F2D_z_r.rows());
    const int Nr_I = static_cast<int>(rho_F2D_z_r.cols());
    const int Nbsp_I = static_cast<int>(block_.labels_S1D_bsp.size());
    const int twoOmega_I = block_.twoOmega_I;
    const double Omega_F = 0.5 * twoOmega_I;
    const double LambdaUp_F = 0.5 * twoOmega_I - 0.5;
    const double LambdaDown_F = 0.5 * twoOmega_I + 0.5;
    const double kramersDegeneracy_F = 2.0; // {+Ω,-Ω} → 2.

    const int spStart_I = block_.indices_I1D_bsp.front();

    // (ρ_{bsp,bsp'},κ_{bsp,bsp'},φ_bsp) → local densities.
    for (int r_I = 0; r_I < Nr_I; ++r_I) {
        for (int z_I = 0; z_I < Nz_I; ++z_I) {
            const double rInv_F = 1.0 / global_basis_.r_F1D_r(r_I);
            const double rInv2_F = rInv_F * rInv_F;

            // Φ(sp,z,r) → local Eigen views.
            Eigen::Map<const Eigen::VectorXd> phi_F1D_bsp(&global_basis_.phi_F3D_sp_z_r(spStart_I, z_I, r_I), Nbsp_I);
            Eigen::Map<const Eigen::VectorXd> dphidr_F1D_bsp(&global_basis_.dphidr_F3D_sp_z_r(spStart_I, z_I, r_I), Nbsp_I);
            Eigen::Map<const Eigen::VectorXd> dphidz_F1D_bsp(&global_basis_.dphidz_F3D_sp_z_r(spStart_I, z_I, r_I), Nbsp_I);
            Eigen::Map<const Eigen::VectorXd> ddphidr_F1D_bsp(&global_basis_.ddphidr_F3D_sp_z_r(spStart_I, z_I, r_I), Nbsp_I);
            Eigen::Map<const Eigen::VectorXd> ddphidz_F1D_bsp(&global_basis_.ddphidz_F3D_sp_z_r(spStart_I, z_I, r_I), Nbsp_I);

            double rho_F = 0.0;
            double tau_F = 0.0;
            double rhoD2_F = 0.0;
            double rhoDr_F = 0.0;
            double rhoDz_F = 0.0;
            double kappa_F = 0.0;
            double dJ_F = 0.0;
            double Jzphi_F = 0.0;
            double Jphiz_F = 0.0;
            double Jphir_F = 0.0;
            double Jrphi_F = 0.0;

            // Σ_1=Σ_2 → {ρ,τ,κ,Δρ,∇ρ,J_{φz},∇·J}.
            const auto add_same_spin_Func = [&](const std::vector<int>& indices_I1D_bsp_, double Lambda_F, double spinSign_F) {
                const double Lambda2_F = Lambda_F * Lambda_F;
                const int Nspin_I = static_cast<int>(indices_I1D_bsp_.size());
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
                        rho_F += rhoWeight_F * rhoKernel_F;
                        tau_F += rhoWeight_F * tauKernel_F;
                        rhoD2_F += rhoWeight_F * rhoD2Kernel_F;
                        rhoDr_F += rhoWeight_F * rhoDrKernel_F;
                        rhoDz_F += rhoWeight_F * rhoDzKernel_F;
                        Jphiz_F += rhoWeight_F * JphizKernel_F;
                        dJ_F += rhoWeight_F * dJKernel_F;

                        kappa_F += 0.5 * symmetryFactor_F * kappa12_F * rhoKernel_F;
                    }
                }
            };

            add_same_spin_Func(block_.indices_I1D_bup, LambdaUp_F, 1.0);
            add_same_spin_Func(block_.indices_I1D_bdn, LambdaDown_F, -1.0);

            // Σ_1≠Σ_2 → {J_{rφ},J_{zφ},J_{φr},∇·J}.
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

                    Jrphi_F += rhoUpDown_F * (phiDown_F * dphiUpdr_F - phiUp_F * dphiDowndr_F);
                    Jzphi_F += rhoUpDown_F * (phiDown_F * dphiUpdz_F - phiUp_F * dphiDowndz_F);
                    Jphir_F += rhoUpDown_F * 2.0 * Omega_F * rInv_F * phiUp_F * phiDown_F;
                    dJ_F += rhoUpDown_F * (2.0 * (-LambdaUp_F * rInv_F * phiUp_F * dphiDowndz_F - LambdaDown_F * rInv_F * phiDown_F * dphiUpdz_F) + 2.0 * (dphiDowndz_F * dphiUpdr_F - dphiUpdz_F * dphiDowndr_F));
                }
            }

            // D_{Ω}(z,r) → D(z,r).
            rho_F2D_z_r(z_I, r_I) += kramersDegeneracy_F * rho_F;
            tau_F2D_z_r(z_I, r_I) += kramersDegeneracy_F * tau_F;
            kappa_F2D_z_r(z_I, r_I) += kramersDegeneracy_F * kappa_F;
            rhoD2_F2D_z_r(z_I, r_I) += kramersDegeneracy_F * rhoD2_F;
            rhoDr_F2D_z_r(z_I, r_I) += kramersDegeneracy_F * rhoDr_F;
            rhoDz_F2D_z_r(z_I, r_I) += kramersDegeneracy_F * rhoDz_F;
            dJ_F2D_z_r(z_I, r_I) += kramersDegeneracy_F * dJ_F;
            Jzphi_F2D_z_r(z_I, r_I) += kramersDegeneracy_F * Jzphi_F;
            Jphiz_F2D_z_r(z_I, r_I) += kramersDegeneracy_F * Jphiz_F;
            Jphir_F2D_z_r(z_I, r_I) += kramersDegeneracy_F * Jphir_F;
            Jrphi_F2D_z_r(z_I, r_I) += kramersDegeneracy_F * Jrphi_F;
        }
    }
}

void AxialHFBDensity::update_density(const AxialBasis& global_basis_, const AxialHFBBlockList& blocklist_) {
    set_zero();
    const int Nblock_I = static_cast<int>(blocklist_.blocks_S1D_block.size());
    const int Nthread_I = std::max(1, omp_get_max_threads());
    std::vector<AxialHFBDensity> partials_S1D_thread;
    partials_S1D_thread.reserve(Nthread_I);
    for (int thread_I = 0; thread_I < Nthread_I; ++thread_I) {
        partials_S1D_thread.push_back(*this);
    }

    #pragma omp parallel num_threads(Nthread_I)
    {
        const int thread_I = omp_get_thread_num();
        #pragma omp for schedule(dynamic)
        for (int block_I = 0; block_I < Nblock_I; ++block_I) {
            const AxialHFBBlock& block_ = blocklist_.blocks_S1D_block[block_I];
            partials_S1D_thread[thread_I].add_density(global_basis_, block_);
        }
    }

    // Σ_thread D_thread → D.
    for (const AxialHFBDensity& partial_ : partials_S1D_thread) {
        rho_F2D_z_r += partial_.rho_F2D_z_r;
        tau_F2D_z_r += partial_.tau_F2D_z_r;
        kappa_F2D_z_r += partial_.kappa_F2D_z_r;
        rhoD2_F2D_z_r += partial_.rhoD2_F2D_z_r;
        rhoDr_F2D_z_r += partial_.rhoDr_F2D_z_r;
        rhoDz_F2D_z_r += partial_.rhoDz_F2D_z_r;
        dJ_F2D_z_r += partial_.dJ_F2D_z_r;
        Jzphi_F2D_z_r += partial_.Jzphi_F2D_z_r;
        Jphiz_F2D_z_r += partial_.Jphiz_F2D_z_r;
        Jphir_F2D_z_r += partial_.Jphir_F2D_z_r;
        Jrphi_F2D_z_r += partial_.Jrphi_F2D_z_r;
    }
}
