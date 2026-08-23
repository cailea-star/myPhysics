/**
 * @file    test_axial_branchless.cpp
 * @author  cailea
 * @date    2026-08-23
 * @brief   Validate axial configuration, basis, and Gaussian kernels.
 */

#include <Eigen/Core>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>
#include "axial_basis.hpp"
#include "axial_gaussian_kernel.hpp"

/**
 * @brief  Test branchless axial numerical modules.
 * @math   config → basis → {G^z,G^r}
 * @output Labeled references, computed values, and errors.
 */
int main() {
    // (b_z,b_r,N_shell,reflection,μ) → encoded configuration.
    const double bz_F = 1.0;
    const double br_F = 1.0;
    const int Nshell_I = 1;
    const bool useReflection_B = false;
    const std::array<double, 1> mu_F1D_g{1.0};
    const double tolBasis_F = 1.0e-11;
    const double tolKernel_F = 1.0e-12;
    AxialConfig config_(bz_F, br_F, Nshell_I, useReflection_B);
    Eigen::Vector4i labelCodeRef_I1D_sp;
    Eigen::Vector4i labelCode_I1D_sp = Eigen::Vector4i::Zero();
    labelCodeRef_I1D_sp << 121000, 1100, 21001, 23100;
    for (int sp_I = 0; sp_I < std::min<int>(4, config_.labels_S1D_sp.size()); ++sp_I) {
        const AxialSPLabel& label_ = config_.labels_S1D_sp[sp_I];
        labelCode_I1D_sp(sp_I) = label_.nz_I + 10 * label_.nr_I + 100 * label_.Lambda_I + 1000 * label_.twoOmega_I + 10000 * (label_.twoSigma_I + 1) + 100000 * static_cast<int>(label_.isParityPositive_B);
    }
    // symmetry blocks → size/index codes.
    Eigen::Vector2i blockSizeRef_I1D_block;
    Eigen::Vector2i blockSize_I1D_block = Eigen::Vector2i::Zero();
    Eigen::Vector2i indexCodeRef_I1D_block;
    Eigen::Vector2i indexCode_I1D_block = Eigen::Vector2i::Zero();
    blockSizeRef_I1D_block << 3, 1;
    indexCodeRef_I1D_block << 12, 3;
    for (int block_I = 0; block_I < std::min<int>(2, config_.labels_S2D_block_spb.size()); ++block_I) {
        blockSize_I1D_block(block_I) = config_.labels_S2D_block_spb[block_I].size();
        for (int spb_I = 0; spb_I < config_.indices_I2D_block_spb[block_I].size(); ++spb_I) {indexCode_I1D_block(block_I) = 10 * indexCode_I1D_block(block_I) + config_.indices_I2D_block_spb[block_I][spb_I];}
    }
    int configError_I = std::max({(labelCode_I1D_sp - labelCodeRef_I1D_sp).cwiseAbs().maxCoeff(), (blockSize_I1D_block - blockSizeRef_I1D_block).cwiseAbs().maxCoeff(), (indexCode_I1D_block - indexCodeRef_I1D_block).cwiseAbs().maxCoeff(), std::abs(static_cast<int>(config_.labels_S1D_sp.size()) - 4), std::abs(static_cast<int>(config_.labels_S2D_block_spb.size()) - 2), std::abs(config_.Nz_I - 10), std::abs(config_.Nr_I - 9)});

    // φ_{n_z}(z) → axial Hermite references.
    std::vector<AxialSPLabel> labelsBasis_S1D_sp{AxialSPLabel(0, 0, 0, 1), AxialSPLabel(1, 1, 0, 1)};
    AxialHermiteBasis zBasis_(bz_F, config_.Nz_I, labelsBasis_S1D_sp, false);
    AxialLaguerreBasis rBasis_(br_F, config_.Nr_I, labelsBasis_S1D_sp);
    const double pi_F = std::acos(-1.0);
    const double sqrt2_F = std::sqrt(2.0);
    Eigen::MatrixXd phiZRef_F2D_sp_z(2, zBasis_.z_F1D_z.size());
    Eigen::MatrixXd dphiZRef_F2D_sp_z(2, zBasis_.z_F1D_z.size());
    Eigen::MatrixXd ddphiZRef_F2D_sp_z(2, zBasis_.z_F1D_z.size());
    for (int z_I = 0; z_I < zBasis_.z_F1D_z.size(); ++z_I) {
        const double z_F = zBasis_.z_F1D_z(z_I);
        const double expHalfZ_F = std::exp(-0.5 * z_F * z_F);
        const double phi0_F = std::pow(pi_F, -0.25) * expHalfZ_F;
        phiZRef_F2D_sp_z(0, z_I) = phi0_F;
        dphiZRef_F2D_sp_z(0, z_I) = -z_F * phi0_F;
        ddphiZRef_F2D_sp_z(0, z_I) = (z_F * z_F - 1.0) * phi0_F;
        phiZRef_F2D_sp_z(1, z_I) = sqrt2_F * z_F * phi0_F;
        dphiZRef_F2D_sp_z(1, z_I) = sqrt2_F * (1.0 - z_F * z_F) * phi0_F;
        ddphiZRef_F2D_sp_z(1, z_I) = sqrt2_F * (z_F * z_F * z_F - 3.0 * z_F) * phi0_F;
    }
    double zBasisError_F = std::max({(zBasis_.phi_F2D_sp_z - phiZRef_F2D_sp_z).cwiseAbs().maxCoeff(), (zBasis_.dphi_F2D_sp_z - dphiZRef_F2D_sp_z).cwiseAbs().maxCoeff(), (zBasis_.ddphi_F2D_sp_z - ddphiZRef_F2D_sp_z).cwiseAbs().maxCoeff()});
    // φ_{n_r}^{Λ}(r) → radial Laguerre references.
    Eigen::MatrixXd phiRRef_F2D_sp_r(2, rBasis_.r_F1D_r.size());
    Eigen::MatrixXd dphiRRef_F2D_sp_r(2, rBasis_.r_F1D_r.size());
    Eigen::MatrixXd ddphiRRef_F2D_sp_r(2, rBasis_.r_F1D_r.size());
    for (int r_I = 0; r_I < rBasis_.r_F1D_r.size(); ++r_I) {
        const double r_F = rBasis_.r_F1D_r(r_I);
        const double expHalfR_F = std::exp(-0.5 * r_F * r_F);
        phiRRef_F2D_sp_r(0, r_I) = sqrt2_F * expHalfR_F;
        dphiRRef_F2D_sp_r(0, r_I) = -sqrt2_F * r_F * expHalfR_F;
        ddphiRRef_F2D_sp_r(0, r_I) = sqrt2_F * (r_F * r_F - 1.0) * expHalfR_F;
        phiRRef_F2D_sp_r(1, r_I) = sqrt2_F * (1.0 - r_F * r_F) * expHalfR_F;
        dphiRRef_F2D_sp_r(1, r_I) = sqrt2_F * (r_F * r_F * r_F - 3.0 * r_F) * expHalfR_F;
        ddphiRRef_F2D_sp_r(1, r_I) = sqrt2_F * (-std::pow(r_F, 4) + 6.0 * r_F * r_F - 3.0) * expHalfR_F;
    }
    double rBasisError_F = std::max({(rBasis_.phi_F2D_sp_r - phiRRef_F2D_sp_r).cwiseAbs().maxCoeff(), (rBasis_.dphi_F2D_sp_r - dphiRRef_F2D_sp_r).cwiseAbs().maxCoeff(), (rBasis_.ddphi_F2D_sp_r - ddphiRRef_F2D_sp_r).cwiseAbs().maxCoeff()});

    // (config,μ) → (G^z,G^r) references.
    AxialGaussianKernel<1> kernel_(config_, mu_F1D_g);
    kernel_.build_tables();
    Eigen::Vector4d kernelRef_F1D_value;
    Eigen::Vector4d kernel_F1D_value;
    kernelRef_F1D_value << 1.0 / std::sqrt(3.0), 2.0 / (3.0 * std::sqrt(3.0)), 1.0 / 3.0, 5.0 / 27.0;
    kernel_F1D_value << kernel_.read_Gz(0, 0, 0, 0)[0], kernel_.read_Gz(1, 0, 1, 0)[0], kernel_.read_Gr(0, 0, 0, 0, 0, 0, 0, 0)[0], kernel_.read_Gr(0, 1, 0, -1, 0, 1, 0, -1)[0];
    double kernelError_F = (kernel_F1D_value - kernelRef_F1D_value).cwiseAbs().maxCoeff();

    // (reference,computed,error) → stdout.
    std::cout << std::scientific << std::setprecision(12);
    std::cout << "[Input] b_z, b_r, N_shell, reflection, μ = " << bz_F << ", " << br_F << ", " << Nshell_I << ", " << useReflection_B << ", " << mu_F1D_g[0] << "\n";
    std::cout << "[Reference] label codes = " << labelCodeRef_I1D_sp.transpose() << "\n[Reference] block sizes = " << blockSizeRef_I1D_block.transpose() << "\n[Reference] block indices = " << indexCodeRef_I1D_block.transpose() << "\n";
    std::cout << "[Reference] phi_z:\n" << phiZRef_F2D_sp_z << "\n[Reference] dphi_z:\n" << dphiZRef_F2D_sp_z << "\n[Reference] ddphi_z:\n" << ddphiZRef_F2D_sp_z << "\n";
    std::cout << "[Reference] phi_r:\n" << phiRRef_F2D_sp_r << "\n[Reference] dphi_r:\n" << dphiRRef_F2D_sp_r << "\n[Reference] ddphi_r:\n" << ddphiRRef_F2D_sp_r << "\n[Reference] Gaussian = " << kernelRef_F1D_value.transpose() << "\n";
    std::cout << "[Computed] label codes = " << labelCode_I1D_sp.transpose() << "\n[Computed] block sizes = " << blockSize_I1D_block.transpose() << "\n[Computed] block indices = " << indexCode_I1D_block.transpose() << "\n";
    std::cout << "[Computed] phi_z:\n" << zBasis_.phi_F2D_sp_z << "\n[Computed] dphi_z:\n" << zBasis_.dphi_F2D_sp_z << "\n[Computed] ddphi_z:\n" << zBasis_.ddphi_F2D_sp_z << "\n";
    std::cout << "[Computed] phi_r:\n" << rBasis_.phi_F2D_sp_r << "\n[Computed] dphi_r:\n" << rBasis_.dphi_F2D_sp_r << "\n[Computed] ddphi_r:\n" << rBasis_.ddphi_F2D_sp_r << "\n[Computed] Gaussian = " << kernel_F1D_value.transpose() << "\n";
    std::cout << "[Error] config = " << configError_I << "\n[Error] z basis = " << zBasisError_F << "\n[Error] r basis = " << rBasisError_F << "\n[Error] Gaussian = " << kernelError_F << "\n";
    // (error_config,error_z,error_r,error_G) < tolerances.
    assert(configError_I == 0);
    assert(zBasisError_F <= tolBasis_F);
    assert(rBasisError_F <= tolBasis_F);
    assert(kernelError_F <= tolKernel_F);
    return 0;
}
