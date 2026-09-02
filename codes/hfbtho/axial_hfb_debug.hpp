/**
 * @file    axial_hfb_debug.hpp
 * @author  cailea
 * @date    2026-05-05
 * @brief   Print axial HFB diagnostic values.
 */

#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Eigenvalues>

#include "axial_basis.hpp"
#include "axial_hfb.hpp"

/**
 * @brief Print axial quadrature nodes.
 * @math (ζ_z,η_r) → stdout
 * @output Selected mesh endpoints.
 */
inline void debug_meshes(const AxialBasis& basis_) {
    const auto& zeta_F1D_z = basis_.zeta_F1D_z;
    const auto& eta_F1D_r = basis_.eta_F1D_r;

    // (ζ_z,η_r) → stdout.
    std::cout << std::scientific << std::setprecision(5) << std::right;
    std::cout << "\n【debug begin】gauss_nodes::C++  (nGH=" << zeta_F1D_z.size() << ", nGL=" << eta_F1D_r.size() << ")" << std::endl;
    std::cout << std::setw(10) << "GH_zeta: " << zeta_F1D_z.head(5).transpose() << " ... " << zeta_F1D_z.tail(5).transpose() << std::endl;
    std::cout << std::setw(10) << "GH_eta : " << eta_F1D_r.head(5).transpose() << " ... " << eta_F1D_r.tail(5).transpose() << std::endl;
    std::cout << "【debug end】gauss_nodes::C++" << std::endl;
}

/**
 * @brief Print axial single-particle labels.
 * @math {α_{sp}} → stdout
 * @output Selected basis labels.
 */
inline void debug_labels(const AxialBasis& basis_) {
    const int Nsp_I = static_cast<int>(basis_.labels_S1D_sp.size());
    const int Nshow_I = std::min(5, Nsp_I);

    // {α_{sp}}_{first} → stdout.
    std::cout << "\n【debug begin】spLabels::C++  (Omega2=1, sp_num =" << Nsp_I << ")" << std::endl;
    for (int sp_I = 0; sp_I < Nshow_I; ++sp_I) {
        const AxialSPLabel& label_ = basis_.labels_S1D_sp[sp_I];
        std::cout << "  [" << sp_I << "] (nz=" << label_.nz_I << ",nr=" << label_.nr_I << ",L=" << label_.Lambda_I << ",N=" << label_.N_I << "," << (label_.twoSigma_I > 0 ? "up" : "dn") << ")" << std::endl;
    }

    // {α_{sp}}_{last} → stdout.
    if (Nsp_I > 5) {std::cout << "  ..." << std::endl;}
    for (int sp_I = std::max(0, Nsp_I - 5); sp_I < Nsp_I; ++sp_I) {
        const AxialSPLabel& label_ = basis_.labels_S1D_sp[sp_I];
        std::cout << "  [" << sp_I << "] (nz=" << label_.nz_I << ",nr=" << label_.nr_I << ",L=" << label_.Lambda_I << ",N=" << label_.N_I << "," << (label_.twoSigma_I > 0 ? "up" : "dn") << ")" << std::endl;
    }
    std::cout << "【debug end】spLabels::C++" << std::endl;
}

/**
 * @brief Print coordinate-space HFB fields.
 * @math F_q(z,r) → stdout
 * @output Selected field submatrices.
 */
inline void debug_field(const AxialHFBField& field_) {
    // F_q(z,r) → stdout.
    std::cout << "\n【debug begin】field::C++  (nz=" << field_.vcent_F2D_z_r.rows() << ", nr=" << field_.vcent_F2D_z_r.cols() << ")" << std::endl;
    std::cout << std::scientific << std::setprecision(5) << std::right;
    std::cout << std::setw(20) << "vcent [top, left]:\n" << field_.vcent_F2D_z_r.topLeftCorner(5, 5) << "\n";
    std::cout << std::setw(20) << "vmass [top, left]:\n" << field_.vmass_F2D_z_r.topLeftCorner(5, 5) << "\n";
    std::cout << std::setw(20) << "vpair [top, left]:\n" << field_.vpair_F2D_z_r.topLeftCorner(5, 5) << "\n";
    std::cout << std::setw(20) << "vD2   [top, left]:\n" << field_.vD2_F2D_z_r.topLeftCorner(5, 5) << "\n";
    std::cout << std::setw(20) << "vDr   [top, left]:\n" << field_.vDr_F2D_z_r.topLeftCorner(5, 5) << "\n";
    std::cout << std::setw(20) << "vDz   [top, left]:\n" << field_.vDz_F2D_z_r.topLeftCorner(5, 5) << "\n";
    std::cout << std::setw(20) << "vJzphi[top, left]:\n" << field_.vJzphi_F2D_z_r.topLeftCorner(5, 5) << "\n";
    std::cout << std::setw(20) << "vJphiz[top, left]:\n" << field_.vJphiz_F2D_z_r.topLeftCorner(5, 5) << "\n";
    std::cout << std::setw(20) << "vJphir[top, left]:\n" << field_.vJphir_F2D_z_r.topLeftCorner(5, 5) << "\n";
    std::cout << std::setw(20) << "vJrphi[top, left]:\n" << field_.vJrphi_F2D_z_r.topLeftCorner(5, 5) << "\n";
    std::cout << std::setw(20) << "vdJ   [top, left]:\n" << field_.vdJ_F2D_z_r.topLeftCorner(5, 5) << "\n";
    std::cout << "【debug end】field::C++" << std::endl;
}

/**
 * @brief Print block particle-hole and pairing fields.
 * @math (Γ_b,Δ_b) → stdout
 * @output Selected block-field submatrices.
 */
inline void debug_Gamma_Delta(const AxialHFBBlock& block_) {
    // (Γ_b,Δ_b) → stdout.
    std::cout << "\n【debug begin】gamdel::C++  (Omega2=" << block_.twoOmega_I << ", sp_num =" << block_.labels_S1D_bsp.size() << ")" << std::endl;
    std::cout << std::scientific << std::setprecision(5) << std::right;
    std::cout << std::setw(20) << "Gamma[top, left]:\n" << block_.Gamma_F2D_bsp_bsp.topLeftCorner(5, 5) << "\n";
    std::cout << std::setw(20) << "Gamma[bot, left]:\n" << block_.Gamma_F2D_bsp_bsp.bottomLeftCorner(5, 5) << "\n\n";
    std::cout << std::setw(20) << "Delta[top, left]:\n" << block_.Delta_F2D_bsp_bsp.topLeftCorner(5, 5) << "\n";
    std::cout << std::setw(20) << "Delta[bot, left]:\n" << block_.Delta_F2D_bsp_bsp.bottomLeftCorner(5, 5) << "\n";
    std::cout << "【debug end】gamdel::C++" << std::endl;
}

/**
 * @brief Print block quasiparticle solutions.
 * @math (Γ_b,Δ_b,λ) → (E,U,V)
 * @output Selected quasiparticle values.
 */
inline void debug_UV(const AxialHFBBlock& blockInput_, double lambda_F) {
    AxialHFBBlock block_ = blockInput_;
    const int Nsp_I = static_cast<int>(block_.labels_S1D_bsp.size());

    // λ ∈ ℝ → diagonalize ℋ.
    std::cout << "\n【debug begin】lambda::C++  (Omega2=" << block_.twoOmega_I << ", sp_num =" << Nsp_I << ", lambda0=" << lambda_F << ")" << std::endl;
    if (std::isfinite(lambda_F)) {
        Eigen::MatrixXd H_F2D_nambu_nambu(2 * Nsp_I, 2 * Nsp_I);
        const Eigen::MatrixXd I_F2D_bsp_bsp = Eigen::MatrixXd::Identity(Nsp_I, Nsp_I);
        H_F2D_nambu_nambu.topLeftCorner(Nsp_I, Nsp_I) = block_.Gamma_F2D_bsp_bsp - lambda_F * I_F2D_bsp_bsp;
        H_F2D_nambu_nambu.topRightCorner(Nsp_I, Nsp_I) = block_.Delta_F2D_bsp_bsp;
        H_F2D_nambu_nambu.bottomLeftCorner(Nsp_I, Nsp_I) = block_.Delta_F2D_bsp_bsp;
        H_F2D_nambu_nambu.bottomRightCorner(Nsp_I, Nsp_I) = -block_.Gamma_F2D_bsp_bsp + lambda_F * I_F2D_bsp_bsp;
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver_(H_F2D_nambu_nambu);
        assert(eigensolver_.info() == Eigen::Success);
        const Eigen::VectorXd& Eqp_F1D_nambu = eigensolver_.eigenvalues();
        const Eigen::MatrixXd& UV_F2D_nambu_nambu = eigensolver_.eigenvectors();
        int bqpCount_I = 0;

        // E_μ>0 → (E_μ,U_μ,V_μ).
        for (int nambu_I = 0; nambu_I < 2 * Nsp_I; ++nambu_I) {
            if (Eqp_F1D_nambu(nambu_I) > 0.0) {
                block_.Eqp_F1D_bqp(bqpCount_I) = Eqp_F1D_nambu(nambu_I);
                block_.U_F2D_bsp_bqp.col(bqpCount_I) = UV_F2D_nambu_nambu.col(nambu_I).head(Nsp_I);
                block_.V_F2D_bsp_bqp.col(bqpCount_I) = UV_F2D_nambu_nambu.col(nambu_I).tail(Nsp_I);
                ++bqpCount_I;
            }
        }
    }

    // {(E_μ,U_μ,V_μ)}_{μ=0}^4 → stdout.
    std::cout << std::scientific << std::setprecision(5) << std::right;
    std::cout << std::setw(10) << "Eqp[0] =" << std::setw(10) << block_.Eqp_F1D_bqp(0) << "\n";
    std::cout << std::setw(10) << "U[0] =" << block_.U_F2D_bsp_bqp.col(0).head(5).transpose() << " ... " << block_.U_F2D_bsp_bqp.col(0).tail(5).transpose() << "\n";
    std::cout << std::setw(10) << "V[0] =" << block_.V_F2D_bsp_bqp.col(0).head(5).transpose() << " ... " << block_.V_F2D_bsp_bqp.col(0).tail(5).transpose() << "\n";
    std::cout << std::setw(10) << "Eqp[1] =" << std::setw(10) << block_.Eqp_F1D_bqp(1) << "\n";
    std::cout << std::setw(10) << "U[1] =" << block_.U_F2D_bsp_bqp.col(1).head(5).transpose() << " ... " << block_.U_F2D_bsp_bqp.col(1).tail(5).transpose() << "\n";
    std::cout << std::setw(10) << "V[1] =" << block_.V_F2D_bsp_bqp.col(1).head(5).transpose() << " ... " << block_.V_F2D_bsp_bqp.col(1).tail(5).transpose() << "\n";
    std::cout << std::setw(10) << "Eqp[2] =" << std::setw(10) << block_.Eqp_F1D_bqp(2) << "\n";
    std::cout << std::setw(10) << "U[2] =" << block_.U_F2D_bsp_bqp.col(2).head(5).transpose() << " ... " << block_.U_F2D_bsp_bqp.col(2).tail(5).transpose() << "\n";
    std::cout << std::setw(10) << "V[2] =" << block_.V_F2D_bsp_bqp.col(2).head(5).transpose() << " ... " << block_.V_F2D_bsp_bqp.col(2).tail(5).transpose() << "\n";
    std::cout << std::setw(10) << "Eqp[3] =" << std::setw(10) << block_.Eqp_F1D_bqp(3) << "\n";
    std::cout << std::setw(10) << "U[3] =" << block_.U_F2D_bsp_bqp.col(3).head(5).transpose() << " ... " << block_.U_F2D_bsp_bqp.col(3).tail(5).transpose() << "\n";
    std::cout << std::setw(10) << "V[3] =" << block_.V_F2D_bsp_bqp.col(3).head(5).transpose() << " ... " << block_.V_F2D_bsp_bqp.col(3).tail(5).transpose() << "\n";
    std::cout << std::setw(10) << "Eqp[4] =" << std::setw(10) << block_.Eqp_F1D_bqp(4) << "\n";
    std::cout << std::setw(10) << "U[4] =" << block_.U_F2D_bsp_bqp.col(4).head(5).transpose() << " ... " << block_.U_F2D_bsp_bqp.col(4).tail(5).transpose() << "\n";
    std::cout << std::setw(10) << "V[4] =" << block_.V_F2D_bsp_bqp.col(4).head(5).transpose() << " ... " << block_.V_F2D_bsp_bqp.col(4).tail(5).transpose() << "\n";
    std::cout << "【debug end】lambda::C++" << std::endl;
}

/**
 * @brief Override block quasiparticle solutions.
 * @math ({Γ_b,Δ_b},λ) → {E_b,U_b,V_b}
 * @output Updated block solutions.
 */
inline void override_UV(std::vector<AxialHFBBlock>& blocks_S1D_block_, double lambda_F) {
    const auto override_block_Func = [&](AxialHFBBlock& block_) {
        const int Nsp_I = static_cast<int>(block_.labels_S1D_bsp.size());
        Eigen::MatrixXd H_F2D_nambu_nambu(2 * Nsp_I, 2 * Nsp_I);
        const Eigen::MatrixXd I_F2D_bsp_bsp = Eigen::MatrixXd::Identity(Nsp_I, Nsp_I);
        H_F2D_nambu_nambu.topLeftCorner(Nsp_I, Nsp_I) = block_.Gamma_F2D_bsp_bsp - lambda_F * I_F2D_bsp_bsp;
        H_F2D_nambu_nambu.topRightCorner(Nsp_I, Nsp_I) = block_.Delta_F2D_bsp_bsp;
        H_F2D_nambu_nambu.bottomLeftCorner(Nsp_I, Nsp_I) = block_.Delta_F2D_bsp_bsp;
        H_F2D_nambu_nambu.bottomRightCorner(Nsp_I, Nsp_I) = -block_.Gamma_F2D_bsp_bsp + lambda_F * I_F2D_bsp_bsp;
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver_(H_F2D_nambu_nambu);
        assert(eigensolver_.info() == Eigen::Success);
        const Eigen::VectorXd& Eqp_F1D_nambu = eigensolver_.eigenvalues();
        const Eigen::MatrixXd& UV_F2D_nambu_nambu = eigensolver_.eigenvectors();
        int bqpCount_I = 0;

        // E_μ>0 → (E_μ,U_μ,V_μ).
        for (int nambu_I = 0; nambu_I < 2 * Nsp_I; ++nambu_I) {
            if (Eqp_F1D_nambu(nambu_I) > 0.0) {
                block_.Eqp_F1D_bqp(bqpCount_I) = Eqp_F1D_nambu(nambu_I);
                block_.U_F2D_bsp_bqp.col(bqpCount_I) = UV_F2D_nambu_nambu.col(nambu_I).head(Nsp_I);
                block_.V_F2D_bsp_bqp.col(bqpCount_I) = UV_F2D_nambu_nambu.col(nambu_I).tail(Nsp_I);
                ++bqpCount_I;
            }
        }
    };

    // {B_b} → {(E_b,U_b,V_b)}.
    std::cout << "\n【debug begin】override lambda::C++  (lambda0=" << lambda_F << ")" << std::endl;
    for (AxialHFBBlock& block_ : blocks_S1D_block_) {
        override_block_Func(block_);
    }
    std::cout << "【debug end】override lambda::C++" << std::endl;
}

/**
 * @brief Print coordinate-space HFB densities.
 * @math D_q(z,r) → stdout
 * @output Density normalization and submatrices.
 */
inline void debug_density(const AxialBasis& basis_, const AxialHFBDensity& density_) {
    const auto weight_Func = [&](int r_I, int z_I) {
        return basis_.w_F2D_z_r(z_I, r_I);
    };
    const auto calc_particle_number_Func = [&](const AxialHFBDensity& densityInput_) {
        double Nparticle_F = 0.0;
        const int Nz_I = static_cast<int>(densityInput_.rho_F2D_z_r.rows());
        const int Nr_I = static_cast<int>(densityInput_.rho_F2D_z_r.cols());

        // N_q = Σ_{z,r}w_{zr}ρ_q(z,r).
        for (int z_I = 0; z_I < Nz_I; ++z_I) {
            for (int r_I = 0; r_I < Nr_I; ++r_I) {
                Nparticle_F += densityInput_.rho_F2D_z_r(z_I, r_I) * weight_Func(r_I, z_I);
            }
        }
        return Nparticle_F;
    };

    // (N_q,D_q) → stdout.
    std::cout << std::scientific << std::setprecision(10);
    std::cout << "\n【debug begin】density::C++  Σρ·whl=" << calc_particle_number_Func(density_) << std::endl;
    std::cout << std::scientific << std::setprecision(5) << std::right;
    std::cout << std::setw(20) << "rho   [top, left]:\n" << density_.rho_F2D_z_r.topLeftCorner(5, 5) << "\n";
    std::cout << std::setw(20) << "tau   [top, left]:\n" << density_.tau_F2D_z_r.topLeftCorner(5, 5) << "\n";
    std::cout << std::setw(20) << "rhoDr [top, left]:\n" << density_.rhoDr_F2D_z_r.topLeftCorner(5, 5) << "\n";
    std::cout << std::setw(20) << "rhoDz [top, left]:\n" << density_.rhoDz_F2D_z_r.topLeftCorner(5, 5) << "\n";
    std::cout << std::setw(20) << "rhoD2 [top, left]:\n" << density_.rhoD2_F2D_z_r.topLeftCorner(5, 5) << "\n";
    std::cout << std::setw(20) << "kappa [top, left]:\n" << density_.kappa_F2D_z_r.topLeftCorner(5, 5) << "\n";
    std::cout << std::setw(20) << "Jzphi [top, left]:\n" << density_.Jzphi_F2D_z_r.topLeftCorner(5, 5) << "\n";
    std::cout << std::setw(20) << "Jphiz [top, left]:\n" << density_.Jphiz_F2D_z_r.topLeftCorner(5, 5) << "\n";
    std::cout << std::setw(20) << "Jrphi [top, left]:\n" << density_.Jrphi_F2D_z_r.topLeftCorner(5, 5) << "\n";
    std::cout << std::setw(20) << "Jphir [top, left]:\n" << density_.Jphir_F2D_z_r.topLeftCorner(5, 5) << "\n";
    std::cout << std::setw(20) << "dJ    [top, left]:\n" << density_.dJ_F2D_z_r.topLeftCorner(5, 5) << "\n";
    std::cout << "【debug end】density::C++" << std::endl;
}
