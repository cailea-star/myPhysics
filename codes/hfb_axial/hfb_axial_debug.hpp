/**
 * @file    hfb_axial_debug.hpp
 * @author  cailea
 * @date    2026-05-05
 * @brief   Print axial HFB diagnostic values.
 */

#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <numbers>
#include <iomanip>
#include <iostream>
#include <vector>

#include <Eigen/Core>

#include "cylindrical_basis.hpp"
#include "hfb_axial.hpp"

/**
 * @brief Print a vector or its endpoint samples.
 * @math x → stdout
 * @output Full vector or first and last five entries.
 */
inline void debug_print_vector(const Eigen::VectorXd& vector_F1D_index_) {
    if (vector_F1D_index_.size() <= 10) {
        std::cout << vector_F1D_index_.transpose() << '\n';
    } else {
        std::cout << vector_F1D_index_.head(5).transpose() << " ... " << vector_F1D_index_.tail(5).transpose() << '\n';
    }
}

/**
 * @brief Print a matrix or its four corners.
 * @math A → stdout
 * @output Full matrix or bounded corner samples.
 */
inline void debug_print_matrix(const char* name_Str_, const Eigen::MatrixXd& matrix_F2D_row_column_) {
    std::cout << name_Str_ << " (" << matrix_F2D_row_column_.rows() << ',' << matrix_F2D_row_column_.cols() << ")\n";
    if (matrix_F2D_row_column_.rows() <= 10 && matrix_F2D_row_column_.cols() <= 10) {
        std::cout << matrix_F2D_row_column_ << '\n';
        return;
    }
    const int Nrow_I = std::min(5, static_cast<int>(matrix_F2D_row_column_.rows()));
    const int Ncolumn_I = std::min(5, static_cast<int>(matrix_F2D_row_column_.cols()));
    // A → {top-left,top-right,bottom-left,bottom-right}.
    std::cout << "[top, left]\n" << matrix_F2D_row_column_.topLeftCorner(Nrow_I, Ncolumn_I) << '\n';
    std::cout << "[top, right]\n" << matrix_F2D_row_column_.topRightCorner(Nrow_I, Ncolumn_I) << '\n';
    std::cout << "[bottom, left]\n" << matrix_F2D_row_column_.bottomLeftCorner(Nrow_I, Ncolumn_I) << '\n';
    std::cout << "[bottom, right]\n" << matrix_F2D_row_column_.bottomRightCorner(Nrow_I, Ncolumn_I) << '\n';
}
/**
 * @brief Print axial quadrature nodes.
 * @math (ζ_z,η_r) → stdout
 * @output Selected mesh endpoints.
 */
inline void debug_meshes(const CylindricalBasis2D& basis_) {
    const auto& zeta_F1D_z = basis_.zeta_F1D_z;
    const auto& eta_F1D_r = basis_.eta_F1D_r;

    // (ζ_z,η_r) → stdout.
    std::cout << std::scientific << std::setprecision(5) << std::right;
    std::cout << "\n【debug begin】gauss_nodes::C++  (nGH=" << zeta_F1D_z.size() << ", nGL=" << eta_F1D_r.size() << ")" << std::endl;
    std::cout << "GH_zeta: ";
    debug_print_vector(zeta_F1D_z);
    std::cout << "GH_eta: ";
    debug_print_vector(eta_F1D_r);
    std::cout << "【debug end】gauss_nodes::C++" << std::endl;
}

/**
 * @brief Print axial single-particle labels.
 * @math {α_{sp}} → stdout
 * @output Selected basis labels.
 */
inline void debug_labels(const CylindricalBasis2D& basis_) {
    const int Nsp_I = static_cast<int>(basis_.labels_S1D_sp.size());
    const int Nshow_I = Nsp_I <= 10 ? Nsp_I : 5;

    // {α_{sp}}_{first} → stdout.
    std::cout << "\n【debug begin】spLabels::C++  (sp_num =" << Nsp_I << ")" << std::endl;
    for (int sp_I = 0; sp_I < Nshow_I; ++sp_I) {
        const CylindricalSPLabel& label_ = basis_.labels_S1D_sp[sp_I];
        std::cout << "  [" << sp_I << "] (nz=" << label_.nz_I << ",nr=" << label_.nr_I << ",L=" << label_.Lambda_I << ",N=" << label_.N_I << "," << (label_.twoSigma_I > 0 ? "up" : "dn") << ")" << std::endl;
    }

    // {α_{sp}}_{last} → stdout.
    if (Nsp_I > 10) {std::cout << "  ..." << std::endl;}
    for (int sp_I = std::max(Nshow_I, Nsp_I - 5); sp_I < Nsp_I; ++sp_I) {
        const CylindricalSPLabel& label_ = basis_.labels_S1D_sp[sp_I];
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
    debug_print_matrix("vcent", field_.vcent_F2D_z_r);
    debug_print_matrix("vmass", field_.vmass_F2D_z_r);
    debug_print_matrix("vpair", field_.vpair_F2D_z_r);
    debug_print_matrix("vD2", field_.vD2_F2D_z_r);
    debug_print_matrix("vDr", field_.vDr_F2D_z_r);
    debug_print_matrix("vDz", field_.vDz_F2D_z_r);
    debug_print_matrix("vJzphi", field_.vJzphi_F2D_z_r);
    debug_print_matrix("vJphiz", field_.vJphiz_F2D_z_r);
    debug_print_matrix("vJphir", field_.vJphir_F2D_z_r);
    debug_print_matrix("vJrphi", field_.vJrphi_F2D_z_r);
    debug_print_matrix("vdJ", field_.vdJ_F2D_z_r);
    std::cout << "【debug end】field::C++" << std::endl;
}

/**
 * @brief Print block particle-hole and pairing fields.
 * @math (Γ_b⁺⁺,Γ_b⁻⁻,Δ_b⁺⁻,Δ_b⁻⁺) → stdout
 * @output Selected block-field submatrices.
 */
inline void debug_Gamma_Delta(const AxialHFBBlock& block_) {
    // (Γ_b⁺⁺,Γ_b⁻⁻,Δ_b⁺⁻,Δ_b⁻⁺) → stdout.
    std::cout << "\n【debug begin】gamdel::C++  (Omega2=" << block_.twoOmega_I << ", sp_num =" << block_.labels_S1D_bsp.size() << ")" << std::endl;
    std::cout << std::scientific << std::setprecision(5) << std::right;
    debug_print_matrix("GammaPosPos", block_.GammaPosPos_F2D_bsp_bsp);
    debug_print_matrix("GammaNegNeg", block_.GammaNegNeg_F2D_bsp_bsp);
    debug_print_matrix("DeltaPosNeg", block_.DeltaPosNeg_F2D_bsp_bsp);
    debug_print_matrix("DeltaNegPos", block_.DeltaNegPos_F2D_bsp_bsp);
    std::cout << "【debug end】gamdel::C++" << std::endl;
}

/**
 * @brief Print representative quasiparticle solutions.
 * @math (Γ⁺⁺,Γ⁻⁻,Δ⁺⁻,Δ⁻⁺,λ) → ℋ⁺ → (E,U⁺,V⁻)
 * @output Selected quasiparticle values; input block remains unchanged.
 */
inline void debug_UV(const AxialHFBBlock& blockInput_, double lambda_F, const HFBSettings& hfbsettings_) {
    AxialHFBBlock block_ = blockInput_;
    if (std::isfinite(lambda_F)) {block_.update_UV_E_rho_kappa(hfbsettings_, lambda_F);}

    // (E,U,V) → bounded quasiparticle samples.
    const auto print_sector_Func = [&](const char* sector_Str_, const Eigen::VectorXd& Eqp_F1D_bqp_, const Eigen::MatrixXd& U_F2D_bsp_bqp_, const Eigen::MatrixXd& V_F2D_bsp_bqp_) {
        std::cout << sector_Str_ << " (Nqp=" << Eqp_F1D_bqp_.size() << ")\nEqp: ";
        debug_print_vector(Eqp_F1D_bqp_);
        for (int bqp_I = 0; bqp_I < std::min(5, static_cast<int>(Eqp_F1D_bqp_.size())); ++bqp_I) {
            std::cout << "U[" << bqp_I << "]: ";
            debug_print_vector(U_F2D_bsp_bqp_.col(bqp_I));
            std::cout << "V[" << bqp_I << "]: ";
            debug_print_vector(V_F2D_bsp_bqp_.col(bqp_I));
        }
    };

    // {Eqp,UPos,VNeg} → stdout.
    std::cout << "\n【debug begin】lambda::C++ (Omega2=" << block_.twoOmega_I << ", lambda=" << lambda_F << ")\n";
    std::cout << std::scientific << std::setprecision(5) << std::right;
    print_sector_Func("Pos: Eqp, UPos, VNeg", block_.Eqp_F1D_bqp, block_.UPos_F2D_bsp_bqp, block_.VNeg_F2D_bsp_bqp);
    std::cout << "【debug end】lambda::C++" << std::endl;
}

/**
 * @brief Recompute block solutions at a specified chemical potential.
 * @math ({Γ⁺⁺,Γ⁻⁻,Δ⁺⁻,Δ⁻⁺},λ,T) → ℋ⁺ → {E,U⁺,V⁻,ρ⁺⁺,ρ⁻⁻,κ⁺⁻,κ⁻⁺}
 * @output Updated block amplitudes and densities.
 * @note Blocking corrections must be reapplied after recomputation.
 */
inline void override_UV(std::vector<AxialHFBBlock>& blocks_X1D_block_, double lambda_F, const HFBSettings& hfbsettings_) {
    assert(std::isfinite(lambda_F));
    std::cout << "\n【debug begin】override lambda::C++ (lambda=" << lambda_F << ")\n";
    for (AxialHFBBlock& block_ : blocks_X1D_block_) {
        block_.update_UV_E_rho_kappa(hfbsettings_, lambda_F);
    }
    std::cout << "【debug end】override lambda::C++" << std::endl;
}
/**
 * @brief Print coordinate-space HFB densities.
 * @math D_q(z,r) → stdout
 * @output Density normalization and submatrices.
 */
inline void debug_density(const CylindricalBasis2D& basis_, const AxialHFBDensity& density_) {
    const auto weight_Func = [&](int r_I, int z_I) {
        return 2.0 * std::numbers::pi * basis_.w_F2D_z_r(z_I, r_I);
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
    debug_print_matrix("rho", density_.rho_F2D_z_r);
    debug_print_matrix("tau", density_.tau_F2D_z_r);
    debug_print_matrix("rhoDr", density_.rhoDr_F2D_z_r);
    debug_print_matrix("rhoDz", density_.rhoDz_F2D_z_r);
    debug_print_matrix("rhoD2", density_.rhoD2_F2D_z_r);
    debug_print_matrix("kappa", density_.kappa_F2D_z_r);
    debug_print_matrix("Jzphi", density_.Jzphi_F2D_z_r);
    debug_print_matrix("Jphiz", density_.Jphiz_F2D_z_r);
    debug_print_matrix("Jrphi", density_.Jrphi_F2D_z_r);
    debug_print_matrix("Jphir", density_.Jphir_F2D_z_r);
    debug_print_matrix("dJ", density_.dJ_F2D_z_r);
    std::cout << "【debug end】density::C++" << std::endl;
}
