/**
 * @file    psm_spherical_debug.hpp
 * @author  cailea
 * @date    2026-09-28
 * @brief   PSM diagnostic printing.
 */

#pragma once

#include <iomanip>
#include <iostream>
#include "psm_spherical.hpp"

namespace psm_debug {

/**
 * @brief Print neutron and proton Nilsson-BCS solutions.
 * @math f → u → v → Ehf → Eqp → Δ.
 * @output Block labels, eigenvectors, amplitudes, energies, and gaps.
 */
inline void print_hfbcs(const PSMSpherical& psm) {
    std::cout << std::fixed << std::setprecision(8);

    // x: first five, ..., last five; size ≤ 10: all.
    const auto print_vector = [](const char* name, const Eigen::VectorXd& data_F1D_i) {
        std::cout << name << " = ";
        if (data_F1D_i.size() <= 10) {
            std::cout << data_F1D_i.transpose();
        } else {
            std::cout << data_F1D_i.head(5).transpose() << " ... " << data_F1D_i.tail(5).transpose();
        }
        std::cout << '\n';
    };

    // M: first/last five rows and columns; dimensions ≤ 10: all.
    const auto print_matrix = [](const char* name, const Eigen::MatrixXd& data_F2D_i_j) {
        std::cout << name << " =\n";
        for (Eigen::Index row_I = 0; row_I < data_F2D_i_j.rows(); ++row_I) {
            if (data_F2D_i_j.rows() > 10 && row_I == 5) {
                std::cout << "...\n";
                row_I = data_F2D_i_j.rows() - 5;
            }
            for (Eigen::Index column_I = 0; column_I < data_F2D_i_j.cols(); ++column_I) {
                if (data_F2D_i_j.cols() > 10 && column_I == 5) {
                    std::cout << " ... ";
                    column_I = data_F2D_i_j.cols() - 5;
                }
                std::cout << std::setw(16) << data_F2D_i_j(row_I, column_I);
            }
            std::cout << '\n';
        }
    };

    // q ∈ {n,p}; block = (N,2K).
    const auto print_species = [&](const char* species, const SphericalSetting& setting, const HFBCS& hfb) {
        std::cout << '\n' << species << '\n';
        for (int block_I = 0; block_I < hfb.Nblock_I; ++block_I) {
            const auto& label = setting.labels_S2D_block_bsp[block_I].front();
            const auto& solution = hfb.solutions[block_I];
            std::cout << "\nN = " << label.N_I << ", 2K = " << label.twom_I << '\n';

            // f(bsp,bhf); u,v,Ehf,Eqp(bhf).
            print_matrix("f", solution.fhf_F2D_sp_bhf);
            print_vector("u", solution.uPos_F1D_bhf);
            print_vector("v", solution.vNeg_F1D_bhf);
            print_vector("Ehf [MeV]", solution.Ehf_F1D_bhf);
            print_vector("Eqp [MeV]", solution.Eqp_F1D_bhf);
        }
        std::cout << "\nDelta [MeV] = " << hfb.Delta_F << '\n';
    };

    print_species("Neutron", psm.sphericalsetting_neutron, psm.hfb_neutron);
    print_species("Proton", psm.sphericalsetting_proton, psm.hfb_proton);
}

} // namespace psm_debug
