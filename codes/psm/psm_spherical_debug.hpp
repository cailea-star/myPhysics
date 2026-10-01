/**
 * @file    psm_spherical_debug.hpp
 * @author  cailea
 * @date    2026-09-28
 * @brief   PSM diagnostic printing.
 */

#pragma once

#include <cassert>
#include <iomanip>
#include <iostream>
#include <initializer_list>
#include <string_view>
#include "psm_spherical.hpp"

namespace {

/**
 * @brief Print vector endpoints horizontally.
 * @math x = (x₀,…,x₄,…,xₙ₋₅,…,xₙ₋₁).
 * @output First/last five entries; n ≤ 10: all.
 */
template<typename T>
inline void print_vector(std::string_view name, const Eigen::MatrixBase<T>& data_T1D_i) {
    assert(data_T1D_i.rows() == 1 || data_T1D_i.cols() == 1 || data_T1D_i.size() == 0);
    std::cout << name << " = ";
    for (Eigen::Index i_I = 0; i_I < data_T1D_i.size(); ++i_I) {
        if (data_T1D_i.size() > 10 && i_I == 5) {
            std::cout << " ...";
            i_I = data_T1D_i.size() - 5;
        }
        if (i_I > 0) {std::cout << ' ';}
        std::cout << data_T1D_i.derived().coeff(i_I);
    }
    std::cout << '\n';
}

/**
 * @brief Print matrix corner blocks.
 * @math M = [Mᵢⱼ].
 * @output First/last five rows and columns; ≤ 10: all.
 */
template<typename T>
inline void print_matrix(std::string_view name, const Eigen::MatrixBase<T>& data_T2D_i_j) {
    constexpr int width_I = Eigen::NumTraits<typename T::Scalar>::IsComplex ? 28 : 16;
    std::cout << name << " =\n";
    for (Eigen::Index row_I = 0; row_I < data_T2D_i_j.rows(); ++row_I) {
        if (data_T2D_i_j.rows() > 10 && row_I == 5) {
            std::cout << "...\n";
            row_I = data_T2D_i_j.rows() - 5;
        }
        for (Eigen::Index column_I = 0; column_I < data_T2D_i_j.cols(); ++column_I) {
            if (data_T2D_i_j.cols() > 10 && column_I == 5) {
                std::cout << " ... ";
                column_I = data_T2D_i_j.cols() - 5;
            }
            std::cout << std::setw(width_I) << data_T2D_i_j(row_I, column_I) << ' ';
        }
        std::cout << '\n';
    }
}

/**
 * @brief Print labeled table endpoints.
 * @math T = [Tᵢⱼ]; headers ↔ columns.
 * @output First/last five rows, all columns; ≤ 10: all.
 */
template<typename T>
inline void print_table(std::initializer_list<std::string_view> headers, const Eigen::MatrixBase<T>& data_T2D_i_j) {
    assert(static_cast<Eigen::Index>(headers.size()) == data_T2D_i_j.cols());
    constexpr int width_I = Eigen::NumTraits<typename T::Scalar>::IsComplex ? 28 : 16;
    for (std::string_view header : headers) {std::cout << std::setw(width_I) << header << ' ';}
    std::cout << '\n';
    for (Eigen::Index row_I = 0; row_I < data_T2D_i_j.rows(); ++row_I) {
        if (data_T2D_i_j.rows() > 10 && row_I == 5) {
            std::cout << "...\n";
            row_I = data_T2D_i_j.rows() - 5;
        }
        for (Eigen::Index column_I = 0; column_I < data_T2D_i_j.cols(); ++column_I) {
            std::cout << std::setw(width_I) << data_T2D_i_j(row_I, column_I) << ' ';
        }
        std::cout << '\n';
    }
}

} // namespace

/**
 * @brief Print neutron and proton Nilsson-BCS solutions.
 * @math f → u → v → Ehf → Eqp → Δ.
 * @output Block labels, eigenvectors, amplitudes, energies, and gaps.
 */
inline void print_hfbcs(const PSMSpherical& psm) {
    std::cout << std::scientific << std::setprecision(8);

    // neutron: block = (N,2K).
    std::cout << "\nHFBCS-Neutron\n";
    for (int block_I = 0; block_I < psm.hfb_neutron.Nblock_I; ++block_I) {
        const auto& label = psm.psm_setting.sphericalsetting_neutron.labels_S2D_block_bsp[block_I].front();
        const auto& solution = psm.hfb_neutron.solutions[block_I];
        std::cout << "\nN = " << label.N_I << ", 2K = " << label.twom_I << '\n';

        // f(bsp,bhf); u,v,Ehf,Eqp(bhf).
        print_matrix("f(hf)", solution.fhf_F2D_sp_bhf);
        print_vector("u(hf)", solution.uPos_F1D_bhf);
        print_vector("v(hf)", solution.vNeg_F1D_bhf);
        print_vector("E(hf) [MeV]", solution.Ehf_F1D_bhf);
        print_vector("E(qp) [MeV]", solution.Eqp_F1D_bhf);
    }
    std::cout << "\nDelta [MeV] = " << psm.hfb_neutron.Delta_F << '\n';

    // proton: block = (N,2K).
    std::cout << "\nHFBCS-Proton\n";
    for (int block_I = 0; block_I < psm.hfb_proton.Nblock_I; ++block_I) {
        const auto& label = psm.psm_setting.sphericalsetting_proton.labels_S2D_block_bsp[block_I].front();
        const auto& solution = psm.hfb_proton.solutions[block_I];
        std::cout << "\nN = " << label.N_I << ", 2K = " << label.twom_I << '\n';

        // f(bsp,bhf); u,v,Ehf,Eqp(bhf).
        print_matrix("f(hf)", solution.fhf_F2D_sp_bhf);
        print_vector("u(hf)", solution.uPos_F1D_bhf);
        print_vector("v(hf)", solution.vNeg_F1D_bhf);
        print_vector("E(hf) [MeV]", solution.Ehf_F1D_bhf);
        print_vector("E(qp) [MeV]", solution.Eqp_F1D_bhf);
    }
    std::cout << "\nDelta [MeV] = " << psm.hfb_proton.Delta_F << '\n';
}

/**
 * @brief Print CI energies and natural-basis coefficients.
 * @math Xci(r,ν); r = eigenN, ν = eigenH; Xci†Xci = 1.
 * @output E: first/last five; Xci: four corner blocks.
 */
inline void print_ci(const PSMSpherical& psm, int TargetTwoI_I) {
    assert(psm.Eci_F1D_eigenH.size() > 0);

    std::cout << std::scientific << std::setprecision(8);
    std::cout << "CI-(I=" << TargetTwoI_I * 0.5 << ")\n";
    print_vector("Eci [MeV]", psm.Eci_F1D_eigenH);
    print_matrix("Xci [eigenN, eigenH]", psm.Xci_C2D_eigenN_eigenH);
}
