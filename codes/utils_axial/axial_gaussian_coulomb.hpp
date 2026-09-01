/**
 * @file    axial_gaussian_coulomb.hpp
 * @author  cailea
 * @date    2026-05-27
 * @brief   Gaussian Coulomb interaction in an axial HO basis.
 */

#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <exception>
#include <filesystem>
#include <iostream>
#include <string>

#include "axial_gaussian_kernel.hpp"
#include "integration_gauss.hpp"

/**
 * @brief Represent Gaussian Coulomb matrix elements.
 */
class AxialGaussianCoulomb {
public:
    static constexpr int Ng_I = 9;
    using KernelBase = AxialGaussianKernel<Ng_I>;
    using GaussianValues = std::array<double, Ng_I>;

    struct GammaElements {
        double vPosPosPosPos_F = 0.0;
        double vPosNegPosNeg_F = 0.0;
        double vNegPosNegPos_F = 0.0;
        double vNegNegNegNeg_F = 0.0;
    };

private:
    struct CoulombExpansion {
        GaussianValues mu_F1D_g{};
        GaussianValues W_F1D_g{};
    };

    AxialConfig axialconfig;
    CoulombExpansion expansion;
    KernelBase kernel;

public:
    /**
     * @brief  Construct the Gaussian Coulomb interaction.
     * @math   r^{-1} \approx \sum_g W_g\exp(-r^2/\mu_g^2)
     * @output Initialized Coulomb interaction.
     */
    AxialGaussianCoulomb(const AxialConfig& axialconfig_, double e2_F_ = 1.439978408596513)
    : axialconfig(axialconfig_), expansion(calc_coulomb_expansion(axialconfig_, e2_F_)), kernel(axialconfig, expansion.mu_F1D_g) {}

    /**
     * @brief  Build or load Gaussian kernel tables.
     * @math   \{\mu_g\} \rightarrow \{G_g^z,G_g^r\}
     * @output Built Coulomb kernel.
     */
    void build_tables() {
        if (kernel.isBuilt_B) {return;}

        // (config, μ_g, cache) → kernel.
        const std::string filepath_Str = kernel_cache_path(axialconfig);
        try {
            kernel = KernelBase::from_cache(filepath_Str, axialconfig, expansion.mu_F1D_g);
            std::cout << "[AxialGaussianCoulomb]: Loaded Gaussian kernel cache: " << filepath_Str << std::endl;
            return;
        } catch (const std::exception& error_) {
            std::cout << "[AxialGaussianCoulomb]: Writing Gaussian kernel cache: " << filepath_Str << "\n (" << error_.what() << ")" << std::endl;
        }
        kernel.build_tables();
        kernel.to_cache(filepath_Str);
    }

    /**
     * @brief  Read four Kramers-sector Coulomb elements.
     * @math   \bar v_{12;34}=\delta_{13}\delta_{24}V_{12,34}-\delta_{14}\delta_{23}V_{12,43}
     * @output Antisymmetrized Γ-sector elements.
     */
    GammaElements read_v(int sp1_I, int sp2_I, int sp3_I, int sp4_I) const {
        assert(kernel.isBuilt_B);
        assert(sp1_I >= 0 && sp1_I < static_cast<int>(axialconfig.labels_S1D_sp.size()));
        assert(sp2_I >= 0 && sp2_I < static_cast<int>(axialconfig.labels_S1D_sp.size()));
        assert(sp3_I >= 0 && sp3_I < static_cast<int>(axialconfig.labels_S1D_sp.size()));
        assert(sp4_I >= 0 && sp4_I < static_cast<int>(axialconfig.labels_S1D_sp.size()));

        // sp_a → α_a.
        const AxialSPLabel& label1_ = axialconfig.labels_S1D_sp[sp1_I];
        const AxialSPLabel& label2_ = axialconfig.labels_S1D_sp[sp2_I];
        const AxialSPLabel& label3_ = axialconfig.labels_S1D_sp[sp3_I];
        const AxialSPLabel& label4_ = axialconfig.labels_S1D_sp[sp4_I];

        // (Λ, 2Σ) → (-Λ, -2Σ).
        const auto apply_reversal_sign_Func = [](int quantum_I, bool isReversed_B) {
            return (1 - 2 * static_cast<int>(isReversed_B)) * quantum_I;
        };
        const auto calc_vbar_Func = [&](bool isReversed1_B, bool isReversed2_B, bool isReversed3_B, bool isReversed4_B) {
            const int twoSigma1_I = apply_reversal_sign_Func(label1_.twoSigma_I, isReversed1_B);
            const int twoSigma2_I = apply_reversal_sign_Func(label2_.twoSigma_I, isReversed2_B);
            const int twoSigma3_I = apply_reversal_sign_Func(label3_.twoSigma_I, isReversed3_B);
            const int twoSigma4_I = apply_reversal_sign_Func(label4_.twoSigma_I, isReversed4_B);
            const int deltaS_13_24_I = static_cast<int>(twoSigma1_I == twoSigma3_I) * static_cast<int>(twoSigma2_I == twoSigma4_I);
            const int deltaS_14_23_I = static_cast<int>(twoSigma1_I == twoSigma4_I) * static_cast<int>(twoSigma2_I == twoSigma3_I);
            const double V1234_F = calc_spatial_v(sp1_I, isReversed1_B, sp2_I, isReversed2_B, sp3_I, isReversed3_B, sp4_I, isReversed4_B);
            const double V1243_F = calc_spatial_v(sp1_I, isReversed1_B, sp2_I, isReversed2_B, sp4_I, isReversed4_B, sp3_I, isReversed3_B);
            return deltaS_13_24_I * V1234_F - deltaS_14_23_I * V1243_F;
        };

        // (v̄⁺⁺, v̄⁺⁻) → (v̄⁻⁻, v̄⁻⁺).
        GammaElements elements_;
        elements_.vPosPosPosPos_F = calc_vbar_Func(false, false, false, false);
        elements_.vPosNegPosNeg_F = calc_vbar_Func(false, true, false, true);
        elements_.vNegPosNegPos_F = elements_.vPosNegPosNeg_F;
        elements_.vNegNegNegNeg_F = elements_.vPosPosPosPos_F;
        return elements_;
    }

private:
    /**
     * @brief  Calculate one spatial Coulomb element.
     * @math   V_{12,34}=\sum_g W_gG_g^zG_g^r
     * @output Spatial Coulomb matrix element.
     */
    double calc_spatial_v(int sp1_I, bool isReversed1_B, int sp2_I, bool isReversed2_B, int sp3_I, bool isReversed3_B, int sp4_I, bool isReversed4_B) const {
        // sp_a → α_a.
        const AxialSPLabel& label1_ = axialconfig.labels_S1D_sp[sp1_I];
        const AxialSPLabel& label2_ = axialconfig.labels_S1D_sp[sp2_I];
        const AxialSPLabel& label3_ = axialconfig.labels_S1D_sp[sp3_I];
        const AxialSPLabel& label4_ = axialconfig.labels_S1D_sp[sp4_I];

        // Λ_a → (-1)^{bar_a} Λ_a.
        const auto apply_reversal_sign_Func = [](int quantum_I, bool isReversed_B) {
            return (1 - 2 * static_cast<int>(isReversed_B)) * quantum_I;
        };
        const int Lambda1_I = apply_reversal_sign_Func(label1_.Lambda_I, isReversed1_B);
        const int Lambda2_I = apply_reversal_sign_Func(label2_.Lambda_I, isReversed2_B);
        const int Lambda3_I = apply_reversal_sign_Func(label3_.Lambda_I, isReversed3_B);
        const int Lambda4_I = apply_reversal_sign_Func(label4_.Lambda_I, isReversed4_B);

        if ((label1_.nz_I + label2_.nz_I + label3_.nz_I + label4_.nz_I) % 2 != 0) {return 0.0;}
        if (-Lambda1_I - Lambda2_I + Lambda3_I + Lambda4_I != 0) {return 0.0;}

        const auto& Gr_F1D_g = kernel.read_Gr(label1_.nr_I, Lambda1_I, label2_.nr_I, Lambda2_I, label3_.nr_I, Lambda3_I, label4_.nr_I, Lambda4_I);
        const auto& Gz_F1D_g = kernel.read_Gz(label1_.nz_I, label2_.nz_I, label3_.nz_I, label4_.nz_I);

        // V = Σ_g W_g G_g^z G_g^r.
        double V_F = 0.0;
        for (int g_I = 0; g_I < Ng_I; ++g_I) {
            V_F += expansion.W_F1D_g[g_I] * Gz_F1D_g[g_I] * Gr_F1D_g[g_I];
        }
        return V_F;
    }

    /**
     * @brief  Generate the Gauss-Legendre Coulomb expansion.
     * @math   r^{-1}\approx\sum_g W_g\exp(-r^2/\mu_g^2)
     * @output Gaussian ranges and weights.
     */
    static CoulombExpansion calc_coulomb_expansion(const AxialConfig& axialconfig_, double e2_F_) {
        assert(std::isfinite(e2_F_) && e2_F_ >= 0.0);

        // (t_g, w_g, b_max, e²) → (μ_g, W_g).
        GaussianValues mu_F1D_g{};
        GaussianValues W_F1D_g{};
        const GaussLegendreMeshes legendre_meshes(Ng_I);
        const double bmax_F = std::max(axialconfig_.br_F, axialconfig_.bz_F);
        const double pi_F = std::acos(-1.0);
        for (int g_I = 0; g_I < Ng_I; ++g_I) {
            const double t_F = 0.5 * (legendre_meshes.x_F1D_x(g_I) + 1.0);
            const double weight_F = 0.5 * legendre_meshes.w_F1D_x(g_I);
            mu_F1D_g[g_I] = bmax_F * std::sqrt(1.0 - t_F * t_F) / t_F;
            W_F1D_g[g_I] = 2.0 * e2_F_ * weight_F / (std::sqrt(pi_F) * bmax_F * std::pow(1.0 - t_F * t_F, 1.5));
        }
        return {mu_F1D_g, W_F1D_g};
    }

    /**
     * @brief  Build the Coulomb-kernel cache path.
     * @math   (N_g,N_{shell}) \rightarrow path
     * @output Cache-file path.
     */
    static std::string kernel_cache_path(const AxialConfig& axialconfig_) {
        // (N_g, N_shell, cwd) → cache path.
        const std::string cacheName_Str = "Kernel-Coulomb" + std::to_string(Ng_I) + "-Nshell" + std::to_string(axialconfig_.Nshell_I) + ".cache";
        const std::filesystem::path cwd_ = std::filesystem::current_path();
        const std::filesystem::path cacheDir_ = (cwd_.filename() == "build") ? cwd_ : cwd_ / "build";
        std::filesystem::create_directories(cacheDir_);
        return (cacheDir_ / cacheName_Str).string();
    }
};
