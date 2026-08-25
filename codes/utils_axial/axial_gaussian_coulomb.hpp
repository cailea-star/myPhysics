/**
 * @file    axial_gaussian_coulomb.hpp
 * @author  cailea
 * @date    2026-05-27
 * @brief   Gaussian Coulomb interaction in an axial HO basis.
 */

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <filesystem>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

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

    struct VElement {
        double vcffff_F = 0.0;
        double vcftft_F = 0.0;

        /**
         * @brief  Sum Coulomb channels.
         * @math   v_C = v_C^{ffff} + v_C^{f\tilde f f\tilde f}
         * @output Coulomb element sum.
         */
        double sum() const {
            return vcffff_F + vcftft_F;
        }
    };

private:
    KernelBase kernel;
    AxialConfig axialconfig;
    GaussianValues mu_F1D_g{};
    GaussianValues W_F1D_g{};
    std::vector<AxialSPLabel> labels_S1D_sp;

public:
    /**
     * @brief  Construct a Gaussian Coulomb interaction.
     * @math   r^{-1} \approx \sum_g W_g\exp(-r^2/\mu_g^2)
     * @output Initialized Coulomb interaction.
     */
    AxialGaussianCoulomb(const AxialConfig& axialconfig_, double e2_F_ = 1.439978408596513)
    : AxialGaussianCoulomb(axialconfig_, calc_coulomb_expansion(axialconfig_, e2_F_)) {}

    /**
     * @brief  Construct from a supplied Gaussian expansion.
     * @math   \{\mu_g,W_g\} \rightarrow V_C
     * @output Initialized Coulomb interaction.
     */
    AxialGaussianCoulomb(const AxialConfig& axialconfig_, const std::pair<GaussianValues, GaussianValues>& expansion_F2D_component_g_)
    : kernel(KernelBase::from_config(axialconfig_, expansion_F2D_component_g_.first)), axialconfig(axialconfig_), mu_F1D_g(expansion_F2D_component_g_.first), W_F1D_g(expansion_F2D_component_g_.second), labels_S1D_sp(axialconfig_.globalSpLabels) {}

    /**
     * @brief  Build or load Gaussian kernel tables.
     * @math   \{\mu_g\} \rightarrow \{G_g^z,G_g^r\}
     * @output Built Coulomb kernel.
     * @note   Reuses a compatible cache when available.
     */
    void build_tables() {
        // cache → validated kernel; fallback → rebuilt kernel.
        if (kernel.built_B) {return;}
        const std::string filepath_Str = kernel_cache_path(axialconfig);
        try {
            KernelBase cachedKernel_ = KernelBase::from_cache(filepath_Str);
            cachedKernel_.load_check(axialconfig, mu_F1D_g, filepath_Str);
            kernel = std::move(cachedKernel_);
            std::cout << "[AxialGaussianCoulomb]: Loaded Gaussian kernel cache: " << filepath_Str << std::endl;
        } catch (const std::exception& error_) {
            std::cout << "[AxialGaussianCoulomb]: Writing Gaussian kernel cache: " << filepath_Str << "\n (" << error_.what() << ")" << std::endl;
            kernel.build_Gaussian_tables();
            kernel.write_to(filepath_Str);
        }
    }

    /**
     * @brief  Read antisymmetrized Coulomb matrix elements.
     * @math   v_{12,34} = \delta_{13}\delta_{24}V_{12,34} - \delta_{14}\delta_{23}V_{12,43}
     * @output Direct spin-channel elements.
     */
    VElement read_v(int sp1_I, int sp2_I, int sp3_I, int sp4_I) const {
        // {sp_a} → {α_a}.
        const AxialSPLabel& label1_ = labels_S1D_sp[sp1_I];
        const AxialSPLabel& label2_ = labels_S1D_sp[sp2_I];
        const AxialSPLabel& label3_ = labels_S1D_sp[sp3_I];
        const AxialSPLabel& label4_ = labels_S1D_sp[sp4_I];

        // (Σ_a,reversal_a) → spin signs.
        const auto calcSpinSign_Func = [](bool isSpinUp_B, bool isReversed_B) {
            return 1 - 2 * static_cast<int>(isSpinUp_B and not isReversed_B);
        };
        const auto calcSignedSpin_Func = [](bool isSpinUp_B, bool isReversed_B) {
            return isSpinUp_B not_eq isReversed_B;
        };
        const auto calcValue_Func = [&](bool isReversed1_B, bool isReversed2_B, bool isReversed3_B, bool isReversed4_B) {
            // (spin deltas,spatial terms) → antisymmetrized value.
            const bool isSpin1SignedUp_B = calcSignedSpin_Func(label1_.isSpinUp_B, isReversed1_B);
            const bool isSpin2SignedUp_B = calcSignedSpin_Func(label2_.isSpinUp_B, isReversed2_B);
            const bool isSpin3SignedUp_B = calcSignedSpin_Func(label3_.isSpinUp_B, isReversed3_B);
            const bool isSpin4SignedUp_B = calcSignedSpin_Func(label4_.isSpinUp_B, isReversed4_B);
            const int deltaS_13_24_I = static_cast<int>(isSpin1SignedUp_B == isSpin3SignedUp_B) * static_cast<int>(isSpin2SignedUp_B == isSpin4SignedUp_B);
            const int deltaS_14_23_I = static_cast<int>(isSpin1SignedUp_B == isSpin4SignedUp_B) * static_cast<int>(isSpin2SignedUp_B == isSpin3SignedUp_B);
            const int spinSign_I = calcSpinSign_Func(label1_.isSpinUp_B, isReversed1_B) * calcSpinSign_Func(label2_.isSpinUp_B, isReversed2_B) * calcSpinSign_Func(label3_.isSpinUp_B, isReversed3_B) * calcSpinSign_Func(label4_.isSpinUp_B, isReversed4_B);
            const double V1234_F = calc_spatial_v(sp1_I, isReversed1_B, sp2_I, isReversed2_B, sp3_I, isReversed3_B, sp4_I, isReversed4_B);
            const double V1243_F = calc_spatial_v(sp1_I, isReversed1_B, sp2_I, isReversed2_B, sp4_I, isReversed4_B, sp3_I, isReversed3_B);
            return spinSign_I * (deltaS_13_24_I * V1234_F - deltaS_14_23_I * V1243_F);
        };

        VElement element_;
        element_.vcffff_F = calcValue_Func(false, false, false, false);
        element_.vcftft_F = calcValue_Func(false, true, false, true);
        return element_;
    }

    /**
     * @brief  Calculate one spatial Coulomb matrix element.
     * @math   V_{12,34} = \sum_g W_gG_g^zG_g^r
     * @output Spatial Coulomb matrix element.
     */
    double calc_spatial_v(int sp1_I, bool isReversed1_B, int sp2_I, bool isReversed2_B, int sp3_I, bool isReversed3_B, int sp4_I, bool isReversed4_B) const {
        // {sp_a} → {α_a}.
        const AxialSPLabel& label1_ = labels_S1D_sp[sp1_I];
        const AxialSPLabel& label2_ = labels_S1D_sp[sp2_I];
        const AxialSPLabel& label3_ = labels_S1D_sp[sp3_I];
        const AxialSPLabel& label4_ = labels_S1D_sp[sp4_I];

        // (Λ_a,reversal_a) → signed Λ_a.
        const auto calcLambda_Func = [](int Lambda_I, bool isReversed_B) {
            return (1 - 2 * static_cast<int>(isReversed_B)) * Lambda_I;
        };
        const int Lambda1_I = calcLambda_Func(label1_.Lambda_I, isReversed1_B);
        const int Lambda2_I = calcLambda_Func(label2_.Lambda_I, isReversed2_B);
        const int Lambda3_I = calcLambda_Func(label3_.Lambda_I, isReversed3_B);
        const int Lambda4_I = calcLambda_Func(label4_.Lambda_I, isReversed4_B);

        // n_z parity and Λ conservation.
        if ((label1_.nz_I + label2_.nz_I + label3_.nz_I + label4_.nz_I) % 2 != 0) {return 0.0;}
        if (-Lambda1_I - Lambda2_I + Lambda3_I + Lambda4_I != 0) {return 0.0;}

        // {α_a} → {G_g^r,G_g^z}.
        const auto& Gr_F1D_g = kernel.read_Gaussian_r_values(label1_.nr_I, Lambda1_I, label2_.nr_I, Lambda2_I, label3_.nr_I, Lambda3_I, label4_.nr_I, Lambda4_I);
        const auto& Gz_F1D_g = kernel.read_Gaussian_z_values(label1_.nz_I, label2_.nz_I, label3_.nz_I, label4_.nz_I);

        // Σ_g W_g G_g^z G_g^r.
        double V_F = 0.0;
        for (int g_I = 0; g_I < Ng_I; ++g_I) {
            V_F += W_F1D_g[g_I] * Gz_F1D_g[g_I] * Gr_F1D_g[g_I];
        }
        return V_F;
    }

private:
    /**
     * @brief  Generate a Gauss-Legendre Coulomb expansion.
     * @math   r^{-1} \approx \sum_g W_g\exp(-r^2/\mu_g^2)
     * @output Gaussian ranges and weights.
     */
    static std::pair<GaussianValues, GaussianValues> calc_coulomb_expansion(const AxialConfig& axialconfig_, double e2_F_) {
        // (x_g,w_g,b_max,e²) → (μ_g,W_g).
        GaussianValues mu_F1D_g{};
        GaussianValues W_F1D_g{};
        const GaussLegendreMeshes legendre_meshes(Ng_I);
        const double bmax_F = std::max(axialconfig_.bperp_F, axialconfig_.bz_F);
        const double pi_F = std::acos(-1.0);
        for (int g_I = 0; g_I < Ng_I; ++g_I) {
            const double t_F = 0.5 * (legendre_meshes.x_F1D_x(g_I) + 1.0);
            const double weight_F = 0.5 * legendre_meshes.w_F1D_x(g_I);
            const double omtx2_F = 1.0 - t_F * t_F;
            mu_F1D_g[g_I] = bmax_F * std::sqrt(omtx2_F) / t_F;
            const double Vcoulg_F = 2.0 / std::sqrt(pi_F) * weight_F / bmax_F / std::pow(omtx2_F, 1.5);
            W_F1D_g[g_I] = e2_F_ * Vcoulg_F;
        }
        return {mu_F1D_g, W_F1D_g};
    }

    /**
     * @brief  Build the Coulomb-kernel cache path.
     * @math   (N_g,N_{shell}) \rightarrow path
     * @output Cache-file path.
     */
    static std::string kernel_cache_path(const AxialConfig& axialconfig_) {
        // (N_g,N_shell,cwd) → cache path.
        const std::string cacheName_Str = "Kernel-Coulomb" + std::to_string(Ng_I) + "-Nshell" + std::to_string(axialconfig_.Nshell_I) + ".cache";
        const std::filesystem::path cwd_ = std::filesystem::current_path();
        const std::filesystem::path cacheDir_ = (cwd_.filename() == "build") ? cwd_ : cwd_ / "build";
        std::filesystem::create_directories(cacheDir_);
        return (cacheDir_ / cacheName_Str).string();
    }
};
