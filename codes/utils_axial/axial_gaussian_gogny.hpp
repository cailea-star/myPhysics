/**
 * @file    axial_gaussian_gogny.hpp
 * @author  cailea
 * @date    2026-05-23
 * @brief   Gaussian Gogny interaction in an axial HO basis.
 */

#pragma once

#include <array>
#include <cassert>
#include <exception>
#include <filesystem>
#include <iostream>
#include <string>
#include <tuple>
#include <utility>

#include "axial_gaussian_kernel.hpp"

/**
 * @brief Represent Gaussian Gogny matrix elements.
 */
class AxialGaussianGogny {
public:
    using GaussianValues = std::array<double, 2>;

    struct GognyChannels {
        double W_F = 0.0;
        double B_F = 0.0;
        double H_F = 0.0;
        double M_F = 0.0;
    };

    struct GognyElements {
        double vSamePosPosPosPos_F = 0.0;
        double vSameNegNegNegNeg_F = 0.0;
        double vSameNegPosNegPos_F = 0.0;
        double vSamePosNegPosNeg_F = 0.0;
        double vSamePosNegNegPos_F = 0.0;
        double vSameNegPosPosNeg_F = 0.0;
        double vCrossPosPosPosPos_F = 0.0;
        double vCrossNegNegNegNeg_F = 0.0;
        double vCrossNegPosNegPos_F = 0.0;
        double vCrossPosNegPosNeg_F = 0.0;
        double vCrossPosNegNegPos_F = 0.0;
        double vCrossNegPosPosNeg_F = 0.0;
    };

private:
    AxialConfig axialconfig;
    std::string forceName_Str;
    GaussianValues mu_F1D_g{};
    GaussianValues W_F1D_g{};
    GaussianValues B_F1D_g{};
    GaussianValues H_F1D_g{};
    GaussianValues M_F1D_g{};
    AxialGaussianKernel<2> kernel;

public:
    /**
     * @brief  Construct the Gaussian Gogny interaction.
     * @math   P_G \oplus C_{axial} \rightarrow V_G
     * @output Initialized Gogny interaction.
     */
    AxialGaussianGogny(const AxialConfig& axialconfig_, const std::string& forceName_Str_, const GaussianValues& mu_F1D_g_, const GaussianValues& W_F1D_g_, const GaussianValues& B_F1D_g_, const GaussianValues& H_F1D_g_, const GaussianValues& M_F1D_g_)
    : axialconfig(axialconfig_), forceName_Str(forceName_Str_), mu_F1D_g(mu_F1D_g_), W_F1D_g(W_F1D_g_), B_F1D_g(B_F1D_g_), H_F1D_g(H_F1D_g_), M_F1D_g(M_F1D_g_), kernel(axialconfig, mu_F1D_g) {}

    /**
     * @brief  Build or load Gaussian kernel tables.
     * @math   \{\mu_g\} \rightarrow \{G_g^z,G_g^r\}
     * @output Built Gogny kernel.
     */
    void build_tables() {
        if (kernel.isBuilt_B) {return;}

        // (P_G, config, cache) → kernel.
        const std::string filepath_Str = kernel_cache_path(forceName_Str, axialconfig);
        try {
            kernel = AxialGaussianKernel<2>::from_cache(filepath_Str, axialconfig, mu_F1D_g);
            std::cout << "[AxialGaussianGogny]: Loaded Gaussian kernel cache: " << filepath_Str << std::endl;
            return;
        } catch (const std::exception& error_) {
            std::cout << "[AxialGaussianGogny]: Writing Gaussian kernel cache: " << filepath_Str << "\n (" << error_.what() << ")" << std::endl;
        }
        kernel.build_tables();
        kernel.to_cache(filepath_Str);
    }

    /**
     * @brief  Read twelve antisymmetrized Gogny elements.
     * @math   \bar v_{12;34}=\langle12|V(1-P_{12})|34\rangle
     * @output Six Kramers sectors per isospin channel.
     */
    GognyElements read_v(int sp1_I, int sp2_I, int sp3_I, int sp4_I) const {
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

        // (V_{12,34}, V_{12,43}) → (v̄_same,v̄_cross).
        const auto calc_vbar_Func = [&](bool isReversed1_B, bool isReversed2_B, bool isReversed3_B, bool isReversed4_B) {
            // (1234,1243) → (V_{12,34},V_{12,43}).
            const GognyChannels V1234_ = calc_spatial_v(sp1_I, isReversed1_B, sp2_I, isReversed2_B, sp3_I, isReversed3_B, sp4_I, isReversed4_B);
            const GognyChannels V1243_ = calc_spatial_v(sp1_I, isReversed1_B, sp2_I, isReversed2_B, sp4_I, isReversed4_B, sp3_I, isReversed3_B);

            const int twoSigma1_I = apply_reversal_sign_Func(label1_.twoSigma_I, isReversed1_B);
            const int twoSigma2_I = apply_reversal_sign_Func(label2_.twoSigma_I, isReversed2_B);
            const int twoSigma3_I = apply_reversal_sign_Func(label3_.twoSigma_I, isReversed3_B);
            const int twoSigma4_I = apply_reversal_sign_Func(label4_.twoSigma_I, isReversed4_B);
            const int deltaS_13_24_I = static_cast<int>(twoSigma1_I == twoSigma3_I) * static_cast<int>(twoSigma2_I == twoSigma4_I);
            const int deltaS_14_23_I = static_cast<int>(twoSigma1_I == twoSigma4_I) * static_cast<int>(twoSigma2_I == twoSigma3_I);

            // δ_T = (1,1) → v̄_same.
            const double vW_F = deltaS_13_24_I * V1234_.W_F - deltaS_14_23_I * V1243_.W_F;
            const double vB_F = deltaS_14_23_I * V1234_.B_F - deltaS_13_24_I * V1243_.B_F;
            const double vH_F = -deltaS_13_24_I * V1234_.H_F + deltaS_14_23_I * V1243_.H_F;
            const double vM_F = -deltaS_14_23_I * V1234_.M_F + deltaS_13_24_I * V1243_.M_F;
            const double vSame_F = vW_F + vB_F + vH_F + vM_F;

            // δ_T = (1,0) → v̄_cross.
            const double vCrossW_F = deltaS_13_24_I * V1234_.W_F;
            const double vCrossB_F = deltaS_14_23_I * V1234_.B_F;
            const double vCrossH_F = deltaS_14_23_I * V1243_.H_F;
            const double vCrossM_F = deltaS_13_24_I * V1243_.M_F;
            const double vCross_F = vCrossW_F + vCrossB_F + vCrossH_F + vCrossM_F;
            return std::pair<double, double>{vSame_F, vCross_F};
        };

        // {++++, +-+-, +--+} → {(same,cross)}.
        GognyElements elements_;
        std::tie(elements_.vSamePosPosPosPos_F, elements_.vCrossPosPosPosPos_F) = calc_vbar_Func(false, false, false, false);
        std::tie(elements_.vSamePosNegPosNeg_F, elements_.vCrossPosNegPosNeg_F) = calc_vbar_Func(false, true, false, true);
        std::tie(elements_.vSamePosNegNegPos_F, elements_.vCrossPosNegNegPos_F) = calc_vbar_Func(false, true, true, false);

        // (----, -+-+, -++-) ← (++++, +-+-, +--+).
        elements_.vSameNegNegNegNeg_F = elements_.vSamePosPosPosPos_F;
        elements_.vSameNegPosNegPos_F = elements_.vSamePosNegPosNeg_F;
        elements_.vSameNegPosPosNeg_F = elements_.vSamePosNegNegPos_F;
        elements_.vCrossNegNegNegNeg_F = elements_.vCrossPosPosPosPos_F;
        elements_.vCrossNegPosNegPos_F = elements_.vCrossPosNegPosNeg_F;
        elements_.vCrossNegPosPosNeg_F = elements_.vCrossPosNegNegPos_F;
        return elements_;
    }

private:
    /**
     * @brief  Calculate one spatial Gogny element.
     * @math   (W,B,H,M)_{12,34}=\sum_g(W_g,B_g,H_g,M_g)G_g^zG_g^r
     * @output Four interaction-channel elements.
     */
    GognyChannels calc_spatial_v(int sp1_I, bool isReversed1_B, int sp2_I, bool isReversed2_B, int sp3_I, bool isReversed3_B, int sp4_I, bool isReversed4_B) const {
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

        GognyChannels channels_;
        if ((label1_.nz_I + label2_.nz_I + label3_.nz_I + label4_.nz_I) % 2 != 0) {return channels_;}
        if (-Lambda1_I - Lambda2_I + Lambda3_I + Lambda4_I != 0) {return channels_;}

        const auto& Gr_F1D_g = kernel.read_Gr(label1_.nr_I, Lambda1_I, label2_.nr_I, Lambda2_I, label3_.nr_I, Lambda3_I, label4_.nr_I, Lambda4_I);
        const auto& Gz_F1D_g = kernel.read_Gz(label1_.nz_I, label2_.nz_I, label3_.nz_I, label4_.nz_I);

        // (W_g,B_g,H_g,M_g)G_g^zG_g^r → (W,B,H,M).
        for (int g_I = 0; g_I < static_cast<int>(mu_F1D_g.size()); ++g_I) {
            const double G_F = Gz_F1D_g[g_I] * Gr_F1D_g[g_I];
            channels_.W_F += W_F1D_g[g_I] * G_F;
            channels_.B_F += B_F1D_g[g_I] * G_F;
            channels_.H_F += H_F1D_g[g_I] * G_F;
            channels_.M_F += M_F1D_g[g_I] * G_F;
        }
        return channels_;
    }

    /**
     * @brief  Build the Gogny-kernel cache path.
     * @math   (P_G,N_{shell}) \rightarrow path
     * @output Cache-file path.
     */
    static std::string kernel_cache_path(const std::string& forceName_Str_, const AxialConfig& axialconfig_) {
        // (P_G, N_shell, cwd) → cache path.
        const std::string cacheName_Str = "Kernel-" + forceName_Str_ + "-Nshell" + std::to_string(axialconfig_.Nshell_I) + ".cache";
        const std::filesystem::path cwd_ = std::filesystem::current_path();
        const std::filesystem::path cacheDir_ = (cwd_.filename() == "build") ? cwd_ : cwd_ / "build";
        std::filesystem::create_directories(cacheDir_);
        return (cacheDir_ / cacheName_Str).string();
    }
};
