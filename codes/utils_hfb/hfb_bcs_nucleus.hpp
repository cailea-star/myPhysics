/**
 * @file    hfb_bcs_nucleus.hpp
 * @author  cailea
 * @date    2026-09-17
 * @brief   Joint neutron-proton HF-BCS iteration.
 */

#pragma once

#include <functional>
#include <limits>

#include "hfb_bcs.hpp"
#include "hfb_setting.hpp"
#include "root_broyden.hpp"

class HFBCSNucleus {
public:
    struct Element {
        double vSamePosPosPosPos_F = 0.0;
        double vSamePosNegPosNeg_F = 0.0;
        double vCrossPosPosPosPos_F = 0.0;
        double vCrossPosNegPosNeg_F = 0.0;
    };

    using GammaElementFunc = std::function<Element(int block13_I, int block24_I, int bsp1_I, int bsp2_I, int bsp3_I, int bsp4_I)>;

    HFBCS hfb_neutron;
    HFBCS hfb_proton;
    HFBSetting hfbsetting;

public:
    /**
     * @brief Initialize species dimensions and workspaces.
     * @math {η_b} → {Nbsp_b}; q ∈ {n,p}.
     * @output Allocated solutions, fields, and workspaces.
     */
    HFBCSNucleus(const std::vector<Eigen::VectorXd>& eta_F2D_block_bsp_, const HFBSetting& hfbsetting_)
    : hfb_neutron(eta_F2D_block_bsp_), hfb_proton(eta_F2D_block_bsp_), hfbsetting(hfbsetting_) {}

    /**
     * @brief Accumulate particle-hole fields using direct contraction.
     * @math Γ_n += v_same ρ_n + v_cross ρ_p; n ↔ p.
     * @output Accumulated neutron and proton Gamma matrices.
     * @note Requires matching species bases and thread-safe callbacks.
     */
    void add_Gamma_from_Element(const GammaElementFunc& read_element_Func);

    /**
     * @brief Initialize one-body fields in the derived model.
     * @math (N,Z) → (h₀,n,h₀,p).
     * @output Initialized neutron and proton one-body fields.
     */
    virtual void initialize_h0() = 0;

    /**
     * @brief Initialize particle-hole fields in the derived model.
     * @math (N,Z) → (Γ_n,Γ_p)_initial.
     * @output Initialized neutron and proton Gamma matrices.
     */
    virtual void initialize_Gamma() = 0;

    /**
     * @brief Update both species using representative densities.
     * @math {ρ_q,b⁺⁺,η_q,b} → {Γ_q,b⁺⁺}.
     * @output Overwritten neutron and proton Gamma matrices.
     */
    virtual void update_Gamma() = 0;

    /**
     * @brief Print the current iteration summary.
     * @math (i,ε,α) → stdout.
     * @output Iteration diagnostics.
     */
    virtual void print_abstract(int iteration_I, double error_F, double mixing_F) const {}

    /**
     * @brief Iterate HF-BCS using modified Broyden mixing.
     * @math h_q = h₀,q + Γ_q; (N,Z) → HFBCS_converged.
     * @output Updated neutron and proton fields and solutions.
     * @note Requires initialized h₀,Γ and species pairing strengths G.
     * @note Continuation starts fresh Broyden history; temperature is zero.
     */
    void iterate(bool useCurrentFields_B = false);
};

/**
 * @brief Accumulate particle-hole fields using parallel direct contraction.
 * @math Γ_n += v_same ρ_n + v_cross ρ_p; n ↔ p.
 * @output Accumulated neutron and proton Gamma matrices.
 * @note Requires matching species bases and thread-safe callbacks.
 */
inline void HFBCSNucleus::add_Gamma_from_Element(const GammaElementFunc& read_element_Func) {
    assert(read_element_Func);
    assert(hfb_neutron.Nbsp_I1D_block == hfb_proton.Nbsp_I1D_block);

    const auto add_Gamma_Func = [&](int block13_I, int bsp1_I, int bsp3_I) {
        double GammaSameN_F = 0.0;
        double GammaCrossN_F = 0.0;
        double GammaSameP_F = 0.0;
        double GammaCrossP_F = 0.0;
        for (int block24_I = 0; block24_I < hfb_neutron.Nblock_I; ++block24_I) {
            const auto& solutionN = hfb_neutron.solutions[block24_I];
            const auto& solutionP = hfb_proton.solutions[block24_I];
            const int Nbsp_I = hfb_neutron.Nbsp_I1D_block[block24_I];
            for (int bsp2_I = 0; bsp2_I < Nbsp_I; ++bsp2_I) {
                for (int bsp4_I = 0; bsp4_I < Nbsp_I; ++bsp4_I) {
                    const Element element = read_element_Func(block13_I, block24_I, bsp1_I, bsp2_I, bsp3_I, bsp4_I);
                    const double rhoNPosPos_F = solutionN.rhoPosPos_F2D_bsp_bsp(bsp4_I, bsp2_I);
                    const double rhoPPosPos_F = solutionP.rhoPosPos_F2D_bsp_bsp(bsp4_I, bsp2_I);
                    const double rhoNNegNeg_F = solutionN.eta_F1D_bsp(bsp4_I) * solutionN.eta_F1D_bsp(bsp2_I) * rhoNPosPos_F;
                    const double rhoPNegNeg_F = solutionP.eta_F1D_bsp(bsp4_I) * solutionP.eta_F1D_bsp(bsp2_I) * rhoPPosPos_F;
                    GammaSameN_F += element.vSamePosPosPosPos_F * rhoNPosPos_F + element.vSamePosNegPosNeg_F * rhoNNegNeg_F;
                    GammaCrossN_F += element.vCrossPosPosPosPos_F * rhoPPosPos_F + element.vCrossPosNegPosNeg_F * rhoPNegNeg_F;
                    GammaSameP_F += element.vSamePosPosPosPos_F * rhoPPosPos_F + element.vSamePosNegPosNeg_F * rhoPNegNeg_F;
                    GammaCrossP_F += element.vCrossPosPosPosPos_F * rhoNPosPos_F + element.vCrossPosNegPosNeg_F * rhoNNegNeg_F;
                }
            }
        }
        hfb_neutron.fields[block13_I].GammaPosPos_F2D_bsp_bsp(bsp1_I, bsp3_I) += GammaSameN_F;
        hfb_neutron.fields[block13_I].GammaPosPos_F2D_bsp_bsp(bsp1_I, bsp3_I) += GammaCrossN_F;
        hfb_proton.fields[block13_I].GammaPosPos_F2D_bsp_bsp(bsp1_I, bsp3_I) += GammaSameP_F;
        hfb_proton.fields[block13_I].GammaPosPos_F2D_bsp_bsp(bsp1_I, bsp3_I) += GammaCrossP_F;
    };

    #pragma omp parallel
    {
        #pragma omp single
        for (int block13_I = 0; block13_I < hfb_neutron.Nblock_I; ++block13_I) {
            const int Nbsp_I = hfb_neutron.Nbsp_I1D_block[block13_I];
            for (int bsp1_I = 0; bsp1_I < Nbsp_I; ++bsp1_I) {
                for (int bsp3_I = 0; bsp3_I < Nbsp_I; ++bsp3_I) {
                    #pragma omp task firstprivate(block13_I, bsp1_I, bsp3_I)
                    {
                        add_Gamma_Func(block13_I, bsp1_I, bsp3_I);
                    }
                }
            }
        }
    }
}

/**
 * @brief Iterate neutron-proton HF-BCS using modified Broyden mixing.
 * @math G:h → (λ,Δ,ρ,κ) → Γ → h₀+Γ.
 * @output Updated fields and solutions at the final iterate.
 * @note Requires initialized fields; each call resets mixing history.
 */
inline void HFBCSNucleus::iterate(bool useCurrentFields_B) {
    assert(hfbsetting.temperature_F == 0.0);
    assert(!hfbsetting.useLipkinNogami_B);
    const double EspCut_F = hfbsetting.useEspCut_B ? hfbsetting.EspCut_F : std::numeric_limits<double>::infinity();
    assert(std::isfinite(hfbsetting.accuracy_F) && hfbsetting.accuracy_F > 0.0);
    assert(hfbsetting.NiterationsMax_I > 0);
    assert(std::isfinite(hfbsetting.mixingMin_F) && std::isfinite(hfbsetting.mixingMax_F) && hfbsetting.mixingMin_F > 0.0 && hfbsetting.mixingMin_F <= hfbsetting.mixingMax_F && hfbsetting.mixingMax_F <= 1.0);

    // Each block packs h = h₀ + Γ.
    int Npacked_I = 0;
    for (const HFBCS* hfb_Ptr : {&hfb_neutron, &hfb_proton}) {
        for (int Nbsp_I : hfb_Ptr->Nbsp_I1D_block) {Npacked_I += Nbsp_I * Nbsp_I;}
    }
    assert(Npacked_I > 0);

    Eigen::VectorXd x_F1D_packed{};
    Eigen::VectorXd Gx_F1D_packed{};
    x_F1D_packed.resize(Npacked_I);
    Gx_F1D_packed.resize(Npacked_I);

    // pack: species → blocks → full h matrices.
    const auto pack_h_Func = [&](Eigen::VectorXd& data_F1D_packed) {
        int packed_I = 0;
        for (const HFBCS* hfb_Ptr : {&hfb_neutron, &hfb_proton}) {
            for (int block_I = 0; block_I < hfb_Ptr->Nblock_I; ++block_I) {
                const int Nbsp_I = hfb_Ptr->Nbsp_I1D_block[block_I];
                const HFBCSBlockField& field = hfb_Ptr->fields[block_I];
                for (int column_I = 0; column_I < Nbsp_I; ++column_I) {
                    for (int row_I = 0; row_I < Nbsp_I; ++row_I) {
                        data_F1D_packed(packed_I++) = field.h0PosPos_F2D_bsp_bsp(row_I, column_I) + field.GammaPosPos_F2D_bsp_bsp(row_I, column_I);
                    }
                }
            }
        }
        assert(packed_I == Npacked_I);
    };

    // unpack: species → blocks → full h matrices.
    const auto unpack_h_Func = [&](const Eigen::VectorXd& data_F1D_packed) {
        int packed_I = 0;
        for (HFBCS* hfb_Ptr : {&hfb_neutron, &hfb_proton}) {
            for (int block_I = 0; block_I < hfb_Ptr->Nblock_I; ++block_I) {
                const int Nbsp_I = hfb_Ptr->Nbsp_I1D_block[block_I];
                HFBCSBlockField& field = hfb_Ptr->fields[block_I];
                for (int column_I = 0; column_I < Nbsp_I; ++column_I) {
                    for (int row_I = 0; row_I < Nbsp_I; ++row_I) {
                        field.GammaPosPos_F2D_bsp_bsp(row_I, column_I) = data_F1D_packed(packed_I++) - field.h0PosPos_F2D_bsp_bsp(row_I, column_I);
                    }
                }
            }
        }
        assert(packed_I == Npacked_I);
    };

    // G:h → (λ,Δ,E,u,v,ρ,κ) → Γ → h₀+Γ.
    const double lambdaAccuracyMin_F = hfbsetting.accuracy_F * 1.0e-6;
    double lambdaAccuracy_F = hfbsetting.accuracy_F;
    const auto calc_Gx_Func = [&](const Eigen::VectorXd& x_F1D_packed_, Eigen::VectorXd& Gx_F1D_packed_) {
        unpack_h_Func(x_F1D_packed_);
        hfb_neutron.search_lambda(EspCut_F, lambdaAccuracy_F);
        hfb_proton.search_lambda(EspCut_F, lambdaAccuracy_F);
        update_Gamma();
        pack_h_Func(Gx_F1D_packed_);
    };

    // Initial fields → G(x₀); fresh x₀ = 0.
    pack_h_Func(Gx_F1D_packed);
    x_F1D_packed.setZero();
    if (useCurrentFields_B) {
        x_F1D_packed = Gx_F1D_packed;
        calc_Gx_Func(x_F1D_packed, Gx_F1D_packed);
    }

    // (x₀,G(x₀)) → Broyden history.
    double alpha_F = hfbsetting.mixingMin_F;
    BroydenIterator broyden_(std::min(7, Npacked_I), calc_Gx_Func, hfbsetting.mixingMin_F, x_F1D_packed, Gx_F1D_packed);
    print_abstract(0, 0.0, hfbsetting.mixingMin_F);

    // ||G(x_i)-x_i||∞ → ε_i; adaptive α and λ tolerance.
    double errorPrevious_F = 1.0;
    for (int iteration_I = 1; iteration_I <= hfbsetting.NiterationsMax_I; ++iteration_I) {
        const double error_F = broyden_.iterate(calc_Gx_Func, alpha_F);
        print_abstract(iteration_I, error_F, alpha_F);
        if (std::isfinite(error_F) && error_F <= hfbsetting.accuracy_F) {break;}
        if (std::isfinite(error_F) && error_F < errorPrevious_F) {
            alpha_F = std::min(hfbsetting.mixingMax_F, alpha_F * 1.10);
            errorPrevious_F = error_F;
            continue;
        }
        alpha_F = hfbsetting.mixingMin_F;
        const bool tightenLambdaAccuracy_B = lambdaAccuracy_F > lambdaAccuracyMin_F * (1.0 + 1.0e-12);
        if (iteration_I > 1 && tightenLambdaAccuracy_B) {
            lambdaAccuracy_F = std::max(lambdaAccuracyMin_F, lambdaAccuracy_F * 0.1);
        }
        errorPrevious_F = error_F;
    }

    // Densities correspond to x_final, not G(x_final).
    unpack_h_Func(broyden_.xnext_F1D_i);
}
