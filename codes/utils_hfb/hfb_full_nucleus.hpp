/**
 * @file    hfb_full_nucleus.hpp
 * @author  cailea
 * @date    2026-09-12
 * @brief   Joint neutron-proton HFB iteration.
 */

#pragma once

#include "hfb_full.hpp"
#include "root_broyden.hpp"

class HFBNucleus {
public:
    HFB hfb_neutron;
    HFB hfb_proton;

    double accuracy_F = 1.0e-5;
    double mixingMin_F = 0.20;
    double mixingMax_F = 0.90;
    int NiterationsMax_I = 100;

public:
    /**
     * @brief  Initialize species dimensions and workspaces.
     * @math   Nsp_n, Nsp_p.
     * @output Allocated solutions, fields, and workspaces.
     */
    HFBNucleus(int NspN_I_, int NspP_I_)
    : hfb_neutron(NspN_I_), hfb_proton(NspP_I_) {}

    /**
     * @brief  Initialize one-body fields in the derived model.
     * @math   (N,Z) → (h₀,n,h₀,p).
     * @output Initialized neutron and proton one-body fields.
     */
    virtual void initialize_h0() = 0;

    /**
     * @brief  Initialize HFB fields in the derived model.
     * @math   (N,Z) → (Γ_n,Δ_n,Γ_p,Δ_p)_initial.
     * @output Initialized Gamma and Delta for both species.
     * @note   Include model-specific LN corrections when enabled.
     */
    virtual void initialize_GammaDelta() = 0;

    /**
     * @brief  Update both species using their joint densities.
     * @math   (ρ_n,κ_n,ρ_p,κ_p) → (Γ_n,Δ_n,Γ_p,Δ_p).
     * @output Overwritten neutron and proton fields.
     * @note   Rebuild bare Gamma, then add LN once when enabled.
     */
    virtual void update_Gamma_Delta() = 0;

    /**
     * @brief  Print the current iteration summary.
     * @math   (i,ε,α) → stdout.
     * @output Iteration diagnostics.
     */
    virtual void print_abstract(int iteration_I, double error_F, double mixing_F) const {}

    /**
     * @brief  Iterate HFB using modified Broyden mixing.
     * @note   Requires initialized fields; continuation starts fresh Broyden history.
     * @math   (N,Z) → HFB_converged.
     * @output Updated neutron and proton fields and solutions.
     */
    void iterate(bool useCurrentFields_B = false);
};

inline void HFBNucleus::iterate(bool useCurrentFields_B) {
    assert(std::isfinite(accuracy_F) && accuracy_F > 0.0);
    assert(NiterationsMax_I > 0);
    assert(std::isfinite(mixingMin_F) && std::isfinite(mixingMax_F) && mixingMin_F > 0.0 && mixingMin_F <= mixingMax_F && mixingMax_F <= 1.0);

    const int Npacked_I = 2 * (hfb_neutron.Nsp_I * hfb_neutron.Nsp_I + hfb_proton.Nsp_I * hfb_proton.Nsp_I);
    assert(Npacked_I >= 7);

    Eigen::VectorXd x_F1D_packed{};
    Eigen::VectorXd Gx_F1D_packed{};
    x_F1D_packed.resize(Npacked_I);
    Gx_F1D_packed.resize(Npacked_I);

    // x ← (h₀,n+Γ_n) ⊕ Δ_n ⊕ (h₀,p+Γ_p) ⊕ Δ_p.
    const auto pack_h_Delta_Func = [&](Eigen::VectorXd& data_F1D_packed) {
        int packed_I = 0;
        for (const HFB* hfb_Ptr : {&hfb_neutron, &hfb_proton}) {
            for (int column_I = 0; column_I < hfb_Ptr->Nsp_I; ++column_I) {
                for (int row_I = 0; row_I < hfb_Ptr->Nsp_I; ++row_I) {
                    data_F1D_packed(packed_I++) = hfb_Ptr->hfb_field.h0_F2D_sp_sp(row_I, column_I) + hfb_Ptr->hfb_field.Gamma_F2D_sp_sp(row_I, column_I);
                }
            }
            for (int column_I = 0; column_I < hfb_Ptr->Nsp_I; ++column_I) {
                for (int row_I = 0; row_I < hfb_Ptr->Nsp_I; ++row_I) {
                    data_F1D_packed(packed_I++) = hfb_Ptr->hfb_field.Delta_F2D_sp_sp(row_I, column_I);
                }
            }
        }
        assert(packed_I == Npacked_I);
    };

    // x → Γ = h-h₀; Δ.
    const auto unpack_h_Delta_Func = [&](const Eigen::VectorXd& data_F1D_packed) {
        int packed_I = 0;
        for (HFB* hfb_Ptr : {&hfb_neutron, &hfb_proton}) {
            for (int column_I = 0; column_I < hfb_Ptr->Nsp_I; ++column_I) {
                for (int row_I = 0; row_I < hfb_Ptr->Nsp_I; ++row_I) {
                    hfb_Ptr->hfb_field.Gamma_F2D_sp_sp(row_I, column_I) = data_F1D_packed(packed_I++) - hfb_Ptr->hfb_field.h0_F2D_sp_sp(row_I, column_I);
                }
            }
            for (int column_I = 0; column_I < hfb_Ptr->Nsp_I; ++column_I) {
                for (int row_I = 0; row_I < hfb_Ptr->Nsp_I; ++row_I) {
                    hfb_Ptr->hfb_field.Delta_F2D_sp_sp(row_I, column_I) = data_F1D_packed(packed_I++);
                }
            }
        }
        assert(packed_I == Npacked_I);
    };

    // G:x → (λ,U,V,E,ρ,κ) → (Γ,Δ) → (h₀+Γ,Δ).
    const double lambdaToleranceMin_F = accuracy_F * 1.0e-6;
    double lambdaTolerance_F = accuracy_F;
    const auto calc_Gx_Func = [&](const Eigen::VectorXd& x_F1D_packed_, Eigen::VectorXd& Gx_F1D_packed_) {
        unpack_h_Delta_Func(x_F1D_packed_);
        hfb_neutron.search_lambda(lambdaTolerance_F);
        hfb_proton.search_lambda(lambdaTolerance_F);
        update_Gamma_Delta();
        pack_h_Delta_Func(Gx_F1D_packed_);
    };

    // Initial fields → G(x₀); fresh x₀ = 0.
    pack_h_Delta_Func(Gx_F1D_packed);
    x_F1D_packed.setZero();
    if (useCurrentFields_B) {
        x_F1D_packed = Gx_F1D_packed;
        calc_Gx_Func(x_F1D_packed, Gx_F1D_packed);
    }

    // (x₀,G(x₀)) → Broyden history.
    double alpha_F = mixingMin_F;
    BroydenIterator broyden_(7, calc_Gx_Func, mixingMin_F, x_F1D_packed, Gx_F1D_packed);
    print_abstract(0, 0.0, mixingMin_F);

    // ||G(x_i)-x_i||∞ → ε_i; adaptive α and λ tolerance.
    double errorPrevious_F = 1.0;
    for (int iteration_I = 1; iteration_I <= NiterationsMax_I; ++iteration_I) {
        const double error_F = broyden_.iterate(calc_Gx_Func, alpha_F);
        print_abstract(iteration_I, error_F, alpha_F);
        if (std::isfinite(error_F) && error_F <= accuracy_F) {break;}
        if (std::isfinite(error_F) && error_F < errorPrevious_F) {
            alpha_F = std::min(mixingMax_F, alpha_F * 1.10);
            errorPrevious_F = error_F;
            continue;
        }
        alpha_F = mixingMin_F;
        const bool tightenLambdaTolerance_B = lambdaTolerance_F > lambdaToleranceMin_F * (1.0 + 1.0e-12);
        if (iteration_I > 1 && tightenLambdaTolerance_B) {
            lambdaTolerance_F = std::max(lambdaToleranceMin_F, lambdaTolerance_F * 0.1);
        }
        errorPrevious_F = error_F;
    }
}
