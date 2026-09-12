/**
 * @file    hfb_full.hpp
 * @author  cailea
 * @date    2026-09-10
 * @brief   Real finite-temperature HFB matrix calculations.
 */

#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <initializer_list>
#include <map>

#include <Eigen/Core>
#include <Eigen/Eigenvalues>

#include "root.hpp"
#include "root_broyden.hpp"

struct HFBSolution {
    Eigen::VectorXd Eqp_F1D_qp{};
    Eigen::VectorXd f_F1D_qp{};
    Eigen::MatrixXd U_F2D_sp_qp{};
    Eigen::MatrixXd V_F2D_sp_qp{};
    Eigen::MatrixXd rho_F2D_sp_sp{};
    Eigen::MatrixXd kappa_F2D_sp_sp{};
};

struct HFBField {
    Eigen::MatrixXd h0_F2D_sp_sp{};
    Eigen::MatrixXd Gamma_F2D_sp_sp{};
    Eigen::MatrixXd Delta_F2D_sp_sp{};
};

class HFB {
public:
    using BlockingFunc = std::function<double(HFBSolution& solution, bool updateTracking_B)>;
    // Antisymmetrized matrix elements v̄_{12;34}.
    using GammaElementFunc = std::function<double(int sp1_I, int sp2_I, int sp3_I, int sp4_I)>;
    using DeltaElementFunc = std::function<double(int sp1_I, int sp2_I, int sp3_I, int sp4_I)>;

    int TargetN_I = 0; // Target particle number.
    double lambda_F = -7.0; // Fermi energy [MeV].
    double lambda2_F = 0.0; // Lipkin-Nogami λ₂ [MeV].
    double temperature_F = 0.0;
    double EspCut_F = 60.0; // Equivalent single-particle energy cutoff [MeV].

    int Nsp_I = 0;

    HFBField hfb_field{};
    HFBSolution hfb_solution{};
    BlockingFunc blocking_Func{}; // updateTracking_B: commit blocking trackers.

    Eigen::MatrixXd H_F2D_2sp_2sp{};
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> H_eigensolver{};

public:
    /**
     * @brief  Initialize the single-species HFB workspace.
     * @math   h = h₀ + Γ.
     * @output Allocated matrices and eigensolver workspace.
     */
    HFB(int Nsp_I_) {
        assert(Nsp_I_ > 0);

        Nsp_I = Nsp_I_;

        // h₀, Γ, Δ ∈ ℝ^{Nsp×Nsp}.
        hfb_field.h0_F2D_sp_sp.resize(Nsp_I, Nsp_I);
        hfb_field.Gamma_F2D_sp_sp.resize(Nsp_I, Nsp_I);
        hfb_field.Delta_F2D_sp_sp.resize(Nsp_I, Nsp_I);

        // E, f ∈ ℝ^{Nsp}; U, V ∈ ℝ^{Nsp×Nsp}.
        hfb_solution.Eqp_F1D_qp.resize(Nsp_I);
        hfb_solution.f_F1D_qp.resize(Nsp_I);
        hfb_solution.U_F2D_sp_qp.resize(Nsp_I, Nsp_I);
        hfb_solution.V_F2D_sp_qp.resize(Nsp_I, Nsp_I);

        // ρ, κ ∈ ℝ^{Nsp×Nsp}.
        hfb_solution.rho_F2D_sp_sp.resize(Nsp_I, Nsp_I);
        hfb_solution.kappa_F2D_sp_sp.resize(Nsp_I, Nsp_I);

        // H ∈ ℝ^{2Nsp×2Nsp}.
        H_F2D_2sp_2sp.resize(2 * Nsp_I, 2 * Nsp_I);
        H_eigensolver = Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>(2 * Nsp_I);

    }

    /**
     * @brief  Solve thermal HFB using self-adjoint eigendecomposition.
     * @math   H = [h₀+Γ-λI, Δ; -Δ, -h₀-Γ+λI].
     * @math   H[U;V] = [U;V]E; f = 1/(1+exp(E/T)).
     * @math   ρ = V(1-f)Vᵀ + UfUᵀ.
     * @math   κ = V(1-f)Uᵀ + UfVᵀ.
     * @output Updated H, U, V, Eqp, f, rho, and kappa.
     * @output Returned mean particle number N = Tr(ρ).
     * @note   ε = λ+E(1-2‖V‖²) cuts ρ,κ; Ecut = ∞ disables.
     */
    double update_UV_E_rho_kappa();

    /**
     * @brief  Search chemical potential using bracket expansion and Brent.
     * @math   Tr(ρ_blocked(λ)) = TargetN.
     * @output Updated chemical potential, solution, and densities.
     */
    void search_lambda(double lambdaTolerance_F);

    /**
     * @brief Accumulate particle-hole fields by direct matrix-element contraction.
     * @math Γ_{13} += Σ_{24}v̄_{12;34}ρ_source,42.
     * @output Accumulated Gamma; source supplies the contracted species.
     */
    void add_Gamma_from_Element(const HFB& source_, const GammaElementFunc& read_element_Func);

    /**
     * @brief Accumulate pairing fields by direct matrix-element contraction.
     * @math Δ_{12} += ½Σ_{34}v̄_{12;34}κ_{34}.
     * @output Accumulated Delta from this species' pairing tensor.
     */
    void add_Delta_from_Element(const DeltaElementFunc& read_element_Func);
};


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
    virtual void print_abstract(int iteration_I, double error_F, double mixing_F) {}

    /**
     * @brief  Iterate HFB using modified Broyden mixing.
     * @note   Requires initialized fields; continuation starts fresh Broyden history.
     * @math   (N,Z) → HFB_converged.
     * @output Updated neutron and proton fields and solutions.
     */
    void iterate(bool useCurrentFields_B = false);
};

inline void HFB::add_Gamma_from_Element(const HFB& source_, const GammaElementFunc& read_element_Func) {
    assert(read_element_Func);
    assert(source_.hfb_solution.rho_F2D_sp_sp.rows() == source_.Nsp_I && source_.hfb_solution.rho_F2D_sp_sp.cols() == source_.Nsp_I);

    // Γ_{13} += Σ_{24}v̄_{12;34}ρ_source,42.
    for (int sp1_I = 0; sp1_I < Nsp_I; ++sp1_I) {
        for (int sp3_I = 0; sp3_I < Nsp_I; ++sp3_I) {
            double Gamma13_F = 0.0;
            for (int sp2_I = 0; sp2_I < source_.Nsp_I; ++sp2_I) {
                for (int sp4_I = 0; sp4_I < source_.Nsp_I; ++sp4_I) {
                    Gamma13_F += read_element_Func(sp1_I, sp2_I, sp3_I, sp4_I) * source_.hfb_solution.rho_F2D_sp_sp(sp4_I, sp2_I);
                }
            }
            hfb_field.Gamma_F2D_sp_sp(sp1_I, sp3_I) += Gamma13_F;
        }
    }
}

inline void HFB::add_Delta_from_Element(const DeltaElementFunc& read_element_Func) {
    assert(read_element_Func);
    assert(hfb_solution.kappa_F2D_sp_sp.rows() == Nsp_I && hfb_solution.kappa_F2D_sp_sp.cols() == Nsp_I);

    // Δ_{12} += ½Σ_{34}v̄_{12;34}κ_{34}.
    for (int sp1_I = 0; sp1_I < Nsp_I; ++sp1_I) {
        for (int sp2_I = 0; sp2_I < Nsp_I; ++sp2_I) {
            double Delta12_F = 0.0;
            for (int sp3_I = 0; sp3_I < Nsp_I; ++sp3_I) {
                for (int sp4_I = 0; sp4_I < Nsp_I; ++sp4_I) {
                    Delta12_F += 0.5 * read_element_Func(sp1_I, sp2_I, sp3_I, sp4_I) * hfb_solution.kappa_F2D_sp_sp(sp3_I, sp4_I);
                }
            }
            hfb_field.Delta_F2D_sp_sp(sp1_I, sp2_I) += Delta12_F;
        }
    }
}

inline double HFB::update_UV_E_rho_kappa() {
    assert(temperature_F >= 0.0);
    assert(hfb_field.h0_F2D_sp_sp.rows() == Nsp_I && hfb_field.h0_F2D_sp_sp.cols() == Nsp_I);
    assert(hfb_field.Gamma_F2D_sp_sp.rows() == Nsp_I && hfb_field.Gamma_F2D_sp_sp.cols() == Nsp_I);
    assert(hfb_field.Delta_F2D_sp_sp.rows() == Nsp_I && hfb_field.Delta_F2D_sp_sp.cols() == Nsp_I);
    assert(hfb_field.h0_F2D_sp_sp.isApprox(hfb_field.h0_F2D_sp_sp.transpose(), 1.0e-12));
    assert(hfb_field.Gamma_F2D_sp_sp.isApprox(hfb_field.Gamma_F2D_sp_sp.transpose(), 1.0e-12));
    assert(hfb_field.Delta_F2D_sp_sp.isApprox(-hfb_field.Delta_F2D_sp_sp.transpose(), 1.0e-12));

    // H₁₁ = h₀ + Γ - λI; Γ includes LN.
    H_F2D_2sp_2sp.topLeftCorner(Nsp_I, Nsp_I) = hfb_field.h0_F2D_sp_sp + hfb_field.Gamma_F2D_sp_sp;
    H_F2D_2sp_2sp.topLeftCorner(Nsp_I, Nsp_I).diagonal().array() -= lambda_F;

    // H₁₂ = Δ; H₂₁ = -Δ; H₂₂ = -H₁₁.
    H_F2D_2sp_2sp.topRightCorner(Nsp_I, Nsp_I) = hfb_field.Delta_F2D_sp_sp;
    H_F2D_2sp_2sp.bottomLeftCorner(Nsp_I, Nsp_I) = -hfb_field.Delta_F2D_sp_sp;
    H_F2D_2sp_2sp.bottomRightCorner(Nsp_I, Nsp_I) = -H_F2D_2sp_2sp.topLeftCorner(Nsp_I, Nsp_I);

    // H[U;V] = [U;V]E.
    H_eigensolver.compute(H_F2D_2sp_2sp);
    assert(H_eigensolver.info() == Eigen::Success);
    const Eigen::VectorXd& eigenvalues_F1D_state = H_eigensolver.eigenvalues();
    const Eigen::MatrixXd& eigenvectors_F2D_2sp_state = H_eigensolver.eigenvectors();

    // E₁ ≤ ⋯ ≤ E₂Nsp; retain Nsp positive-energy modes.
    assert(eigenvalues_F1D_state(Nsp_I - 1) < 0.0);
    assert(eigenvalues_F1D_state(Nsp_I) > 0.0);
    hfb_solution.Eqp_F1D_qp = eigenvalues_F1D_state.tail(Nsp_I);
    hfb_solution.U_F2D_sp_qp = eigenvectors_F2D_2sp_state.topRightCorner(Nsp_I, Nsp_I);
    hfb_solution.V_F2D_sp_qp = eigenvectors_F2D_2sp_state.bottomRightCorner(Nsp_I, Nsp_I);

    // ε = λ + E(1-2‖V‖²); ε > Ecut + tail → factors = 0.
    const double EspCutTail_F = std::log(1.0 / 1.0e-6 - 1.0) / 100.0;
    Eigen::VectorXd factorU_F1D_qp{};
    Eigen::VectorXd factorV_F1D_qp{};
    factorU_F1D_qp.resize(Nsp_I);
    factorV_F1D_qp.resize(Nsp_I);
    factorU_F1D_qp.setZero();
    factorV_F1D_qp.setZero();

    // f = (1-tanh(E/(2T)))/2; active factors = (f,1-f).
    for (int qp_I = 0; qp_I < Nsp_I; ++qp_I) {
        const double Eqp_F = hfb_solution.Eqp_F1D_qp(qp_I);
        hfb_solution.f_F1D_qp(qp_I) = temperature_F > 1.0e-12 ? 0.5 * (1.0 - std::tanh(0.5 * Eqp_F / temperature_F)) : 0.0;
        const double Esp_F = lambda_F + Eqp_F * (1.0 - 2.0 * hfb_solution.V_F2D_sp_qp.col(qp_I).squaredNorm());
        if (Esp_F > EspCut_F + EspCutTail_F) {continue;}
        factorU_F1D_qp(qp_I) = hfb_solution.f_F1D_qp(qp_I);
        factorV_F1D_qp(qp_I) = 1.0 - hfb_solution.f_F1D_qp(qp_I);
    }

    // ρ = V(1-f)Vᵀ + UfUᵀ.
    hfb_solution.rho_F2D_sp_sp.noalias() = hfb_solution.V_F2D_sp_qp * factorV_F1D_qp.asDiagonal() * hfb_solution.V_F2D_sp_qp.transpose();
    hfb_solution.rho_F2D_sp_sp.noalias() += hfb_solution.U_F2D_sp_qp * factorU_F1D_qp.asDiagonal() * hfb_solution.U_F2D_sp_qp.transpose();

    // κ = V(1-f)Uᵀ + UfVᵀ.
    hfb_solution.kappa_F2D_sp_sp.noalias() = hfb_solution.V_F2D_sp_qp * factorV_F1D_qp.asDiagonal() * hfb_solution.U_F2D_sp_qp.transpose();
    hfb_solution.kappa_F2D_sp_sp.noalias() += hfb_solution.U_F2D_sp_qp * factorU_F1D_qp.asDiagonal() * hfb_solution.V_F2D_sp_qp.transpose();

    return hfb_solution.rho_F2D_sp_sp.trace();
}

inline void HFB::search_lambda(double lambdaTolerance_F) {
    assert(TargetN_I >= 0 && TargetN_I <= Nsp_I);
    assert(std::isfinite(lambda_F));
    assert(std::isfinite(lambdaTolerance_F) && lambdaTolerance_F > 0.0);
    const double Ntolerance_F = std::max(1.0e-8, 1.0e-10 * std::max(1, TargetN_I));
    std::map<double, double> NerrorByLambda_Map{};

    // Fixed fields; trial blocking preserves external trackers.
    const auto calc_N_Func = [&](double lambdaTrial_F, bool updateTracking_B) {
        lambda_F = lambdaTrial_F;
        const double Ncalc_F = update_UV_E_rho_kappa();
        if (blocking_Func) {return blocking_Func(hfb_solution, updateTracking_B);}
        return Ncalc_F;
    };

    // λ → N(λ)-TargetN; cache trial residuals.
    const auto calc_Nerror_Func = [&](double lambdaTrial_F) {
        const auto NerrorIt_ = NerrorByLambda_Map.find(lambdaTrial_F);
        if (NerrorIt_ != NerrorByLambda_Map.end()) {return NerrorIt_->second;}
        const double Nerror_F = calc_N_Func(lambdaTrial_F, false) - static_cast<double>(TargetN_I);
        assert(std::isfinite(Nerror_F));
        NerrorByLambda_Map.emplace(lambdaTrial_F, Nerror_F);
        return Nerror_F;
    };

    // λ_initial → N(λ_initial); commit the accepted blocking state.
    const double lambdaInitial_F = lambda_F;
    const double NerrorInitial_F = calc_Nerror_Func(lambdaInitial_F);
    if (std::abs(NerrorInitial_F) <= Ntolerance_F) {
        calc_N_Func(lambdaInitial_F, true);
        return;
    }

    // [λ_min,λ_max] = [λ_initial-2,λ_initial+2].
    double lambdaMax_F = lambdaInitial_F + 2.0;
    double lambdaMin_F = lambdaInitial_F - 2.0;
    double NerrorMin_F = calc_Nerror_Func(lambdaMin_F);
    double NerrorMax_F = calc_Nerror_Func(lambdaMax_F);
    const double NerrorSlope_F = (NerrorMax_F - NerrorMin_F) / (lambdaMax_F - lambdaMin_F);

    // f(λ_min)f(λ_max)>0 → expand one boundary.
    const bool shouldShiftLambdaMin_B = NerrorMin_F * NerrorSlope_F > 0.0;
    double& lambdaBound_F = shouldShiftLambdaMin_B ? lambdaMin_F : lambdaMax_F;
    double& NerrorBound_F = shouldShiftLambdaMin_B ? NerrorMin_F : NerrorMax_F;
    const double lambdaStep_F = shouldShiftLambdaMin_B ? -5.0 : 5.0;
    constexpr int NexpandMax_I = 100;
    int Nexpand_I = 0;
    while (NerrorMin_F * NerrorMax_F > 0.0 && NerrorSlope_F != 0.0 && Nexpand_I < NexpandMax_I) {
        lambdaBound_F += lambdaStep_F;
        NerrorBound_F = calc_Nerror_Func(lambdaBound_F);
        ++Nexpand_I;
    }
    assert(NerrorMin_F * NerrorMax_F <= 0.0);

    // Recompute at λ_root and commit blocking trackers.
    const double lambdaRoot_F = root_brent(calc_Nerror_Func, lambdaMin_F, lambdaMax_F, lambdaTolerance_F);
    assert(std::isfinite(lambdaRoot_F));
    calc_N_Func(lambdaRoot_F, true);
}

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
