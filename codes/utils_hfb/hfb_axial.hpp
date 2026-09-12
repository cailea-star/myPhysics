/**
 * @file    hfb_axial.hpp
 * @author  cailea
 * @date    2026-09-10
 * @brief   Axial HFB without time-reversal symmetry.
 */

#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <functional>
#include <initializer_list>
#include <map>
#include <numeric>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Eigenvalues>

#include "root.hpp"
#include "root_broyden.hpp"

struct HFBAxialBlockSolution {
    Eigen::VectorXd EqpPos_F1D_bqp{};
    Eigen::VectorXd EqpNeg_F1D_bqp{};
    Eigen::VectorXd fPos_F1D_bqp{};
    Eigen::VectorXd fNeg_F1D_bqp{};

    // X⁺ = [U⁺; V⁻].
    Eigen::MatrixXd UPos_F2D_bsp_bqp{};
    Eigen::MatrixXd VNeg_F2D_bsp_bqp{};

    // X⁻ = [U⁻; V⁺].
    Eigen::MatrixXd UNeg_F2D_bsp_bqp{};
    Eigen::MatrixXd VPos_F2D_bsp_bqp{};

    Eigen::MatrixXd rhoPosPos_F2D_bsp_bsp{};
    Eigen::MatrixXd rhoNegNeg_F2D_bsp_bsp{};
    Eigen::MatrixXd kappaPosNeg_F2D_bsp_bsp{};
};

struct HFBAxialBlockField {
    Eigen::MatrixXd h0PosPos_F2D_bsp_bsp{};
    Eigen::MatrixXd h0NegNeg_F2D_bsp_bsp{};

    Eigen::MatrixXd GammaPosPos_F2D_bsp_bsp{};
    Eigen::MatrixXd GammaNegNeg_F2D_bsp_bsp{};
    Eigen::MatrixXd DeltaPosNeg_F2D_bsp_bsp{};
};

class HFBAxial {
public:
    using BlockingFunc = std::function<double(std::vector<HFBAxialBlockSolution>& solutions, bool updateTracking_B)>;
    // v̄_{12;34}: Γ uses block13, block24; Δ uses block12, block34.
    // bsp indices follow 1,2,3,4.
    // Γ entries: ++++, +-+-, -+-+, ----; Δ entry: +-+-.
    using GammaElementFunc = std::function<std::array<double, 4>(int block13_I, int block24_I, int bsp1_I, int bsp2_I, int bsp3_I, int bsp4_I)>;
    using DeltaElementFunc = std::function<double(int block12_I, int block34_I, int bsp1_I, int bsp2_I, int bsp3_I, int bsp4_I)>;

    int TargetN_I = 0; // Target particle number.
    double lambda_F = -7.0; // Fermi energy [MeV].
    double lambda2_F = 0.0; // Lipkin-Nogami λ₂ [MeV].
    double temperature_F = 0.0;
    double EspCut_F = 60.0; // Equivalent single-particle energy cutoff [MeV].

    int Nblock_I = 0;
    std::vector<int> Nbsp_I1D_block{};

    std::vector<HFBAxialBlockSolution> hfb_axial_solutions{};
    BlockingFunc blocking_Func{}; // updateTracking_B: commit blocking trackers.
    std::vector<HFBAxialBlockField> hfb_axial_fields{};

    std::vector<Eigen::MatrixXd> HPos_F3D_block_2bsp_2bsp{};
    std::vector<Eigen::MatrixXd> HNeg_F3D_block_2bsp_2bsp{};
    std::vector<Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>> HPos_eigensolvers{};
    std::vector<Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>> HNeg_eigensolvers{};

public:
    /**
     * @brief  Allocate single-species axial HFB blocks.
     * @math   dim(bsp_b⁺) = dim(bsp_b⁻) = Nbsp_b.
     * @output Allocated fields, solutions, and eigensolver workspaces.
     */
    HFBAxial(const std::vector<int>& Nbsp_I1D_block_) {
        Nbsp_I1D_block = Nbsp_I1D_block_;
        Nblock_I = static_cast<int>(Nbsp_I1D_block.size());

        // {Nbsp_b} → {field_b, solution_b, workspace_b}.
        hfb_axial_fields.resize(Nblock_I);
        hfb_axial_solutions.resize(Nblock_I);
        HPos_F3D_block_2bsp_2bsp.resize(Nblock_I);
        HNeg_F3D_block_2bsp_2bsp.resize(Nblock_I);
        HPos_eigensolvers.resize(Nblock_I);
        HNeg_eigensolvers.resize(Nblock_I);

        for (int block_I = 0; block_I < Nblock_I; ++block_I) {
            const int Nbsp_I = Nbsp_I1D_block[block_I];
            assert(Nbsp_I > 0);
            HFBAxialBlockField& field = hfb_axial_fields[block_I];
            HFBAxialBlockSolution& solution = hfb_axial_solutions[block_I];

            // h₀, Γ, Δ ∈ ℝ^{Nbsp×Nbsp}.
            field.h0PosPos_F2D_bsp_bsp.resize(Nbsp_I, Nbsp_I);
            field.h0NegNeg_F2D_bsp_bsp.resize(Nbsp_I, Nbsp_I);
            field.GammaPosPos_F2D_bsp_bsp.resize(Nbsp_I, Nbsp_I);
            field.GammaNegNeg_F2D_bsp_bsp.resize(Nbsp_I, Nbsp_I);
            field.DeltaPosNeg_F2D_bsp_bsp.resize(Nbsp_I, Nbsp_I);

            // ρ, κ ∈ ℝ^{Nbsp×Nbsp}.
            solution.rhoPosPos_F2D_bsp_bsp.resize(Nbsp_I, Nbsp_I);
            solution.rhoNegNeg_F2D_bsp_bsp.resize(Nbsp_I, Nbsp_I);
            solution.kappaPosNeg_F2D_bsp_bsp.resize(Nbsp_I, Nbsp_I);

            // U±, V± ∈ ℝ^{Nbsp×Nbsp} initially.
            solution.UPos_F2D_bsp_bqp.resize(Nbsp_I, Nbsp_I);
            solution.VNeg_F2D_bsp_bqp.resize(Nbsp_I, Nbsp_I);
            solution.UNeg_F2D_bsp_bqp.resize(Nbsp_I, Nbsp_I);
            solution.VPos_F2D_bsp_bqp.resize(Nbsp_I, Nbsp_I);

            // H± ∈ ℝ^{2Nbsp×2Nbsp}.
            HPos_F3D_block_2bsp_2bsp[block_I].resize(2 * Nbsp_I, 2 * Nbsp_I);
            HNeg_F3D_block_2bsp_2bsp[block_I].resize(2 * Nbsp_I, 2 * Nbsp_I);
            HPos_eigensolvers[block_I] = Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>(2 * Nbsp_I);
            HNeg_eigensolvers[block_I] = Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>(2 * Nbsp_I);
        }
    }

    /**
     * @brief  Solve thermal branches using real symmetric eigendecomposition.
     * @math   H⁺ = [h⁺⁺-λI, Δ⁺⁻; (Δ⁺⁻)ᵀ, -h⁻⁻+λI].
     * @math   H⁻ = [h⁻⁻-λI, -(Δ⁺⁻)ᵀ; -Δ⁺⁻, -h⁺⁺+λI].
     * @output Updated quasiparticle solutions and densities.
     * @output Mean particle number N = Σ_b Tr(ρ_b⁺⁺+ρ_b⁻⁻).
     * @note   ε = λ+E(1-2‖V‖²) cuts ρ,κ; Ecut = ∞ disables.
     */
    double update_UV_E_rho_kappa();

    /**
     * @brief  Search chemical potential using bracket expansion and Brent.
     * @math   N_blocked(λ) = TargetN.
     * @output Updated chemical potential, blocked solutions, and densities.
     */
    void search_lambda(double lambdaTolerance_F);

    /**
     * @brief Accumulate particle-hole fields by direct matrix-element contraction.
     * @math Γ_{13} += Σ_{24}v̄_{12;34}ρ_source,42.
     * @output Accumulated Gamma; source supplies the contracted species.
     */
    void add_Gamma_from_Element(const HFBAxial& source_, const GammaElementFunc& read_element_Func);

    /**
     * @brief Accumulate pairing fields by direct matrix-element contraction.
     * @math Δ⁺⁻_{12} += Σ_{34}v̄⁺⁻⁺⁻_{12;34}κ⁺⁻_{34}.
     * @output Accumulated Delta from this species' pairing tensor.
     */
    void add_Delta_from_Element(const DeltaElementFunc& read_element_Func);
};

class HFBAxialNucleus {
public:
    HFBAxial hfb_axial_neutron;
    HFBAxial hfb_axial_proton;

    double accuracy_F = 1.0e-5;
    double mixingMin_F = 0.20;
    double mixingMax_F = 0.90;
    int NiterationsMax_I = 100;

public:
    /**
     * @brief  Initialize species dimensions and workspaces.
     * @math   {Nbsp_n,b, Nbsp_p,b}.
     * @output Allocated solutions, fields, and workspaces.
     */
    HFBAxialNucleus(const std::vector<int>& NbspN_I1D_block_, const std::vector<int>& NbspP_I1D_block_)
    : hfb_axial_neutron(NbspN_I1D_block_), hfb_axial_proton(NbspP_I1D_block_) {}

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
     * @brief  Update both species using all block densities.
     * @math   {ρ_n,b,κ_n,b,ρ_p,b,κ_p,b} → {Γ_n,b,Δ_n,b,Γ_p,b,Δ_p,b}.
     * @output Overwritten neutron and proton Gamma and Delta.
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

inline void HFBAxial::add_Gamma_from_Element(const HFBAxial& source_, const GammaElementFunc& read_element_Func) {
    assert(read_element_Func);

    // Γ⁺⁺_{13} += Σ_{24}(v̄⁺⁺⁺⁺ρ⁺⁺_{42}+v̄⁺⁻⁺⁻ρ⁻⁻_{42}).
    for (int block13_I = 0; block13_I < Nblock_I; ++block13_I) {
        for (int bsp1_I = 0; bsp1_I < Nbsp_I1D_block[block13_I]; ++bsp1_I) {
            for (int bsp3_I = 0; bsp3_I < Nbsp_I1D_block[block13_I]; ++bsp3_I) {
                double Gamma13PosPos_F = 0.0;
                double Gamma13NegNeg_F = 0.0;
                for (int block24_I = 0; block24_I < source_.Nblock_I; ++block24_I) {
                    const auto& solution = source_.hfb_axial_solutions[block24_I];
                    for (int bsp2_I = 0; bsp2_I < source_.Nbsp_I1D_block[block24_I]; ++bsp2_I) {
                        for (int bsp4_I = 0; bsp4_I < source_.Nbsp_I1D_block[block24_I]; ++bsp4_I) {
                            const auto v_F1D_branch = read_element_Func(block13_I, block24_I, bsp1_I, bsp2_I, bsp3_I, bsp4_I);
                            const double rho42PosPos_F = solution.rhoPosPos_F2D_bsp_bsp(bsp4_I, bsp2_I);
                            const double rho42NegNeg_F = solution.rhoNegNeg_F2D_bsp_bsp(bsp4_I, bsp2_I);
                            Gamma13PosPos_F += v_F1D_branch[0] * rho42PosPos_F + v_F1D_branch[1] * rho42NegNeg_F;
                            Gamma13NegNeg_F += v_F1D_branch[2] * rho42PosPos_F + v_F1D_branch[3] * rho42NegNeg_F;
                        }
                    }
                }
                hfb_axial_fields[block13_I].GammaPosPos_F2D_bsp_bsp(bsp1_I, bsp3_I) += Gamma13PosPos_F;
                hfb_axial_fields[block13_I].GammaNegNeg_F2D_bsp_bsp(bsp1_I, bsp3_I) += Gamma13NegNeg_F;
            }
        }
    }
}

inline void HFBAxial::add_Delta_from_Element(const DeltaElementFunc& read_element_Func) {
    assert(read_element_Func);

    // κ⁻⁺ = −(κ⁺⁻)ᵀ and v̄_{12;43} = −v̄_{12;34} cancel ½.
    for (int block12_I = 0; block12_I < Nblock_I; ++block12_I) {
        for (int bsp1_I = 0; bsp1_I < Nbsp_I1D_block[block12_I]; ++bsp1_I) {
            for (int bsp2_I = 0; bsp2_I < Nbsp_I1D_block[block12_I]; ++bsp2_I) {
                double Delta12PosNeg_F = 0.0;
                for (int block34_I = 0; block34_I < Nblock_I; ++block34_I) {
                    const auto& solution = hfb_axial_solutions[block34_I];
                    for (int bsp3_I = 0; bsp3_I < Nbsp_I1D_block[block34_I]; ++bsp3_I) {
                        for (int bsp4_I = 0; bsp4_I < Nbsp_I1D_block[block34_I]; ++bsp4_I) {
                            Delta12PosNeg_F += read_element_Func(block12_I, block34_I, bsp1_I, bsp2_I, bsp3_I, bsp4_I) * solution.kappaPosNeg_F2D_bsp_bsp(bsp3_I, bsp4_I);
                        }
                    }
                }
                hfb_axial_fields[block12_I].DeltaPosNeg_F2D_bsp_bsp(bsp1_I, bsp2_I) += Delta12PosNeg_F;
            }
        }
    }
}

inline double HFBAxial::update_UV_E_rho_kappa() {
    assert(temperature_F >= 0.0);
    assert(static_cast<int>(hfb_axial_fields.size()) == Nblock_I);
    assert(static_cast<int>(hfb_axial_solutions.size()) == Nblock_I);
    double N_F = 0.0;

    // block_b → (E_b, U_b, V_b, f_b, ρ_b, κ_b).
    for (int block_I = 0; block_I < Nblock_I; ++block_I) {
        const int Nbsp_I = Nbsp_I1D_block[block_I];
        const HFBAxialBlockField& field = hfb_axial_fields[block_I];
        HFBAxialBlockSolution& solution = hfb_axial_solutions[block_I];
        Eigen::MatrixXd& HPos_F2D_2bsp_2bsp = HPos_F3D_block_2bsp_2bsp[block_I];
        Eigen::MatrixXd& HNeg_F2D_2bsp_2bsp = HNeg_F3D_block_2bsp_2bsp[block_I];
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>& HPos_eigensolver = HPos_eigensolvers[block_I];
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>& HNeg_eigensolver = HNeg_eigensolvers[block_I];

        assert(field.h0PosPos_F2D_bsp_bsp.isApprox(field.h0PosPos_F2D_bsp_bsp.transpose(), 1.0e-12));
        assert(field.h0NegNeg_F2D_bsp_bsp.isApprox(field.h0NegNeg_F2D_bsp_bsp.transpose(), 1.0e-12));
        assert(field.GammaPosPos_F2D_bsp_bsp.isApprox(field.GammaPosPos_F2D_bsp_bsp.transpose(), 1.0e-12));
        assert(field.GammaNegNeg_F2D_bsp_bsp.isApprox(field.GammaNegNeg_F2D_bsp_bsp.transpose(), 1.0e-12));

        // H±₁₁ = h₀±± + Γ±± - λI; Γ±± includes LN.
        HPos_F2D_2bsp_2bsp.topLeftCorner(Nbsp_I, Nbsp_I) = field.h0PosPos_F2D_bsp_bsp + field.GammaPosPos_F2D_bsp_bsp;
        HNeg_F2D_2bsp_2bsp.topLeftCorner(Nbsp_I, Nbsp_I) = field.h0NegNeg_F2D_bsp_bsp + field.GammaNegNeg_F2D_bsp_bsp;
        HPos_F2D_2bsp_2bsp.topLeftCorner(Nbsp_I, Nbsp_I).diagonal().array() -= lambda_F;
        HNeg_F2D_2bsp_2bsp.topLeftCorner(Nbsp_I, Nbsp_I).diagonal().array() -= lambda_F;

        // H⁺ = [h⁺⁺-λI, Δ; Δᵀ, -h⁻⁻+λI].
        HPos_F2D_2bsp_2bsp.topRightCorner(Nbsp_I, Nbsp_I) = field.DeltaPosNeg_F2D_bsp_bsp;
        HPos_F2D_2bsp_2bsp.bottomLeftCorner(Nbsp_I, Nbsp_I) = field.DeltaPosNeg_F2D_bsp_bsp.transpose();
        HPos_F2D_2bsp_2bsp.bottomRightCorner(Nbsp_I, Nbsp_I) = -HNeg_F2D_2bsp_2bsp.topLeftCorner(Nbsp_I, Nbsp_I);

        // H⁻ = [h⁻⁻-λI, -Δᵀ; -Δ, -h⁺⁺+λI].
        HNeg_F2D_2bsp_2bsp.topRightCorner(Nbsp_I, Nbsp_I) = -field.DeltaPosNeg_F2D_bsp_bsp.transpose();
        HNeg_F2D_2bsp_2bsp.bottomLeftCorner(Nbsp_I, Nbsp_I) = -field.DeltaPosNeg_F2D_bsp_bsp;
        HNeg_F2D_2bsp_2bsp.bottomRightCorner(Nbsp_I, Nbsp_I) = -HPos_F2D_2bsp_2bsp.topLeftCorner(Nbsp_I, Nbsp_I);

        // H± → (E±, X±).
        HPos_eigensolver.compute(HPos_F2D_2bsp_2bsp);
        HNeg_eigensolver.compute(HNeg_F2D_2bsp_2bsp);
        assert(HPos_eigensolver.info() == Eigen::Success);
        assert(HNeg_eigensolver.info() == Eigen::Success);
        const Eigen::VectorXd& eigenvaluesPos_F1D_state = HPos_eigensolver.eigenvalues();
        const Eigen::VectorXd& eigenvaluesNeg_F1D_state = HNeg_eigensolver.eigenvalues();
        const Eigen::MatrixXd& eigenvectorsPos_F2D_2bsp_state = HPos_eigensolver.eigenvectors();
        const Eigen::MatrixXd& eigenvectorsNeg_F2D_2bsp_state = HNeg_eigensolver.eigenvectors();
        assert((eigenvaluesPos_F1D_state.array() != 0.0).all());
        assert((eigenvaluesNeg_F1D_state.array() != 0.0).all());

        // E± > 0; Nqp⁺ + Nqp⁻ = 2Nbsp.
        solution.EqpPos_F1D_bqp.resize((eigenvaluesPos_F1D_state.array() > 0.0).count());
        solution.EqpNeg_F1D_bqp.resize((eigenvaluesNeg_F1D_state.array() > 0.0).count());
        assert(solution.EqpPos_F1D_bqp.size() + solution.EqpNeg_F1D_bqp.size() == 2 * Nbsp_I);
        solution.EqpPos_F1D_bqp = eigenvaluesPos_F1D_state.tail(solution.EqpPos_F1D_bqp.size());
        solution.EqpNeg_F1D_bqp = eigenvaluesNeg_F1D_state.tail(solution.EqpNeg_F1D_bqp.size());

        // X⁺ = [U⁺; V⁻]; X⁻ = [U⁻; V⁺].
        solution.UPos_F2D_bsp_bqp = eigenvectorsPos_F2D_2bsp_state.topRightCorner(Nbsp_I, solution.EqpPos_F1D_bqp.size());
        solution.VNeg_F2D_bsp_bqp = eigenvectorsPos_F2D_2bsp_state.bottomRightCorner(Nbsp_I, solution.EqpPos_F1D_bqp.size());
        solution.UNeg_F2D_bsp_bqp = eigenvectorsNeg_F2D_2bsp_state.topRightCorner(Nbsp_I, solution.EqpNeg_F1D_bqp.size());
        solution.VPos_F2D_bsp_bqp = eigenvectorsNeg_F2D_2bsp_state.bottomRightCorner(Nbsp_I, solution.EqpNeg_F1D_bqp.size());

        // ε± = λ + E±(1-2‖V∓‖²); ε± > Ecut + tail → factors = 0.
        const double EspCutTail_F = std::log(1.0 / 1.0e-6 - 1.0) / 100.0;
        solution.fPos_F1D_bqp.resize(solution.EqpPos_F1D_bqp.size());
        Eigen::VectorXd factorUPos_F1D_bqp{};
        factorUPos_F1D_bqp.resize(solution.EqpPos_F1D_bqp.size());
        factorUPos_F1D_bqp.setZero();
        Eigen::VectorXd factorVPos_F1D_bqp{};
        factorVPos_F1D_bqp.resize(solution.EqpPos_F1D_bqp.size());
        factorVPos_F1D_bqp.setZero();

        // f = (1-tanh(E/(2T)))/2; active factors = (f,1-f).
        for (int bqp_I = 0; bqp_I < solution.EqpPos_F1D_bqp.size(); ++bqp_I) {
            const double Eqp_F = solution.EqpPos_F1D_bqp(bqp_I);
            solution.fPos_F1D_bqp(bqp_I) = temperature_F > 1.0e-12 ? 0.5 * (1.0 - std::tanh(0.5 * Eqp_F / temperature_F)) : 0.0;
            const double Esp_F = lambda_F + Eqp_F * (1.0 - 2.0 * solution.VNeg_F2D_bsp_bqp.col(bqp_I).squaredNorm());
            if (Esp_F > EspCut_F + EspCutTail_F) {continue;}
            factorUPos_F1D_bqp(bqp_I) = solution.fPos_F1D_bqp(bqp_I);
            factorVPos_F1D_bqp(bqp_I) = 1.0 - solution.fPos_F1D_bqp(bqp_I);
        }

        solution.fNeg_F1D_bqp.resize(solution.EqpNeg_F1D_bqp.size());
        Eigen::VectorXd factorUNeg_F1D_bqp{};
        factorUNeg_F1D_bqp.resize(solution.EqpNeg_F1D_bqp.size());
        factorUNeg_F1D_bqp.setZero();
        Eigen::VectorXd factorVNeg_F1D_bqp{};
        factorVNeg_F1D_bqp.resize(solution.EqpNeg_F1D_bqp.size());
        factorVNeg_F1D_bqp.setZero();

        // f = (1-tanh(E/(2T)))/2; active factors = (f,1-f).
        for (int bqp_I = 0; bqp_I < solution.EqpNeg_F1D_bqp.size(); ++bqp_I) {
            const double Eqp_F = solution.EqpNeg_F1D_bqp(bqp_I);
            solution.fNeg_F1D_bqp(bqp_I) = temperature_F > 1.0e-12 ? 0.5 * (1.0 - std::tanh(0.5 * Eqp_F / temperature_F)) : 0.0;
            const double Esp_F = lambda_F + Eqp_F * (1.0 - 2.0 * solution.VPos_F2D_bsp_bqp.col(bqp_I).squaredNorm());
            if (Esp_F > EspCut_F + EspCutTail_F) {continue;}
            factorUNeg_F1D_bqp(bqp_I) = solution.fNeg_F1D_bqp(bqp_I);
            factorVNeg_F1D_bqp(bqp_I) = 1.0 - solution.fNeg_F1D_bqp(bqp_I);
        }

        // ρ⁺⁺ = V⁺(1-f⁻)(V⁺)ᵀ + U⁺f⁺(U⁺)ᵀ.
        solution.rhoPosPos_F2D_bsp_bsp.noalias() = solution.VPos_F2D_bsp_bqp * factorVNeg_F1D_bqp.asDiagonal() * solution.VPos_F2D_bsp_bqp.transpose();
        solution.rhoPosPos_F2D_bsp_bsp.noalias() += solution.UPos_F2D_bsp_bqp * factorUPos_F1D_bqp.asDiagonal() * solution.UPos_F2D_bsp_bqp.transpose();

        // ρ⁻⁻ = V⁻(1-f⁺)(V⁻)ᵀ + U⁻f⁻(U⁻)ᵀ.
        solution.rhoNegNeg_F2D_bsp_bsp.noalias() = solution.VNeg_F2D_bsp_bqp * factorVPos_F1D_bqp.asDiagonal() * solution.VNeg_F2D_bsp_bqp.transpose();
        solution.rhoNegNeg_F2D_bsp_bsp.noalias() += solution.UNeg_F2D_bsp_bqp * factorUNeg_F1D_bqp.asDiagonal() * solution.UNeg_F2D_bsp_bqp.transpose();

        // κ⁺⁻ = V⁺(1-f⁻)(U⁻)ᵀ + U⁺f⁺(V⁻)ᵀ.
        solution.kappaPosNeg_F2D_bsp_bsp.noalias() = solution.VPos_F2D_bsp_bqp * factorVNeg_F1D_bqp.asDiagonal() * solution.UNeg_F2D_bsp_bqp.transpose();
        solution.kappaPosNeg_F2D_bsp_bsp.noalias() += solution.UPos_F2D_bsp_bqp * factorUPos_F1D_bqp.asDiagonal() * solution.VNeg_F2D_bsp_bqp.transpose();

        // N = Σ_b Tr(ρ_b⁺⁺+ρ_b⁻⁻).
        N_F += solution.rhoPosPos_F2D_bsp_bsp.trace() + solution.rhoNegNeg_F2D_bsp_bsp.trace();
    }
    return N_F;
}

inline void HFBAxial::search_lambda(double lambdaTolerance_F) {
    assert(TargetN_I >= 0 && TargetN_I <= 2 * std::accumulate(Nbsp_I1D_block.begin(), Nbsp_I1D_block.end(), 0));
    assert(std::isfinite(lambda_F));
    assert(std::isfinite(lambdaTolerance_F) && lambdaTolerance_F > 0.0);
    const double Ntolerance_F = std::max(1.0e-8, 1.0e-10 * std::max(1, TargetN_I));
    std::map<double, double> NerrorByLambda_Map{};

    // Fixed fields; trial blocking preserves external trackers.
    const auto calc_N_Func = [&](double lambdaTrial_F, bool updateTracking_B) {
        lambda_F = lambdaTrial_F;
        const double Ncalc_F = update_UV_E_rho_kappa();
        if (blocking_Func) {return blocking_Func(hfb_axial_solutions, updateTracking_B);}
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

inline void HFBAxialNucleus::iterate(bool useCurrentFields_B) {
    assert(std::isfinite(accuracy_F) && accuracy_F > 0.0);
    assert(NiterationsMax_I > 0);
    assert(std::isfinite(mixingMin_F) && std::isfinite(mixingMax_F) && mixingMin_F > 0.0 && mixingMin_F <= mixingMax_F && mixingMax_F <= 1.0);

    // Each block packs 3 full matrices.
    int Npacked_I = 0;
    for (const HFBAxial* hfb_Ptr : {&hfb_axial_neutron, &hfb_axial_proton}) {
        for (int Nbsp_I : hfb_Ptr->Nbsp_I1D_block) {Npacked_I += 3 * Nbsp_I * Nbsp_I;}
    }
    assert(Npacked_I >= 7);

    Eigen::VectorXd x_F1D_packed{};
    Eigen::VectorXd Gx_F1D_packed{};
    x_F1D_packed.resize(Npacked_I);
    Gx_F1D_packed.resize(Npacked_I);

    // pack: species → blocks → full h and Delta matrices.
    const auto pack_h_Delta_Func = [&](Eigen::VectorXd& data_F1D_packed) {
        int packed_I = 0;
        for (const HFBAxial* hfb_Ptr : {&hfb_axial_neutron, &hfb_axial_proton}) {
            for (int block_I = 0; block_I < hfb_Ptr->Nblock_I; ++block_I) {
                const int Nbsp_I = hfb_Ptr->Nbsp_I1D_block[block_I];
                const HFBAxialBlockField& field = hfb_Ptr->hfb_axial_fields[block_I];
                for (int column_I = 0; column_I < Nbsp_I; ++column_I) {
                    for (int row_I = 0; row_I < Nbsp_I; ++row_I) {
                        data_F1D_packed(packed_I++) = field.h0PosPos_F2D_bsp_bsp(row_I, column_I) + field.GammaPosPos_F2D_bsp_bsp(row_I, column_I);
                    }
                }
                for (int column_I = 0; column_I < Nbsp_I; ++column_I) {
                    for (int row_I = 0; row_I < Nbsp_I; ++row_I) {
                        data_F1D_packed(packed_I++) = field.h0NegNeg_F2D_bsp_bsp(row_I, column_I) + field.GammaNegNeg_F2D_bsp_bsp(row_I, column_I);
                    }
                }
                for (int column_I = 0; column_I < Nbsp_I; ++column_I) {
                    for (int row_I = 0; row_I < Nbsp_I; ++row_I) {
                        data_F1D_packed(packed_I++) = field.DeltaPosNeg_F2D_bsp_bsp(row_I, column_I);
                    }
                }
            }
        }
        assert(packed_I == Npacked_I);
    };

    // unpack: species → blocks → full h and Delta matrices.
    const auto unpack_h_Delta_Func = [&](const Eigen::VectorXd& data_F1D_packed) {
        int packed_I = 0;
        for (HFBAxial* hfb_Ptr : {&hfb_axial_neutron, &hfb_axial_proton}) {
            for (int block_I = 0; block_I < hfb_Ptr->Nblock_I; ++block_I) {
                const int Nbsp_I = hfb_Ptr->Nbsp_I1D_block[block_I];
                HFBAxialBlockField& field = hfb_Ptr->hfb_axial_fields[block_I];
                for (int column_I = 0; column_I < Nbsp_I; ++column_I) {
                    for (int row_I = 0; row_I < Nbsp_I; ++row_I) {
                        field.GammaPosPos_F2D_bsp_bsp(row_I, column_I) = data_F1D_packed(packed_I++) - field.h0PosPos_F2D_bsp_bsp(row_I, column_I);
                    }
                }
                for (int column_I = 0; column_I < Nbsp_I; ++column_I) {
                    for (int row_I = 0; row_I < Nbsp_I; ++row_I) {
                        field.GammaNegNeg_F2D_bsp_bsp(row_I, column_I) = data_F1D_packed(packed_I++) - field.h0NegNeg_F2D_bsp_bsp(row_I, column_I);
                    }
                }
                for (int column_I = 0; column_I < Nbsp_I; ++column_I) {
                    for (int row_I = 0; row_I < Nbsp_I; ++row_I) {
                        field.DeltaPosNeg_F2D_bsp_bsp(row_I, column_I) = data_F1D_packed(packed_I++);
                    }
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
        hfb_axial_neutron.search_lambda(lambdaTolerance_F);
        hfb_axial_proton.search_lambda(lambdaTolerance_F);
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
