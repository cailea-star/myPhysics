/**
 * @file    hfb_axial.hpp
 * @author  cailea
 * @date    2026-09-10
 * @brief   Axial HFB without time-reversal symmetry.
 */

#pragma once

#include <cassert>
#include <cmath>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Eigenvalues>

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
    double lambda_F = -7.0; // Fermi energy [MeV].
    double lambda2_F = 0.0; // Lipkin-Nogami λ₂ [MeV].
    double temperature_F = 0.0;
    double EspCut_F = 60.0; // Equivalent single-particle energy cutoff [MeV].

    int Nblock_I = 0;
    std::vector<int> Nbsp_I1D_block{};

    std::vector<HFBAxialBlockSolution> hfb_axial_solutions{};
    std::vector<HFBAxialBlockField> hfb_axial_fields{};

    std::vector<Eigen::MatrixXd> HPos_F3D_block_2bsp_2bsp{};
    std::vector<Eigen::MatrixXd> HNeg_F3D_block_2bsp_2bsp{};
    std::vector<Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>> HPos_eigensolvers{};
    std::vector<Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>> HNeg_eigensolvers{};

public:
    /**
     * @brief  Allocate single-species axial HFB blocks.
     * @math   dim(bsp_b⁺) = dim(bsp_b⁻) = Nbsp_b.
     * @output Initialized fields, solutions, and eigensolver workspaces.
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
            field.h0PosPos_F2D_bsp_bsp.setZero();
            field.h0NegNeg_F2D_bsp_bsp.resize(Nbsp_I, Nbsp_I);
            field.h0NegNeg_F2D_bsp_bsp.setZero();
            field.GammaPosPos_F2D_bsp_bsp.resize(Nbsp_I, Nbsp_I);
            field.GammaPosPos_F2D_bsp_bsp.setZero();
            field.GammaNegNeg_F2D_bsp_bsp.resize(Nbsp_I, Nbsp_I);
            field.GammaNegNeg_F2D_bsp_bsp.setZero();
            field.DeltaPosNeg_F2D_bsp_bsp.resize(Nbsp_I, Nbsp_I);
            field.DeltaPosNeg_F2D_bsp_bsp.setZero();

            // ρ, κ ∈ ℝ^{Nbsp×Nbsp}.
            solution.rhoPosPos_F2D_bsp_bsp.resize(Nbsp_I, Nbsp_I);
            solution.rhoPosPos_F2D_bsp_bsp.setZero();
            solution.rhoNegNeg_F2D_bsp_bsp.resize(Nbsp_I, Nbsp_I);
            solution.rhoNegNeg_F2D_bsp_bsp.setZero();
            solution.kappaPosNeg_F2D_bsp_bsp.resize(Nbsp_I, Nbsp_I);
            solution.kappaPosNeg_F2D_bsp_bsp.setZero();

            // U±, V± ∈ ℝ^{Nbsp×Nbsp} initially.
            solution.UPos_F2D_bsp_bqp.resize(Nbsp_I, Nbsp_I);
            solution.VNeg_F2D_bsp_bqp.resize(Nbsp_I, Nbsp_I);
            solution.UNeg_F2D_bsp_bqp.resize(Nbsp_I, Nbsp_I);
            solution.VPos_F2D_bsp_bqp.resize(Nbsp_I, Nbsp_I);
            solution.UPos_F2D_bsp_bqp.setZero();
            solution.VNeg_F2D_bsp_bqp.setZero();
            solution.UNeg_F2D_bsp_bqp.setZero();
            solution.VPos_F2D_bsp_bqp.setZero();

            // H± ∈ ℝ^{2Nbsp×2Nbsp}.
            HPos_F3D_block_2bsp_2bsp[block_I].resize(2 * Nbsp_I, 2 * Nbsp_I);
            HNeg_F3D_block_2bsp_2bsp[block_I].resize(2 * Nbsp_I, 2 * Nbsp_I);
            HPos_F3D_block_2bsp_2bsp[block_I].setZero();
            HNeg_F3D_block_2bsp_2bsp[block_I].setZero();
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
     * @note   Real fields; no zero modes; T = k_B T_phys ≥ 0.
     * @note   ε = λ + E(1-2‖V‖²); cutoff affects ρ and κ.
     * @note   EspCut_F = ∞ disables the cutoff.
     * @note   h_LN = h₀ + Γ + 4λ₂ρ - 2λ₂I; pre-solve ρ.
     */
    double update_UV_E_rho_kappa();
};

class HFBAxialNucleus {
public:
    HFBAxial hfb_axial_neutron;
    HFBAxial hfb_axial_proton;

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
     * @math   h₀,n, h₀,p.
     * @output Initialized neutron and proton one-body fields.
     */
    virtual void initialize_h0() = 0;

    /**
     * @brief  Initialize HFB fields in the derived model.
     * @math   (N,Z) → (Γ_n,Δ_n,Γ_p,Δ_p)_initial.
     * @output Initialized Gamma and Delta for both species.
     */
    virtual void initialize_GammaDelta(int TargetN_I, int TargetZ_I) = 0;

    /**
     * @brief  Update both species using all block densities.
     * @math   {ρ_n,b,κ_n,b,ρ_p,b,κ_p,b} → {Γ_n,b,Δ_n,b,Γ_p,b,Δ_p,b}.
     * @output Overwritten neutron and proton Gamma and Delta.
     */
    virtual void update_Gamma_Delta() = 0;

    /**
     * @brief  Iterate the unblocked HFB equations.
     * @math   (N,Z) → HFB_converged.
     * @output Updated converged solver state.
     */
    void iterate(int TargetN_I, int TargetZ_I);
};

/**
 * @brief  Solve thermal branches using real symmetric eigendecomposition.
 * @math   H±X± = X±E±; f± = 1/(1+exp(E±/T)).
 * @output Updated quasiparticle solutions and densities.
 * @output Mean particle number N = Σ_b Tr(ρ_b⁺⁺+ρ_b⁻⁻).
 */
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

        // H±₁₁ = h₀±± + Γ±± + 4λ₂ρ±± - (λ+2λ₂)I.
        HPos_F2D_2bsp_2bsp.topLeftCorner(Nbsp_I, Nbsp_I) = field.h0PosPos_F2D_bsp_bsp + field.GammaPosPos_F2D_bsp_bsp + 4.0 * lambda2_F * solution.rhoPosPos_F2D_bsp_bsp;
        HNeg_F2D_2bsp_2bsp.topLeftCorner(Nbsp_I, Nbsp_I) = field.h0NegNeg_F2D_bsp_bsp + field.GammaNegNeg_F2D_bsp_bsp + 4.0 * lambda2_F * solution.rhoNegNeg_F2D_bsp_bsp;
        HPos_F2D_2bsp_2bsp.topLeftCorner(Nbsp_I, Nbsp_I).diagonal().array() -= lambda_F + 2.0 * lambda2_F;
        HNeg_F2D_2bsp_2bsp.topLeftCorner(Nbsp_I, Nbsp_I).diagonal().array() -= lambda_F + 2.0 * lambda2_F;

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
