/**
 * @file    hfb_axial.hpp
 * @author  cailea
 * @date    2026-09-10
 * @brief   Axial HFB without time-reversal symmetry.
 */

#pragma once

#include <cassert>
#include <cmath>
#include <functional>
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
     */
    double update_UV_E_rho_kappa(double lambda_F, double temperature_F);
};

class HFBAxialNucleus {
public:
    using OneBodyFunc = std::function<void(int block_I, Eigen::MatrixXd& h0PosPos_F2D_bsp_bsp, Eigen::MatrixXd& h0NegNeg_F2D_bsp_bsp)>;

    using TwoBodyFunc = std::function<void(const std::vector<HFBAxialBlockSolution>& solutions_n, const std::vector<HFBAxialBlockSolution>& solutions_p, std::vector<HFBAxialBlockField>& fields_n, std::vector<HFBAxialBlockField>& fields_p)>;

    HFBAxial hfb_axial_neutron;
    HFBAxial hfb_axial_proton;

    OneBodyFunc build_onebody_neutron{};
    OneBodyFunc build_onebody_proton{};
    TwoBodyFunc build_twobody{};

public:
    /**
     * @brief  Initialize species and build their one-body matrices.
     * @math   h_q,b = h₀,q,b + Γ_q,b; q ∈ {n,p}.
     * @output Initialized species and stored field callbacks.
     */
    HFBAxialNucleus(const std::vector<int>& NbspN_I1D_block_, const std::vector<int>& NbspP_I1D_block_, const OneBodyFunc& build_onebodyN_, const OneBodyFunc& build_onebodyP_, const TwoBodyFunc& build_twobody_)
    : hfb_axial_neutron(NbspN_I1D_block_), hfb_axial_proton(NbspP_I1D_block_) {
        assert(build_onebodyN_);
        assert(build_onebodyP_);
        assert(build_twobody_);
        build_onebody_neutron = build_onebodyN_;
        build_onebody_proton = build_onebodyP_;
        build_twobody = build_twobody_;

        // block_n → (h₀,n⁺⁺, h₀,n⁻⁻).
        for (int block_I = 0; block_I < hfb_axial_neutron.Nblock_I; ++block_I) {
            build_onebody_neutron(block_I, hfb_axial_neutron.hfb_axial_fields[block_I].h0PosPos_F2D_bsp_bsp, hfb_axial_neutron.hfb_axial_fields[block_I].h0NegNeg_F2D_bsp_bsp);
        }

        // block_p → (h₀,p⁺⁺, h₀,p⁻⁻).
        for (int block_I = 0; block_I < hfb_axial_proton.Nblock_I; ++block_I) {
            build_onebody_proton(block_I, hfb_axial_proton.hfb_axial_fields[block_I].h0PosPos_F2D_bsp_bsp, hfb_axial_proton.hfb_axial_fields[block_I].h0NegNeg_F2D_bsp_bsp);
        }
    }

    /**
     * @brief  Update both species using all block densities.
     * @math   {ρ_n,b,κ_n,b,ρ_p,b,κ_p,b} → {Γ_n,b,Δ_n,b,Γ_p,b,Δ_p,b}.
     * @output Overwritten neutron and proton Gamma and Delta.
     */
    void update_Gamma_Delta();
};

/**
 * @brief  Solve thermal branches using real symmetric eigendecomposition.
 * @math   H±X± = X±E±; f± = 1/(1+exp(E±/T)).
 * @output Updated quasiparticle solutions and densities.
 * @output Mean particle number N = Σ_b Tr(ρ_b⁺⁺+ρ_b⁻⁻).
 */
inline double HFBAxial::update_UV_E_rho_kappa(double lambda_F, double temperature_F) {
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

        // H⁺₁₁ = h⁺⁺ - λI; H⁻₁₁ = h⁻⁻ - λI.
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

        // T ≤ 10⁻¹² → f± = 0.
        solution.fPos_F1D_bqp.resize(solution.EqpPos_F1D_bqp.size());
        solution.fNeg_F1D_bqp.resize(solution.EqpNeg_F1D_bqp.size());
        solution.fPos_F1D_bqp.setZero();
        solution.fNeg_F1D_bqp.setZero();
        if (temperature_F > 1.0e-12) {
            // f± = e⁻ᴱ±ᐟᵀ/(1+e⁻ᴱ±ᐟᵀ).
            for (int bqp_I = 0; bqp_I < solution.EqpPos_F1D_bqp.size(); ++bqp_I) {
                const double expMinusEOverT_F = std::exp(-solution.EqpPos_F1D_bqp(bqp_I) / temperature_F);
                solution.fPos_F1D_bqp(bqp_I) = expMinusEOverT_F / (1.0 + expMinusEOverT_F);
            }
            for (int bqp_I = 0; bqp_I < solution.EqpNeg_F1D_bqp.size(); ++bqp_I) {
                const double expMinusEOverT_F = std::exp(-solution.EqpNeg_F1D_bqp(bqp_I) / temperature_F);
                solution.fNeg_F1D_bqp(bqp_I) = expMinusEOverT_F / (1.0 + expMinusEOverT_F);
            }
        }

        // ρ⁺⁺ = V⁺(1-f⁻)(V⁺)ᵀ + U⁺f⁺(U⁺)ᵀ.
        solution.rhoPosPos_F2D_bsp_bsp.noalias() = solution.VPos_F2D_bsp_bqp * (1.0 - solution.fNeg_F1D_bqp.array()).matrix().asDiagonal() * solution.VPos_F2D_bsp_bqp.transpose();
        solution.rhoPosPos_F2D_bsp_bsp.noalias() += solution.UPos_F2D_bsp_bqp * solution.fPos_F1D_bqp.asDiagonal() * solution.UPos_F2D_bsp_bqp.transpose();

        // ρ⁻⁻ = V⁻(1-f⁺)(V⁻)ᵀ + U⁻f⁻(U⁻)ᵀ.
        solution.rhoNegNeg_F2D_bsp_bsp.noalias() = solution.VNeg_F2D_bsp_bqp * (1.0 - solution.fPos_F1D_bqp.array()).matrix().asDiagonal() * solution.VNeg_F2D_bsp_bqp.transpose();
        solution.rhoNegNeg_F2D_bsp_bsp.noalias() += solution.UNeg_F2D_bsp_bqp * solution.fNeg_F1D_bqp.asDiagonal() * solution.UNeg_F2D_bsp_bqp.transpose();

        // κ⁺⁻ = V⁺(1-f⁻)(U⁻)ᵀ + U⁺f⁺(V⁻)ᵀ.
        solution.kappaPosNeg_F2D_bsp_bsp.noalias() = solution.VPos_F2D_bsp_bqp * (1.0 - solution.fNeg_F1D_bqp.array()).matrix().asDiagonal() * solution.UNeg_F2D_bsp_bqp.transpose();
        solution.kappaPosNeg_F2D_bsp_bsp.noalias() += solution.UPos_F2D_bsp_bqp * solution.fPos_F1D_bqp.asDiagonal() * solution.VNeg_F2D_bsp_bqp.transpose();

        // N = Σ_b Tr(ρ_b⁺⁺+ρ_b⁻⁻).
        N_F += solution.rhoPosPos_F2D_bsp_bsp.trace() + solution.rhoNegNeg_F2D_bsp_bsp.trace();
    }
    return N_F;
}

/**
 * @brief  Update both species using all block densities.
 * @math   {ρ_n,b,κ_n,b,ρ_p,b,κ_p,b} → {Γ_n,b,Δ_n,b,Γ_p,b,Δ_p,b}.
 * @output Overwritten neutron and proton Gamma and Delta.
 */
inline void HFBAxialNucleus::update_Gamma_Delta() {
    build_twobody(hfb_axial_neutron.hfb_axial_solutions, hfb_axial_proton.hfb_axial_solutions, hfb_axial_neutron.hfb_axial_fields, hfb_axial_proton.hfb_axial_fields);
}
