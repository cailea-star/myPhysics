/**
 * @file    hfb_pfaffian.hpp
 * @author  cailea
 * @date    2026-09-10
 * @brief   HFB contractions and multiquasiparticle Pfaffian kernels.
 */

#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <complex>
#include <vector>

#include <Eigen/Core>
#include <Eigen/LU>
#include <unsupported/Eigen/CXX11/Tensor>

using doubleC = std::complex<double>;

/**
 * @brief  Evaluate HFB configuration kernels using Pfaffians.
 * @note   Updated contractions; overlap ≠ 0; unconjugated x; references alias outputs.
 */
class HFBPfaffian {
public:
    int Nsp_I = 0;


    Eigen::MatrixXcd A_C2D_qp1_qp2{};

    doubleC overlap_C = {0.0, 0.0};

    Eigen::MatrixXcd SpDagSp_C2D_sp_sp{};
    Eigen::MatrixXcd SpSpDag_C2D_sp_sp{};
    Eigen::MatrixXcd SpDagSpDag_C2D_sp_sp{};
    Eigen::MatrixXcd SpSp_C2D_sp_sp{};

    Eigen::MatrixXcd Qp1Qp1_C2D_qp1_qp1{};
    Eigen::MatrixXcd Qp2DagQp2Dag_C2D_qp2_qp2{};
    Eigen::MatrixXcd Qp1Qp2Dag_C2D_qp1_qp2{};

    Eigen::MatrixXcd Qp1Sp_C2D_qp1_sp{};
    Eigen::MatrixXcd Qp1SpDag_C2D_qp1_sp{};
    Eigen::MatrixXcd SpQp2Dag_C2D_sp_qp2{};
    Eigen::MatrixXcd SpDagQp2Dag_C2D_sp_qp2{};

    std::vector<std::vector<int>> config1_I2D_cfg1_cqp1{};
    std::vector<std::vector<int>> config2_I2D_cfg2_cqp2{};

    // Signed Pfaffian minors; excludes vacuum overlap.
    std::vector<std::vector<doubleC>> F0_C2D_cfg1_cfg2{};
    std::vector<std::vector<Eigen::VectorXcd>> F1_C3D_cfg1_cfg2_cqp12{};
    std::vector<std::vector<Eigen::MatrixXcd>> F2_C4D_cfg1_cfg2_cqp12_cqp12{};
    // F⁴: fully antisymmetric; repeated indices vanish.
    std::vector<std::vector<Eigen::Tensor<doubleC, 4, Eigen::ColMajor>>> F4_C6D_cfg1_cfg2_cqp12_cqp12_cqp12_cqp12{};

    Eigen::MatrixXcd Overlap_C2D_cfg1_cfg2{};
    Eigen::MatrixXcd OneFermion_C2D_cfg1_cfg2{};
    Eigen::MatrixXcd OneBody_C2D_cfg1_cfg2{};
    Eigen::MatrixXcd TwoBody_C2D_cfg1_cfg2{};

private:
    Eigen::FullPivLU<Eigen::MatrixXcd> U1_lu{};
    Eigen::FullPivLU<Eigen::MatrixXcd> U2_lu{};
    Eigen::FullPivLU<Eigen::MatrixXcd> A_lu{};
    Eigen::VectorXcd Sworkspace_C1D_element{};

public:
    /**
     * @brief  Construct empty Pfaffian workspaces.
     * @math   Nsp = 0; C₁ = C₂ = ∅.
     * @output Empty workspaces.
     */
    HFBPfaffian() = default;

    /**
     * @brief  Store explicit configurations and allocate HFB workspaces.
     * @math   U₁,V₁,U₂,V₂ ∈ ℂ^{Nsp×Nsp}.
     * @output Stored both configuration tables; allocated matrices and workspaces.
     * @note   0 ≤ q₁ < ⋯ < qᵣ < Nsp; ∅ = vacuum.
     */
    HFBPfaffian(int Nsp_I_, const std::vector<std::vector<int>>& config1_I2D_cfg1_cqp1_, const std::vector<std::vector<int>>& config2_I2D_cfg2_cqp2_) {
        assert(Nsp_I_ > 0);

        Nsp_I = Nsp_I_;
        config1_I2D_cfg1_cqp1 = config1_I2D_cfg1_cqp1_;
        config2_I2D_cfg2_cqp2 = config2_I2D_cfg2_cqp2_;

        // Z₁ = (V₁U₁⁻¹)*; Z₂ = (V₂U₂⁻¹)*.
        U1_lu = Eigen::FullPivLU<Eigen::MatrixXcd>(Nsp_I, Nsp_I);
        U2_lu = Eigen::FullPivLU<Eigen::MatrixXcd>(Nsp_I, Nsp_I);

        // A ∈ ℂ^{Nsp×Nsp}.
        A_C2D_qp1_qp2.resize(Nsp_I, Nsp_I);
        A_lu = Eigen::FullPivLU<Eigen::MatrixXcd>(Nsp_I, Nsp_I);

        // ⟨c†c⟩,⟨cc†⟩,⟨c†c†⟩,⟨cc⟩ ∈ ℂ^{Nsp×Nsp}.
        SpDagSp_C2D_sp_sp.resize(Nsp_I, Nsp_I);
        SpSpDag_C2D_sp_sp.resize(Nsp_I, Nsp_I);
        SpDagSpDag_C2D_sp_sp.resize(Nsp_I, Nsp_I);
        SpSp_C2D_sp_sp.resize(Nsp_I, Nsp_I);

        // ⟨β₁β₁⟩,⟨β₂†β₂†⟩,⟨β₁β₂†⟩ ∈ ℂ^{Nsp×Nsp}.
        Qp1Qp1_C2D_qp1_qp1.resize(Nsp_I, Nsp_I);
        Qp2DagQp2Dag_C2D_qp2_qp2.resize(Nsp_I, Nsp_I);
        Qp1Qp2Dag_C2D_qp1_qp2.resize(Nsp_I, Nsp_I);

        // ⟨β₁c⟩,⟨β₁c†⟩,⟨cβ₂†⟩,⟨c†β₂†⟩ ∈ ℂ^{Nsp×Nsp}.
        Qp1Sp_C2D_qp1_sp.resize(Nsp_I, Nsp_I);
        Qp1SpDag_C2D_qp1_sp.resize(Nsp_I, Nsp_I);
        SpQp2Dag_C2D_sp_qp2.resize(Nsp_I, Nsp_I);
        SpDagQp2Dag_C2D_sp_qp2.resize(Nsp_I, Nsp_I);

        // κ₁ = (μ₁,…,μ_{r₁}); 0 ≤ μ₁ < ⋯ < Nsp.
        for (const auto& config_I1D_cqp : config1_I2D_cfg1_cqp1) {
            for (int cqp_I = 0; cqp_I < static_cast<int>(config_I1D_cqp.size()); ++cqp_I) {
                assert(config_I1D_cqp[cqp_I] >= 0 && config_I1D_cqp[cqp_I] < Nsp_I);
                assert(cqp_I == 0 || config_I1D_cqp[cqp_I - 1] < config_I1D_cqp[cqp_I]);
            }
        }

        // κ₂ = (ν₁,…,ν_{r₂}); 0 ≤ ν₁ < ⋯ < Nsp.
        for (const auto& config_I1D_cqp : config2_I2D_cfg2_cqp2) {
            for (int cqp_I = 0; cqp_I < static_cast<int>(config_I1D_cqp.size()); ++cqp_I) {
                assert(config_I1D_cqp[cqp_I] >= 0 && config_I1D_cqp[cqp_I] < Nsp_I);
                assert(cqp_I == 0 || config_I1D_cqp[cqp_I - 1] < config_I1D_cqp[cqp_I]);
            }
        }

        // N, OBTD, TBTD ∈ ℂ^{Ncfg1×Ncfg2}.
        Overlap_C2D_cfg1_cfg2.resize(config1_I2D_cfg1_cqp1.size(), config2_I2D_cfg2_cqp2.size());
        OneFermion_C2D_cfg1_cfg2.resize(config1_I2D_cfg1_cqp1.size(), config2_I2D_cfg2_cqp2.size());
        OneBody_C2D_cfg1_cfg2.resize(config1_I2D_cfg1_cqp1.size(), config2_I2D_cfg2_cqp2.size());
        TwoBody_C2D_cfg1_cfg2.resize(config1_I2D_cfg1_cqp1.size(), config2_I2D_cfg2_cqp2.size());

        // Lmax = 2Nsp+4; size(Sworkspace) = Lmax².
        const int Lmax_I = 2 * Nsp_I + 4;
        Sworkspace_C1D_element.resize(Lmax_I * Lmax_I);
    }

    /**
     * @brief  Update overlap and contractions using Pfaffians and full-pivot LU.
     * @math   A = U₁ᵀU₂* + V₁ᵀV₂*; Qp1Qp2Dag = A⁻ᵀ.
     * @math   νₐ = √|det Uₐ|; overlap ← exp[i(θ₂−θ₁)] overlap.
     * @output Updated overlap, A, and eleven contractions.
     * @note   Canonical U,V; numerically invertible U₁,U₂,A.
     */
    void update_contractions(const Eigen::MatrixXcd& U1_C2D_sp_qp1_, const Eigen::MatrixXcd& V1_C2D_sp_qp1_, const Eigen::MatrixXcd& U2_C2D_sp_qp2_, const Eigen::MatrixXcd& V2_C2D_sp_qp2_, doubleC phase_C);

    /**
     * @brief Update signed Pfaffian minors for configuration pairs.
     * @math F⁰ = pf(S); F¹,F²,F⁴ = signed deletion minors.
     * @output Updated F0, F1, F2, F4 configuration caches.
     * @note Requires current contractions; excludes vacuum overlap.
     */
    void build_pfaffian_minors();

    /**
     * @brief  Calculate configuration overlaps using cached Pfaffian minors.
     * @math   N_{κ₁κ₂} = ⟨Φ₁;κ₁|Φ₂;κ₂⟩.
     * @output Updated Overlap_C2D_cfg1_cfg2 and its const reference.
     */
    const Eigen::MatrixXcd& calc_overlap();

    /**
     * @brief  Calculate single-fermion kernels using cached Pfaffian minors.
     * @math   X = Σα xα c†α (true), Σα xα cα (false).
     * @output Updated OneFermion_C2D_cfg1_cfg2 and its const reference.
     * @note   Shared output overwritten per call; coefficients remain unconjugated.
     */
    const Eigen::MatrixXcd& calc_one_fermion(const Eigen::VectorXcd& x_C1D_sp, bool creation_B);

    /**
     * @brief  Calculate one-body kernels using cached Pfaffian minors.
     * @math   O = Σαβ O_αβ cα†cβ; K = F⁰T⁰ + Σi<j F²ij T²ij.
     * @output Updated OneBody_C2D_cfg1_cfg2 and its const reference.
     */
    const Eigen::MatrixXcd& calc_one_body(const Eigen::MatrixXcd& OneBody_C2D_sp_sp);

    /**
     * @brief  Calculate two-body kernels using cached Pfaffian minors.
     * @math   O₂ = Σαβγδ O_αβγδ cα†cβ†cδcγ.
     * @output Updated TwoBody_C2D_cfg1_cfg2 and its const reference.
     * @note   Input coefficients include all numerical prefactors.
     */
    const Eigen::MatrixXcd& calc_two_body(const Eigen::Tensor<doubleC, 4, Eigen::ColMajor>& TwoBody_C4D_sp_sp_sp_sp);

    /**
     * @brief  Calculate Pfaffians using pivoted skew-symmetric elimination.
     * @math   pf(S)² = det(S); pf(∅) = 1.
     * @output Pfaffian value; input matrix overwritten.
     * @note   Finite S = −Sᵀ ∈ ℂ²ⁿˣ²ⁿ; zero pivot → 0; no cutoff.
     */
    static doubleC calc_pfaffian(Eigen::Ref<Eigen::MatrixXcd> S_C2D_chain_chain);

};


inline void HFBPfaffian::update_contractions(const Eigen::MatrixXcd& U1_C2D_sp_qp1_, const Eigen::MatrixXcd& V1_C2D_sp_qp1_, const Eigen::MatrixXcd& U2_C2D_sp_qp2_, const Eigen::MatrixXcd& V2_C2D_sp_qp2_, doubleC phase_C) {
    // U₁,V₁,U₂,V₂ ∈ ℂ^{Nsp×Nsp}.
    assert(U1_C2D_sp_qp1_.rows() == Nsp_I && U1_C2D_sp_qp1_.cols() == Nsp_I);
    assert(V1_C2D_sp_qp1_.rows() == Nsp_I && V1_C2D_sp_qp1_.cols() == Nsp_I);
    assert(U2_C2D_sp_qp2_.rows() == Nsp_I && U2_C2D_sp_qp2_.cols() == Nsp_I);
    assert(V2_C2D_sp_qp2_.rows() == Nsp_I && V2_C2D_sp_qp2_.cols() == Nsp_I);

    // UᵀZ† = Vᵀ; ν = √|det U| > 0.
    assert(U1_C2D_sp_qp1_.allFinite() && V1_C2D_sp_qp1_.allFinite());
    assert(U2_C2D_sp_qp2_.allFinite() && V2_C2D_sp_qp2_.allFinite());
    U1_lu.compute(U1_C2D_sp_qp1_);
    U2_lu.compute(U2_C2D_sp_qp2_);
    assert(U1_lu.isInvertible() && U2_lu.isInvertible());

    // ν₁ = √|det U₁|; ν₂ = √|det U₂|.
    const double nu1_F = std::sqrt(std::abs(U1_lu.determinant()));
    const double nu2_F = std::sqrt(std::abs(U2_lu.determinant()));

    Eigen::Map<Eigen::MatrixXcd> X_C2D_2sp_2sp(Sworkspace_C1D_element.data(), 2 * Nsp_I, 2 * Nsp_I);

    // X₁₁ = Z₂ = (V₂U₂⁻¹)*.
    X_C2D_2sp_2sp.topLeftCorner(Nsp_I, Nsp_I) = U2_lu.transpose().solve(V2_C2D_sp_qp2_.transpose());
    X_C2D_2sp_2sp.topLeftCorner(Nsp_I, Nsp_I).adjointInPlace();

    // X₂₂ = −Z₁* = −V₁U₁⁻¹.
    X_C2D_2sp_2sp.bottomRightCorner(Nsp_I, Nsp_I) = U1_lu.transpose().solve(V1_C2D_sp_qp1_.transpose());
    X_C2D_2sp_2sp.bottomRightCorner(Nsp_I, Nsp_I).transposeInPlace();
    X_C2D_2sp_2sp.bottomRightCorner(Nsp_I, Nsp_I) *= -1.0;

    // X₁₂ = −I; X₂₁ = I.
    X_C2D_2sp_2sp.topRightCorner(Nsp_I, Nsp_I).setIdentity();
    X_C2D_2sp_2sp.topRightCorner(Nsp_I, Nsp_I) *= -1.0;
    X_C2D_2sp_2sp.bottomLeftCorner(Nsp_I, Nsp_I).setIdentity();

    // overlap = phase × ν₁ν₂s_Nsp pf(X).
    const int sign_I = (Nsp_I % 4 == 0 || Nsp_I % 4 == 3) ? 1 : -1;
    overlap_C = phase_C * nu1_F * nu2_F * static_cast<double>(sign_I) * calc_pfaffian(X_C2D_2sp_2sp);

    // A = U₁ᵀU₂* + V₁ᵀV₂*.
    A_C2D_qp1_qp2.noalias() = U1_C2D_sp_qp1_.transpose() * U2_C2D_sp_qp2_.conjugate();
    A_C2D_qp1_qp2.noalias() += V1_C2D_sp_qp1_.transpose() * V2_C2D_sp_qp2_.conjugate();
    assert(A_C2D_qp1_qp2.allFinite());
    A_lu.compute(A_C2D_qp1_qp2);
    assert(A_lu.isInvertible());
    Qp1Qp2Dag_C2D_qp1_qp2 = A_lu.inverse().transpose();

    // ⟨β₁c⟩ = A⁻ᵀV₂†; ⟨β₁c†⟩ = A⁻ᵀU₂†.
    Qp1Sp_C2D_qp1_sp.noalias() = Qp1Qp2Dag_C2D_qp1_qp2 * V2_C2D_sp_qp2_.adjoint();
    Qp1SpDag_C2D_qp1_sp.noalias() = Qp1Qp2Dag_C2D_qp1_qp2 * U2_C2D_sp_qp2_.adjoint();

    // ⟨cβ₂†⟩ = U₁A⁻ᵀ; ⟨c†β₂†⟩ = V₁A⁻ᵀ.
    SpQp2Dag_C2D_sp_qp2.noalias() = U1_C2D_sp_qp1_ * Qp1Qp2Dag_C2D_qp1_qp2;
    SpDagQp2Dag_C2D_sp_qp2.noalias() = V1_C2D_sp_qp1_ * Qp1Qp2Dag_C2D_qp1_qp2;

    // ⟨c†c⟩ = V₁⟨β₁c⟩; ⟨cc†⟩ = U₁⟨β₁c†⟩.
    SpDagSp_C2D_sp_sp.noalias() = V1_C2D_sp_qp1_ * Qp1Sp_C2D_qp1_sp;
    SpSpDag_C2D_sp_sp.noalias() = U1_C2D_sp_qp1_ * Qp1SpDag_C2D_qp1_sp;

    // ⟨c†c†⟩ = V₁⟨β₁c†⟩; ⟨cc⟩ = U₁⟨β₁c⟩.
    SpDagSpDag_C2D_sp_sp.noalias() = V1_C2D_sp_qp1_ * Qp1SpDag_C2D_qp1_sp;
    SpSp_C2D_sp_sp.noalias() = U1_C2D_sp_qp1_ * Qp1Sp_C2D_qp1_sp;

    // ⟨β₁β₁⟩ = ⟨β₁c⟩U₁* + ⟨β₁c†⟩V₁*.
    Qp1Qp1_C2D_qp1_qp1.noalias() = Qp1Sp_C2D_qp1_sp * U1_C2D_sp_qp1_.conjugate();
    Qp1Qp1_C2D_qp1_qp1.noalias() += Qp1SpDag_C2D_qp1_sp * V1_C2D_sp_qp1_.conjugate();

    // ⟨β₂†β₂†⟩ = V₂ᵀ⟨cβ₂†⟩ + U₂ᵀ⟨c†β₂†⟩.
    Qp2DagQp2Dag_C2D_qp2_qp2.noalias() = V2_C2D_sp_qp2_.transpose() * SpQp2Dag_C2D_sp_qp2;
    Qp2DagQp2Dag_C2D_qp2_qp2.noalias() += U2_C2D_sp_qp2_.transpose() * SpDagQp2Dag_C2D_sp_qp2;

    build_pfaffian_minors();
}

inline void HFBPfaffian::build_pfaffian_minors() {
    // F⁰,F¹,F²,F⁴: signed Pfaffian minors.
    const int Ncfg1_I = static_cast<int>(config1_I2D_cfg1_cqp1.size());
    const int Ncfg2_I = static_cast<int>(config2_I2D_cfg2_cqp2.size());
    // Columns: permutation of (i,j,k,l), antisymmetric sign.
    static constexpr std::array<std::array<int, 5>, 24> permutationSign_I2D_permutation_entry{{
        {0, 1, 2, 3, +1},
        {0, 1, 3, 2, -1},
        {0, 2, 1, 3, -1},
        {0, 2, 3, 1, +1},
        {0, 3, 1, 2, +1},
        {0, 3, 2, 1, -1},
        {1, 0, 2, 3, -1},
        {1, 0, 3, 2, +1},
        {1, 2, 0, 3, +1},
        {1, 2, 3, 0, -1},
        {1, 3, 0, 2, -1},
        {1, 3, 2, 0, +1},
        {2, 0, 1, 3, +1},
        {2, 0, 3, 1, -1},
        {2, 1, 0, 3, -1},
        {2, 1, 3, 0, +1},
        {2, 3, 0, 1, +1},
        {2, 3, 1, 0, -1},
        {3, 0, 1, 2, -1},
        {3, 0, 2, 1, +1},
        {3, 1, 0, 2, +1},
        {3, 1, 2, 0, -1},
        {3, 2, 0, 1, -1},
        {3, 2, 1, 0, +1},
    }};

    // Ncqp12 ≤ 2Nsp; reuse minor storage across configurations.
    Eigen::VectorXcd MinorWorkspace_C1D_element(4 * Nsp_I * Nsp_I);
    std::vector<int> remaining_cqp12{};
    remaining_cqp12.reserve(2 * Nsp_I);

    F0_C2D_cfg1_cfg2.resize(Ncfg1_I);
    F1_C3D_cfg1_cfg2_cqp12.resize(Ncfg1_I);
    F2_C4D_cfg1_cfg2_cqp12_cqp12.resize(Ncfg1_I);
    F4_C6D_cfg1_cfg2_cqp12_cqp12_cqp12_cqp12.resize(Ncfg1_I);
    for (int cfg1_I = 0; cfg1_I < Ncfg1_I; ++cfg1_I) {
        F0_C2D_cfg1_cfg2[cfg1_I].assign(Ncfg2_I, 0.0);
        F1_C3D_cfg1_cfg2_cqp12[cfg1_I].resize(Ncfg2_I);
        F2_C4D_cfg1_cfg2_cqp12_cqp12[cfg1_I].resize(Ncfg2_I);
        F4_C6D_cfg1_cfg2_cqp12_cqp12_cqp12_cqp12[cfg1_I].resize(Ncfg2_I);
        for (int cfg2_I = 0; cfg2_I < Ncfg2_I; ++cfg2_I) {
            // z = (reversed left β, forward right β†).
            const auto& config1_I1D_cqp1 = config1_I2D_cfg1_cqp1[cfg1_I];
            const auto& config2_I1D_cqp2 = config2_I2D_cfg2_cqp2[cfg2_I];

            const int Ncqp1_I = static_cast<int>(config1_I1D_cqp1.size());
            const int Ncqp2_I = static_cast<int>(config2_I1D_cqp2.size());
            const int Ncqp12_I = Ncqp1_I + Ncqp2_I;

            // Left indices reversed; right indices forward.
            Eigen::Map<Eigen::MatrixXcd> S_C2D_cqp12_cqp12(Sworkspace_C1D_element.data(), Ncqp12_I, Ncqp12_I);
            S_C2D_cqp12_cqp12.topLeftCorner(Ncqp1_I, Ncqp1_I) = Qp1Qp1_C2D_qp1_qp1(config1_I1D_cqp1, config1_I1D_cqp1).reverse();
            S_C2D_cqp12_cqp12.topRightCorner(Ncqp1_I, Ncqp2_I) = Qp1Qp2Dag_C2D_qp1_qp2(config1_I1D_cqp1, config2_I1D_cqp2).colwise().reverse();
            S_C2D_cqp12_cqp12.bottomRightCorner(Ncqp2_I, Ncqp2_I) = Qp2DagQp2Dag_C2D_qp2_qp2(config2_I1D_cqp2, config2_I1D_cqp2);

            // Sᵀ = −S.
            S_C2D_cqp12_cqp12.diagonal().setZero();
            S_C2D_cqp12_cqp12.triangularView<Eigen::StrictlyLower>() = -S_C2D_cqp12_cqp12.transpose();

            // Phase = sum of removed one-based positions.
            const auto calc_signed_pfaffian_minor = [&](int chain1_I = -1, int chain2_I = -1, int chain3_I = -1, int chain4_I = -1) {
                const int Nremove_I = (chain1_I >= 0) + (chain2_I >= 0) + (chain3_I >= 0) + (chain4_I >= 0);
                int phase_I = Nremove_I;
                phase_I += std::max(chain1_I, 0) + std::max(chain2_I, 0) + std::max(chain3_I, 0) + std::max(chain4_I, 0);
                double sign_F = 1.0 - 2.0 * (phase_I % 2);

                remaining_cqp12.clear();
                for (int chain_I = 0; chain_I < Ncqp12_I; ++chain_I) {
                    if (chain_I != chain1_I && chain_I != chain2_I && chain_I != chain3_I && chain_I != chain4_I) { remaining_cqp12.push_back(chain_I); }
                }

                const int Nremaining_I = Ncqp12_I - Nremove_I;
                Eigen::Map<Eigen::MatrixXcd> minor_C2D_remaincqp12_remaincqp12(MinorWorkspace_C1D_element.data(), Nremaining_I, Nremaining_I);
                minor_C2D_remaincqp12_remaincqp12 = S_C2D_cqp12_cqp12(remaining_cqp12, remaining_cqp12);
                return sign_F * HFBPfaffian::calc_pfaffian(minor_C2D_remaincqp12_remaincqp12);
            };

            auto& F1_C1D_cqp12 = F1_C3D_cfg1_cfg2_cqp12[cfg1_I][cfg2_I];
            auto& F2_C2D_cqp12_cqp12 = F2_C4D_cfg1_cfg2_cqp12_cqp12[cfg1_I][cfg2_I];
            auto& F4_C4D_cqp12_cqp12_cqp12_cqp12 = F4_C6D_cfg1_cfg2_cqp12_cqp12_cqp12_cqp12[cfg1_I][cfg2_I];

            // F⁰ = pf(S).
            F0_C2D_cfg1_cfg2[cfg1_I][cfg2_I] = calc_signed_pfaffian_minor();

            // Zero-based i: F¹ᵢ = (−1)ⁱ⁺¹ pf(S without i).
            F1_C1D_cqp12.setZero(Ncqp12_I);
            if (Ncqp12_I % 2 != 0) {
                for (int i_I = 0; i_I < Ncqp12_I; ++i_I) {F1_C1D_cqp12(i_I) = calc_signed_pfaffian_minor(i_I); }
                continue;
            }

            // F²ᵢ₁ᵢ₂ = −F²ᵢ₂ᵢ₁.
            F2_C2D_cqp12_cqp12.setZero(Ncqp12_I, Ncqp12_I);
            for (int i1_I = 0; i1_I < Ncqp12_I; ++i1_I) {
                for (int i2_I = i1_I + 1; i2_I < Ncqp12_I; ++i2_I) {
                    F2_C2D_cqp12_cqp12(i1_I, i2_I) = calc_signed_pfaffian_minor(i1_I, i2_I);
                    F2_C2D_cqp12_cqp12(i2_I, i1_I) = -F2_C2D_cqp12_cqp12(i1_I, i2_I);
                }
            }

            // F⁴: i₁ < i₂ < i₃ < i₄ → all permutations.
            F4_C4D_cqp12_cqp12_cqp12_cqp12 = Eigen::Tensor<doubleC, 4, Eigen::ColMajor>(Ncqp12_I, Ncqp12_I, Ncqp12_I, Ncqp12_I).setZero();
            for (int i1_I = 0; i1_I < Ncqp12_I; ++i1_I) {
                for (int i2_I = i1_I + 1; i2_I < Ncqp12_I; ++i2_I) {
                    for (int i3_I = i2_I + 1; i3_I < Ncqp12_I; ++i3_I) {
                        for (int i4_I = i3_I + 1; i4_I < Ncqp12_I; ++i4_I) {
                            const doubleC F4_C = calc_signed_pfaffian_minor(i1_I, i2_I, i3_I, i4_I);
                            const std::array<int, 4> chain_I1D_cqp{i1_I, i2_I, i3_I, i4_I};
                            // F⁴[p(i₁,i₂,i₃,i₄)] = sign(p) F⁴[i₁,i₂,i₃,i₄].
                            for (const auto& [a_I, b_I, c_I, d_I, sign_I] : permutationSign_I2D_permutation_entry) {
                                F4_C4D_cqp12_cqp12_cqp12_cqp12(chain_I1D_cqp[a_I], chain_I1D_cqp[b_I], chain_I1D_cqp[c_I], chain_I1D_cqp[d_I]) = static_cast<double>(sign_I) * F4_C;
                            }
                        }
                    }
                }
            }
        }
    }
}

inline const Eigen::MatrixXcd& HFBPfaffian::calc_overlap() {
    // Nκ₁κ₂ = ⟨Φ₁|Φ₂⟩ F⁰κ₁κ₂.
    for (int cfg2_I = 0; cfg2_I < static_cast<int>(config2_I2D_cfg2_cqp2.size()); ++cfg2_I) {
        for (int cfg1_I = 0; cfg1_I < static_cast<int>(config1_I2D_cfg1_cqp1.size()); ++cfg1_I) {
            Overlap_C2D_cfg1_cfg2(cfg1_I, cfg2_I) = overlap_C * F0_C2D_cfg1_cfg2[cfg1_I][cfg2_I];
        }
    }
    return Overlap_C2D_cfg1_cfg2;
}

inline const Eigen::MatrixXcd& HFBPfaffian::calc_one_fermion(const Eigen::VectorXcd& x_C1D_sp, bool creation_B) {
    assert(x_C1D_sp.size() == Nsp_I && x_C1D_sp.allFinite());

    // ⟨β₁X⟩; ⟨Xβ₂†⟩.
    const Eigen::VectorXcd Qp1X_C1D_qp1 = (creation_B ? Qp1SpDag_C2D_qp1_sp : Qp1Sp_C2D_qp1_sp) * x_C1D_sp;
    const Eigen::RowVectorXcd XQp2Dag_C1D_qp2 = x_C1D_sp.transpose() * (creation_B ? SpDagQp2Dag_C2D_sp_qp2 : SpQp2Dag_C2D_sp_qp2);

    OneFermion_C2D_cfg1_cfg2.setZero();
    for (int cfg2_I = 0; cfg2_I < static_cast<int>(config2_I2D_cfg2_cqp2.size()); ++cfg2_I) {
        const auto& config2_I1D_cqp2 = config2_I2D_cfg2_cqp2[cfg2_I];
        const Eigen::Index Ncqp2_I = config2_I1D_cqp2.size();
        for (int cfg1_I = 0; cfg1_I < static_cast<int>(config1_I2D_cfg1_cqp1.size()); ++cfg1_I) {
            const auto& config1_I1D_cqp1 = config1_I2D_cfg1_cqp1[cfg1_I];
            const Eigen::Index Ncqp1_I = config1_I1D_cqp1.size();
            if ((Ncqp1_I + Ncqp2_I) % 2 == 0) { continue; }

            const auto& F1_C1D_cqp12 = F1_C3D_cfg1_cfg2_cqp12[cfg1_I][cfg2_I];
            const double sign_F = 1.0 - 2.0 * (Ncqp1_I % 2);

            // Left configuration reversed; transpose without conjugation.
            OneFermion_C2D_cfg1_cfg2(cfg1_I, cfg2_I) = overlap_C * sign_F * (F1_C1D_cqp12.head(Ncqp1_I).reverse().transpose() * Qp1X_C1D_qp1(config1_I1D_cqp1) - XQp2Dag_C1D_qp2(config2_I1D_cqp2) * F1_C1D_cqp12.tail(Ncqp2_I)).value();
        }
    }
    return OneFermion_C2D_cfg1_cfg2;
}

inline const Eigen::MatrixXcd& HFBPfaffian::calc_one_body(const Eigen::MatrixXcd& OneBody_C2D_sp_sp) {
    assert(OneBody_C2D_sp_sp.rows() == Nsp_I && OneBody_C2D_sp_sp.cols() == Nsp_I && OneBody_C2D_sp_sp.allFinite());

    // S₁ = [⟨β₁x₁⟩; −⟨x₁β₂†⟩ᵀ], x₁ = c†,  S₂ = [⟨β₁x₂⟩; −⟨x₂β₂†⟩ᵀ], x₂ = c.
    Eigen::MatrixXcd S1_C2D_qp12_sp(2 * Nsp_I, Nsp_I);
    Eigen::MatrixXcd S2_C2D_qp12_sp(2 * Nsp_I, Nsp_I);
    Eigen::MatrixXcd T2_C2D_qp12_qp12(2 * Nsp_I, 2 * Nsp_I);
    std::vector<int> qp_I1D_cqp12{};
    qp_I1D_cqp12.reserve(2 * Nsp_I);

    OneBody_C2D_cfg1_cfg2.setZero();
    // T⁰ = Σαβ hαβ ⟨cα†cβ⟩; T² = S₁hS₂ᵀ.
    S1_C2D_qp12_sp.topRows(Nsp_I) = Qp1SpDag_C2D_qp1_sp;
    S1_C2D_qp12_sp.bottomRows(Nsp_I) = -SpDagQp2Dag_C2D_sp_qp2.transpose();
    S2_C2D_qp12_sp.topRows(Nsp_I) = Qp1Sp_C2D_qp1_sp;
    S2_C2D_qp12_sp.bottomRows(Nsp_I) = -SpQp2Dag_C2D_sp_qp2.transpose();
    const doubleC T0_C = OneBody_C2D_sp_sp.cwiseProduct(SpDagSp_C2D_sp_sp).sum();
    T2_C2D_qp12_qp12.noalias() = S1_C2D_qp12_sp * OneBody_C2D_sp_sp * S2_C2D_qp12_sp.transpose();
    T2_C2D_qp12_qp12 = (T2_C2D_qp12_qp12 - T2_C2D_qp12_qp12.transpose()).eval();

    // T² ← T² − T²ᵀ; K = F⁰T⁰ + ½Σij F²ij T²ij.
    for (int cfg2_I = 0; cfg2_I < static_cast<int>(config2_I2D_cfg2_cqp2.size()); ++cfg2_I) {
        const auto& config2_I1D_cqp = config2_I2D_cfg2_cqp2[cfg2_I];
        for (int cfg1_I = 0; cfg1_I < static_cast<int>(config1_I2D_cfg1_cqp1.size()); ++cfg1_I) {
            const auto& config1_I1D_cqp = config1_I2D_cfg1_cqp1[cfg1_I];
            const int Ncqp1_I = static_cast<int>(config1_I1D_cqp.size());
            const int Nchain_I = Ncqp1_I + static_cast<int>(config2_I1D_cqp.size());
            if (Nchain_I % 2 != 0) { continue; }
            // Left indices reversed; right indices shifted by Nsp.
            qp_I1D_cqp12.clear();
            for (int cqp1_I = Ncqp1_I - 1; cqp1_I >= 0; --cqp1_I) { qp_I1D_cqp12.push_back(config1_I1D_cqp[cqp1_I]); }
            for (const int qp_I : config2_I1D_cqp) { qp_I1D_cqp12.push_back(Nsp_I + qp_I); }
            const auto& F2_C2D_cqp12_cqp12 = F2_C4D_cfg1_cfg2_cqp12_cqp12[cfg1_I][cfg2_I];
            const doubleC kernel_C = F0_C2D_cfg1_cfg2[cfg1_I][cfg2_I] * T0_C + 0.5 * F2_C2D_cqp12_cqp12.cwiseProduct(T2_C2D_qp12_qp12(qp_I1D_cqp12, qp_I1D_cqp12)).sum();
            OneBody_C2D_cfg1_cfg2(cfg1_I, cfg2_I) = kernel_C;
        }
    }
    OneBody_C2D_cfg1_cfg2 *= overlap_C;
    return OneBody_C2D_cfg1_cfg2;
}

inline const Eigen::MatrixXcd& HFBPfaffian::calc_two_body(const Eigen::Tensor<doubleC, 4, Eigen::ColMajor>& TwoBody_C4D_sp_sp_sp_sp) {
    assert(TwoBody_C4D_sp_sp_sp_sp.dimension(0) == Nsp_I && TwoBody_C4D_sp_sp_sp_sp.dimension(1) == Nsp_I && TwoBody_C4D_sp_sp_sp_sp.dimension(2) == Nsp_I && TwoBody_C4D_sp_sp_sp_sp.dimension(3) == Nsp_I);

    // X = (x₁†,x₂†,x₄,x₃); workspaces shared across β,γ,δ.
    // Sₐ = [⟨β₁Xₐ⟩; −⟨Xₐβ₂†⟩ᵀ], a = 1,2,3,4.
    Eigen::MatrixXcd S_C2D_qp12_X(2 * Nsp_I, 4);
    auto S1_C1D_qp12 = S_C2D_qp12_X.col(0);
    auto S2_C1D_qp12 = S_C2D_qp12_X.col(1);
    auto S3_C1D_qp12 = S_C2D_qp12_X.col(2);
    auto S4_C1D_qp12 = S_C2D_qp12_X.col(3);
    Eigen::Matrix4cd XX_C2D_X_X;
    Eigen::Matrix4cd SF2S_C2D_X_X;
    // Ncqp12 ≤ 2Nsp; capacity = (2Nsp)².
    Eigen::VectorXcd S1S2workspace_C1D_element(4 * Nsp_I * Nsp_I);
    Eigen::VectorXcd S3S4workspace_C1D_element(4 * Nsp_I * Nsp_I);
    const auto accumulate_kernel = [&](Eigen::Ref<const Eigen::VectorXcd> x1_C1D_sp, int sp2_I, int sp3_I, int sp4_I) {
        S1_C1D_qp12.head(Nsp_I).noalias() = Qp1SpDag_C2D_qp1_sp * x1_C1D_sp;
        S1_C1D_qp12.tail(Nsp_I).noalias() = -SpDagQp2Dag_C2D_sp_qp2.transpose() * x1_C1D_sp;
        S2_C1D_qp12.head(Nsp_I) = Qp1SpDag_C2D_qp1_sp.col(sp2_I);
        S2_C1D_qp12.tail(Nsp_I) = -SpDagQp2Dag_C2D_sp_qp2.row(sp2_I).transpose();
        S3_C1D_qp12.head(Nsp_I) = Qp1Sp_C2D_qp1_sp.col(sp4_I);
        S3_C1D_qp12.tail(Nsp_I) = -SpQp2Dag_C2D_sp_qp2.row(sp4_I).transpose();
        S4_C1D_qp12.head(Nsp_I) = Qp1Sp_C2D_qp1_sp.col(sp3_I);
        S4_C1D_qp12.tail(Nsp_I) = -SpQp2Dag_C2D_sp_qp2.row(sp3_I).transpose();

        // Six upper-triangle contractions in the same insertion order.
        XX_C2D_X_X.diagonal().setZero();
        XX_C2D_X_X(0, 1) = (x1_C1D_sp.transpose() * SpDagSpDag_C2D_sp_sp.col(sp2_I)).value();
        XX_C2D_X_X(0, 2) = (x1_C1D_sp.transpose() * SpDagSp_C2D_sp_sp.col(sp4_I)).value();
        XX_C2D_X_X(0, 3) = (x1_C1D_sp.transpose() * SpDagSp_C2D_sp_sp.col(sp3_I)).value();
        XX_C2D_X_X(1, 2) = SpDagSp_C2D_sp_sp(sp2_I, sp4_I);
        XX_C2D_X_X(1, 3) = SpDagSp_C2D_sp_sp(sp2_I, sp3_I);
        XX_C2D_X_X(2, 3) = SpSp_C2D_sp_sp(sp4_I, sp3_I);
        // XX_ji = −XX_ij; no complex conjugation.
        XX_C2D_X_X(1, 0) = -XX_C2D_X_X(0, 1);
        XX_C2D_X_X(2, 0) = -XX_C2D_X_X(0, 2);
        XX_C2D_X_X(3, 0) = -XX_C2D_X_X(0, 3);
        XX_C2D_X_X(2, 1) = -XX_C2D_X_X(1, 2);
        XX_C2D_X_X(3, 1) = -XX_C2D_X_X(1, 3);
        XX_C2D_X_X(3, 2) = -XX_C2D_X_X(2, 3);
        const doubleC T0_C = XX_C2D_X_X(0, 1) * XX_C2D_X_X(2, 3)
                          - XX_C2D_X_X(0, 2) * XX_C2D_X_X(1, 3)
                          + XX_C2D_X_X(0, 3) * XX_C2D_X_X(1, 2);

        for (int cfg2_I = 0; cfg2_I < static_cast<int>(config2_I2D_cfg2_cqp2.size()); ++cfg2_I) {
            const auto& config2_I1D_cqp2 = config2_I2D_cfg2_cqp2[cfg2_I];
            const Eigen::Index Ncqp2_I = config2_I1D_cqp2.size();
            for (int cfg1_I = 0; cfg1_I < static_cast<int>(config1_I2D_cfg1_cqp1.size()); ++cfg1_I) {
                const auto& config1_I1D_cqp1 = config1_I2D_cfg1_cqp1[cfg1_I];
                const Eigen::Index Ncqp1_I = config1_I1D_cqp1.size();
                const Eigen::Index Ncqp12_I = Ncqp1_I + Ncqp2_I;
                if (Ncqp12_I % 2 != 0) { continue; }

                // S: selected configurations; left order reversed.
                Eigen::Map<Eigen::MatrixXcd> S_C2D_cqp12_X(Sworkspace_C1D_element.data(), Ncqp12_I, 4);
                S_C2D_cqp12_X.topRows(Ncqp1_I) = S_C2D_qp12_X.topRows(Nsp_I)(config1_I1D_cqp1, Eigen::placeholders::all).colwise().reverse();
                S_C2D_cqp12_X.bottomRows(Ncqp2_I) = S_C2D_qp12_X.bottomRows(Nsp_I)(config2_I1D_cqp2, Eigen::placeholders::all);

                // SF²S = SᵀF²S; six complementary contractions.
                SF2S_C2D_X_X.noalias() = S_C2D_cqp12_X.transpose()
                    * F2_C4D_cfg1_cfg2_cqp12_cqp12[cfg1_I][cfg2_I] * S_C2D_cqp12_X;
                doubleC kernel_C = F0_C2D_cfg1_cfg2[cfg1_I][cfg2_I] * T0_C
                    + SF2S_C2D_X_X(0, 1) * XX_C2D_X_X(2, 3)
                    - SF2S_C2D_X_X(0, 2) * XX_C2D_X_X(1, 3)
                    + SF2S_C2D_X_X(0, 3) * XX_C2D_X_X(1, 2)
                    + SF2S_C2D_X_X(1, 2) * XX_C2D_X_X(0, 3)
                    - SF2S_C2D_X_X(1, 3) * XX_C2D_X_X(0, 2)
                    + SF2S_C2D_X_X(2, 3) * XX_C2D_X_X(0, 1);

                // K⁴ = vec(S₁S₂ᵀ)ᵀ F⁴_(ij,kl) vec(S₃S₄ᵀ).
                if (Ncqp12_I >= 4) {
                    Eigen::Map<Eigen::MatrixXcd> S1S2_C2D_cqp12_cqp12(S1S2workspace_C1D_element.data(), Ncqp12_I, Ncqp12_I);
                    Eigen::Map<Eigen::MatrixXcd> S3S4_C2D_cqp12_cqp12(S3S4workspace_C1D_element.data(), Ncqp12_I, Ncqp12_I);
                    S1S2_C2D_cqp12_cqp12.noalias() = S_C2D_cqp12_X.col(0) * S_C2D_cqp12_X.col(1).transpose();
                    S3S4_C2D_cqp12_cqp12.noalias() = S_C2D_cqp12_X.col(2) * S_C2D_cqp12_X.col(3).transpose();
                    // ColMajor: (i,j) → i+nj; (k,l) → k+nl.
                    const Eigen::Map<const Eigen::MatrixXcd> F4_C2D_cqp12cqp12_cqp12cqp12(F4_C6D_cfg1_cfg2_cqp12_cqp12_cqp12_cqp12[cfg1_I][cfg2_I].data(), Ncqp12_I * Ncqp12_I, Ncqp12_I * Ncqp12_I);
                    const Eigen::Map<const Eigen::VectorXcd> S1S2_C1D_cqp12cqp12(S1S2_C2D_cqp12_cqp12.data(), S1S2_C2D_cqp12_cqp12.size());
                    const Eigen::Map<const Eigen::VectorXcd> S3S4_C1D_cqp12cqp12(S3S4_C2D_cqp12_cqp12.data(), S3S4_C2D_cqp12_cqp12.size());
                    kernel_C += (S1S2_C1D_cqp12cqp12.transpose() * F4_C2D_cqp12cqp12_cqp12cqp12 * S3S4_C1D_cqp12cqp12).value();
                }
                TwoBody_C2D_cfg1_cfg2(cfg1_I, cfg2_I) += kernel_C;
            }
        }
    };

    // H₂ = Σβγδ (Σα TwoBody_αβγδ cα†) cβ†cδcγ.
    TwoBody_C2D_cfg1_cfg2.setZero();
    for (int sp4_I = 0; sp4_I < Nsp_I; ++sp4_I) {
        for (int sp3_I = 0; sp3_I < Nsp_I; ++sp3_I) {
            for (int sp2_I = 0; sp2_I < Nsp_I; ++sp2_I) {
                const Eigen::Index offset_I = static_cast<Eigen::Index>(Nsp_I) * (sp2_I + static_cast<Eigen::Index>(Nsp_I) * (sp3_I + static_cast<Eigen::Index>(Nsp_I) * sp4_I));
                const Eigen::Map<const Eigen::VectorXcd> x1_C1D_sp(TwoBody_C4D_sp_sp_sp_sp.data() + offset_I, Nsp_I);
                assert(x1_C1D_sp.allFinite());
                if (x1_C1D_sp.isZero(0.0)) { continue; }
                accumulate_kernel(x1_C1D_sp, sp2_I, sp3_I, sp4_I);
            }
        }
    }
    TwoBody_C2D_cfg1_cfg2 *= overlap_C;
    return TwoBody_C2D_cfg1_cfg2;
}


inline doubleC HFBPfaffian::calc_pfaffian(Eigen::Ref<Eigen::MatrixXcd> S_C2D_chain_chain) {
    const Eigen::Index Nchain_I = S_C2D_chain_chain.rows();
    assert(S_C2D_chain_chain.cols() == Nchain_I && Nchain_I % 2 == 0);
    assert(S_C2D_chain_chain.allFinite());
    assert(S_C2D_chain_chain.isApprox(-S_C2D_chain_chain.transpose(), 1.0e-12));
    doubleC pf_C(1.0, 0.0);

    // S → PSPᵀ → 2×2 pivot ⊕ Schur complement.
    for (Eigen::Index k_I = 0; k_I < Nchain_I; k_I += 2) {
        // Pivot: argmax_{j>k}|S_kj|.
        Eigen::Index pivotOffset_I = 0;
        const double pivotAbs_F = S_C2D_chain_chain.row(k_I).segment(k_I + 1, Nchain_I - k_I - 1).cwiseAbs().maxCoeff(&pivotOffset_I);
        const Eigen::Index pivot_I = k_I + 1 + pivotOffset_I;
        if (pivotAbs_F == 0.0) {return doubleC(0.0, 0.0);}

        // One simultaneous row-column swap contributes -1.
        if (pivot_I != k_I + 1) {
            S_C2D_chain_chain.row(k_I + 1).swap(S_C2D_chain_chain.row(pivot_I));
            S_C2D_chain_chain.col(k_I + 1).swap(S_C2D_chain_chain.col(pivot_I));
            pf_C = -pf_C;
        }
        const doubleC pivot_C = S_C2D_chain_chain(k_I, k_I + 1);
        pf_C *= pivot_C;

        // S_remain ← S_remain + (vᵀu − uᵀv)/pivot.
        const Eigen::Index Nremain_I = Nchain_I - k_I - 2;
        const auto u_C1D_remain = S_C2D_chain_chain.row(k_I).tail(Nremain_I);
        const auto v_C1D_remain = S_C2D_chain_chain.row(k_I + 1).tail(Nremain_I);
        auto S_C2D_remain_remain = S_C2D_chain_chain.bottomRightCorner(Nremain_I, Nremain_I);
        S_C2D_remain_remain.triangularView<Eigen::StrictlyUpper>() += ((v_C1D_remain.transpose() * u_C1D_remain - u_C1D_remain.transpose() * v_C1D_remain) / pivot_C).eval();
        S_C2D_remain_remain.triangularView<Eigen::StrictlyLower>() = (-S_C2D_remain_remain.transpose()).eval();
    }
    return pf_C;
}
