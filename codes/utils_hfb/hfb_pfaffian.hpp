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
 * @note   Updated contractions; overlap ≠ 0; caller-owned outputs.
 */
class HFBPfaffian {
public:
    int Nsp_I = 0;
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

    // Sorted unique quasiparticles used by each configuration table.
    std::vector<int> qp1_I1D_compactqp1{};
    std::vector<int> qp2_I1D_compactqp2{};
    std::vector<std::vector<int>> compactqp1_I2D_cfg1_cqp1{};
    std::vector<std::vector<int>> compactqp2_I2D_cfg2_cqp2{};

    // Signed Pfaffian minors; excludes vacuum overlap.
    std::vector<std::vector<doubleC>> F0_C2D_cfg1_cfg2{};
    std::vector<std::vector<Eigen::VectorXcd>> F1_C3D_cfg1_cfg2_cqp12{};
    std::vector<std::vector<Eigen::MatrixXcd>> F2_C4D_cfg1_cfg2_cqp12_cqp12{};
    std::vector<std::vector<Eigen::Tensor<doubleC, 4, Eigen::ColMajor>>> F4_C6D_cfg1_cfg2_cqp12_cqp12_cqp12_cqp12{};


private:
    Eigen::MatrixXcd U1_C2D_sp_qp1{};
    Eigen::MatrixXcd V1_C2D_sp_qp1{};
    Eigen::MatrixXcd Z1_C2D_sp_sp{};
    double nu1_F = 0.0;
    Eigen::MatrixXcd U2_C2D_sp_qp2{};
    Eigen::MatrixXcd V2_C2D_sp_qp2{};
    Eigen::MatrixXcd Z2_C2D_sp_sp{};
    double nu2_F = 0.0;
    Eigen::MatrixXcd A_C2D_qp1_qp2{};
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

        // U₁,V₁,U₂,V₂,Z₁,Z₂ ∈ ℂ^{Nsp×Nsp}.
        U1_C2D_sp_qp1.resize(Nsp_I, Nsp_I);
        V1_C2D_sp_qp1.resize(Nsp_I, Nsp_I);
        U2_C2D_sp_qp2.resize(Nsp_I, Nsp_I);
        V2_C2D_sp_qp2.resize(Nsp_I, Nsp_I);
        Z1_C2D_sp_sp.resize(Nsp_I, Nsp_I);
        Z2_C2D_sp_sp.resize(Nsp_I, Nsp_I);

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


        // compactqp₁ = ∪κ₁; compactqp₂ = ∪κ₂.
        for (const auto& config_I1D_cqp : config1_I2D_cfg1_cqp1) { qp1_I1D_compactqp1.insert(qp1_I1D_compactqp1.end(), config_I1D_cqp.begin(), config_I1D_cqp.end()); }
        for (const auto& config_I1D_cqp : config2_I2D_cfg2_cqp2) { qp2_I1D_compactqp2.insert(qp2_I1D_compactqp2.end(), config_I1D_cqp.begin(), config_I1D_cqp.end()); }
        std::sort(qp1_I1D_compactqp1.begin(), qp1_I1D_compactqp1.end());
        std::sort(qp2_I1D_compactqp2.begin(), qp2_I1D_compactqp2.end());
        qp1_I1D_compactqp1.erase(std::unique(qp1_I1D_compactqp1.begin(), qp1_I1D_compactqp1.end()), qp1_I1D_compactqp1.end());
        qp2_I1D_compactqp2.erase(std::unique(qp2_I1D_compactqp2.begin(), qp2_I1D_compactqp2.end()), qp2_I1D_compactqp2.end());

        // Original qp → compactqp; unused entries = −1.
        std::vector<int> compactqp1_I1D_qp1(Nsp_I, -1);
        std::vector<int> compactqp2_I1D_qp2(Nsp_I, -1);
        for (int compactqp1_I = 0; compactqp1_I < qp1_I1D_compactqp1.size(); ++compactqp1_I) { compactqp1_I1D_qp1[qp1_I1D_compactqp1[compactqp1_I]] = compactqp1_I; }
        for (int compactqp2_I = 0; compactqp2_I < qp2_I1D_compactqp2.size(); ++compactqp2_I) { compactqp2_I1D_qp2[qp2_I1D_compactqp2[compactqp2_I]] = compactqp2_I; }

        // Configurations: original qp → compactqp; order unchanged.
        compactqp1_I2D_cfg1_cqp1 = config1_I2D_cfg1_cqp1;
        compactqp2_I2D_cfg2_cqp2 = config2_I2D_cfg2_cqp2;
        for (auto& compactqp1_I1D_cqp1 : compactqp1_I2D_cfg1_cqp1) {
            for (int& qp1_I : compactqp1_I1D_cqp1) { qp1_I = compactqp1_I1D_qp1[qp1_I]; }
        }
        for (auto& compactqp2_I1D_cqp2 : compactqp2_I2D_cfg2_cqp2) {
            for (int& qp2_I : compactqp2_I1D_cqp2) { qp2_I = compactqp2_I1D_qp2[qp2_I]; }
        }

        // Lmax = 2Nsp+4; size(Sworkspace) = Lmax².
        const int Lmax_I = 2 * Nsp_I + 4;
        Sworkspace_C1D_element.resize(Lmax_I * Lmax_I);
    }

    /**
     * @brief Cache state 1 using full-pivot LU.
     * @math Z1 = (V1U1⁻¹)*; ν1 = √|det U1|.
     * @output Stored U1, V1, Z1, ν1.
     */
    void update_U1V1(const Eigen::MatrixXcd& U1_C2D_sp_qp1_, const Eigen::MatrixXcd& V1_C2D_sp_qp1_);

    /**
     * @brief Cache state 2 using full-pivot LU.
     * @math Z2 = (V2U2⁻¹)*; ν2 = √|det U2|.
     * @output Stored U2, V2, Z2, ν2.
     */
    void update_U2V2(const Eigen::MatrixXcd& U2_C2D_sp_qp2_, const Eigen::MatrixXcd& V2_C2D_sp_qp2_);

    /**
     * @brief Update contractions and signed minors using LU/Pfaffians.
     * @math A = U₁ᵀU₂* + V₁ᵀV₂*; Qp1Qp2Dag = A⁻ᵀ.
     * @math F⁰ = pf(S); F¹,F²,F⁴ = signed deletion minors.
     * @output Updated overlap, A, eleven contractions, F0/F1/F2/F4.
     * @note Requires update_U1V1(), update_U2V2(); numerically invertible A.
     */
    void update_contractions(doubleC phase_C);

    /**
     * @brief  Calculate configuration overlaps using cached Pfaffian minors.
     * @math   N_{κ₁κ₂} = ⟨Φ₁;κ₁|Φ₂;κ₂⟩.
     * @output Overwritten caller-provided result_C2D_cfg1_cfg2.
     */
    void calc_overlap(Eigen::Ref<Eigen::MatrixXcd> result_C2D_cfg1_cfg2);

    /**
     * @brief  Calculate single-fermion kernels using cached Pfaffian minors.
     * @math   K_abα = ⟨Φ₁;a|cα† (true), cα (false)|Φ₂;b⟩.
     * @output Overwritten caller-provided result_C3D_cfg1_cfg2_sp.
     * @note   Caller allocates Ncfg1 × Ncfg2 × Nsp.
     */
    void calc_one_fermion(Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& result_C3D_cfg1_cfg2_sp, bool creation_B);

    /**
     * @brief  Calculate one-body kernels using cached Pfaffian minors.
     * @math   O = Σαβ O_αβ cα†cβ; K = F⁰T⁰ − Σi<j F²ij T²ij.
     * @output Overwritten caller-provided result_C2D_cfg1_cfg2.
     */
    void calc_one_body(Eigen::Ref<Eigen::MatrixXcd> result_C2D_cfg1_cfg2, const Eigen::MatrixXcd& OneBody_C2D_sp_sp);

    /**
     * @brief  Calculate two-body kernels using cached Pfaffian minors.
     * @math   O₂ = Σαβγδ O_αβγδ cα†cβ†cδcγ.
     * @output Overwritten caller-provided result_C2D_cfg1_cfg2.
     * @note   Input coefficients include all numerical prefactors.
     */
    void calc_two_body(Eigen::Ref<Eigen::MatrixXcd> result_C2D_cfg1_cfg2, const Eigen::Tensor<doubleC, 4, Eigen::ColMajor>& TwoBody_C4D_sp_sp_sp_sp);

    /**
     * @brief  Calculate separable Q/P kernels using cached Pfaffian minors.
     * @math   K = Σq gQq(Qq†Qq − one-body contraction) + Σp gPp Pp†Pp.
     * @math   P† = ½Σab Pab c†a c†b; Pᵀ = −P.
     * @output Overwritten caller-provided result_C2D_cfg1_cfg2.
     * @note   gQ,gP include signs; pairing ¼ is internal.
     */
    void calc_two_body(Eigen::Ref<Eigen::MatrixXcd> result_C2D_cfg1_cfg2, const Eigen::VectorXd& gQ_F1D_q, const Eigen::VectorXd& gP_F1D_p, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& Q_C3D_sp_sp_q, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& P_C3D_sp_sp_p);

    /**
     * @brief  Calculate Pfaffians using pivoted skew-symmetric elimination.
     * @math   pf(S)² = det(S); pf(∅) = 1.
     * @output Pfaffian value; input matrix overwritten.
     * @note   Finite S = −Sᵀ ∈ ℂ²ⁿˣ²ⁿ; zero pivot → 0; no cutoff.
     */
    static doubleC calc_pfaffian(Eigen::Ref<Eigen::MatrixXcd> S_C2D_chain_chain);

};


inline void HFBPfaffian::update_U1V1(const Eigen::MatrixXcd& U1_C2D_sp_qp1_, const Eigen::MatrixXcd& V1_C2D_sp_qp1_) {
    // U,V ∈ ℂ^{Nsp×Nsp}; UᵀZ† = Vᵀ.
    assert(U1_C2D_sp_qp1_.rows() == Nsp_I && U1_C2D_sp_qp1_.cols() == Nsp_I && U1_C2D_sp_qp1_.allFinite());
    assert(V1_C2D_sp_qp1_.rows() == Nsp_I && V1_C2D_sp_qp1_.cols() == Nsp_I && V1_C2D_sp_qp1_.allFinite());
    U1_C2D_sp_qp1 = U1_C2D_sp_qp1_;
    V1_C2D_sp_qp1 = V1_C2D_sp_qp1_;
    U1_lu.compute(U1_C2D_sp_qp1);
    assert(U1_lu.isInvertible());
    Z1_C2D_sp_sp = U1_lu.transpose().solve(V1_C2D_sp_qp1.transpose());
    Z1_C2D_sp_sp.adjointInPlace();
    nu1_F = std::sqrt(std::abs(U1_lu.determinant()));
}

inline void HFBPfaffian::update_U2V2(const Eigen::MatrixXcd& U2_C2D_sp_qp2_, const Eigen::MatrixXcd& V2_C2D_sp_qp2_) {
    // U,V ∈ ℂ^{Nsp×Nsp}; UᵀZ† = Vᵀ.
    assert(U2_C2D_sp_qp2_.rows() == Nsp_I && U2_C2D_sp_qp2_.cols() == Nsp_I && U2_C2D_sp_qp2_.allFinite());
    assert(V2_C2D_sp_qp2_.rows() == Nsp_I && V2_C2D_sp_qp2_.cols() == Nsp_I && V2_C2D_sp_qp2_.allFinite());
    U2_C2D_sp_qp2 = U2_C2D_sp_qp2_;
    V2_C2D_sp_qp2 = V2_C2D_sp_qp2_;
    U2_lu.compute(U2_C2D_sp_qp2);
    assert(U2_lu.isInvertible());
    Z2_C2D_sp_sp = U2_lu.transpose().solve(V2_C2D_sp_qp2.transpose());
    Z2_C2D_sp_sp.adjointInPlace();
    nu2_F = std::sqrt(std::abs(U2_lu.determinant()));
}

inline void HFBPfaffian::update_contractions(doubleC phase_C) {
    assert(U1_C2D_sp_qp1.rows() == Nsp_I && U2_C2D_sp_qp2.rows() == Nsp_I && Nsp_I > 0);
    // X = [Z₂,−I; I,−Z₁*].
    Eigen::Map<Eigen::MatrixXcd> X_C2D_2sp_2sp(Sworkspace_C1D_element.data(), 2 * Nsp_I, 2 * Nsp_I);
    X_C2D_2sp_2sp.topLeftCorner(Nsp_I, Nsp_I) = Z2_C2D_sp_sp;
    X_C2D_2sp_2sp.bottomRightCorner(Nsp_I, Nsp_I) = -Z1_C2D_sp_sp.conjugate();

    // X₁₂ = −I; X₂₁ = I.
    X_C2D_2sp_2sp.topRightCorner(Nsp_I, Nsp_I).setIdentity();
    X_C2D_2sp_2sp.topRightCorner(Nsp_I, Nsp_I) *= -1.0;
    X_C2D_2sp_2sp.bottomLeftCorner(Nsp_I, Nsp_I).setIdentity();

    // overlap = phase × ν₁ν₂s_Nsp pf(X).
    const int sign_I = (Nsp_I % 4 == 0 || Nsp_I % 4 == 3) ? 1 : -1;
    overlap_C = phase_C * nu1_F * nu2_F * static_cast<double>(sign_I) * calc_pfaffian(X_C2D_2sp_2sp);

    // A = U₁ᵀU₂* + V₁ᵀV₂*.
    A_C2D_qp1_qp2.noalias() = U1_C2D_sp_qp1.transpose() * U2_C2D_sp_qp2.conjugate();
    A_C2D_qp1_qp2.noalias() += V1_C2D_sp_qp1.transpose() * V2_C2D_sp_qp2.conjugate();
    assert(A_C2D_qp1_qp2.allFinite());
    A_lu.compute(A_C2D_qp1_qp2);
    assert(A_lu.isInvertible());
    Qp1Qp2Dag_C2D_qp1_qp2 = A_lu.inverse().transpose();

    // ⟨β₁c⟩ = A⁻ᵀV₂†; ⟨β₁c†⟩ = A⁻ᵀU₂†.
    Qp1Sp_C2D_qp1_sp.noalias() = Qp1Qp2Dag_C2D_qp1_qp2 * V2_C2D_sp_qp2.adjoint();
    Qp1SpDag_C2D_qp1_sp.noalias() = Qp1Qp2Dag_C2D_qp1_qp2 * U2_C2D_sp_qp2.adjoint();

    // ⟨cβ₂†⟩ = U₁A⁻ᵀ; ⟨c†β₂†⟩ = V₁A⁻ᵀ.
    SpQp2Dag_C2D_sp_qp2.noalias() = U1_C2D_sp_qp1 * Qp1Qp2Dag_C2D_qp1_qp2;
    SpDagQp2Dag_C2D_sp_qp2.noalias() = V1_C2D_sp_qp1 * Qp1Qp2Dag_C2D_qp1_qp2;

    // ⟨c†c⟩ = V₁⟨β₁c⟩; ⟨cc†⟩ = U₁⟨β₁c†⟩.
    SpDagSp_C2D_sp_sp.noalias() = V1_C2D_sp_qp1 * Qp1Sp_C2D_qp1_sp;
    SpSpDag_C2D_sp_sp.noalias() = U1_C2D_sp_qp1 * Qp1SpDag_C2D_qp1_sp;

    // ⟨c†c†⟩ = V₁⟨β₁c†⟩; ⟨cc⟩ = U₁⟨β₁c⟩.
    SpDagSpDag_C2D_sp_sp.noalias() = V1_C2D_sp_qp1 * Qp1SpDag_C2D_qp1_sp;
    SpSp_C2D_sp_sp.noalias() = U1_C2D_sp_qp1 * Qp1Sp_C2D_qp1_sp;

    // ⟨β₁β₁⟩ = ⟨β₁c⟩U₁* + ⟨β₁c†⟩V₁*.
    Qp1Qp1_C2D_qp1_qp1.noalias() = Qp1Sp_C2D_qp1_sp * U1_C2D_sp_qp1.conjugate();
    Qp1Qp1_C2D_qp1_qp1.noalias() += Qp1SpDag_C2D_qp1_sp * V1_C2D_sp_qp1.conjugate();

    // ⟨β₂†β₂†⟩ = V₂ᵀ⟨cβ₂†⟩ + U₂ᵀ⟨c†β₂†⟩.
    Qp2DagQp2Dag_C2D_qp2_qp2.noalias() = V2_C2D_sp_qp2.transpose() * SpQp2Dag_C2D_sp_qp2;
    Qp2DagQp2Dag_C2D_qp2_qp2.noalias() += U2_C2D_sp_qp2.transpose() * SpDagQp2Dag_C2D_sp_qp2;

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

inline void HFBPfaffian::calc_overlap(Eigen::Ref<Eigen::MatrixXcd> result_C2D_cfg1_cfg2) {
    assert(result_C2D_cfg1_cfg2.rows() == static_cast<Eigen::Index>(config1_I2D_cfg1_cqp1.size()) && result_C2D_cfg1_cfg2.cols() == static_cast<Eigen::Index>(config2_I2D_cfg2_cqp2.size()));
    // Nκ₁κ₂ = ⟨Φ₁|Φ₂⟩ F⁰κ₁κ₂.
    for (int cfg2_I = 0; cfg2_I < static_cast<int>(config2_I2D_cfg2_cqp2.size()); ++cfg2_I) {
        for (int cfg1_I = 0; cfg1_I < static_cast<int>(config1_I2D_cfg1_cqp1.size()); ++cfg1_I) {
            result_C2D_cfg1_cfg2(cfg1_I, cfg2_I) = overlap_C * F0_C2D_cfg1_cfg2[cfg1_I][cfg2_I];
        }
    }
}

inline void HFBPfaffian::calc_one_fermion(Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& result_C3D_cfg1_cfg2_sp, bool creation_B) {
    assert(result_C3D_cfg1_cfg2_sp.dimension(0) == static_cast<Eigen::Index>(compactqp1_I2D_cfg1_cqp1.size()) && result_C3D_cfg1_cfg2_sp.dimension(1) == static_cast<Eigen::Index>(compactqp2_I2D_cfg2_cqp2.size()) && result_C3D_cfg1_cfg2_sp.dimension(2) == Nsp_I);
    const Eigen::Index Ncfg1_I = compactqp1_I2D_cfg1_cqp1.size();
    const Eigen::Index Ncfg2_I = compactqp2_I2D_cfg2_cqp2.size();
    // ⟨β₁cα†⟩ or ⟨β₁cα⟩; ⟨cα†β₂†⟩ or ⟨cαβ₂†⟩.
    const auto Qp1X_C2D_compactqp1_sp = (creation_B ? Qp1SpDag_C2D_qp1_sp : Qp1Sp_C2D_qp1_sp)(qp1_I1D_compactqp1, Eigen::placeholders::all);
    const auto XQp2Dag_C2D_sp_compactqp2 = (creation_B ? SpDagQp2Dag_C2D_sp_qp2 : SpQp2Dag_C2D_sp_qp2)(Eigen::placeholders::all, qp2_I1D_compactqp2);

    result_C3D_cfg1_cfg2_sp.setZero();
    for (Eigen::Index cfg2_I = 0; cfg2_I < Ncfg2_I; ++cfg2_I) {
        const auto& compactqp2_I1D_cqp2 = compactqp2_I2D_cfg2_cqp2[cfg2_I];
        const Eigen::Index Ncqp2_I = compactqp2_I1D_cqp2.size();
        for (Eigen::Index cfg1_I = 0; cfg1_I < Ncfg1_I; ++cfg1_I) {
            const auto& compactqp1_I1D_cqp1 = compactqp1_I2D_cfg1_cqp1[cfg1_I];
            const Eigen::Index Ncqp1_I = compactqp1_I1D_cqp1.size();
            if ((Ncqp1_I + Ncqp2_I) % 2 == 0) { continue; }

            const auto& F1_C1D_cqp12 = F1_C3D_cfg1_cfg2_cqp12[cfg1_I][cfg2_I];
            const double sign_F = 1.0 - 2.0 * (Ncqp1_I % 2);
            // K_: = overlap × (−1)^Ncqp1 [(F¹L)revᵀ C_L − (F¹R)ᵀ C_Rᵀ].
            // offset = cfg1 + Ncfg1 cfg2 + Ncfg1 Ncfg2 sp.
            Eigen::Map<Eigen::RowVectorXcd, 0, Eigen::InnerStride<Eigen::Dynamic>> OneFermion_C1D_sp(result_C3D_cfg1_cfg2_sp.data() + cfg1_I + Ncfg1_I * cfg2_I, Nsp_I, Eigen::InnerStride<Eigen::Dynamic>(Ncfg1_I * Ncfg2_I));
            OneFermion_C1D_sp.noalias() = F1_C1D_cqp12.head(Ncqp1_I).reverse().transpose() * Qp1X_C2D_compactqp1_sp(compactqp1_I1D_cqp1, Eigen::placeholders::all);
            OneFermion_C1D_sp.noalias() -= F1_C1D_cqp12.tail(Ncqp2_I).transpose() * XQp2Dag_C2D_sp_compactqp2(Eigen::placeholders::all, compactqp2_I1D_cqp2).transpose();
            OneFermion_C1D_sp *= overlap_C * sign_F;
        }
    }
}

inline void HFBPfaffian::calc_one_body(Eigen::Ref<Eigen::MatrixXcd> result_C2D_cfg1_cfg2, const Eigen::MatrixXcd& OneBody_C2D_sp_sp) {
    assert(result_C2D_cfg1_cfg2.rows() == static_cast<Eigen::Index>(compactqp1_I2D_cfg1_cqp1.size()) && result_C2D_cfg1_cfg2.cols() == static_cast<Eigen::Index>(compactqp2_I2D_cfg2_cqp2.size()));
    assert(OneBody_C2D_sp_sp.rows() == Nsp_I && OneBody_C2D_sp_sp.cols() == Nsp_I && OneBody_C2D_sp_sp.allFinite());

    const Eigen::Index Ncompactqp1_I = qp1_I1D_compactqp1.size();
    const Eigen::Index Ncompactqp2_I = qp2_I1D_compactqp2.size();
    const Eigen::Index Ncompactqp12_I = Ncompactqp1_I + Ncompactqp2_I;

    // S₁ = [⟨β₁x₁⟩; −⟨x₁β₂†⟩ᵀ], x₁ = c†,  S₂ = [⟨β₁x₂⟩; −⟨x₂β₂†⟩ᵀ], x₂ = c.
    Eigen::MatrixXcd S1_C2D_compactqp12_sp(Ncompactqp12_I, Nsp_I);
    Eigen::MatrixXcd S2_C2D_compactqp12_sp(Ncompactqp12_I, Nsp_I);
    Eigen::MatrixXcd T2_C2D_compactqp12_compactqp12(Ncompactqp12_I, Ncompactqp12_I);
    std::vector<int> compactqp12_I1D_cqp12{};
    compactqp12_I1D_cqp12.reserve(Ncompactqp12_I);

    result_C2D_cfg1_cfg2.setZero();
    S1_C2D_compactqp12_sp.topRows(Ncompactqp1_I) = Qp1SpDag_C2D_qp1_sp(qp1_I1D_compactqp1, Eigen::placeholders::all);
    S2_C2D_compactqp12_sp.topRows(Ncompactqp1_I) = Qp1Sp_C2D_qp1_sp(qp1_I1D_compactqp1, Eigen::placeholders::all);
    S1_C2D_compactqp12_sp.bottomRows(Ncompactqp2_I) = -SpDagQp2Dag_C2D_sp_qp2(Eigen::placeholders::all, qp2_I1D_compactqp2).transpose();
    S2_C2D_compactqp12_sp.bottomRows(Ncompactqp2_I) = -SpQp2Dag_C2D_sp_qp2(Eigen::placeholders::all, qp2_I1D_compactqp2).transpose();
    // T⁰ = Σαβ hαβ pf([0,t;−t,0]) = Σαβ hαβ t, t = ⟨cα†cβ⟩.
    const doubleC T0_C = OneBody_C2D_sp_sp.cwiseProduct(SpDagSp_C2D_sp_sp).sum();
    // S_{ij;αβ} = [S₁iα,S₂iβ; S₁jα,S₂jβ].
    // T²ij = Σαβ hαβ pf([0,S_{ij;αβ};−S_{ij;αβ}ᵀ,0]) = −Σαβ hαβ det(S_{ij;αβ}).
    // T² = (S₁hS₂ᵀ)ᵀ − S₁hS₂ᵀ.
    T2_C2D_compactqp12_compactqp12.noalias() = S1_C2D_compactqp12_sp * OneBody_C2D_sp_sp * S2_C2D_compactqp12_sp.transpose();
    T2_C2D_compactqp12_compactqp12 = (T2_C2D_compactqp12_compactqp12.transpose() - T2_C2D_compactqp12_compactqp12).eval();

    // K = F⁰T⁰ − ½Σij F²ij T²ij.
    for (int cfg2_I = 0; cfg2_I < static_cast<int>(compactqp2_I2D_cfg2_cqp2.size()); ++cfg2_I) {
        const auto& compactqp2_I1D_cqp2 = compactqp2_I2D_cfg2_cqp2[cfg2_I];
        for (int cfg1_I = 0; cfg1_I < static_cast<int>(compactqp1_I2D_cfg1_cqp1.size()); ++cfg1_I) {
            const auto& compactqp1_I1D_cqp1 = compactqp1_I2D_cfg1_cqp1[cfg1_I];
            const int Ncqp1_I = static_cast<int>(compactqp1_I1D_cqp1.size());
            const int Ncqp12_I = Ncqp1_I + static_cast<int>(compactqp2_I1D_cqp2.size());
            if (Ncqp12_I % 2 != 0) { continue; }
            // Left indices reversed; right indices shifted by Ncompactqp1.

            compactqp12_I1D_cqp12.clear();
            for (int cqp1_I = Ncqp1_I - 1; cqp1_I >= 0; --cqp1_I) { compactqp12_I1D_cqp12.push_back(compactqp1_I1D_cqp1[cqp1_I]); }
            for (const int compactqp2_I : compactqp2_I1D_cqp2) { compactqp12_I1D_cqp12.push_back(Ncompactqp1_I + compactqp2_I); }
            const auto T2_C2D_cqp12_cqp12 = T2_C2D_compactqp12_compactqp12(compactqp12_I1D_cqp12, compactqp12_I1D_cqp12);
            const auto& F2_C2D_cqp12_cqp12 = F2_C4D_cfg1_cfg2_cqp12_cqp12[cfg1_I][cfg2_I];
            const doubleC kernel_C = F0_C2D_cfg1_cfg2[cfg1_I][cfg2_I] * T0_C - 0.5 * F2_C2D_cqp12_cqp12.cwiseProduct(T2_C2D_cqp12_cqp12).sum();
            result_C2D_cfg1_cfg2(cfg1_I, cfg2_I) = kernel_C;
        }
    }
    result_C2D_cfg1_cfg2 *= overlap_C;
}

inline void HFBPfaffian::calc_two_body(Eigen::Ref<Eigen::MatrixXcd> result_C2D_cfg1_cfg2, const Eigen::Tensor<doubleC, 4, Eigen::ColMajor>& TwoBody_C4D_sp_sp_sp_sp) {
    assert(result_C2D_cfg1_cfg2.rows() == static_cast<Eigen::Index>(config1_I2D_cfg1_cqp1.size()) && result_C2D_cfg1_cfg2.cols() == static_cast<Eigen::Index>(config2_I2D_cfg2_cqp2.size()));
    assert(TwoBody_C4D_sp_sp_sp_sp.dimension(0) == Nsp_I && TwoBody_C4D_sp_sp_sp_sp.dimension(1) == Nsp_I && TwoBody_C4D_sp_sp_sp_sp.dimension(2) == Nsp_I && TwoBody_C4D_sp_sp_sp_sp.dimension(3) == Nsp_I);

    const Eigen::Index Ncompactqp1_I = qp1_I1D_compactqp1.size();
    const Eigen::Index Ncompactqp2_I = qp2_I1D_compactqp2.size();
    const Eigen::Index Ncompactqp12_I = Ncompactqp1_I + Ncompactqp2_I;

    // X = (x₁†,x₂†,x₄,x₃); workspaces shared across β,γ,δ.
    // Sₐ = [⟨β₁Xₐ⟩; −⟨Xₐβ₂†⟩ᵀ], a = 1,2,3,4.
    Eigen::MatrixXcd S_C2D_compactqp12_X(Ncompactqp12_I, 4);
    auto S1_C1D_compactqp12 = S_C2D_compactqp12_X.col(0);
    auto S2_C1D_compactqp12 = S_C2D_compactqp12_X.col(1);
    auto S3_C1D_compactqp12 = S_C2D_compactqp12_X.col(2);
    auto S4_C1D_compactqp12 = S_C2D_compactqp12_X.col(3);
    Eigen::Matrix4cd SF0_C2D_X_X;
    Eigen::Matrix4cd SF2S_C2D_X_X;
    // Ncqp12 ≤ 2Nsp; capacity = (2Nsp)².
    Eigen::VectorXcd S1S2workspace_C1D_element(4 * Nsp_I * Nsp_I);
    Eigen::VectorXcd S3S4workspace_C1D_element(4 * Nsp_I * Nsp_I);
    const auto accumulate_kernel = [&](Eigen::Ref<const Eigen::VectorXcd> x1_C1D_sp, int sp2_I, int sp3_I, int sp4_I) {
        S1_C1D_compactqp12.head(Ncompactqp1_I).noalias() = Qp1SpDag_C2D_qp1_sp(qp1_I1D_compactqp1, Eigen::placeholders::all) * x1_C1D_sp;
        S1_C1D_compactqp12.tail(Ncompactqp2_I).noalias() = -SpDagQp2Dag_C2D_sp_qp2(Eigen::placeholders::all, qp2_I1D_compactqp2).transpose() * x1_C1D_sp;
        S2_C1D_compactqp12.head(Ncompactqp1_I) = Qp1SpDag_C2D_qp1_sp(qp1_I1D_compactqp1, sp2_I);
        S2_C1D_compactqp12.tail(Ncompactqp2_I) = -SpDagQp2Dag_C2D_sp_qp2(sp2_I, qp2_I1D_compactqp2).transpose();
        S3_C1D_compactqp12.head(Ncompactqp1_I) = Qp1Sp_C2D_qp1_sp(qp1_I1D_compactqp1, sp4_I);
        S3_C1D_compactqp12.tail(Ncompactqp2_I) = -SpQp2Dag_C2D_sp_qp2(sp4_I, qp2_I1D_compactqp2).transpose();
        S4_C1D_compactqp12.head(Ncompactqp1_I) = Qp1Sp_C2D_qp1_sp(qp1_I1D_compactqp1, sp3_I);
        S4_C1D_compactqp12.tail(Ncompactqp2_I) = -SpQp2Dag_C2D_sp_qp2(sp3_I, qp2_I1D_compactqp2).transpose();

        // Six upper-triangle contractions in the same insertion order.
        SF0_C2D_X_X.diagonal().setZero();
        SF0_C2D_X_X(0, 1) = (x1_C1D_sp.transpose() * SpDagSpDag_C2D_sp_sp.col(sp2_I)).value();
        SF0_C2D_X_X(0, 2) = (x1_C1D_sp.transpose() * SpDagSp_C2D_sp_sp.col(sp4_I)).value();
        SF0_C2D_X_X(0, 3) = (x1_C1D_sp.transpose() * SpDagSp_C2D_sp_sp.col(sp3_I)).value();
        SF0_C2D_X_X(1, 2) = SpDagSp_C2D_sp_sp(sp2_I, sp4_I);
        SF0_C2D_X_X(1, 3) = SpDagSp_C2D_sp_sp(sp2_I, sp3_I);
        SF0_C2D_X_X(2, 3) = SpSp_C2D_sp_sp(sp4_I, sp3_I);
        SF0_C2D_X_X.triangularView<Eigen::StrictlyLower>() = -SF0_C2D_X_X.transpose();
        // T⁰ = pf(S_F⁰) = Sₚ₁ₚ₂Sₚ₃ₚ₄ − Sₚ₁ₚ₃Sₚ₂ₚ₄ + Sₚ₁ₚ₄Sₚ₂ₚ₃.
        const doubleC T0_C = calc_pfaffian(SF0_C2D_X_X);

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
                S_C2D_cqp12_X.topRows(Ncqp1_I) = S_C2D_compactqp12_X.topRows(Ncompactqp1_I)(compactqp1_I2D_cfg1_cqp1[cfg1_I], Eigen::placeholders::all).colwise().reverse();
                S_C2D_cqp12_X.bottomRows(Ncqp2_I) = S_C2D_compactqp12_X.bottomRows(Ncompactqp2_I)(compactqp2_I2D_cfg2_cqp2[cfg2_I], Eigen::placeholders::all);

                // pf([0,S;−Sᵀ,0]) = −det(S) = S₁₂S₂₁ − S₁₁S₂₂, S ∈ ℂ²ˣ².
                // −Σi<j F²ij pf([0,S;−Sᵀ,0]) = Σi<j F²ij det(S) = Σi<j F²ij(SᵢaSⱼb−SᵢbSⱼa) = (SᵀF²S)ab.
                const auto& F2_C2D_cqp12_cqp12 = F2_C4D_cfg1_cfg2_cqp12_cqp12[cfg1_I][cfg2_I];
                SF2S_C2D_X_X.noalias() = S_C2D_cqp12_X.transpose() * F2_C2D_cqp12_cqp12 * S_C2D_cqp12_X;
                doubleC TwoBody_C = F0_C2D_cfg1_cfg2[cfg1_I][cfg2_I] * T0_C
                    + SF2S_C2D_X_X(0, 1) * SF0_C2D_X_X(2, 3)
                    - SF2S_C2D_X_X(0, 2) * SF0_C2D_X_X(1, 3)
                    + SF2S_C2D_X_X(0, 3) * SF0_C2D_X_X(1, 2)
                    + SF2S_C2D_X_X(1, 2) * SF0_C2D_X_X(0, 3)
                    - SF2S_C2D_X_X(1, 3) * SF0_C2D_X_X(0, 2)
                    + SF2S_C2D_X_X(2, 3) * SF0_C2D_X_X(0, 1);

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
                    TwoBody_C += (S1S2_C1D_cqp12cqp12.transpose() * F4_C2D_cqp12cqp12_cqp12cqp12 * S3S4_C1D_cqp12cqp12).value();
                }
                result_C2D_cfg1_cfg2(cfg1_I, cfg2_I) += TwoBody_C;
            }
        }
    };

    // H₂ = Σβγδ (Σα TwoBody_αβγδ cα†) cβ†cδcγ.
    result_C2D_cfg1_cfg2.setZero();
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
    result_C2D_cfg1_cfg2 *= overlap_C;
}


inline void HFBPfaffian::calc_two_body(Eigen::Ref<Eigen::MatrixXcd> result_C2D_cfg1_cfg2, const Eigen::VectorXd& gQ_F1D_q, const Eigen::VectorXd& gP_F1D_p, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& Q_C3D_sp_sp_q, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& P_C3D_sp_sp_p) {
    assert(result_C2D_cfg1_cfg2.rows() == static_cast<Eigen::Index>(config1_I2D_cfg1_cqp1.size()) && result_C2D_cfg1_cfg2.cols() == static_cast<Eigen::Index>(config2_I2D_cfg2_cqp2.size()));
    assert(Q_C3D_sp_sp_q.dimension(2) == gQ_F1D_q.size());
    assert(P_C3D_sp_sp_p.dimension(2) == gP_F1D_p.size());
    assert(gQ_F1D_q.size() == 0 || (Q_C3D_sp_sp_q.dimension(0) == Nsp_I && Q_C3D_sp_sp_q.dimension(1) == Nsp_I));
    assert(gP_F1D_p.size() == 0 || (P_C3D_sp_sp_p.dimension(0) == Nsp_I && P_C3D_sp_sp_p.dimension(1) == Nsp_I));

    const Eigen::Index Ncompactqp1_I = qp1_I1D_compactqp1.size();
    const Eigen::Index Ncompactqp2_I = qp2_I1D_compactqp2.size();
    const Eigen::Index Ncompactqp12_I = Ncompactqp1_I + Ncompactqp2_I;

    // Internal contractions and antisymmetric pair factors.
    Eigen::MatrixXcd x1x2_C2D_sp_sp(Nsp_I, Nsp_I);
    Eigen::MatrixXcd x1x3_C2D_sp_sp(Nsp_I, Nsp_I);
    Eigen::MatrixXcd x1x4_C2D_sp_sp(Nsp_I, Nsp_I);
    Eigen::MatrixXcd x2x3_C2D_sp_sp(Nsp_I, Nsp_I);
    Eigen::MatrixXcd x2x4_C2D_sp_sp(Nsp_I, Nsp_I);
    Eigen::MatrixXcd x3x4_C2D_sp_sp(Nsp_I, Nsp_I);

    Eigen::MatrixXcd eta1_C2D_compactqp12_sp(Ncompactqp12_I, Nsp_I);
    Eigen::MatrixXcd eta2_C2D_compactqp12_sp(Ncompactqp12_I, Nsp_I);
    Eigen::MatrixXcd eta3_C2D_compactqp12_sp(Ncompactqp12_I, Nsp_I);
    Eigen::MatrixXcd eta4_C2D_compactqp12_sp(Ncompactqp12_I, Nsp_I);

    doubleC T0_C = 0.0;
    Eigen::MatrixXcd T2_C2D_compactqp12_compactqp12(Ncompactqp12_I, Ncompactqp12_I);
    Eigen::MatrixXcd PI14_C2D_compactqp12_compactqp12(Ncompactqp12_I, Ncompactqp12_I);
    Eigen::MatrixXcd PI23_C2D_compactqp12_compactqp12(Ncompactqp12_I, Ncompactqp12_I);
    Eigen::MatrixXcd PI12_C2D_compactqp12_compactqp12(Ncompactqp12_I, Ncompactqp12_I);
    Eigen::MatrixXcd PI34_C2D_compactqp12_compactqp12(Ncompactqp12_I, Ncompactqp12_I);

    const auto calc_constract_Q = [&](const Eigen::Ref<const Eigen::MatrixXcd>& Q_C2D_sp_sp) {
        // x₁x₂ = Q*κ̄Q; x₁x₃ = Q*ρ; x₂x₃ = Qᵀρ.
        x1x2_C2D_sp_sp.noalias() = Q_C2D_sp_sp.conjugate() * SpDagSpDag_C2D_sp_sp * Q_C2D_sp_sp;
        x1x3_C2D_sp_sp.noalias() = Q_C2D_sp_sp.conjugate() * SpDagSp_C2D_sp_sp;
        x2x3_C2D_sp_sp.noalias() = Q_C2D_sp_sp.transpose() * SpDagSp_C2D_sp_sp;

        // ηⁱ: left ⟨β₁xᵢ⟩; right −⟨xᵢβ₂†⟩ᵀ.
        eta1_C2D_compactqp12_sp.topRows(Ncompactqp1_I).noalias() = Qp1SpDag_C2D_qp1_sp(qp1_I1D_compactqp1, Eigen::placeholders::all) * Q_C2D_sp_sp.adjoint();
        eta1_C2D_compactqp12_sp.bottomRows(Ncompactqp2_I).noalias() = -SpDagQp2Dag_C2D_sp_qp2(Eigen::placeholders::all, qp2_I1D_compactqp2).transpose() * Q_C2D_sp_sp.adjoint();
        eta2_C2D_compactqp12_sp.topRows(Ncompactqp1_I).noalias() = Qp1SpDag_C2D_qp1_sp(qp1_I1D_compactqp1, Eigen::placeholders::all) * Q_C2D_sp_sp;
        eta2_C2D_compactqp12_sp.bottomRows(Ncompactqp2_I).noalias() = -SpDagQp2Dag_C2D_sp_qp2(Eigen::placeholders::all, qp2_I1D_compactqp2).transpose() * Q_C2D_sp_sp;
        eta3_C2D_compactqp12_sp.topRows(Ncompactqp1_I) = Qp1Sp_C2D_qp1_sp(qp1_I1D_compactqp1, Eigen::placeholders::all);
        eta3_C2D_compactqp12_sp.bottomRows(Ncompactqp2_I) = -SpQp2Dag_C2D_sp_qp2(Eigen::placeholders::all, qp2_I1D_compactqp2).transpose();
        eta4_C2D_compactqp12_sp = eta3_C2D_compactqp12_sp;

        // Reference PSM: direct contractions only.
        T0_C = 0.0;
        // T0_C += (x1x2_C2D_sp_sp.array() * SpSp_C2D_sp_sp.transpose().array()).sum();
        // T0_C -= (x1x3_C2D_sp_sp.array() * x2x3_C2D_sp_sp.transpose().array()).sum();
        T0_C += x1x3_C2D_sp_sp.trace() * x2x3_C2D_sp_sp.trace();

        // T²_Q = antisym[−⟨23⟩η¹(η⁴)ᵀ − ⟨14⟩η²(η³)ᵀ].
        T2_C2D_compactqp12_compactqp12.setZero();
        // T2_C2D_compactqp12_compactqp12.noalias() -= eta1_C2D_compactqp12_sp * SpSp_C2D_sp_sp.transpose() * eta2_C2D_compactqp12_sp.transpose();
        // T2_C2D_compactqp12_compactqp12.noalias() += eta1_C2D_compactqp12_sp * x2x3_C2D_sp_sp.transpose() * eta3_C2D_compactqp12_sp.transpose();
        T2_C2D_compactqp12_compactqp12.noalias() -= x2x3_C2D_sp_sp.trace() * eta1_C2D_compactqp12_sp * eta4_C2D_compactqp12_sp.transpose();
        T2_C2D_compactqp12_compactqp12.noalias() -= x1x3_C2D_sp_sp.trace() * eta2_C2D_compactqp12_sp * eta3_C2D_compactqp12_sp.transpose();
        // T2_C2D_compactqp12_compactqp12.noalias() += eta2_C2D_compactqp12_sp * x1x3_C2D_sp_sp.transpose() * eta4_C2D_compactqp12_sp.transpose();
        // T2_C2D_compactqp12_compactqp12.noalias() -= eta3_C2D_compactqp12_sp * x1x2_C2D_sp_sp.transpose() * eta4_C2D_compactqp12_sp.transpose();

        // T⁴ = Π¹⁴ ∧ Π²³.
        PI14_C2D_compactqp12_compactqp12.noalias() = eta1_C2D_compactqp12_sp * eta4_C2D_compactqp12_sp.transpose();
        PI23_C2D_compactqp12_compactqp12.noalias() = eta2_C2D_compactqp12_sp * eta3_C2D_compactqp12_sp.transpose();
        PI14_C2D_compactqp12_compactqp12 = (PI14_C2D_compactqp12_compactqp12 - PI14_C2D_compactqp12_compactqp12.transpose()).eval();
        PI23_C2D_compactqp12_compactqp12 = (PI23_C2D_compactqp12_compactqp12 - PI23_C2D_compactqp12_compactqp12.transpose()).eval();
        T2_C2D_compactqp12_compactqp12 = (T2_C2D_compactqp12_compactqp12 - T2_C2D_compactqp12_compactqp12.transpose()).eval();
    };

    const auto calc_constract_P = [&](const Eigen::Ref<const Eigen::MatrixXcd>& P_C2D_sp_sp) {
        // tr(Pᵀκ̄) = Σαβ Pαβ κ̄αβ; tr(κP*) = Σαβ καβ P*βα.
        const doubleC trace12_C = P_C2D_sp_sp.cwiseProduct(SpDagSpDag_C2D_sp_sp).sum();
        const doubleC trace34_C = SpSp_C2D_sp_sp.cwiseProduct(P_C2D_sp_sp.adjoint()).sum();

        // x₁x₃ = Pᵀρ; x₁x₄ = PᵀρP*.
        x1x3_C2D_sp_sp.noalias() = P_C2D_sp_sp.transpose() * SpDagSp_C2D_sp_sp;
        x1x4_C2D_sp_sp.noalias() = x1x3_C2D_sp_sp * P_C2D_sp_sp.conjugate();
        x2x3_C2D_sp_sp = SpDagSp_C2D_sp_sp;
        x2x4_C2D_sp_sp.noalias() = SpDagSp_C2D_sp_sp * P_C2D_sp_sp.conjugate();

        // ηⁱ: left ⟨β₁xᵢ⟩; right −⟨xᵢβ₂†⟩ᵀ.
        eta1_C2D_compactqp12_sp.topRows(Ncompactqp1_I).noalias() = Qp1SpDag_C2D_qp1_sp(qp1_I1D_compactqp1, Eigen::placeholders::all) * P_C2D_sp_sp;
        eta1_C2D_compactqp12_sp.bottomRows(Ncompactqp2_I).noalias() = -SpDagQp2Dag_C2D_sp_qp2(Eigen::placeholders::all, qp2_I1D_compactqp2).transpose() * P_C2D_sp_sp;
        eta2_C2D_compactqp12_sp.topRows(Ncompactqp1_I) = Qp1SpDag_C2D_qp1_sp(qp1_I1D_compactqp1, Eigen::placeholders::all);
        eta2_C2D_compactqp12_sp.bottomRows(Ncompactqp2_I) = -SpDagQp2Dag_C2D_sp_qp2(Eigen::placeholders::all, qp2_I1D_compactqp2).transpose();
        eta3_C2D_compactqp12_sp.topRows(Ncompactqp1_I) = Qp1Sp_C2D_qp1_sp(qp1_I1D_compactqp1, Eigen::placeholders::all);
        eta3_C2D_compactqp12_sp.bottomRows(Ncompactqp2_I) = -SpQp2Dag_C2D_sp_qp2(Eigen::placeholders::all, qp2_I1D_compactqp2).transpose();
        eta4_C2D_compactqp12_sp.topRows(Ncompactqp1_I).noalias() = Qp1Sp_C2D_qp1_sp(qp1_I1D_compactqp1, Eigen::placeholders::all) * P_C2D_sp_sp.conjugate();
        eta4_C2D_compactqp12_sp.bottomRows(Ncompactqp2_I).noalias() = -SpQp2Dag_C2D_sp_qp2(Eigen::placeholders::all, qp2_I1D_compactqp2).transpose() * P_C2D_sp_sp.conjugate();

        // Reference PSM: direct contractions only.
        T0_C = 0.0;
        T0_C += trace12_C * trace34_C;
        // T0_C -= (x1x3_C2D_sp_sp.array() * x2x4_C2D_sp_sp.array()).sum();
        // T0_C += (x1x4_C2D_sp_sp.array() * x2x3_C2D_sp_sp.array()).sum();

        // T²_P = antisym[−⟨34⟩η¹(η²)ᵀ − ⟨12⟩η³(η⁴)ᵀ].
        T2_C2D_compactqp12_compactqp12.setZero();
        T2_C2D_compactqp12_compactqp12.noalias() -= trace34_C * eta1_C2D_compactqp12_sp * eta2_C2D_compactqp12_sp.transpose();
        // T2_C2D_compactqp12_compactqp12.noalias() += eta1_C2D_compactqp12_sp * x2x4_C2D_sp_sp * eta3_C2D_compactqp12_sp.transpose();
        // T2_C2D_compactqp12_compactqp12.noalias() -= eta1_C2D_compactqp12_sp * x2x3_C2D_sp_sp * eta4_C2D_compactqp12_sp.transpose();
        // T2_C2D_compactqp12_compactqp12.noalias() -= eta2_C2D_compactqp12_sp * x1x4_C2D_sp_sp * eta3_C2D_compactqp12_sp.transpose();
        // T2_C2D_compactqp12_compactqp12.noalias() += eta2_C2D_compactqp12_sp * x1x3_C2D_sp_sp * eta4_C2D_compactqp12_sp.transpose();
        T2_C2D_compactqp12_compactqp12.noalias() -= trace12_C * eta3_C2D_compactqp12_sp * eta4_C2D_compactqp12_sp.transpose();

        T2_C2D_compactqp12_compactqp12 = (T2_C2D_compactqp12_compactqp12 - T2_C2D_compactqp12_compactqp12.transpose()).eval();

        // T⁴ = Π¹² ∧ Π³⁴.
        PI12_C2D_compactqp12_compactqp12.noalias() = eta1_C2D_compactqp12_sp * eta2_C2D_compactqp12_sp.transpose();
        PI34_C2D_compactqp12_compactqp12.noalias() = eta3_C2D_compactqp12_sp * eta4_C2D_compactqp12_sp.transpose();
        PI12_C2D_compactqp12_compactqp12 = (PI12_C2D_compactqp12_compactqp12 - PI12_C2D_compactqp12_compactqp12.transpose()).eval();
        PI34_C2D_compactqp12_compactqp12 = (PI34_C2D_compactqp12_compactqp12 - PI34_C2D_compactqp12_compactqp12.transpose()).eval();
    };

    // Ncqp12 ≤ Ncompactqp12.
    std::vector<int> compactqp12_I1D_cqp12{};
    compactqp12_I1D_cqp12.reserve(Ncompactqp12_I);
    Eigen::VectorXcd PI1workspace_C1D_element(4 * Nsp_I * Nsp_I);
    Eigen::VectorXcd PI2workspace_C1D_element(4 * Nsp_I * Nsp_I);

    const auto calc_kernel024 = [&](double g_F, const Eigen::MatrixXcd& PI1_C2D_compactqp12_compactqp12, const Eigen::MatrixXcd& PI2_C2D_compactqp12_compactqp12) {
        const doubleC factor_C = overlap_C * g_F;
        for (int cfg2_I = 0; cfg2_I < static_cast<int>(compactqp2_I2D_cfg2_cqp2.size()); ++cfg2_I) {
            const auto& compactqp2_I1D_cqp2 = compactqp2_I2D_cfg2_cqp2[cfg2_I];
            for (int cfg1_I = 0; cfg1_I < static_cast<int>(compactqp1_I2D_cfg1_cqp1.size()); ++cfg1_I) {
                const auto& compactqp1_I1D_cqp1 = compactqp1_I2D_cfg1_cqp1[cfg1_I];
                const int Ncqp1_I = static_cast<int>(compactqp1_I1D_cqp1.size());
                const Eigen::Index Ncqp12_I = Ncqp1_I + static_cast<Eigen::Index>(compactqp2_I1D_cqp2.size());
                if (Ncqp12_I % 2 != 0) { continue; }
                // Left reversed; right shifted by Ncompactqp1.
                compactqp12_I1D_cqp12.clear();
                for (int cqp1_I = Ncqp1_I - 1; cqp1_I >= 0; --cqp1_I) { compactqp12_I1D_cqp12.push_back(compactqp1_I1D_cqp1[cqp1_I]); }
                for (const int compactqp2_I : compactqp2_I1D_cqp2) { compactqp12_I1D_cqp12.push_back(Ncompactqp1_I + compactqp2_I); }
                const auto& F2_C2D_cqp12_cqp12 = F2_C4D_cfg1_cfg2_cqp12_cqp12[cfg1_I][cfg2_I];
                // K⁰ + K² = F⁰T⁰ − ½Σij F²ij T²ij.
                doubleC kernel_C = F0_C2D_cfg1_cfg2[cfg1_I][cfg2_I] * T0_C
                    - 0.5 * F2_C2D_cqp12_cqp12.cwiseProduct(T2_C2D_compactqp12_compactqp12(compactqp12_I1D_cqp12, compactqp12_I1D_cqp12)).sum();

                // K⁴ = ¼ vec(Π₁)ᵀ F⁴_(ij,kl) vec(Π₂).
                if (Ncqp12_I >= 4) {
                    Eigen::Map<Eigen::MatrixXcd> PI1_C2D_cqp12_cqp12(PI1workspace_C1D_element.data(), Ncqp12_I, Ncqp12_I);
                    Eigen::Map<Eigen::MatrixXcd> PI2_C2D_cqp12_cqp12(PI2workspace_C1D_element.data(), Ncqp12_I, Ncqp12_I);
                    PI1_C2D_cqp12_cqp12 = PI1_C2D_compactqp12_compactqp12(compactqp12_I1D_cqp12, compactqp12_I1D_cqp12);
                    PI2_C2D_cqp12_cqp12 = PI2_C2D_compactqp12_compactqp12(compactqp12_I1D_cqp12, compactqp12_I1D_cqp12);
                    // ColMajor: (i,j) → i+nj; (k,l) → k+nl.
                    const Eigen::Map<const Eigen::MatrixXcd> F4_C2D_cqp12cqp12_cqp12cqp12(F4_C6D_cfg1_cfg2_cqp12_cqp12_cqp12_cqp12[cfg1_I][cfg2_I].data(), Ncqp12_I * Ncqp12_I, Ncqp12_I * Ncqp12_I);
                    const Eigen::Map<const Eigen::VectorXcd> PI1_C1D_cqp12cqp12(PI1_C2D_cqp12_cqp12.data(), PI1_C2D_cqp12_cqp12.size());
                    const Eigen::Map<const Eigen::VectorXcd> PI2_C1D_cqp12cqp12(PI2_C2D_cqp12_cqp12.data(), PI2_C2D_cqp12_cqp12.size());
                    kernel_C += 0.25 * (PI1_C1D_cqp12cqp12.transpose() * F4_C2D_cqp12cqp12_cqp12cqp12 * PI2_C1D_cqp12cqp12).value();
                }
                result_C2D_cfg1_cfg2(cfg1_I, cfg2_I) += factor_C * kernel_C;
            }
        }
    };

    result_C2D_cfg1_cfg2.setZero();
    // Q: gQ (F⁰T⁰ − ½F²T² + ¼F⁴Π¹⁴Π²³).
    for (Eigen::Index iQ_I = 0; iQ_I < gQ_F1D_q.size(); ++iQ_I) {
        if (gQ_F1D_q(iQ_I) == 0.0) { continue; }
        const Eigen::Index offset_I = iQ_I * Nsp_I * Nsp_I;
        const Eigen::Map<const Eigen::MatrixXcd> Q_C2D_sp_sp(Q_C3D_sp_sp_q.data() + offset_I, Nsp_I, Nsp_I);
        calc_constract_Q(Q_C2D_sp_sp);
        calc_kernel024(gQ_F1D_q(iQ_I), PI14_C2D_compactqp12_compactqp12, PI23_C2D_compactqp12_compactqp12);
    }

    // P: ¼gP (F⁰T⁰ − ½F²T² + ¼F⁴Π¹²Π³⁴).
    for (Eigen::Index iP_I = 0; iP_I < gP_F1D_p.size(); ++iP_I) {
        if (gP_F1D_p(iP_I) == 0.0) { continue; }
        const Eigen::Index offset_I = iP_I * Nsp_I * Nsp_I;
        const Eigen::Map<const Eigen::MatrixXcd> P_C2D_sp_sp(P_C3D_sp_sp_p.data() + offset_I, Nsp_I, Nsp_I);
        calc_constract_P(P_C2D_sp_sp);
        calc_kernel024(0.25 * gP_F1D_p(iP_I), PI12_C2D_compactqp12_compactqp12, PI34_C2D_compactqp12_compactqp12);
    }

    // h₀Q = Σᵢ gQi Qi†Qi.
    // Eigen::MatrixXcd h0Q_C2D_sp_sp = Eigen::MatrixXcd::Zero(Nsp_I, Nsp_I);
    // for (Eigen::Index iQ_I = 0; iQ_I < gQ_F1D_q.size(); ++iQ_I) {
    //     const Eigen::Index offset_I = iQ_I * Nsp_I * Nsp_I;
    //     const Eigen::Map<const Eigen::MatrixXcd> Q_C2D_sp_sp(Q_C3D_sp_sp_q.data() + offset_I, Nsp_I, Nsp_I);
    //     h0Q_C2D_sp_sp.noalias() += gQ_F1D_q(iQ_I) * Q_C2D_sp_sp.adjoint() * Q_C2D_sp_sp;
    // }
    // result ← result + ⟨Φ₁;cfg1|h₀Q|Φ₂;cfg2⟩.
    // Eigen::MatrixXcd h0Q_C2D_cfg1_cfg2(result_C2D_cfg1_cfg2.rows(), result_C2D_cfg1_cfg2.cols());
    // calc_one_body(h0Q_C2D_cfg1_cfg2, h0Q_C2D_sp_sp);
    // result_C2D_cfg1_cfg2 += h0Q_C2D_cfg1_cfg2;
}

inline doubleC HFBPfaffian::calc_pfaffian(Eigen::Ref<Eigen::MatrixXcd> S_C2D_chain_chain) {
    const Eigen::Index Nchain_I = S_C2D_chain_chain.rows();
    assert(S_C2D_chain_chain.cols() == Nchain_I && Nchain_I % 2 == 0);
    assert(S_C2D_chain_chain.allFinite());
    assert(S_C2D_chain_chain.isApprox(-S_C2D_chain_chain.transpose(), 1.0e-12));
    // pf(∅) = 1; pf(S₂) = S₀₁; pf(S₄) = S₀₁S₂₃ − S₀₂S₁₃ + S₀₃S₁₂.
    if (Nchain_I == 0) { return doubleC(1.0, 0.0); }
    if (Nchain_I == 2) { return S_C2D_chain_chain(0, 1); }
    if (Nchain_I == 4) {
        return S_C2D_chain_chain(0, 1) * S_C2D_chain_chain(2, 3)
             - S_C2D_chain_chain(0, 2) * S_C2D_chain_chain(1, 3)
             + S_C2D_chain_chain(0, 3) * S_C2D_chain_chain(1, 2);
    }
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

        // Sᵢⱼ ← Sᵢⱼ + (Sₖ₊₁,ᵢSₖⱼ − SₖᵢSₖ₊₁,ⱼ)/pivot; Sⱼᵢ = −Sᵢⱼ.
        for (Eigen::Index j_I = k_I + 3; j_I < Nchain_I; ++j_I) {
            for (Eigen::Index i_I = k_I + 2; i_I < j_I; ++i_I) {
                S_C2D_chain_chain(i_I, j_I) += (S_C2D_chain_chain(k_I + 1, i_I) * S_C2D_chain_chain(k_I, j_I) - S_C2D_chain_chain(k_I, i_I) * S_C2D_chain_chain(k_I + 1, j_I)) / pivot_C;
                S_C2D_chain_chain(j_I, i_I) = -S_C2D_chain_chain(i_I, j_I);
            }
        }
    }
    return pf_C;
}
