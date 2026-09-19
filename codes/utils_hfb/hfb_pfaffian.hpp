/**
 * @file    hfb_pfaffian.hpp
 * @author  cailea
 * @date    2026-09-10
 * @brief   HFB contractions and multiquasiparticle Pfaffian kernels.
 */

#pragma once

#include <cassert>
#include <cmath>
#include <complex>
#include <vector>

#include <Eigen/Core>
#include <Eigen/LU>

using doubleC = std::complex<double>;

/**
 * @brief  Evaluate HFB configuration kernels using Pfaffians.
 * @note   Configuration kernels require current contractions and nonzero vacuum overlap.
 * @note   Kernels include x coefficients and vacuum overlap; exclude prefactors.
 * @note   Returned references alias reusable outputs; indices are zero-based.
 */
class HFBPfaffian {
public:
    int Nsp_I = 0;

    Eigen::MatrixXcd U1_C2D_sp_qp1{};
    Eigen::MatrixXcd V1_C2D_sp_qp1{};
    Eigen::MatrixXcd U2_C2D_sp_qp2{};
    Eigen::MatrixXcd V2_C2D_sp_qp2{};
    Eigen::MatrixXcd Z1_C2D_sp_sp{};
    Eigen::MatrixXcd Z2_C2D_sp_sp{};

    Eigen::MatrixXcd A_C2D_qp1_qp2{};
    Eigen::MatrixXcd AInv_C2D_qp2_qp1{};

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

    Eigen::MatrixXcd overlap_C2D_cfg1_cfg2{};
    Eigen::MatrixXcd creator_C2D_cfg1_cfg2{};
    Eigen::MatrixXcd annihilator_C2D_cfg1_cfg2{};
    Eigen::MatrixXcd OBTD_C2D_cfg1_cfg2{};
    Eigen::MatrixXcd TBTD_C2D_cfg1_cfg2{};

private:
    Eigen::FullPivLU<Eigen::MatrixXcd> U1_lu{};
    Eigen::FullPivLU<Eigen::MatrixXcd> U2_lu{};
    Eigen::FullPivLU<Eigen::MatrixXcd> A_lu{};
    Eigen::VectorXcd Sworkspace_C1D_element{};
    Eigen::VectorXcd Qp1Xworkspace_C1D_element{};
    Eigen::VectorXcd XQp2Dagworkspace_C1D_element{};
    Eigen::VectorXcd XXworkspace_C1D_element{};

public:
    /**
     * @brief  Store explicit configurations and allocate HFB workspaces.
     * @math   U₁,V₁,U₂,V₂ ∈ ℂ^{Nsp×Nsp}.
     * @output Stored both configuration tables; allocated matrices and workspaces.
     * @note   Indices increase strictly; each lies in [0,Nsp).
     * @note   An empty inner list denotes the vacuum.
     */
    HFBPfaffian(int Nsp_I_, const std::vector<std::vector<int>>& config1_I2D_cfg1_cqp1_, const std::vector<std::vector<int>>& config2_I2D_cfg2_cqp2_) {
        assert(Nsp_I_ > 0);

        Nsp_I = Nsp_I_;
        config1_I2D_cfg1_cqp1 = config1_I2D_cfg1_cqp1_;
        config2_I2D_cfg2_cqp2 = config2_I2D_cfg2_cqp2_;

        // U₁,V₁,U₂,V₂ ∈ ℂ^{Nsp×Nsp}.
        U1_C2D_sp_qp1.resize(Nsp_I, Nsp_I);
        V1_C2D_sp_qp1.resize(Nsp_I, Nsp_I);
        U2_C2D_sp_qp2.resize(Nsp_I, Nsp_I);
        V2_C2D_sp_qp2.resize(Nsp_I, Nsp_I);

        // Z₁ = (V₁U₁⁻¹)*; Z₂ = (V₂U₂⁻¹)*.
        U1_lu = Eigen::FullPivLU<Eigen::MatrixXcd>(Nsp_I, Nsp_I);
        U2_lu = Eigen::FullPivLU<Eigen::MatrixXcd>(Nsp_I, Nsp_I);
        Z1_C2D_sp_sp.resize(Nsp_I, Nsp_I);
        Z2_C2D_sp_sp.resize(Nsp_I, Nsp_I);

        // A,A⁻¹ ∈ ℂ^{Nsp×Nsp}.
        A_C2D_qp1_qp2.resize(Nsp_I, Nsp_I);
        AInv_C2D_qp2_qp1.resize(Nsp_I, Nsp_I);
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
        overlap_C2D_cfg1_cfg2.resize(config1_I2D_cfg1_cqp1.size(), config2_I2D_cfg2_cqp2.size());
        creator_C2D_cfg1_cfg2.resize(config1_I2D_cfg1_cqp1.size(), config2_I2D_cfg2_cqp2.size());
        annihilator_C2D_cfg1_cfg2.resize(config1_I2D_cfg1_cqp1.size(), config2_I2D_cfg2_cqp2.size());
        OBTD_C2D_cfg1_cfg2.resize(config1_I2D_cfg1_cqp1.size(), config2_I2D_cfg2_cqp2.size());
        TBTD_C2D_cfg1_cfg2.resize(config1_I2D_cfg1_cqp1.size(), config2_I2D_cfg2_cqp2.size());

        // Lmax = 2Nsp+4; size(Sworkspace) = Lmax².
        const int Lmax_I = 2 * Nsp_I + 4;
        Sworkspace_C1D_element.resize(Lmax_I * Lmax_I);
        // NX ≤ 4; buffers are mapped to current insertion dimensions.
        Qp1Xworkspace_C1D_element.resize(4 * Nsp_I);
        XQp2Dagworkspace_C1D_element.resize(4 * Nsp_I);
        XXworkspace_C1D_element.resize(16);
    }

    /**
     * @brief  Update overlap and contractions using Pfaffians and full-pivot LU.
     * @math   A = U₁ᵀU₂* + V₁ᵀV₂*; Qp1Qp2Dag = A⁻ᵀ.
     * @output Updated U,V,Z, overlap, A, AInv, and eleven contractions.
     * @note   Requires canonical U,V and numerically invertible U₁,U₂,A.
     * @note   Positive reference phases: ν_a = √|det U_a|.
     * @note   phase_C = exp[i(θ₂-θ₁)] multiplies the vacuum overlap.
     */
    void update_contractions(const Eigen::MatrixXcd& U1_C2D_sp_qp1_, const Eigen::MatrixXcd& V1_C2D_sp_qp1_, const Eigen::MatrixXcd& U2_C2D_sp_qp2_, const Eigen::MatrixXcd& V2_C2D_sp_qp2_, doubleC phase_C);

    /**
     * @brief  Assemble configuration kernels using Wick Pfaffians.
     * @math   Kκ₁κ₂ = overlap × pf(S[β₁, X, β₂†]).
     * @output Filled preallocated result_C2D_cfg1_cfg2.
     * @note   X follows the physical operator order; NX = 0,1,2,4.
     * @note   XX is fully initialized, with XXᵀ = −XX.
     * @note   Inputs must not alias output or Sworkspace.
     */
    void calc_kernel(Eigen::Ref<const Eigen::MatrixXcd> Qp1X_C2D_qp1_X, Eigen::Ref<const Eigen::MatrixXcd> XQp2Dag_C2D_X_qp2, Eigen::Ref<const Eigen::MatrixXcd> XX_C2D_X_X, Eigen::MatrixXcd& result_C2D_cfg1_cfg2);

    /**
     * @brief  Calculate configuration overlaps using Wick Pfaffians.
     * @math   N_{κ₁κ₂} = ⟨Φ₁;κ₁|Φ₂;κ₂⟩.
     * @output Updated overlap_C2D_cfg1_cfg2 and its const reference.
     */
    const Eigen::MatrixXcd& calc_overlap();

    /**
     * @brief  Calculate creation kernels using Wick Pfaffians.
     * @math   C_{κ₁κ₂}[x] = ⟨Φ₁;κ₁|Σα xα c†α|Φ₂;κ₂⟩.
     * @note   Coefficients enter linearly, without conjugation.
     * @output Updated creator_C2D_cfg1_cfg2 and its const reference.
     * @note   Nonzero only for odd combined quasiparticle counts.
     */
    const Eigen::MatrixXcd& calc_creator(const Eigen::VectorXcd& x_C1D_sp);

    /**
     * @brief  Calculate annihilation kernels using Wick Pfaffians.
     * @math   A_{κ₁κ₂}[x] = ⟨Φ₁;κ₁|Σα xα cα|Φ₂;κ₂⟩.
     * @note   Coefficients enter linearly, without conjugation.
     * @output Updated annihilator_C2D_cfg1_cfg2 and its const reference.
     * @note   Nonzero only for odd combined quasiparticle counts.
     */
    const Eigen::MatrixXcd& calc_annihilator(const Eigen::VectorXcd& x_C1D_sp);

    /**
     * @brief  Calculate one-body kernels using Wick Pfaffians.
     * @math   OBTD_{κ₁κ₂}[x₁,x₂] = ⟨Φ₁;κ₁|(Σα x₁α c†α)(Σβ x₂β cβ)|Φ₂;κ₂⟩.
     * @note   Both coefficient vectors enter without conjugation.
     * @output Updated OBTD_C2D_cfg1_cfg2 and its const reference.
     */
    const Eigen::MatrixXcd& calc_obtd(const Eigen::VectorXcd& x1_C1D_sp, const Eigen::VectorXcd& x2_C1D_sp);

    /**
     * @brief  Calculate two-body kernels using Wick Pfaffians.
     * @math   TBTD_{κ₁κ₂}[x₁,x₂,x₃,x₄] = ⟨Φ₁;κ₁|x₁†x₂†x₄x₃|Φ₂;κ₂⟩.
     * @note   x₁†,x₂† = Σα x₁α c†α,Σα x₂α c†α.
     * @note   x₃,x₄ = Σα x₃α cα,Σα x₄α cα.
     * @note   All coefficient vectors enter without conjugation.
     * @output Updated TBTD_C2D_cfg1_cfg2 and its const reference.
     */
    const Eigen::MatrixXcd& calc_tbtd(const Eigen::VectorXcd& x1_C1D_sp, const Eigen::VectorXcd& x2_C1D_sp, const Eigen::VectorXcd& x3_C1D_sp, const Eigen::VectorXcd& x4_C1D_sp);

    /**
     * @brief  Calculate Pfaffians using pivoted skew-symmetric elimination.
     * @math   pf(S)² = det(S); pf(∅) = 1.
     * @output Pfaffian value; input matrix overwritten.
     * @note   Requires finite, even-order square S with Sᵀ = -S.
     * @note   Exact-zero pivots return zero; no magnitude cutoff.
     */
    static doubleC calc_pfaffian(Eigen::Ref<Eigen::MatrixXcd> S_C2D_chain_chain);

};


inline void HFBPfaffian::update_contractions(const Eigen::MatrixXcd& U1_C2D_sp_qp1_, const Eigen::MatrixXcd& V1_C2D_sp_qp1_, const Eigen::MatrixXcd& U2_C2D_sp_qp2_, const Eigen::MatrixXcd& V2_C2D_sp_qp2_, doubleC phase_C) {
    // U₁,V₁,U₂,V₂ ∈ ℂ^{Nsp×Nsp}.
    assert(U1_C2D_sp_qp1_.rows() == Nsp_I && U1_C2D_sp_qp1_.cols() == Nsp_I);
    assert(V1_C2D_sp_qp1_.rows() == Nsp_I && V1_C2D_sp_qp1_.cols() == Nsp_I);
    assert(U2_C2D_sp_qp2_.rows() == Nsp_I && U2_C2D_sp_qp2_.cols() == Nsp_I);
    assert(V2_C2D_sp_qp2_.rows() == Nsp_I && V2_C2D_sp_qp2_.cols() == Nsp_I);
    U1_C2D_sp_qp1 = U1_C2D_sp_qp1_;
    V1_C2D_sp_qp1 = V1_C2D_sp_qp1_;
    U2_C2D_sp_qp2 = U2_C2D_sp_qp2_;
    V2_C2D_sp_qp2 = V2_C2D_sp_qp2_;

    // UᵀZ† = Vᵀ; ν = √|det U| > 0.
    assert(U1_C2D_sp_qp1.allFinite() && V1_C2D_sp_qp1.allFinite());
    assert(U2_C2D_sp_qp2.allFinite() && V2_C2D_sp_qp2.allFinite());
    U1_lu.compute(U1_C2D_sp_qp1);
    U2_lu.compute(U2_C2D_sp_qp2);
    assert(U1_lu.isInvertible() && U2_lu.isInvertible());
    Z1_C2D_sp_sp = U1_lu.transpose().solve(V1_C2D_sp_qp1.transpose());
    Z2_C2D_sp_sp = U2_lu.transpose().solve(V2_C2D_sp_qp2.transpose());
    Z1_C2D_sp_sp.adjointInPlace();
    Z2_C2D_sp_sp.adjointInPlace();
    const double nu1_F = std::sqrt(std::abs(U1_lu.determinant()));
    const double nu2_F = std::sqrt(std::abs(U2_lu.determinant()));

    // X = [Z₂,-I; I,-Z₁*]; overlap = phase × ν₁ν₂s_Nsp pf(X).
    Eigen::Map<Eigen::MatrixXcd> X_C2D_sp_sp(Sworkspace_C1D_element.data(), 2 * Nsp_I, 2 * Nsp_I);
    X_C2D_sp_sp.topLeftCorner(Nsp_I, Nsp_I) = Z2_C2D_sp_sp;
    X_C2D_sp_sp.bottomRightCorner(Nsp_I, Nsp_I) = -Z1_C2D_sp_sp.conjugate();
    X_C2D_sp_sp.topRightCorner(Nsp_I, Nsp_I).setIdentity();
    X_C2D_sp_sp.topRightCorner(Nsp_I, Nsp_I) *= -1.0;
    X_C2D_sp_sp.bottomLeftCorner(Nsp_I, Nsp_I).setIdentity();
    const int sign_I = (Nsp_I % 4 == 0 || Nsp_I % 4 == 3) ? 1 : -1;
    overlap_C = phase_C * nu1_F * nu2_F * static_cast<double>(sign_I) * calc_pfaffian(X_C2D_sp_sp);

    // A = U₁ᵀU₂* + V₁ᵀV₂*.
    A_C2D_qp1_qp2.noalias() = U1_C2D_sp_qp1.transpose() * U2_C2D_sp_qp2.conjugate();
    A_C2D_qp1_qp2.noalias() += V1_C2D_sp_qp1.transpose() * V2_C2D_sp_qp2.conjugate();
    assert(A_C2D_qp1_qp2.allFinite());
    A_lu.compute(A_C2D_qp1_qp2);
    assert(A_lu.isInvertible());
    AInv_C2D_qp2_qp1 = A_lu.inverse();
    Qp1Qp2Dag_C2D_qp1_qp2 = AInv_C2D_qp2_qp1.transpose();

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
}

inline void HFBPfaffian::calc_kernel(Eigen::Ref<const Eigen::MatrixXcd> Qp1X_C2D_qp1_X, Eigen::Ref<const Eigen::MatrixXcd> XQp2Dag_C2D_X_qp2, Eigen::Ref<const Eigen::MatrixXcd> XX_C2D_X_X, Eigen::MatrixXcd& result_C2D_cfg1_cfg2) {
    assert(XX_C2D_X_X.rows() == 0 || XX_C2D_X_X.rows() == 1 || XX_C2D_X_X.rows() == 2 || XX_C2D_X_X.rows() == 4);
    const int NX_I = static_cast<int>(XX_C2D_X_X.rows());
    assert(XX_C2D_X_X.cols() == NX_I);
    assert(Qp1X_C2D_qp1_X.rows() == Nsp_I && Qp1X_C2D_qp1_X.cols() == NX_I && Qp1X_C2D_qp1_X.allFinite());
    assert(XQp2Dag_C2D_X_qp2.rows() == NX_I && XQp2Dag_C2D_X_qp2.cols() == Nsp_I && XQp2Dag_C2D_X_qp2.allFinite());
    assert(result_C2D_cfg1_cfg2.rows() == static_cast<Eigen::Index>(config1_I2D_cfg1_cqp1.size()) && result_C2D_cfg1_cfg2.cols() == static_cast<Eigen::Index>(config2_I2D_cfg2_cqp2.size()));
    assert(overlap_C != doubleC(0.0, 0.0));

    // (β₁,μr₁,…,β₁,μ₁, X₁,…,X_NX, β₂,ν₁†,…,β₂,νr₂†).
    for (int cfg2_I = 0; cfg2_I < static_cast<int>(config2_I2D_cfg2_cqp2.size()); ++cfg2_I) {
        const int Ncqp2_I = static_cast<int>(config2_I2D_cfg2_cqp2[cfg2_I].size());
        for (int cfg1_I = 0; cfg1_I < static_cast<int>(config1_I2D_cfg1_cqp1.size()); ++cfg1_I) {
            const int Ncqp1_I = static_cast<int>(config1_I2D_cfg1_cqp1[cfg1_I].size());
            const int Nchain_I = Ncqp1_I + NX_I + Ncqp2_I;

            // Odd chains vanish; pf(∅) = 1.
            if (Nchain_I % 2 != 0) {
                result_C2D_cfg1_cfg2(cfg1_I, cfg2_I) = 0.0;
                continue;
            }

            // S ∈ ℂ^{(r₁+NX+r₂)×(r₁+NX+r₂)}; S_ii = 0.
            assert(static_cast<Eigen::Index>(Nchain_I) * Nchain_I <= Sworkspace_C1D_element.size());
            Eigen::Map<Eigen::MatrixXcd> S_C2D_chain_chain(Sworkspace_C1D_element.data(), Nchain_I, Nchain_I);
            for (int chain1_I = 0; chain1_I < Nchain_I; ++chain1_I) {S_C2D_chain_chain(chain1_I, chain1_I) = 0.0;}

            // ⟨β₁β₁⟩: reversed left configuration.
            for (int chain2_I = 1; chain2_I < Ncqp1_I; ++chain2_I) {
                for (int chain1_I = 0; chain1_I < chain2_I; ++chain1_I) {
                    S_C2D_chain_chain(chain1_I, chain2_I) = Qp1Qp1_C2D_qp1_qp1(config1_I2D_cfg1_cqp1[cfg1_I][Ncqp1_I - 1 - chain1_I], config1_I2D_cfg1_cqp1[cfg1_I][Ncqp1_I - 1 - chain2_I]);
                }
            }

            // ⟨β₁X⟩ and ⟨Xβ₂†⟩ preserve insertion order.
            for (int X_I = 0; X_I < NX_I; ++X_I) {
                for (int chain1_I = 0; chain1_I < Ncqp1_I; ++chain1_I) {
                    S_C2D_chain_chain(chain1_I, Ncqp1_I + X_I) = Qp1X_C2D_qp1_X(config1_I2D_cfg1_cqp1[cfg1_I][Ncqp1_I - 1 - chain1_I], X_I);
                }
                for (int chain2_I = 0; chain2_I < Ncqp2_I; ++chain2_I) {
                    S_C2D_chain_chain(Ncqp1_I + X_I, Ncqp1_I + NX_I + chain2_I) = XQp2Dag_C2D_X_qp2(X_I, config2_I2D_cfg2_cqp2[cfg2_I][chain2_I]);
                }
            }

            // ⟨X_i X_j⟩, i<j; lower triangle is unused.
            for (int chain2_I = 1; chain2_I < NX_I; ++chain2_I) {
                for (int chain1_I = 0; chain1_I < chain2_I; ++chain1_I) {
                    S_C2D_chain_chain(Ncqp1_I + chain1_I, Ncqp1_I + chain2_I) = XX_C2D_X_X(chain1_I, chain2_I);
                }
            }

            // ⟨β₁β₂†⟩: upper-right block.
            for (int chain2_I = 0; chain2_I < Ncqp2_I; ++chain2_I) {
                for (int chain1_I = 0; chain1_I < Ncqp1_I; ++chain1_I) {
                    S_C2D_chain_chain(chain1_I, Ncqp1_I + NX_I + chain2_I) = Qp1Qp2Dag_C2D_qp1_qp2(config1_I2D_cfg1_cqp1[cfg1_I][Ncqp1_I - 1 - chain1_I], config2_I2D_cfg2_cqp2[cfg2_I][chain2_I]);
                }
            }

            // ⟨β₂†β₂†⟩: forward right configuration.
            for (int chain2_I = 1; chain2_I < Ncqp2_I; ++chain2_I) {
                for (int chain1_I = 0; chain1_I < chain2_I; ++chain1_I) {
                    S_C2D_chain_chain(Ncqp1_I + NX_I + chain1_I, Ncqp1_I + NX_I + chain2_I) = Qp2DagQp2Dag_C2D_qp2_qp2(config2_I2D_cfg2_cqp2[cfg2_I][chain1_I], config2_I2D_cfg2_cqp2[cfg2_I][chain2_I]);
                }
            }

            // S_ji = −S_ij; vacuum overlap enters exactly once.
            for (int chain2_I = 1; chain2_I < Nchain_I; ++chain2_I) {
                for (int chain1_I = 0; chain1_I < chain2_I; ++chain1_I) {S_C2D_chain_chain(chain2_I, chain1_I) = -S_C2D_chain_chain(chain1_I, chain2_I);}
            }
            result_C2D_cfg1_cfg2(cfg1_I, cfg2_I) = overlap_C * calc_pfaffian(S_C2D_chain_chain);
        }
    }
}

inline const Eigen::MatrixXcd& HFBPfaffian::calc_overlap() {
    // X = ∅; NX = 0.
    Eigen::Map<Eigen::MatrixXcd> Qp1X_C2D_qp1_X(Qp1Xworkspace_C1D_element.data(), Nsp_I, 0);
    Eigen::Map<Eigen::MatrixXcd> XQp2Dag_C2D_X_qp2(XQp2Dagworkspace_C1D_element.data(), 0, Nsp_I);
    Eigen::Map<Eigen::MatrixXcd> XX_C2D_X_X(XXworkspace_C1D_element.data(), 0, 0);
    calc_kernel(Qp1X_C2D_qp1_X, XQp2Dag_C2D_X_qp2, XX_C2D_X_X, overlap_C2D_cfg1_cfg2);
    return overlap_C2D_cfg1_cfg2;
}

inline const Eigen::MatrixXcd& HFBPfaffian::calc_creator(const Eigen::VectorXcd& x_C1D_sp) {
    assert(x_C1D_sp.size() == Nsp_I && x_C1D_sp.allFinite());
    assert(overlap_C != doubleC(0.0, 0.0));

    // NX = 1; workspace views require no copies.
    Eigen::Map<Eigen::MatrixXcd> Qp1X_C2D_qp1_X(Qp1Xworkspace_C1D_element.data(), Nsp_I, 1);
    Eigen::Map<Eigen::MatrixXcd> XQp2Dag_C2D_X_qp2(XQp2Dagworkspace_C1D_element.data(), 1, Nsp_I);
    Eigen::Map<Eigen::MatrixXcd> XX_C2D_X_X(XXworkspace_C1D_element.data(), 1, 1);
    XX_C2D_X_X.setZero();

    // Coefficients enter linearly, without conjugation.
    Qp1X_C2D_qp1_X.col(0).noalias() = Qp1SpDag_C2D_qp1_sp * x_C1D_sp;
    XQp2Dag_C2D_X_qp2.row(0).noalias() = x_C1D_sp.transpose() * SpDagQp2Dag_C2D_sp_qp2;
    calc_kernel(Qp1X_C2D_qp1_X, XQp2Dag_C2D_X_qp2, XX_C2D_X_X, creator_C2D_cfg1_cfg2);
    return creator_C2D_cfg1_cfg2;
}

inline const Eigen::MatrixXcd& HFBPfaffian::calc_annihilator(const Eigen::VectorXcd& x_C1D_sp) {
    assert(x_C1D_sp.size() == Nsp_I && x_C1D_sp.allFinite());
    assert(overlap_C != doubleC(0.0, 0.0));

    // NX = 1; workspace views require no copies.
    Eigen::Map<Eigen::MatrixXcd> Qp1X_C2D_qp1_X(Qp1Xworkspace_C1D_element.data(), Nsp_I, 1);
    Eigen::Map<Eigen::MatrixXcd> XQp2Dag_C2D_X_qp2(XQp2Dagworkspace_C1D_element.data(), 1, Nsp_I);
    Eigen::Map<Eigen::MatrixXcd> XX_C2D_X_X(XXworkspace_C1D_element.data(), 1, 1);
    XX_C2D_X_X.setZero();

    // Coefficients enter linearly, without conjugation.
    Qp1X_C2D_qp1_X.col(0).noalias() = Qp1Sp_C2D_qp1_sp * x_C1D_sp;
    XQp2Dag_C2D_X_qp2.row(0).noalias() = x_C1D_sp.transpose() * SpQp2Dag_C2D_sp_qp2;
    calc_kernel(Qp1X_C2D_qp1_X, XQp2Dag_C2D_X_qp2, XX_C2D_X_X, annihilator_C2D_cfg1_cfg2);
    return annihilator_C2D_cfg1_cfg2;
}

inline const Eigen::MatrixXcd& HFBPfaffian::calc_obtd(const Eigen::VectorXcd& x1_C1D_sp, const Eigen::VectorXcd& x2_C1D_sp) {
    assert(x1_C1D_sp.size() == Nsp_I && x1_C1D_sp.allFinite());
    assert(x2_C1D_sp.size() == Nsp_I && x2_C1D_sp.allFinite());
    assert(overlap_C != doubleC(0.0, 0.0));

    // X = (x₁†,x₂); workspace views require no copies.
    Eigen::Map<Eigen::MatrixXcd> Qp1X_C2D_qp1_X(Qp1Xworkspace_C1D_element.data(), Nsp_I, 2);
    Eigen::Map<Eigen::MatrixXcd> XQp2Dag_C2D_X_qp2(XQp2Dagworkspace_C1D_element.data(), 2, Nsp_I);
    Eigen::Map<Eigen::MatrixXcd> XX_C2D_X_X(XXworkspace_C1D_element.data(), 2, 2);

    // Creation/annihilation coefficients enter without conjugation.
    Qp1X_C2D_qp1_X.col(0).noalias() = Qp1SpDag_C2D_qp1_sp * x1_C1D_sp;
    Qp1X_C2D_qp1_X.col(1).noalias() = Qp1Sp_C2D_qp1_sp * x2_C1D_sp;
    XQp2Dag_C2D_X_qp2.row(0).noalias() = x1_C1D_sp.transpose() * SpDagQp2Dag_C2D_sp_qp2;
    XQp2Dag_C2D_X_qp2.row(1).noalias() = x2_C1D_sp.transpose() * SpQp2Dag_C2D_sp_qp2;

    // XX₀₁ = ⟨x₁†x₂⟩; XXᵀ = −XX.
    XX_C2D_X_X.diagonal().setZero();
    XX_C2D_X_X(0, 1) = (x1_C1D_sp.transpose() * SpDagSp_C2D_sp_sp * x2_C1D_sp)(0, 0);
    XX_C2D_X_X(1, 0) = -XX_C2D_X_X(0, 1);
    calc_kernel(Qp1X_C2D_qp1_X, XQp2Dag_C2D_X_qp2, XX_C2D_X_X, OBTD_C2D_cfg1_cfg2);
    return OBTD_C2D_cfg1_cfg2;
}

inline const Eigen::MatrixXcd& HFBPfaffian::calc_tbtd(const Eigen::VectorXcd& x1_C1D_sp, const Eigen::VectorXcd& x2_C1D_sp, const Eigen::VectorXcd& x3_C1D_sp, const Eigen::VectorXcd& x4_C1D_sp) {
    assert(x1_C1D_sp.size() == Nsp_I && x1_C1D_sp.allFinite());
    assert(x2_C1D_sp.size() == Nsp_I && x2_C1D_sp.allFinite());
    assert(x3_C1D_sp.size() == Nsp_I && x3_C1D_sp.allFinite());
    assert(x4_C1D_sp.size() == Nsp_I && x4_C1D_sp.allFinite());
    assert(overlap_C != doubleC(0.0, 0.0));

    // X = (x₁†,x₂†,x₄,x₃); NX = 4.
    Eigen::Map<Eigen::MatrixXcd> Qp1X_C2D_qp1_X(Qp1Xworkspace_C1D_element.data(), Nsp_I, 4);
    Eigen::Map<Eigen::MatrixXcd> XQp2Dag_C2D_X_qp2(XQp2Dagworkspace_C1D_element.data(), 4, Nsp_I);
    Eigen::Map<Eigen::MatrixXcd> XX_C2D_X_X(XXworkspace_C1D_element.data(), 4, 4);

    // Creation/annihilation coefficients enter without conjugation.
    Qp1X_C2D_qp1_X.col(0).noalias() = Qp1SpDag_C2D_qp1_sp * x1_C1D_sp;
    Qp1X_C2D_qp1_X.col(1).noalias() = Qp1SpDag_C2D_qp1_sp * x2_C1D_sp;
    Qp1X_C2D_qp1_X.col(2).noalias() = Qp1Sp_C2D_qp1_sp * x4_C1D_sp;
    Qp1X_C2D_qp1_X.col(3).noalias() = Qp1Sp_C2D_qp1_sp * x3_C1D_sp;
    XQp2Dag_C2D_X_qp2.row(0).noalias() = x1_C1D_sp.transpose() * SpDagQp2Dag_C2D_sp_qp2;
    XQp2Dag_C2D_X_qp2.row(1).noalias() = x2_C1D_sp.transpose() * SpDagQp2Dag_C2D_sp_qp2;
    XQp2Dag_C2D_X_qp2.row(2).noalias() = x4_C1D_sp.transpose() * SpQp2Dag_C2D_sp_qp2;
    XQp2Dag_C2D_X_qp2.row(3).noalias() = x3_C1D_sp.transpose() * SpQp2Dag_C2D_sp_qp2;

    // Six upper-triangle contractions in the same insertion order.
    XX_C2D_X_X.diagonal().setZero();
    XX_C2D_X_X(0, 1) = (x1_C1D_sp.transpose() * SpDagSpDag_C2D_sp_sp * x2_C1D_sp)(0, 0);
    XX_C2D_X_X(0, 2) = (x1_C1D_sp.transpose() * SpDagSp_C2D_sp_sp * x4_C1D_sp)(0, 0);
    XX_C2D_X_X(0, 3) = (x1_C1D_sp.transpose() * SpDagSp_C2D_sp_sp * x3_C1D_sp)(0, 0);
    XX_C2D_X_X(1, 2) = (x2_C1D_sp.transpose() * SpDagSp_C2D_sp_sp * x4_C1D_sp)(0, 0);
    XX_C2D_X_X(1, 3) = (x2_C1D_sp.transpose() * SpDagSp_C2D_sp_sp * x3_C1D_sp)(0, 0);
    XX_C2D_X_X(2, 3) = (x4_C1D_sp.transpose() * SpSp_C2D_sp_sp * x3_C1D_sp)(0, 0);
    // XX_ji = −XX_ij; no complex conjugation.
    XX_C2D_X_X(1, 0) = -XX_C2D_X_X(0, 1);
    XX_C2D_X_X(2, 0) = -XX_C2D_X_X(0, 2);
    XX_C2D_X_X(3, 0) = -XX_C2D_X_X(0, 3);
    XX_C2D_X_X(2, 1) = -XX_C2D_X_X(1, 2);
    XX_C2D_X_X(3, 1) = -XX_C2D_X_X(1, 3);
    XX_C2D_X_X(3, 2) = -XX_C2D_X_X(2, 3);
    calc_kernel(Qp1X_C2D_qp1_X, XQp2Dag_C2D_X_qp2, XX_C2D_X_X, TBTD_C2D_cfg1_cfg2);
    return TBTD_C2D_cfg1_cfg2;
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

        // S_ij ← S_ij + (S_{k+1,i}S_kj-S_ki S_{k+1,j})/S_{k,k+1}.
        for (Eigen::Index chain2_I = k_I + 3; chain2_I < Nchain_I; ++chain2_I) {
            for (Eigen::Index chain1_I = k_I + 2; chain1_I < chain2_I; ++chain1_I) {
                S_C2D_chain_chain(chain1_I, chain2_I) += S_C2D_chain_chain(k_I + 1, chain1_I) * (S_C2D_chain_chain(k_I, chain2_I) / pivot_C) - S_C2D_chain_chain(k_I, chain1_I) * (S_C2D_chain_chain(k_I + 1, chain2_I) / pivot_C);
                S_C2D_chain_chain(chain2_I, chain1_I) = -S_C2D_chain_chain(chain1_I, chain2_I);
            }
        }
    }
    return pf_C;
}
