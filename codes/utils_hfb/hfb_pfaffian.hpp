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
 * @note   calc_overlap/obtd/tbtd require current contractions and nonzero overlap.
 * @note   Returned references alias reusable outputs; indices are zero-based.
 */
class HFBPfaffian {
public:
    int Nsp_I = 0;
    int NcqpMax1_I = 0;
    int NcqpMax2_I = 0;

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
    Eigen::MatrixXcd OBTD_C2D_cfg1_cfg2{};
    Eigen::MatrixXcd TBTD_C2D_cfg1_cfg2{};

private:
    Eigen::FullPivLU<Eigen::MatrixXcd> U1_lu{};
    Eigen::FullPivLU<Eigen::MatrixXcd> U2_lu{};
    Eigen::FullPivLU<Eigen::MatrixXcd> A_lu{};
    Eigen::VectorXcd Sworkspace_C1D_element{};

public:
    /**
     * @brief  Enumerate lexicographic configurations and allocate HFB workspaces.
     * @math   U₁,V₁,U₂,V₂ ∈ ℂ^{Nsp×Nsp}.
     * @output Generated both configuration tables; allocated matrices and workspaces.
     * @note   Each side uses its maximum's parity; step = 2.
     */
    HFBPfaffian(int Nsp_I_, int NcqpMax1_I_, int NcqpMax2_I_) {
        assert(Nsp_I_ > 0 && NcqpMax1_I_ >= 0 && NcqpMax1_I_ <= Nsp_I_ && NcqpMax2_I_ >= 0 && NcqpMax2_I_ <= Nsp_I_);

        Nsp_I = Nsp_I_;
        NcqpMax1_I = NcqpMax1_I_;
        NcqpMax2_I = NcqpMax2_I_;

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

        // cfg₁: Ncqp = NcqpMax1 mod 2, …, NcqpMax1; step = 2.
        for (int Ncqp_I = NcqpMax1_I % 2; Ncqp_I <= NcqpMax1_I; Ncqp_I += 2) {
            const std::vector<std::vector<int>> config_I2D_cfg_cqp{build_configs(Nsp_I, Ncqp_I)};
            config1_I2D_cfg1_cqp1.insert(config1_I2D_cfg1_cqp1.end(), config_I2D_cfg_cqp.begin(), config_I2D_cfg_cqp.end());
        }

        // cfg₂: Ncqp = NcqpMax2 mod 2, …, NcqpMax2; step = 2.
        for (int Ncqp_I = NcqpMax2_I % 2; Ncqp_I <= NcqpMax2_I; Ncqp_I += 2) {
            const std::vector<std::vector<int>> config_I2D_cfg_cqp{build_configs(Nsp_I, Ncqp_I)};
            config2_I2D_cfg2_cqp2.insert(config2_I2D_cfg2_cqp2.end(), config_I2D_cfg_cqp.begin(), config_I2D_cfg_cqp.end());
        }

        // N, OBTD, TBTD ∈ ℂ^{Ncfg1×Ncfg2}.
        overlap_C2D_cfg1_cfg2.resize(config1_I2D_cfg1_cqp1.size(), config2_I2D_cfg2_cqp2.size());
        OBTD_C2D_cfg1_cfg2.resize(config1_I2D_cfg1_cqp1.size(), config2_I2D_cfg2_cqp2.size());
        TBTD_C2D_cfg1_cfg2.resize(config1_I2D_cfg1_cqp1.size(), config2_I2D_cfg2_cqp2.size());

        // Lmax = 2Nsp+4; size(Sworkspace) = Lmax².
        const int Lmax_I = 2 * Nsp_I + 4;
        Sworkspace_C1D_element.resize(Lmax_I * Lmax_I);
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
     * @brief  Calculate configuration overlaps using Wick Pfaffians.
     * @math   N_ab = ⟨Φ₁;a|Φ₂;b⟩.
     * @output Updated overlap_C2D_cfg1_cfg2 and its const reference.
     */
    const Eigen::MatrixXcd& calc_overlap();

    /**
     * @brief  Calculate one-body kernels using Wick Pfaffians.
     * @math   OBTD_ab(i,j) = ⟨Φ₁;a|c†ᵢcⱼ|Φ₂;b⟩.
     * @output Updated OBTD_C2D_cfg1_cfg2 and its const reference.
     * @note   Includes vacuum overlap; excludes operator coefficients.
     */
    const Eigen::MatrixXcd& calc_obtd(int sp1_I, int sp2_I);

    /**
     * @brief  Calculate two-body kernels using Wick Pfaffians.
     * @math   TBTD_ab(i,j,k,l) = ⟨Φ₁;a|c†ᵢc†ⱼcₗcₖ|Φ₂;b⟩.
     * @output Updated TBTD_C2D_cfg1_cfg2 and its const reference.
     * @note   Includes vacuum overlap; excludes coefficients and symmetry factors.
     */
    const Eigen::MatrixXcd& calc_tbtd(int sp1_I, int sp2_I, int sp3_I, int sp4_I);

    /**
     * @brief  Calculate Pfaffians using pivoted skew-symmetric elimination.
     * @math   pf(X)² = det(X); pf(∅) = 1.
     * @output Pfaffian value; input matrix overwritten.
     * @note   Requires finite, even-order square X with Xᵀ = -X.
     * @note   Exact-zero pivots return zero; no magnitude cutoff.
     */
    static doubleC calc_pfaffian(Eigen::Ref<Eigen::MatrixXcd> X_C2D_chain_chain);

    /**
     * @brief  Enumerate fixed-size configurations in lexicographic order.
     * @math   0 ≤ μ₀ < ⋯ < μ_{Ncqp-1} < Nsp.
     * @output All C(Nsp,Ncqp) configurations; Ncqp = 0 returns {{}}.
     */
    static std::vector<std::vector<int>> build_configs(int Nsp_I, int Ncqp_I);
};

/**
 * @brief  Update overlap and contractions using Pfaffians and full-pivot LU.
 * @math   A = U₁ᵀU₂* + V₁ᵀV₂*; Qp1Qp2Dag = A⁻ᵀ.
 * @output Updated U,V,Z, overlap, A, AInv, and eleven contractions.
 */
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

/**
 * @brief  Calculate configuration overlaps using Wick Pfaffians.
 * @math   N_ab = ⟨Φ₁|Φ₂⟩ pf(S_ab).
 * @output Updated overlap_C2D_cfg1_cfg2 and its const reference.
 */
inline const Eigen::MatrixXcd& HFBPfaffian::calc_overlap() {
    assert(overlap_C != doubleC(0.0, 0.0));

    // (cfg₁,cfg₂) → ⟨Φ₁;cfg₁|Φ₂;cfg₂⟩.
    for (int cfg2_I = 0; cfg2_I < static_cast<int>(config2_I2D_cfg2_cqp2.size()); ++cfg2_I) {
        const int Ncqp2_I = static_cast<int>(config2_I2D_cfg2_cqp2[cfg2_I].size());
        for (int cfg1_I = 0; cfg1_I < static_cast<int>(config1_I2D_cfg1_cqp1.size()); ++cfg1_I) {
            const int Ncqp1_I = static_cast<int>(config1_I2D_cfg1_cqp1[cfg1_I].size());
            const int Nchain_I = Ncqp1_I + Ncqp2_I;

            // Odd chains vanish; pf(∅) = 1.
            if (Nchain_I % 2 != 0) {
                overlap_C2D_cfg1_cfg2(cfg1_I, cfg2_I) = 0.0;
                continue;
            }

            // S ∈ ℂ^{(m+n)×(m+n)}; S_ii = 0.
            assert(static_cast<Eigen::Index>(Nchain_I) * Nchain_I <= Sworkspace_C1D_element.size());
            Eigen::Map<Eigen::MatrixXcd> S_C2D_chain_chain(Sworkspace_C1D_element.data(), Nchain_I, Nchain_I);
            for (int chain1_I = 0; chain1_I < Nchain_I; ++chain1_I) {S_C2D_chain_chain(chain1_I, chain1_I) = 0.0;}

            // (β₁,μₘ,…,β₁,μ₁): upper-left block.
            for (int chain2_I = 1; chain2_I < Ncqp1_I; ++chain2_I) {
                for (int chain1_I = 0; chain1_I < chain2_I; ++chain1_I) {
                    S_C2D_chain_chain(chain1_I, chain2_I) = Qp1Qp1_C2D_qp1_qp1(config1_I2D_cfg1_cqp1[cfg1_I][Ncqp1_I - 1 - chain1_I], config1_I2D_cfg1_cqp1[cfg1_I][Ncqp1_I - 1 - chain2_I]);
                }
            }

            // ⟨β₁,μ β₂,ν†⟩: upper-right block.
            for (int chain2_I = 0; chain2_I < Ncqp2_I; ++chain2_I) {
                for (int chain1_I = 0; chain1_I < Ncqp1_I; ++chain1_I) {
                    S_C2D_chain_chain(chain1_I, Ncqp1_I + chain2_I) = Qp1Qp2Dag_C2D_qp1_qp2(config1_I2D_cfg1_cqp1[cfg1_I][Ncqp1_I - 1 - chain1_I], config2_I2D_cfg2_cqp2[cfg2_I][chain2_I]);
                }
            }

            // (β₂,ν₁†,…,β₂,νₙ†): lower-right block.
            for (int chain2_I = 1; chain2_I < Ncqp2_I; ++chain2_I) {
                for (int chain1_I = 0; chain1_I < chain2_I; ++chain1_I) {
                    S_C2D_chain_chain(Ncqp1_I + chain1_I, Ncqp1_I + chain2_I) = Qp2DagQp2Dag_C2D_qp2_qp2(config2_I2D_cfg2_cqp2[cfg2_I][chain1_I], config2_I2D_cfg2_cqp2[cfg2_I][chain2_I]);
                }
            }

            // S_ji = -S_ij; N_ab = overlap × pf(S).
            for (int chain2_I = 1; chain2_I < Nchain_I; ++chain2_I) {
                for (int chain1_I = 0; chain1_I < chain2_I; ++chain1_I) {S_C2D_chain_chain(chain2_I, chain1_I) = -S_C2D_chain_chain(chain1_I, chain2_I);}
            }
            overlap_C2D_cfg1_cfg2(cfg1_I, cfg2_I) = overlap_C * calc_pfaffian(S_C2D_chain_chain);
        }
    }
    return overlap_C2D_cfg1_cfg2;
}

/**
 * @brief  Calculate one-body kernels using Wick Pfaffians.
 * @math   OBTD_ab(i,j) = ⟨Φ₁;a|c†ᵢcⱼ|Φ₂;b⟩.
 * @output Updated OBTD_C2D_cfg1_cfg2 and its const reference.
 */
inline const Eigen::MatrixXcd& HFBPfaffian::calc_obtd(int sp1_I, int sp2_I) {
    assert(sp1_I >= 0 && sp1_I < Nsp_I && sp2_I >= 0 && sp2_I < Nsp_I);
    assert(overlap_C != doubleC(0.0, 0.0));

    // (cfg₁,cfg₂) → ⟨Φ₁;cfg₁|c†₁c₂|Φ₂;cfg₂⟩.
    for (int cfg2_I = 0; cfg2_I < static_cast<int>(config2_I2D_cfg2_cqp2.size()); ++cfg2_I) {
        const int Ncqp2_I = static_cast<int>(config2_I2D_cfg2_cqp2[cfg2_I].size());
        for (int cfg1_I = 0; cfg1_I < static_cast<int>(config1_I2D_cfg1_cqp1.size()); ++cfg1_I) {
            const int Ncqp1_I = static_cast<int>(config1_I2D_cfg1_cqp1[cfg1_I].size());
            const int Nchain_I = Ncqp1_I + Ncqp2_I + 2;

            // Odd chains vanish.
            if (Nchain_I % 2 != 0) {
                OBTD_C2D_cfg1_cfg2(cfg1_I, cfg2_I) = 0.0;
                continue;
            }

            // S ∈ ℂ^{(m+n+2)×(m+n+2)}; S_ii = 0.
            assert(static_cast<Eigen::Index>(Nchain_I) * Nchain_I <= Sworkspace_C1D_element.size());
            Eigen::Map<Eigen::MatrixXcd> S_C2D_chain_chain(Sworkspace_C1D_element.data(), Nchain_I, Nchain_I);
            for (int chain1_I = 0; chain1_I < Nchain_I; ++chain1_I) {S_C2D_chain_chain(chain1_I, chain1_I) = 0.0;}

            // (β₁,μₘ,…,β₁,μ₁): upper-left block.
            for (int chain2_I = 1; chain2_I < Ncqp1_I; ++chain2_I) {
                for (int chain1_I = 0; chain1_I < chain2_I; ++chain1_I) {
                    S_C2D_chain_chain(chain1_I, chain2_I) = Qp1Qp1_C2D_qp1_qp1(config1_I2D_cfg1_cqp1[cfg1_I][Ncqp1_I - 1 - chain1_I], config1_I2D_cfg1_cqp1[cfg1_I][Ncqp1_I - 1 - chain2_I]);
                }
            }

            // ⟨β₁(c†₁,c₂)⟩.
            for (int chain1_I = 0; chain1_I < Ncqp1_I; ++chain1_I) {
                S_C2D_chain_chain(chain1_I, Ncqp1_I) = Qp1SpDag_C2D_qp1_sp(config1_I2D_cfg1_cqp1[cfg1_I][Ncqp1_I - 1 - chain1_I], sp1_I);
                S_C2D_chain_chain(chain1_I, Ncqp1_I + 1) = Qp1Sp_C2D_qp1_sp(config1_I2D_cfg1_cqp1[cfg1_I][Ncqp1_I - 1 - chain1_I], sp2_I);
            }

            S_C2D_chain_chain(Ncqp1_I, Ncqp1_I + 1) = SpDagSp_C2D_sp_sp(sp1_I, sp2_I);

            // ⟨(c†₁,c₂)β₂†⟩.
            for (int chain2_I = 0; chain2_I < Ncqp2_I; ++chain2_I) {
                S_C2D_chain_chain(Ncqp1_I, Ncqp1_I + 2 + chain2_I) = SpDagQp2Dag_C2D_sp_qp2(sp1_I, config2_I2D_cfg2_cqp2[cfg2_I][chain2_I]);
                S_C2D_chain_chain(Ncqp1_I + 1, Ncqp1_I + 2 + chain2_I) = SpQp2Dag_C2D_sp_qp2(sp2_I, config2_I2D_cfg2_cqp2[cfg2_I][chain2_I]);
            }

            // ⟨β₁,μ β₂,ν†⟩: upper-right block.
            for (int chain2_I = 0; chain2_I < Ncqp2_I; ++chain2_I) {
                for (int chain1_I = 0; chain1_I < Ncqp1_I; ++chain1_I) {
                    S_C2D_chain_chain(chain1_I, Ncqp1_I + 2 + chain2_I) = Qp1Qp2Dag_C2D_qp1_qp2(config1_I2D_cfg1_cqp1[cfg1_I][Ncqp1_I - 1 - chain1_I], config2_I2D_cfg2_cqp2[cfg2_I][chain2_I]);
                }
            }

            // (β₂,ν₁†,…,β₂,νₙ†): lower-right block.
            for (int chain2_I = 1; chain2_I < Ncqp2_I; ++chain2_I) {
                for (int chain1_I = 0; chain1_I < chain2_I; ++chain1_I) {
                    S_C2D_chain_chain(Ncqp1_I + 2 + chain1_I, Ncqp1_I + 2 + chain2_I) = Qp2DagQp2Dag_C2D_qp2_qp2(config2_I2D_cfg2_cqp2[cfg2_I][chain1_I], config2_I2D_cfg2_cqp2[cfg2_I][chain2_I]);
                }
            }

            // S_ji = -S_ij; OBTD_ab = overlap × pf(S).
            for (int chain2_I = 1; chain2_I < Nchain_I; ++chain2_I) {
                for (int chain1_I = 0; chain1_I < chain2_I; ++chain1_I) {S_C2D_chain_chain(chain2_I, chain1_I) = -S_C2D_chain_chain(chain1_I, chain2_I);}
            }
            OBTD_C2D_cfg1_cfg2(cfg1_I, cfg2_I) = overlap_C * calc_pfaffian(S_C2D_chain_chain);
        }
    }
    return OBTD_C2D_cfg1_cfg2;
}

/**
 * @brief  Calculate two-body kernels using Wick Pfaffians.
 * @math   TBTD_ab(i,j,k,l) = ⟨Φ₁;a|c†ᵢc†ⱼcₗcₖ|Φ₂;b⟩.
 * @output Updated TBTD_C2D_cfg1_cfg2 and its const reference.
 */
inline const Eigen::MatrixXcd& HFBPfaffian::calc_tbtd(int sp1_I, int sp2_I, int sp3_I, int sp4_I) {
    assert(sp1_I >= 0 && sp1_I < Nsp_I && sp2_I >= 0 && sp2_I < Nsp_I);
    assert(sp3_I >= 0 && sp3_I < Nsp_I && sp4_I >= 0 && sp4_I < Nsp_I);
    assert(overlap_C != doubleC(0.0, 0.0));

    // (cfg₁,cfg₂) → ⟨Φ₁;cfg₁|c†₁c†₂c₄c₃|Φ₂;cfg₂⟩.
    for (int cfg2_I = 0; cfg2_I < static_cast<int>(config2_I2D_cfg2_cqp2.size()); ++cfg2_I) {
        const int Ncqp2_I = static_cast<int>(config2_I2D_cfg2_cqp2[cfg2_I].size());
        for (int cfg1_I = 0; cfg1_I < static_cast<int>(config1_I2D_cfg1_cqp1.size()); ++cfg1_I) {
            const int Ncqp1_I = static_cast<int>(config1_I2D_cfg1_cqp1[cfg1_I].size());
            const int Nchain_I = Ncqp1_I + Ncqp2_I + 4;

            // Odd chains vanish.
            if (Nchain_I % 2 != 0) {
                TBTD_C2D_cfg1_cfg2(cfg1_I, cfg2_I) = 0.0;
                continue;
            }

            // S ∈ ℂ^{(m+n+4)×(m+n+4)}; S_ii = 0.
            assert(static_cast<Eigen::Index>(Nchain_I) * Nchain_I <= Sworkspace_C1D_element.size());
            Eigen::Map<Eigen::MatrixXcd> S_C2D_chain_chain(Sworkspace_C1D_element.data(), Nchain_I, Nchain_I);
            for (int chain1_I = 0; chain1_I < Nchain_I; ++chain1_I) {S_C2D_chain_chain(chain1_I, chain1_I) = 0.0;}

            // (β₁,μₘ,…,β₁,μ₁): upper-left block.
            for (int chain2_I = 1; chain2_I < Ncqp1_I; ++chain2_I) {
                for (int chain1_I = 0; chain1_I < chain2_I; ++chain1_I) {
                    S_C2D_chain_chain(chain1_I, chain2_I) = Qp1Qp1_C2D_qp1_qp1(config1_I2D_cfg1_cqp1[cfg1_I][Ncqp1_I - 1 - chain1_I], config1_I2D_cfg1_cqp1[cfg1_I][Ncqp1_I - 1 - chain2_I]);
                }
            }

            // ⟨β₁(c†₁,c†₂,c₄,c₃)⟩.
            for (int chain1_I = 0; chain1_I < Ncqp1_I; ++chain1_I) {
                S_C2D_chain_chain(chain1_I, Ncqp1_I) = Qp1SpDag_C2D_qp1_sp(config1_I2D_cfg1_cqp1[cfg1_I][Ncqp1_I - 1 - chain1_I], sp1_I);
                S_C2D_chain_chain(chain1_I, Ncqp1_I + 1) = Qp1SpDag_C2D_qp1_sp(config1_I2D_cfg1_cqp1[cfg1_I][Ncqp1_I - 1 - chain1_I], sp2_I);
                S_C2D_chain_chain(chain1_I, Ncqp1_I + 2) = Qp1Sp_C2D_qp1_sp(config1_I2D_cfg1_cqp1[cfg1_I][Ncqp1_I - 1 - chain1_I], sp4_I);
                S_C2D_chain_chain(chain1_I, Ncqp1_I + 3) = Qp1Sp_C2D_qp1_sp(config1_I2D_cfg1_cqp1[cfg1_I][Ncqp1_I - 1 - chain1_I], sp3_I);
            }

            // (c†₁,c†₂,c₄,c₃): six upper-triangle contractions.
            S_C2D_chain_chain(Ncqp1_I, Ncqp1_I + 1) = SpDagSpDag_C2D_sp_sp(sp1_I, sp2_I);
            S_C2D_chain_chain(Ncqp1_I, Ncqp1_I + 2) = SpDagSp_C2D_sp_sp(sp1_I, sp4_I);
            S_C2D_chain_chain(Ncqp1_I, Ncqp1_I + 3) = SpDagSp_C2D_sp_sp(sp1_I, sp3_I);
            S_C2D_chain_chain(Ncqp1_I + 1, Ncqp1_I + 2) = SpDagSp_C2D_sp_sp(sp2_I, sp4_I);
            S_C2D_chain_chain(Ncqp1_I + 1, Ncqp1_I + 3) = SpDagSp_C2D_sp_sp(sp2_I, sp3_I);
            S_C2D_chain_chain(Ncqp1_I + 2, Ncqp1_I + 3) = SpSp_C2D_sp_sp(sp4_I, sp3_I);

            // ⟨(c†₁,c†₂,c₄,c₃)β₂†⟩.
            for (int chain2_I = 0; chain2_I < Ncqp2_I; ++chain2_I) {
                S_C2D_chain_chain(Ncqp1_I, Ncqp1_I + 4 + chain2_I) = SpDagQp2Dag_C2D_sp_qp2(sp1_I, config2_I2D_cfg2_cqp2[cfg2_I][chain2_I]);
                S_C2D_chain_chain(Ncqp1_I + 1, Ncqp1_I + 4 + chain2_I) = SpDagQp2Dag_C2D_sp_qp2(sp2_I, config2_I2D_cfg2_cqp2[cfg2_I][chain2_I]);
                S_C2D_chain_chain(Ncqp1_I + 2, Ncqp1_I + 4 + chain2_I) = SpQp2Dag_C2D_sp_qp2(sp4_I, config2_I2D_cfg2_cqp2[cfg2_I][chain2_I]);
                S_C2D_chain_chain(Ncqp1_I + 3, Ncqp1_I + 4 + chain2_I) = SpQp2Dag_C2D_sp_qp2(sp3_I, config2_I2D_cfg2_cqp2[cfg2_I][chain2_I]);
            }

            // ⟨β₁,μ β₂,ν†⟩: upper-right block.
            for (int chain2_I = 0; chain2_I < Ncqp2_I; ++chain2_I) {
                for (int chain1_I = 0; chain1_I < Ncqp1_I; ++chain1_I) {
                    S_C2D_chain_chain(chain1_I, Ncqp1_I + 4 + chain2_I) = Qp1Qp2Dag_C2D_qp1_qp2(config1_I2D_cfg1_cqp1[cfg1_I][Ncqp1_I - 1 - chain1_I], config2_I2D_cfg2_cqp2[cfg2_I][chain2_I]);
                }
            }

            // (β₂,ν₁†,…,β₂,νₙ†): lower-right block.
            for (int chain2_I = 1; chain2_I < Ncqp2_I; ++chain2_I) {
                for (int chain1_I = 0; chain1_I < chain2_I; ++chain1_I) {
                    S_C2D_chain_chain(Ncqp1_I + 4 + chain1_I, Ncqp1_I + 4 + chain2_I) = Qp2DagQp2Dag_C2D_qp2_qp2(config2_I2D_cfg2_cqp2[cfg2_I][chain1_I], config2_I2D_cfg2_cqp2[cfg2_I][chain2_I]);
                }
            }

            // S_ji = -S_ij; TBTD_ab = overlap × pf(S).
            for (int chain2_I = 1; chain2_I < Nchain_I; ++chain2_I) {
                for (int chain1_I = 0; chain1_I < chain2_I; ++chain1_I) {S_C2D_chain_chain(chain2_I, chain1_I) = -S_C2D_chain_chain(chain1_I, chain2_I);}
            }
            TBTD_C2D_cfg1_cfg2(cfg1_I, cfg2_I) = overlap_C * calc_pfaffian(S_C2D_chain_chain);
        }
    }
    return TBTD_C2D_cfg1_cfg2;
}

/**
 * @brief  Calculate Pfaffians using pivoted skew-symmetric elimination.
 * @math   pf(PXPᵀ) = det(P)pf(X); pf(∅) = 1.
 * @output Pfaffian value; input matrix overwritten.
 */
inline doubleC HFBPfaffian::calc_pfaffian(Eigen::Ref<Eigen::MatrixXcd> X_C2D_chain_chain) {
    const Eigen::Index Nchain_I = X_C2D_chain_chain.rows();
    assert(X_C2D_chain_chain.cols() == Nchain_I && Nchain_I % 2 == 0);
    assert(X_C2D_chain_chain.allFinite());
    assert(X_C2D_chain_chain.isApprox(-X_C2D_chain_chain.transpose(), 1.0e-12));
    doubleC pf_C(1.0, 0.0);

    // X → PXPᵀ → 2×2 pivot ⊕ Schur complement.
    for (Eigen::Index k_I = 0; k_I < Nchain_I; k_I += 2) {
        // Pivot: argmax_{j>k}|X_kj|.
        Eigen::Index pivotOffset_I = 0;
        const double pivotAbs_F = X_C2D_chain_chain.row(k_I).segment(k_I + 1, Nchain_I - k_I - 1).cwiseAbs().maxCoeff(&pivotOffset_I);
        const Eigen::Index pivot_I = k_I + 1 + pivotOffset_I;
        if (pivotAbs_F == 0.0) {return doubleC(0.0, 0.0);}

        // One simultaneous row-column swap contributes -1.
        if (pivot_I != k_I + 1) {
            X_C2D_chain_chain.row(k_I + 1).swap(X_C2D_chain_chain.row(pivot_I));
            X_C2D_chain_chain.col(k_I + 1).swap(X_C2D_chain_chain.col(pivot_I));
            pf_C = -pf_C;
        }
        const doubleC pivot_C = X_C2D_chain_chain(k_I, k_I + 1);
        pf_C *= pivot_C;

        // X_ij ← X_ij + (X_{k+1,i}X_kj-X_ki X_{k+1,j})/X_{k,k+1}.
        for (Eigen::Index chain2_I = k_I + 3; chain2_I < Nchain_I; ++chain2_I) {
            for (Eigen::Index chain1_I = k_I + 2; chain1_I < chain2_I; ++chain1_I) {
                X_C2D_chain_chain(chain1_I, chain2_I) += X_C2D_chain_chain(k_I + 1, chain1_I) * (X_C2D_chain_chain(k_I, chain2_I) / pivot_C) - X_C2D_chain_chain(k_I, chain1_I) * (X_C2D_chain_chain(k_I + 1, chain2_I) / pivot_C);
                X_C2D_chain_chain(chain2_I, chain1_I) = -X_C2D_chain_chain(chain1_I, chain2_I);
            }
        }
    }
    return pf_C;
}

inline std::vector<std::vector<int>> HFBPfaffian::build_configs(int Nsp_I, int Ncqp_I) {
    assert(Ncqp_I >= 0 && Ncqp_I <= Nsp_I);

    // μ = (0,1,…,Ncqp-1).
    std::vector<std::vector<int>> config_I2D_cfg_cqp{};
    std::vector<int> config_I1D_cqp{};
    config_I1D_cqp.resize(Ncqp_I);
    for (int cqp_I = 0; cqp_I < Ncqp_I; ++cqp_I) {config_I1D_cqp[cqp_I] = cqp_I;}

    // μ_p < Nsp-Ncqp+p → μ_p+1; reset μ_{p+1},….
    while (true) {
        config_I2D_cfg_cqp.push_back(config_I1D_cqp);
        int pivot_I = Ncqp_I - 1;
        while (pivot_I >= 0 && config_I1D_cqp[pivot_I] == Nsp_I - Ncqp_I + pivot_I) {
            --pivot_I;
        }
        if (pivot_I < 0) {break;}
        ++config_I1D_cqp[pivot_I];
        for (int cqp_I = pivot_I + 1; cqp_I < Ncqp_I; ++cqp_I) {config_I1D_cqp[cqp_I] = config_I1D_cqp[cqp_I - 1] + 1;}
    }
    return config_I2D_cfg_cqp;
}
