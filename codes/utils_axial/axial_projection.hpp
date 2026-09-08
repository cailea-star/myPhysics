/**
 * @file    axial_projection.hpp
 * @author  cailea
 * @date    2026-09-07
 * @brief   Axial Bogoliubov overlaps and transition densities.
 */

#pragma once

#include <cassert>
#include <vector>
#include <cmath>
#include <complex>
#include <limits>
#include <utility>
#include <Eigen/Core>
#include <Eigen/LU>
#include <unsupported/Eigen/CXX11/Tensor>

#include "axial_rotation.hpp"

using doubleC = std::complex<double>;

class AxialProjection {
public:
    int Nbeta_I = 1; // Number of β quadrature nodes.
    int Nphi_I = 0; // Number of gauge-angle quadrature nodes.
    Eigen::VectorXd beta_F1D_beta{}; // β_i = acos x_i ∈ (0,π) [rad].
    Eigen::VectorXd weight_F1D_beta{}; // Σ_i w_i f(β_i) ≈ ∫₀^π sinβ f(β)dβ.
    std::vector<AxialSPLabel> labels_S1D_sp{}; // Positive-Ω labels.

    doubleC overlap_C = doubleC(0, 0); // ⟨Φ₁|R_y(β)e^{iφN̂}|Φ₂⟩.
    Eigen::MatrixXcd rho_C2D_2sp_2sp{}; // ρ(g).
    Eigen::MatrixXcd kappa_C2D_2sp_2sp{}; // κ(g).
    Eigen::MatrixXcd kappaBar_C2D_2sp_2sp{}; // κ̄(g).

private:
    Eigen::VectorXcd detRySimplexP_C1D_beta{}; // det d₊(β); det d₋ = (det d₊)*.
    Eigen::Tensor<double, 3, Eigen::ColMajor> Ry_F3D_2sp_2sp_beta{}; // ⟨φ_a|R_y(β_i)|φ_b⟩.

    // s_a = sgn det U_a; positive vacuum phases.
    int signDetU1_I = 0; // Zero until reference initialization succeeds.
    int signDetU2_I = 0; // Zero until reference initialization succeeds.
    int betaCurrent_I = -1; // −1: β cache invalid.

    // Kramers representatives: U⁺ = U; V⁻ = ηV.
    Eigen::MatrixXd U1Kramers_F2D_sp_qp{}; // U₁.
    Eigen::MatrixXd V1Kramers_F2D_sp_qp{}; // V₁.
    Eigen::MatrixXd U2Kramers_F2D_sp_qp{}; // U₂.
    Eigen::MatrixXd V2Kramers_F2D_sp_qp{}; // V₂.

    // U₂β = R_y^{-T}U₂; V₂β = R_yV₂.
    Eigen::MatrixXd U2beta_F2D_2sp_2qp{}; // U₂(β); full (+Ω,-Ω) ordering.
    Eigen::MatrixXd V2beta_F2D_2sp_2qp{}; // V₂(β); full (+Ω,-Ω) ordering.
    Eigen::MatrixXd U1TU2beta_F2D_2qp_2qp{}; // A_Uβ = U₁ᵀU₂β.
    Eigen::MatrixXd V1TV2beta_F2D_2qp_2qp{}; // A_Vβ = V₁ᵀV₂β.
    Eigen::MatrixXcd Ag_C2D_2qp_2qp{}; // A(β,φ).
    Eigen::FullPivLU<Eigen::MatrixXcd> Ag_LU{}; // PAQ = LU.

    Eigen::MatrixXd Tmp_F2D_2sp_2sp{}; // Real β workspace.
    Eigen::MatrixXcd Tmp_C2D_2sp_2sp{}; // AX = V₁ᵀ or U₁ᵀ.

public:
    /**
     * @brief  Build y rotations using harmonic-oscillator generating functions.
     * @math   (C_axial,{α_sp},N_β,N_φ) → {β,w_β,R_y,det d₊}; α = γ = 0.
     * @output β quadrature and rotations; references remain unset.
     * @note   Single species; unblocked K = 0 vacua; ordering (+Ω,-Ω).
     * @note   Weights exclude projection normalization factors.
     */
    AxialProjection(const AxialConfig& axialconfig_, const std::vector<AxialSPLabel>& labels_S1D_sp_, int Nbeta_I_, int Nphi_I_) {
        labels_S1D_sp = labels_S1D_sp_;
        Nbeta_I = Nbeta_I_;
        Nphi_I = Nphi_I_;
        assert(Nbeta_I > 0 && Nphi_I > 0);

        // Retain β geometry; release rotation-generation workspace.
        {
            AxialRotation axialrotation(axialconfig_, labels_S1D_sp, 1, Nbeta_I, 1);
            axialrotation.build_Ry();
            beta_F1D_beta = std::move(axialrotation.beta_F1D_beta);
            weight_F1D_beta = std::move(axialrotation.weight_F1D_beta);
            detRySimplexP_C1D_beta = std::move(axialrotation.detRySimplexP_C1D_beta);
            Ry_F3D_2sp_2sp_beta = std::move(axialrotation.Ry_F3D_2sp_2sp_beta);
        }

        const Eigen::Index Nsp_I = static_cast<Eigen::Index>(labels_S1D_sp.size()); // #{a: Ω_a > 0}.

        // {U₁,V₁,U₂,V₂}: N_sp×N_qp; N_qp = N_sp.
        U1Kramers_F2D_sp_qp.resize(Nsp_I, Nsp_I);
        V1Kramers_F2D_sp_qp.resize(Nsp_I, Nsp_I);
        U2Kramers_F2D_sp_qp.resize(Nsp_I, Nsp_I);
        V2Kramers_F2D_sp_qp.resize(Nsp_I, Nsp_I);

        // {U₂β,V₂β,A_Uβ,A_Vβ,A_g,X,ρ,κ,κ̄}: full (+Ω,-Ω) sectors.
        U2beta_F2D_2sp_2qp.resize(2 * Nsp_I, 2 * Nsp_I);
        V2beta_F2D_2sp_2qp.resize(2 * Nsp_I, 2 * Nsp_I);
        U1TU2beta_F2D_2qp_2qp.resize(2 * Nsp_I, 2 * Nsp_I);
        V1TV2beta_F2D_2qp_2qp.resize(2 * Nsp_I, 2 * Nsp_I);
        Ag_C2D_2qp_2qp.resize(2 * Nsp_I, 2 * Nsp_I);
        rho_C2D_2sp_2sp.resize(2 * Nsp_I, 2 * Nsp_I);
        kappa_C2D_2sp_2sp.resize(2 * Nsp_I, 2 * Nsp_I);
        kappaBar_C2D_2sp_2sp.resize(2 * Nsp_I, 2 * Nsp_I);
        Tmp_F2D_2sp_2sp.resize(2 * Nsp_I, 2 * Nsp_I);
        Tmp_C2D_2sp_2sp.resize(2 * Nsp_I, 2 * Nsp_I);

        // Preallocate LU workspace; factorization remains deferred.
        Ag_LU = Eigen::FullPivLU<Eigen::MatrixXcd>(2 * Nsp_I, 2 * Nsp_I);
    }

    /**
     * @brief  Cache Kramers amplitudes and pivoted-LU determinant signs.
     * @math   U_full = diag(U,ηU); V_full = [0,-V;ηV,0]; η_a = 2Σ_a.
     * @output Left reference and sign(det U₁); β cache invalidated.
     * @note   Real canonical U,V; invertible U; positive vacuum phase.
     * @note   Rows: positive-Ω labels; columns: Ω blocks; V includes η.
     */
    void update_U1V1(const Eigen::MatrixXd& U1Kramers_F2D_sp_qp_, const Eigen::MatrixXd& V1Kramers_F2D_sp_qp_);

    /**
     * @brief  Cache Kramers amplitudes and pivoted-LU determinant signs.
     * @math   U_full = diag(U,ηU); V_full = [0,-V;ηV,0]; η_a = 2Σ_a.
     * @output Right reference and sign(det U₂); β cache invalidated.
     * @note   Real canonical U,V; invertible U; positive vacuum phase.
     * @note   Rows: positive-Ω labels; columns: Ω blocks; V includes η.
     */
    void update_U2V2(const Eigen::MatrixXd& U2Kramers_F2D_sp_qp_, const Eigen::MatrixXd& V2Kramers_F2D_sp_qp_);

    /**
     * @brief  Update transformed matrices using full-pivoting LU.
     * @math   R_yᵀU₂β = U₂; V₂β = R_yV₂; A_Uβ = U₁ᵀU₂β; A_Vβ = V₁ᵀV₂β.
     * @output U₂β, V₂β, A_Uβ, A_Vβ and current β index.
     * @note   Requires update_U1V1() and update_U2V2() first.
     * @note   Requires invertible R_y.
     */
    void update_AUV(int beta_I);

    /**
     * @brief  Update overlap and densities using full-pivoting LU.
     * @math   A = e^{-iφ}A_Uβ+e^{iφ}A_Vβ; n = s₁s₂e^{iN_spφ}det d₊det A₊.
     * @output n, ρ, κ, κ̄ at β[beta_I], φ = 2π phi_I/N_φ.
     * @note   Requires update_U1V1() and update_U2V2() first.
     * @note   Requires invertible A.
     */
    void update_densities(int beta_I, int phi_I);

};

inline void AxialProjection::update_U1V1(const Eigen::MatrixXd& U1Kramers_F2D_sp_qp_, const Eigen::MatrixXd& V1Kramers_F2D_sp_qp_) {
    // N_sp = N_qp; Kramers representatives.
    const Eigen::Index Nsp_I = static_cast<Eigen::Index>(labels_S1D_sp.size());
    assert(U1Kramers_F2D_sp_qp_.rows() == Nsp_I && U1Kramers_F2D_sp_qp_.cols() == Nsp_I);
    assert(V1Kramers_F2D_sp_qp_.rows() == Nsp_I && V1Kramers_F2D_sp_qp_.cols() == Nsp_I);
    assert(U1Kramers_F2D_sp_qp_.allFinite());
    assert(V1Kramers_F2D_sp_qp_.allFinite());

    U1Kramers_F2D_sp_qp = U1Kramers_F2D_sp_qp_;
    V1Kramers_F2D_sp_qp = V1Kramers_F2D_sp_qp_;

    // Reference changes invalidate both overlap and pointwise kernels.
    const double detU1_F = U1Kramers_F2D_sp_qp.determinant();
    signDetU1_I = (detU1_F > 0.0) - (detU1_F < 0.0);
    betaCurrent_I = -1;

    // UᵀU + VᵀV = I; UᵀV = VᵀU.
    assert((U1Kramers_F2D_sp_qp.transpose() * U1Kramers_F2D_sp_qp + V1Kramers_F2D_sp_qp.transpose() * V1Kramers_F2D_sp_qp - Eigen::MatrixXd::Identity(Nsp_I, Nsp_I)).norm() <= 1.0e-10 * (1.0 + Nsp_I));
    assert((U1Kramers_F2D_sp_qp.transpose() * V1Kramers_F2D_sp_qp - V1Kramers_F2D_sp_qp.transpose() * U1Kramers_F2D_sp_qp).norm() <= 1.0e-10 * (1.0 + Nsp_I));
    assert(std::isfinite(detU1_F) && detU1_F != 0.0);
}

inline void AxialProjection::update_U2V2(const Eigen::MatrixXd& U2Kramers_F2D_sp_qp_, const Eigen::MatrixXd& V2Kramers_F2D_sp_qp_) {
    // N_sp = N_qp; Kramers representatives.
    const Eigen::Index Nsp_I = static_cast<Eigen::Index>(labels_S1D_sp.size());
    assert(U2Kramers_F2D_sp_qp_.rows() == Nsp_I && U2Kramers_F2D_sp_qp_.cols() == Nsp_I);
    assert(V2Kramers_F2D_sp_qp_.rows() == Nsp_I && V2Kramers_F2D_sp_qp_.cols() == Nsp_I);
    assert(U2Kramers_F2D_sp_qp_.allFinite());
    assert(V2Kramers_F2D_sp_qp_.allFinite());

    U2Kramers_F2D_sp_qp = U2Kramers_F2D_sp_qp_;
    V2Kramers_F2D_sp_qp = V2Kramers_F2D_sp_qp_;

    // Reference changes invalidate both overlap and pointwise kernels.
    const double detU2_F = U2Kramers_F2D_sp_qp.determinant();
    signDetU2_I = (detU2_F > 0.0) - (detU2_F < 0.0);
    betaCurrent_I = -1;

    // UᵀU + VᵀV = I; UᵀV = VᵀU.
    assert((U2Kramers_F2D_sp_qp.transpose() * U2Kramers_F2D_sp_qp + V2Kramers_F2D_sp_qp.transpose() * V2Kramers_F2D_sp_qp - Eigen::MatrixXd::Identity(Nsp_I, Nsp_I)).norm() <= 1.0e-10 * (1.0 + Nsp_I));
    assert((U2Kramers_F2D_sp_qp.transpose() * V2Kramers_F2D_sp_qp - V2Kramers_F2D_sp_qp.transpose() * U2Kramers_F2D_sp_qp).norm() <= 1.0e-10 * (1.0 + Nsp_I));
    assert(std::isfinite(detU2_F) && detU2_F != 0.0);
}

inline void AxialProjection::update_AUV(int beta_I) {
    betaCurrent_I = -1; // β cache invalid.

    // N_sp = N_qp; fixed β.
    const Eigen::Index Nsp_I = static_cast<Eigen::Index>(labels_S1D_sp.size());
    assert(beta_I >= 0 && beta_I < beta_F1D_beta.size());

    // ColMajor slices: offset = (2N_sp)² i.
    assert(Ry_F3D_2sp_2sp_beta.dimension(0) == 2 * Nsp_I && Ry_F3D_2sp_2sp_beta.dimension(1) == 2 * Nsp_I && beta_I < Ry_F3D_2sp_2sp_beta.dimension(2));
    const Eigen::Map<const Eigen::MatrixXd> Ry_F2D_2sp_2sp(Ry_F3D_2sp_2sp_beta.data() + 4 * Nsp_I * Nsp_I * beta_I, 2 * Nsp_I, 2 * Nsp_I);

    assert(signDetU1_I != 0 && signDetU2_I != 0); // Both references are initialized.
    assert(detRySimplexP_C1D_beta.size() == beta_F1D_beta.size());

    // R_yᵀU₂β = U₂,full.
    Ag_LU.compute(Ry_F2D_2sp_2sp.transpose().cast<doubleC>());
    assert(Ag_LU.isInvertible());

    // U₂,full = diag(U₂,ηU₂).
    Tmp_F2D_2sp_2sp.setZero();
    Tmp_F2D_2sp_2sp.topLeftCorner(Nsp_I, Nsp_I) = U2Kramers_F2D_sp_qp;
    Tmp_F2D_2sp_2sp.bottomRightCorner(Nsp_I, Nsp_I) = U2Kramers_F2D_sp_qp;
    for (Eigen::Index sp_I = 0; sp_I < Nsp_I; ++sp_I) {Tmp_F2D_2sp_2sp.row(Nsp_I + sp_I) *= labels_S1D_sp[sp_I].twoSigma_I;}
    U2beta_F2D_2sp_2qp = Ag_LU.solve(Tmp_F2D_2sp_2sp.cast<doubleC>()).real();

    // V₂,full = [0,-V₂;ηV₂,0]; V₂β = R_yV₂,full.
    Tmp_F2D_2sp_2sp.setZero();
    Tmp_F2D_2sp_2sp.topRightCorner(Nsp_I, Nsp_I) = -V2Kramers_F2D_sp_qp;
    Tmp_F2D_2sp_2sp.bottomLeftCorner(Nsp_I, Nsp_I) = V2Kramers_F2D_sp_qp;
    for (Eigen::Index sp_I = 0; sp_I < Nsp_I; ++sp_I) {Tmp_F2D_2sp_2sp.row(Nsp_I + sp_I) *= labels_S1D_sp[sp_I].twoSigma_I;}
    V2beta_F2D_2sp_2qp.noalias() = Ry_F2D_2sp_2sp * Tmp_F2D_2sp_2sp;

    // U₁,fullᵀ = diag(U₁ᵀ,U₁ᵀη).
    Tmp_F2D_2sp_2sp.setZero();
    Tmp_F2D_2sp_2sp.topLeftCorner(Nsp_I, Nsp_I) = U1Kramers_F2D_sp_qp.transpose();
    Tmp_F2D_2sp_2sp.bottomRightCorner(Nsp_I, Nsp_I) = U1Kramers_F2D_sp_qp.transpose();
    for (Eigen::Index sp_I = 0; sp_I < Nsp_I; ++sp_I) {Tmp_F2D_2sp_2sp.col(Nsp_I + sp_I) *= labels_S1D_sp[sp_I].twoSigma_I;}

    // A_Uβ = U₁,fullᵀU₂β.
    U1TU2beta_F2D_2qp_2qp.noalias() = Tmp_F2D_2sp_2sp * U2beta_F2D_2sp_2qp;

    // V₁,fullᵀ = [0,V₁ᵀη;-V₁ᵀ,0].
    Tmp_F2D_2sp_2sp.setZero();
    Tmp_F2D_2sp_2sp.topRightCorner(Nsp_I, Nsp_I) = V1Kramers_F2D_sp_qp.transpose();
    Tmp_F2D_2sp_2sp.bottomLeftCorner(Nsp_I, Nsp_I) = -V1Kramers_F2D_sp_qp.transpose();
    for (Eigen::Index sp_I = 0; sp_I < Nsp_I; ++sp_I) {Tmp_F2D_2sp_2sp.col(Nsp_I + sp_I) *= labels_S1D_sp[sp_I].twoSigma_I;}

    // A_Vβ = V₁,fullᵀV₂β.
    V1TV2beta_F2D_2qp_2qp.noalias() = Tmp_F2D_2sp_2sp * V2beta_F2D_2sp_2qp;

    // Cache validity follows successful β preparation.
    betaCurrent_I = beta_I;
}

inline void AxialProjection::update_densities(int beta_I, int phi_I) {
    // Invalidate overlap before evaluation.
    overlap_C = doubleC(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());

    const Eigen::Index Nsp_I = static_cast<Eigen::Index>(labels_S1D_sp.size());
    assert(signDetU1_I != 0 && signDetU2_I != 0);
    assert(beta_I >= 0 && beta_I < detRySimplexP_C1D_beta.size());
    assert(Nphi_I > 0 && phi_I >= 0 && phi_I < Nphi_I);

    if (beta_I != betaCurrent_I) {update_AUV(beta_I);}

    // φ_j = 2πj/N_φ.
    const double pi_F = std::acos(-1.0);
    const double phi_F = 2.0 * pi_F * phi_I / Nphi_I;
    const doubleC phase_C = std::polar(1.0, phi_F);

    // A(β,φ) = e^{-iφ}A_Uβ + e^{iφ}A_Vβ.
    Ag_C2D_2qp_2qp = std::conj(phase_C) * U1TU2beta_F2D_2qp_2qp.cast<doubleC>() + phase_C * V1TV2beta_F2D_2qp_2qp.cast<doubleC>();
    assert(Ag_C2D_2qp_2qp.allFinite());
    Ag_LU.compute(Ag_C2D_2qp_2qp);

    // A₊ = ½(A^{++}+A^{--}+iA^{+-}-iA^{-+}).
    Tmp_C2D_2sp_2sp.topLeftCorner(Nsp_I, Nsp_I) = 0.5 * (Ag_C2D_2qp_2qp.topLeftCorner(Nsp_I, Nsp_I) + Ag_C2D_2qp_2qp.bottomRightCorner(Nsp_I, Nsp_I) + doubleC(0, 1) * Ag_C2D_2qp_2qp.topRightCorner(Nsp_I, Nsp_I) - doubleC(0, 1) * Ag_C2D_2qp_2qp.bottomLeftCorner(Nsp_I, Nsp_I));
    assert(Tmp_C2D_2sp_2sp.topLeftCorner(Nsp_I, Nsp_I).allFinite());

    // n = s₁s₂ exp(iN_sp φ) det(d₊) det(A₊).
    overlap_C = static_cast<double>(signDetU1_I * signDetU2_I) * std::polar(1.0, static_cast<double>(Nsp_I) * phi_F) * detRySimplexP_C1D_beta(betaCurrent_I) * Tmp_C2D_2sp_2sp.topLeftCorner(Nsp_I, Nsp_I).determinant();
    assert(std::isfinite(overlap_C.real()) && std::isfinite(overlap_C.imag()));

    assert(Ag_LU.isInvertible());

    // V₁,fullᵀ = [0,V₁ᵀη;-V₁ᵀ,0].
    Tmp_C2D_2sp_2sp.setZero();
    Tmp_C2D_2sp_2sp.topRightCorner(Nsp_I, Nsp_I) = V1Kramers_F2D_sp_qp.transpose().cast<doubleC>();
    Tmp_C2D_2sp_2sp.bottomLeftCorner(Nsp_I, Nsp_I) = -V1Kramers_F2D_sp_qp.transpose().cast<doubleC>();
    for (Eigen::Index sp_I = 0; sp_I < Nsp_I; ++sp_I) {Tmp_C2D_2sp_2sp.col(Nsp_I + sp_I) *= labels_S1D_sp[sp_I].twoSigma_I;}

    // AX = V₁ᵀ; ρ = e^{iφ}V₂βX; κ̄ = e^{-iφ}U₂βX.
    Tmp_C2D_2sp_2sp = Ag_LU.solve(Tmp_C2D_2sp_2sp).eval();
    rho_C2D_2sp_2sp.noalias() = phase_C * (V2beta_F2D_2sp_2qp.cast<doubleC>() * Tmp_C2D_2sp_2sp);
    kappaBar_C2D_2sp_2sp.noalias() = std::conj(phase_C) * (U2beta_F2D_2sp_2qp.cast<doubleC>() * Tmp_C2D_2sp_2sp);

    // U₁,fullᵀ = diag(U₁ᵀ,U₁ᵀη).
    Tmp_C2D_2sp_2sp.setZero();
    Tmp_C2D_2sp_2sp.topLeftCorner(Nsp_I, Nsp_I) = U1Kramers_F2D_sp_qp.transpose().cast<doubleC>();
    Tmp_C2D_2sp_2sp.bottomRightCorner(Nsp_I, Nsp_I) = U1Kramers_F2D_sp_qp.transpose().cast<doubleC>();
    for (Eigen::Index sp_I = 0; sp_I < Nsp_I; ++sp_I) {Tmp_C2D_2sp_2sp.col(Nsp_I + sp_I) *= labels_S1D_sp[sp_I].twoSigma_I;}

    // AX = U₁ᵀ; κ = e^{iφ}V₂βX.
    Tmp_C2D_2sp_2sp = Ag_LU.solve(Tmp_C2D_2sp_2sp).eval();
    kappa_C2D_2sp_2sp.noalias() = phase_C * (V2beta_F2D_2sp_2qp.cast<doubleC>() * Tmp_C2D_2sp_2sp);
}
