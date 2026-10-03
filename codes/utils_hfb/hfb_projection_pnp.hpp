/**
 * @file    hfb_projection_pnp.hpp
 * @author  cailea
 * @date    2026-09-18
 * @brief   Particle-number projection of HFB configurations.
 */

#pragma once

#include <cassert>
#include <complex>
#include <functional>
#include <numbers>
#include <vector>

#include <Eigen/Core>
#include <unsupported/Eigen/CXX11/Tensor>

#include "hfb_pfaffian.hpp"

/**
 * @brief  Project HFB configurations using uniform gauge quadrature.
 * @math   Pᴺ ≈ (1/Nφ) Σφ exp[-iφ(N̂-N)].
 * @note   Shared basis; unitary D; nonzero reference overlaps; reusable outputs.
 */
class HFBProjectionPNP {
public:
    int Nsp_I = 0;
    int Ncfg1_I = 0;
    int Ncfg2_I = 0;

    int TargetN_I = 0;
    int Nphi_I = 0;

    Eigen::MatrixXcd U1_C2D_sp_qp1{};
    Eigen::MatrixXcd V1_C2D_sp_qp1{};
    Eigen::MatrixXcd U2_C2D_sp_qp2{};
    Eigen::MatrixXcd V2_C2D_sp_qp2{};
    Eigen::MatrixXcd D_C2D_sp_sp{};

    HFBPfaffian hfb_pfaffian;

private:
    Eigen::MatrixXcd result_C2D_cfg1_cfg2{};
    Eigen::Tensor<doubleC, 3, Eigen::ColMajor> result_C3D_sp_cfg1_cfg2{};
    Eigen::MatrixXcd DU2_C2D_sp_qp2{};
    Eigen::MatrixXcd DV2_C2D_sp_qp2{};
    Eigen::MatrixXcd gDU2_C2D_sp_qp2{};
    Eigen::MatrixXcd gDV2_C2D_sp_qp2{};


public:
    /**
     * @brief  Construct empty particle-number projection workspaces.
     * @math   Nsp = Ncfg1 = Ncfg2 = Nφ = 0.
     * @output Empty workspaces.
     */
    HFBProjectionPNP() = default;

    /**
     * @brief  Allocate projection and Pfaffian workspaces.
     * @math   φ_k = 2πk/Nφ; k = 0,…,Nφ-1.
     * @output Stored dimensions, target number, and allocated workspaces.
     */
    HFBProjectionPNP(int Nsp_I_, const std::vector<std::vector<int>>& config1_I2D_cfg1_cqp1_, const std::vector<std::vector<int>>& config2_I2D_cfg2_cqp2_, int TargetN_I_, int Nphi_I_)
    : hfb_pfaffian(Nsp_I_, config1_I2D_cfg1_cqp1_, config2_I2D_cfg2_cqp2_) {
        assert(TargetN_I_ >= 0 && TargetN_I_ <= Nsp_I_);
        assert(Nphi_I_ > 0);

        // (Nsp,config1,config2,N,Nφ) → fixed projection parameters.
        Nsp_I = Nsp_I_;
        Ncfg1_I = static_cast<int>(hfb_pfaffian.config1_I2D_cfg1_cqp1.size());
        Ncfg2_I = static_cast<int>(hfb_pfaffian.config2_I2D_cfg2_cqp2.size());
        TargetN_I = TargetN_I_;
        Nphi_I = Nphi_I_;

        // U₁,V₁,U₂,V₂ ∈ ℂ^{Nsp×Nsp}.
        U1_C2D_sp_qp1.resize(Nsp_I, Nsp_I);
        V1_C2D_sp_qp1.resize(Nsp_I, Nsp_I);
        U2_C2D_sp_qp2.resize(Nsp_I, Nsp_I);
        V2_C2D_sp_qp2.resize(Nsp_I, Nsp_I);

        D_C2D_sp_sp.resize(Nsp_I, Nsp_I);
        D_C2D_sp_sp.setIdentity();

        result_C2D_cfg1_cfg2.resize(Ncfg1_I, Ncfg2_I);
        result_C3D_sp_cfg1_cfg2.resize(Nsp_I, Ncfg1_I, Ncfg2_I);
        DU2_C2D_sp_qp2.resize(Nsp_I, Nsp_I);
        DV2_C2D_sp_qp2.resize(Nsp_I, Nsp_I);
        gDU2_C2D_sp_qp2.resize(Nsp_I, Nsp_I);
        gDV2_C2D_sp_qp2.resize(Nsp_I, Nsp_I);

    }

    /**
     * @brief  Store left and right Bogoliubov matrices.
     * @math   U₁,V₁,U₂,V₂ ∈ ℝ^{Nsp×Nsp} → complex caches.
     * @output Updated U₁,V₁,U₂,V₂.
     * @note   Requires canonical pairs and invertible U₁,U₂.
     */
    void update_UV(const Eigen::MatrixXd& U1_F2D_sp_qp1_, const Eigen::MatrixXd& V1_F2D_sp_qp1_, const Eigen::MatrixXd& U2_F2D_sp_qp2_, const Eigen::MatrixXd& V2_F2D_sp_qp2_);

    /**
     * @brief  Integrate overlaps using uniform gauge quadrature.
     * @math   result = (1/Nφ) Σφ eⁱᴺφ ⟨Φ₁|e⁻ⁱφᴺ̂ D̂|Φ₂⟩.
     * @output Const reference to result_C2D_cfg1_cfg2.
     */
    const Eigen::MatrixXcd& calc_overlap_pnp();

    /**
     * @brief  Integrate fermion kernels using uniform gauge quadrature.
     * @math   result_αab = (1/Nφ) Σφ eⁱᴺφ ⟨Φ₁;a|xα e⁻ⁱφᴺ̂ D̂|Φ₂;b⟩.
     * @output Const reference to shared result_C3D_sp_cfg1_cfg2.
     * @note   creation_B: true → xα = cα†; false → xα = cα.
     */
    const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& calc_one_fermion_pnp(bool creation_B);

    /**
     * @brief  Integrate one-body kernels using uniform gauge quadrature.
     * @math   result = (1/Nφ) Σφ eⁱᴺφ Σ_ij OneBody_ij OBTD(i,j;φ,D).
     * @output Const reference to result_C2D_cfg1_cfg2.
     */
    const Eigen::MatrixXcd& calc_one_body_pnp(const Eigen::MatrixXcd& OneBody_C2D_sp_sp);

    /**
     * @brief  Integrate Q/P kernels using Pfaffians and gauge quadrature.
     * @math   K = Σi gQi(Qi†Qi − one-body contraction) + Σi gPi Pi†Pi.
     * @math   P† = ½Σab Pab c†a c†b; Pᵀ = −P.
     * @output Weighted kernel in result_C2D_cfg1_cfg2.
     * @note   gQ,gP include signs; pairing ¼ is internal.
     */
    const Eigen::MatrixXcd& calc_two_body_pnp(const Eigen::VectorXd& gQ_F1D_i, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& Q_C3D_sp_sp_i, const Eigen::VectorXd& gP_F1D_i, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& P_C3D_sp_sp_i);

protected:
    /**
     * @brief  Integrate kernels using uniform gauge quadrature.
     * @math   result = (1/Nφ) Σφ exp(iNφ) K(φ,D).
     * @output Filled caller-provided contiguous result view.
     * @note   accumulate_Func(phi_I) accumulates exp(iNφ) K(φ,D), including vacuum overlap.
     */
    void integrate_pnp(Eigen::Ref<Eigen::VectorXcd> result_C1D_element, const std::function<void(int)>& accumulate_Func);

};

inline void HFBProjectionPNP::integrate_pnp(Eigen::Ref<Eigen::VectorXcd> result_C1D_element, const std::function<void(int)>& accumulate_Func) {
    result_C1D_element.setZero();
    assert(D_C2D_sp_sp.rows() == Nsp_I && D_C2D_sp_sp.cols() == Nsp_I && D_C2D_sp_sp.allFinite());
    // DU₂ = D U₂; DV₂ = D* V₂; independent of φ.
    DU2_C2D_sp_qp2.noalias() = D_C2D_sp_sp * U2_C2D_sp_qp2;
    DV2_C2D_sp_qp2.noalias() = D_C2D_sp_sp.conjugate() * V2_C2D_sp_qp2;
    // Contractions are prepared once per gauge node.
    for (int phi_I = 0; phi_I < Nphi_I; ++phi_I) {
        const double phi_F = 2.0 * std::numbers::pi * phi_I / Nphi_I;
        // U₂φ = exp(-iφ) DU₂; V₂φ = exp(iφ) DV₂.
        const doubleC gauge_C = std::exp(doubleC(0.0, -phi_F));
        gDU2_C2D_sp_qp2 = gauge_C * DU2_C2D_sp_qp2;
        gDV2_C2D_sp_qp2 = std::conj(gauge_C) * DV2_C2D_sp_qp2;
        hfb_pfaffian.update_contractions(U1_C2D_sp_qp1, V1_C2D_sp_qp1, gDU2_C2D_sp_qp2, gDV2_C2D_sp_qp2, doubleC(1.0, 0.0));
        accumulate_Func(phi_I);
    }
    result_C1D_element /= Nphi_I;
}

inline void HFBProjectionPNP::update_UV(const Eigen::MatrixXd& U1_F2D_sp_qp1_, const Eigen::MatrixXd& V1_F2D_sp_qp1_, const Eigen::MatrixXd& U2_F2D_sp_qp2_, const Eigen::MatrixXd& V2_F2D_sp_qp2_) {
    assert(U1_F2D_sp_qp1_.rows() == Nsp_I && U1_F2D_sp_qp1_.cols() == Nsp_I);
    assert(V1_F2D_sp_qp1_.rows() == Nsp_I && V1_F2D_sp_qp1_.cols() == Nsp_I);
    assert(U2_F2D_sp_qp2_.rows() == Nsp_I && U2_F2D_sp_qp2_.cols() == Nsp_I);
    assert(V2_F2D_sp_qp2_.rows() == Nsp_I && V2_F2D_sp_qp2_.cols() == Nsp_I);
    // (U₁,V₁,U₂,V₂) ∈ ℝ → ℂ.
    U1_C2D_sp_qp1 = U1_F2D_sp_qp1_.cast<doubleC>();
    V1_C2D_sp_qp1 = V1_F2D_sp_qp1_.cast<doubleC>();
    U2_C2D_sp_qp2 = U2_F2D_sp_qp2_.cast<doubleC>();
    V2_C2D_sp_qp2 = V2_F2D_sp_qp2_.cast<doubleC>();
}

inline const Eigen::MatrixXcd& HFBProjectionPNP::calc_overlap_pnp() {
    // Σφ exp(iNφ) O(φ,D).
    Eigen::Map<Eigen::VectorXcd> result_C1D_element(result_C2D_cfg1_cfg2.data(), result_C2D_cfg1_cfg2.size());
    integrate_pnp(result_C1D_element, [&](int phi_I) {
        const double phi_F = 2.0 * std::numbers::pi * phi_I / Nphi_I;
        const doubleC factorPhi_C = std::exp(doubleC(0.0, TargetN_I * phi_F));
        result_C2D_cfg1_cfg2 += factorPhi_C * hfb_pfaffian.calc_overlap();
    });
    return result_C2D_cfg1_cfg2;
}

inline const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& HFBProjectionPNP::calc_one_fermion_pnp(bool creation_B) {
    Eigen::Map<Eigen::VectorXcd> result_C1D_element(result_C3D_sp_cfg1_cfg2.data(), result_C3D_sp_cfg1_cfg2.size());
    // (1/Nφ) Σφ exp(iNφ) ⟨xα exp(-iφN̂) D̂⟩; all α together.
    integrate_pnp(result_C1D_element, [&](int phi_I) {
        const double phi_F = 2.0 * std::numbers::pi * phi_I / Nphi_I;
        const doubleC factorPhi_C = std::exp(doubleC(0.0, TargetN_I * phi_F));
        const auto& kernel_C3D_sp_cfg1_cfg2 = hfb_pfaffian.calc_one_fermion(creation_B);
        const Eigen::Map<const Eigen::VectorXcd> kernel_C1D_element(kernel_C3D_sp_cfg1_cfg2.data(), kernel_C3D_sp_cfg1_cfg2.size());
        result_C1D_element += factorPhi_C * kernel_C1D_element;
    });
    return result_C3D_sp_cfg1_cfg2;
}

inline const Eigen::MatrixXcd& HFBProjectionPNP::calc_one_body_pnp(const Eigen::MatrixXcd& OneBody_C2D_sp_sp) {
    Eigen::Map<Eigen::VectorXcd> result_C1D_element(result_C2D_cfg1_cfg2.data(), result_C2D_cfg1_cfg2.size());
    integrate_pnp(result_C1D_element, [&](int phi_I) {
        const double phi_F = 2.0 * std::numbers::pi * phi_I / Nphi_I;
        const doubleC factorPhi_C = std::exp(doubleC(0.0, TargetN_I * phi_F));
        result_C2D_cfg1_cfg2 += factorPhi_C * hfb_pfaffian.calc_one_body(OneBody_C2D_sp_sp);
    });
    return result_C2D_cfg1_cfg2;
}

inline const Eigen::MatrixXcd& HFBProjectionPNP::calc_two_body_pnp(const Eigen::VectorXd& gQ_F1D_i, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& Q_C3D_sp_sp_i, const Eigen::VectorXd& gP_F1D_i, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& P_C3D_sp_sp_i) {
    assert(Q_C3D_sp_sp_i.dimension(2) == gQ_F1D_i.size());
    assert(P_C3D_sp_sp_i.dimension(2) == gP_F1D_i.size());
    assert(gQ_F1D_i.size() == 0 || (Q_C3D_sp_sp_i.dimension(0) == Nsp_I && Q_C3D_sp_sp_i.dimension(1) == Nsp_I));
    assert(gP_F1D_i.size() == 0 || (P_C3D_sp_sp_i.dimension(0) == Nsp_I && P_C3D_sp_sp_i.dimension(1) == Nsp_I));

    // hQ = Σi gQi Qi†Qi.
    Eigen::MatrixXcd hQ_C2D_sp_sp = Eigen::MatrixXcd::Zero(Nsp_I, Nsp_I);
    for (Eigen::Index iQ_I = 0; iQ_I < gQ_F1D_i.size(); ++iQ_I) {
        const Eigen::Index offset_I = iQ_I * Nsp_I * Nsp_I;
        const Eigen::Map<const Eigen::MatrixXcd> Q_C2D_sp_sp(Q_C3D_sp_sp_i.data() + offset_I, Nsp_I, Nsp_I);
        hQ_C2D_sp_sp.noalias() += gQ_F1D_i(iQ_I) * Q_C2D_sp_sp.adjoint() * Q_C2D_sp_sp;
    }

    // Internal contractions and antisymmetric pair factors.
    Eigen::MatrixXcd x1x2_C2D_sp_sp(Nsp_I, Nsp_I);
    Eigen::MatrixXcd x1x3_C2D_sp_sp(Nsp_I, Nsp_I);
    Eigen::MatrixXcd x1x4_C2D_sp_sp(Nsp_I, Nsp_I);
    Eigen::MatrixXcd x2x3_C2D_sp_sp(Nsp_I, Nsp_I);
    Eigen::MatrixXcd x2x4_C2D_sp_sp(Nsp_I, Nsp_I);
    Eigen::MatrixXcd x3x4_C2D_sp_sp(Nsp_I, Nsp_I);

    Eigen::MatrixXcd eta1_C2D_qp_sp(2 * Nsp_I, Nsp_I);
    Eigen::MatrixXcd eta2_C2D_qp_sp(2 * Nsp_I, Nsp_I);
    Eigen::MatrixXcd eta3_C2D_qp_sp(2 * Nsp_I, Nsp_I);
    Eigen::MatrixXcd eta4_C2D_qp_sp(2 * Nsp_I, Nsp_I);

    doubleC T0_C = 0.0;
    Eigen::MatrixXcd T2_C2D_qp_qp(2 * Nsp_I, 2 * Nsp_I);
    Eigen::MatrixXcd PI14_C2D_qp_qp(2 * Nsp_I, 2 * Nsp_I);
    Eigen::MatrixXcd PI23_C2D_qp_qp(2 * Nsp_I, 2 * Nsp_I);
    Eigen::MatrixXcd PI12_C2D_qp_qp(2 * Nsp_I, 2 * Nsp_I);
    Eigen::MatrixXcd PI34_C2D_qp_qp(2 * Nsp_I, 2 * Nsp_I);

    const auto calc_constract_Q = [&](const Eigen::Ref<const Eigen::MatrixXcd>& Q_C2D_sp_sp) {
        // x₁x₂ = Q*κ̄Q; x₁x₃ = Q*ρ; x₂x₃ = Qᵀρ.
        x1x2_C2D_sp_sp.noalias() = Q_C2D_sp_sp.conjugate() * hfb_pfaffian.SpDagSpDag_C2D_sp_sp * Q_C2D_sp_sp;
        x1x3_C2D_sp_sp.noalias() = Q_C2D_sp_sp.conjugate() * hfb_pfaffian.SpDagSp_C2D_sp_sp;
        x2x3_C2D_sp_sp.noalias() = Q_C2D_sp_sp.transpose() * hfb_pfaffian.SpDagSp_C2D_sp_sp;

        // ηⁱ: left ⟨β₁xᵢ⟩; right −⟨xᵢβ₂†⟩ᵀ.
        eta1_C2D_qp_sp.topRows(Nsp_I).noalias() = hfb_pfaffian.Qp1SpDag_C2D_qp1_sp * Q_C2D_sp_sp.adjoint();
        eta1_C2D_qp_sp.bottomRows(Nsp_I).noalias() = -hfb_pfaffian.SpDagQp2Dag_C2D_sp_qp2.transpose() * Q_C2D_sp_sp.adjoint();
        eta2_C2D_qp_sp.topRows(Nsp_I).noalias() = hfb_pfaffian.Qp1SpDag_C2D_qp1_sp * Q_C2D_sp_sp;
        eta2_C2D_qp_sp.bottomRows(Nsp_I).noalias() = -hfb_pfaffian.SpDagQp2Dag_C2D_sp_qp2.transpose() * Q_C2D_sp_sp;
        eta3_C2D_qp_sp.topRows(Nsp_I) = hfb_pfaffian.Qp1Sp_C2D_qp1_sp;
        eta3_C2D_qp_sp.bottomRows(Nsp_I) = -hfb_pfaffian.SpQp2Dag_C2D_sp_qp2.transpose();
        eta4_C2D_qp_sp = eta3_C2D_qp_sp;

        // T⁰ = Σbd (C01 C23 − C02 C13 + C03 C12).
        T0_C = (x1x2_C2D_sp_sp.array() * hfb_pfaffian.SpSp_C2D_sp_sp.transpose().array()).sum()
            - (x1x3_C2D_sp_sp.array() * x2x3_C2D_sp_sp.transpose().array()).sum()
            + x1x3_C2D_sp_sp.trace() * x2x3_C2D_sp_sp.trace();

        // T² = antisym[−η¹ C34 (η²)ᵀ + η¹ C24 (η³)ᵀ − ⋯].
        T2_C2D_qp_qp.noalias() = -eta1_C2D_qp_sp * hfb_pfaffian.SpSp_C2D_sp_sp.transpose() * eta2_C2D_qp_sp.transpose();
        T2_C2D_qp_qp.noalias() += eta1_C2D_qp_sp * x2x3_C2D_sp_sp.transpose() * eta3_C2D_qp_sp.transpose();
        T2_C2D_qp_qp.noalias() -= x2x3_C2D_sp_sp.trace() * eta1_C2D_qp_sp * eta4_C2D_qp_sp.transpose();
        T2_C2D_qp_qp.noalias() -= x1x3_C2D_sp_sp.trace() * eta2_C2D_qp_sp * eta3_C2D_qp_sp.transpose();
        T2_C2D_qp_qp.noalias() += eta2_C2D_qp_sp * x1x3_C2D_sp_sp.transpose() * eta4_C2D_qp_sp.transpose();
        T2_C2D_qp_qp.noalias() -= eta3_C2D_qp_sp * x1x2_C2D_sp_sp.transpose() * eta4_C2D_qp_sp.transpose();

        // T⁴ = Π¹⁴ ∧ Π²³.
        PI14_C2D_qp_qp.noalias() = eta1_C2D_qp_sp * eta4_C2D_qp_sp.transpose();
        PI23_C2D_qp_qp.noalias() = eta2_C2D_qp_sp * eta3_C2D_qp_sp.transpose();
        PI14_C2D_qp_qp = (PI14_C2D_qp_qp - PI14_C2D_qp_qp.transpose()).eval();
        PI23_C2D_qp_qp = (PI23_C2D_qp_qp - PI23_C2D_qp_qp.transpose()).eval();
        T2_C2D_qp_qp = (T2_C2D_qp_qp - T2_C2D_qp_qp.transpose()).eval();
    };

    const auto calc_constract_P = [&](const Eigen::Ref<const Eigen::MatrixXcd>& P_C2D_sp_sp) {
        // x₁x₂ = Pᵀκ̄; x₁x₃ = Pᵀρ; x₁x₄ = PᵀρP*.
        x1x2_C2D_sp_sp.noalias() = P_C2D_sp_sp.transpose() * hfb_pfaffian.SpDagSpDag_C2D_sp_sp;
        x1x3_C2D_sp_sp.noalias() = P_C2D_sp_sp.transpose() * hfb_pfaffian.SpDagSp_C2D_sp_sp;
        x1x4_C2D_sp_sp.noalias() = x1x3_C2D_sp_sp * P_C2D_sp_sp.conjugate();
        x2x3_C2D_sp_sp = hfb_pfaffian.SpDagSp_C2D_sp_sp;
        x2x4_C2D_sp_sp.noalias() = hfb_pfaffian.SpDagSp_C2D_sp_sp * P_C2D_sp_sp.conjugate();
        x3x4_C2D_sp_sp.noalias() = hfb_pfaffian.SpSp_C2D_sp_sp * P_C2D_sp_sp.conjugate();

        // ηⁱ: left ⟨β₁xᵢ⟩; right −⟨xᵢβ₂†⟩ᵀ.
        eta1_C2D_qp_sp.topRows(Nsp_I).noalias() = hfb_pfaffian.Qp1SpDag_C2D_qp1_sp * P_C2D_sp_sp;
        eta1_C2D_qp_sp.bottomRows(Nsp_I).noalias() = -hfb_pfaffian.SpDagQp2Dag_C2D_sp_qp2.transpose() * P_C2D_sp_sp;
        eta2_C2D_qp_sp.topRows(Nsp_I) = hfb_pfaffian.Qp1SpDag_C2D_qp1_sp;
        eta2_C2D_qp_sp.bottomRows(Nsp_I) = -hfb_pfaffian.SpDagQp2Dag_C2D_sp_qp2.transpose();
        eta3_C2D_qp_sp.topRows(Nsp_I) = hfb_pfaffian.Qp1Sp_C2D_qp1_sp;
        eta3_C2D_qp_sp.bottomRows(Nsp_I) = -hfb_pfaffian.SpQp2Dag_C2D_sp_qp2.transpose();
        eta4_C2D_qp_sp.topRows(Nsp_I).noalias() = hfb_pfaffian.Qp1Sp_C2D_qp1_sp * P_C2D_sp_sp.conjugate();
        eta4_C2D_qp_sp.bottomRows(Nsp_I).noalias() = -hfb_pfaffian.SpQp2Dag_C2D_sp_qp2.transpose() * P_C2D_sp_sp.conjugate();

        // T⁰ = Σbd (C01 C23 − C02 C13 + C03 C12).
        T0_C = x1x2_C2D_sp_sp.trace() * x3x4_C2D_sp_sp.trace()
            - (x1x3_C2D_sp_sp.array() * x2x4_C2D_sp_sp.array()).sum()
            + (x1x4_C2D_sp_sp.array() * x2x3_C2D_sp_sp.array()).sum();

        // T²: pair two external z with two insertions.
        T2_C2D_qp_qp.noalias() = -x3x4_C2D_sp_sp.trace() * eta1_C2D_qp_sp * eta2_C2D_qp_sp.transpose();
        T2_C2D_qp_qp.noalias() += eta1_C2D_qp_sp * x2x4_C2D_sp_sp * eta3_C2D_qp_sp.transpose();
        T2_C2D_qp_qp.noalias() -= eta1_C2D_qp_sp * x2x3_C2D_sp_sp * eta4_C2D_qp_sp.transpose();
        T2_C2D_qp_qp.noalias() -= eta2_C2D_qp_sp * x1x4_C2D_sp_sp * eta3_C2D_qp_sp.transpose();
        T2_C2D_qp_qp.noalias() += eta2_C2D_qp_sp * x1x3_C2D_sp_sp * eta4_C2D_qp_sp.transpose();
        T2_C2D_qp_qp.noalias() -= x1x2_C2D_sp_sp.trace() * eta3_C2D_qp_sp * eta4_C2D_qp_sp.transpose();

        T2_C2D_qp_qp = (T2_C2D_qp_qp - T2_C2D_qp_qp.transpose()).eval();

        // T⁴ = Π¹² ∧ Π³⁴.
        PI12_C2D_qp_qp.noalias() = eta1_C2D_qp_sp * eta2_C2D_qp_sp.transpose();
        PI34_C2D_qp_qp.noalias() = eta3_C2D_qp_sp * eta4_C2D_qp_sp.transpose();
        PI12_C2D_qp_qp = (PI12_C2D_qp_qp - PI12_C2D_qp_qp.transpose()).eval();
        PI34_C2D_qp_qp = (PI34_C2D_qp_qp - PI34_C2D_qp_qp.transpose()).eval();
    };

    // Ncqp12 ≤ 2Nsp; reuse selected pair matrices.
    std::vector<int> qp_I1D_cqp12{};
    qp_I1D_cqp12.reserve(2 * Nsp_I);
    Eigen::VectorXcd PI1workspace_C1D_element(4 * Nsp_I * Nsp_I);
    Eigen::VectorXcd PI2workspace_C1D_element(4 * Nsp_I * Nsp_I);

    /**
     * @brief Assemble kernels using Pfaffian minors.
     * @math K = g⟨Φ₁|Φ₂⟩(F⁰T⁰ − ½F²T² + ¼F⁴Π₁Π₂).
     * @output Accumulated result_C2D_cfg1_cfg2.
     */
    const auto calc_kernel024 = [&](double g_F, const Eigen::MatrixXcd& PI1_C2D_qp_qp, const Eigen::MatrixXcd& PI2_C2D_qp_qp, doubleC factorPhi_C) {
        const doubleC factor_C = factorPhi_C * hfb_pfaffian.overlap_C * g_F;
        for (int cfg2_I = 0; cfg2_I < Ncfg2_I; ++cfg2_I) {
            const auto& config2_I1D_cqp = hfb_pfaffian.config2_I2D_cfg2_cqp2[cfg2_I];
            for (int cfg1_I = 0; cfg1_I < Ncfg1_I; ++cfg1_I) {
                const auto& config1_I1D_cqp = hfb_pfaffian.config1_I2D_cfg1_cqp1[cfg1_I];
                const int Ncqp1_I = static_cast<int>(config1_I1D_cqp.size());
                const Eigen::Index Ncqp12_I = Ncqp1_I + static_cast<Eigen::Index>(config2_I1D_cqp.size());
                if (Ncqp12_I % 2 != 0) { continue; }
                // Left indices reversed; right indices shifted by Nsp.
                qp_I1D_cqp12.clear();
                for (int cqp1_I = Ncqp1_I - 1; cqp1_I >= 0; --cqp1_I) { qp_I1D_cqp12.push_back(config1_I1D_cqp[cqp1_I]); }
                for (const int qp_I : config2_I1D_cqp) { qp_I1D_cqp12.push_back(Nsp_I + qp_I); }
                const auto& F2_C2D_cqp12_cqp12 = hfb_pfaffian.F2_C4D_cfg1_cfg2_cqp12_cqp12[cfg1_I][cfg2_I];
                // K⁰ + K² = F⁰T⁰ − ½Σij F²ij T²ij.
                doubleC kernel_C = hfb_pfaffian.F0_C2D_cfg1_cfg2[cfg1_I][cfg2_I] * T0_C
                    - 0.5 * F2_C2D_cqp12_cqp12.cwiseProduct(T2_C2D_qp_qp(qp_I1D_cqp12, qp_I1D_cqp12)).sum();

                // K⁴ = ¼ vec(Π₁)ᵀ F⁴_(ij,kl) vec(Π₂).
                if (Ncqp12_I >= 4) {
                    Eigen::Map<Eigen::MatrixXcd> PI1_C2D_cqp12_cqp12(PI1workspace_C1D_element.data(), Ncqp12_I, Ncqp12_I);
                    Eigen::Map<Eigen::MatrixXcd> PI2_C2D_cqp12_cqp12(PI2workspace_C1D_element.data(), Ncqp12_I, Ncqp12_I);
                    PI1_C2D_cqp12_cqp12 = PI1_C2D_qp_qp(qp_I1D_cqp12, qp_I1D_cqp12);
                    PI2_C2D_cqp12_cqp12 = PI2_C2D_qp_qp(qp_I1D_cqp12, qp_I1D_cqp12);
                    // ColMajor: (i,j) → i+nj; (k,l) → k+nl.
                    const Eigen::Map<const Eigen::MatrixXcd> F4_C2D_cqp12cqp12_cqp12cqp12(hfb_pfaffian.F4_C6D_cfg1_cfg2_cqp12_cqp12_cqp12_cqp12[cfg1_I][cfg2_I].data(), Ncqp12_I * Ncqp12_I, Ncqp12_I * Ncqp12_I);
                    const Eigen::Map<const Eigen::VectorXcd> PI1_C1D_cqp12cqp12(PI1_C2D_cqp12_cqp12.data(), PI1_C2D_cqp12_cqp12.size());
                    const Eigen::Map<const Eigen::VectorXcd> PI2_C1D_cqp12cqp12(PI2_C2D_cqp12_cqp12.data(), PI2_C2D_cqp12_cqp12.size());
                    kernel_C += 0.25 * (PI1_C1D_cqp12cqp12.transpose() * F4_C2D_cqp12cqp12_cqp12cqp12 * PI2_C1D_cqp12cqp12).value();
                }
                result_C2D_cfg1_cfg2(cfg1_I, cfg2_I) += factor_C * kernel_C;
            }
        }
    };

    Eigen::Map<Eigen::VectorXcd> result_C1D_element(result_C2D_cfg1_cfg2.data(), result_C2D_cfg1_cfg2.size());
    integrate_pnp(result_C1D_element, [&](int phi_I) {
        const double phi_F = 2.0 * std::numbers::pi * phi_I / Nphi_I;
        const doubleC factorPhi_C = std::exp(doubleC(0.0, TargetN_I * phi_F));

        // Q: gQ (F⁰T⁰ − ½F²T² + ¼F⁴Π¹⁴Π²³).
        for (Eigen::Index iQ_I = 0; iQ_I < gQ_F1D_i.size(); ++iQ_I) {
            if (gQ_F1D_i(iQ_I) == 0.0) { continue; }
            const Eigen::Index offset_I = iQ_I * Nsp_I * Nsp_I;
            const Eigen::Map<const Eigen::MatrixXcd> Q_C2D_sp_sp(Q_C3D_sp_sp_i.data() + offset_I, Nsp_I, Nsp_I);
            calc_constract_Q(Q_C2D_sp_sp);
            calc_kernel024(gQ_F1D_i(iQ_I), PI14_C2D_qp_qp, PI23_C2D_qp_qp, factorPhi_C);
        }

        // P: ¼gP (F⁰T⁰ − ½F²T² + ¼F⁴Π¹²Π³⁴).
        for (Eigen::Index iP_I = 0; iP_I < gP_F1D_i.size(); ++iP_I) {
            if (gP_F1D_i(iP_I) == 0.0) { continue; }
            const Eigen::Index offset_I = iP_I * Nsp_I * Nsp_I;
            const Eigen::Map<const Eigen::MatrixXcd> P_C2D_sp_sp(P_C3D_sp_sp_i.data() + offset_I, Nsp_I, Nsp_I);
            calc_constract_P(P_C2D_sp_sp);
            calc_kernel024(0.25 * gP_F1D_i(iP_I), PI12_C2D_qp_qp, PI34_C2D_qp_qp, factorPhi_C);
        }
    });
    return result_C2D_cfg1_cfg2;
}
