/**
 * @file    hfb_projection_pnp.hpp
 * @author  cailea
 * @date    2026-09-18
 * @brief   Particle-number projection of HFB configurations.
 */

#pragma once

#include <cassert>
#include <complex>
#include <numbers>
#include <vector>

#include <Eigen/Core>
#include <unsupported/Eigen/CXX11/Tensor>

#include "hfb_pfaffian.hpp"

/**
 * @brief  Project HFB configurations using uniform gauge quadrature.
 * @math   Pᴺ ≈ (1/Nφ) Σφ exp[-iφ(N̂-N)].
 * @note   U,V,D share one ordered single-particle basis.
 * @note   Initialize U,V before kernel calculations.
 * @note   D is unitary; defaults to the identity.
 * @note   Quadrature nodes require nonzero reference overlaps.
 * @note   Returned references alias reusable output workspaces.
 */
class HFBProjectionPNP {
public:
    int TargetN_I = 0;

    int Nsp_I = 0;
    int Ncfg1_I = 0;
    int Ncfg2_I = 0;

    int Nphi_I = 0;

    Eigen::MatrixXcd U1_C2D_sp_qp1{};
    Eigen::MatrixXcd V1_C2D_sp_qp1{};
    Eigen::MatrixXcd U2_C2D_sp_qp2{};
    Eigen::MatrixXcd V2_C2D_sp_qp2{};
    Eigen::MatrixXcd D_C2D_sp_sp{};

    HFBPfaffian hfb_pfaffian;

protected:
    Eigen::VectorXcd x1_C1D_sp{};
    Eigen::VectorXcd x2_C1D_sp{};
    Eigen::VectorXcd x3_C1D_sp{};
    Eigen::VectorXcd x4_C1D_sp{};

private:
    Eigen::MatrixXcd result_C2D_cfg1_cfg2{};
    Eigen::Tensor<doubleC, 3, Eigen::ColMajor> result_C3D_cfg1_cfg2_sp{};
    Eigen::MatrixXcd DU2_C2D_sp_qp2{};
    Eigen::MatrixXcd DV2_C2D_sp_qp2{};
    Eigen::MatrixXcd gDU2_C2D_sp_qp2{};
    Eigen::MatrixXcd gDV2_C2D_sp_qp2{};

    // Q contractions; QDag = Qᵀ for real input.
    Eigen::MatrixXcd Qp1SpDagQDag_C2D_qp1_sp{};
    Eigen::MatrixXcd Qp1SpDagQ_C2D_qp1_sp{};
    Eigen::MatrixXcd QSpDagQp2Dag_C2D_sp_qp2{};
    Eigen::MatrixXcd QDagSpDagQp2Dag_C2D_sp_qp2{};
    Eigen::MatrixXcd QSpDagSpDagQ_C2D_sp_sp{};
    Eigen::MatrixXcd QSpDagSp_C2D_sp_sp{};
    Eigen::MatrixXcd QDagSpDagSp_C2D_sp_sp{};
    Eigen::MatrixXcd QSpDagSpDag_C2D_sp_sp{};

    // P contractions; PDag = Pᵀ for real input.
    Eigen::MatrixXcd Qp1SpDagP_C2D_qp1_sp{};
    Eigen::MatrixXcd Qp1SpP_C2D_qp1_sp{};
    Eigen::MatrixXcd PDagSpDagQp2Dag_C2D_sp_qp2{};
    Eigen::MatrixXcd PDagSpQp2Dag_C2D_sp_qp2{};
    Eigen::MatrixXcd PDagSpDagSpDag_C2D_sp_sp{};
    Eigen::MatrixXcd PDagSpDagSp_C2D_sp_sp{};
    Eigen::MatrixXcd PDagSpDagSpP_C2D_sp_sp{};
    Eigen::MatrixXcd SpDagSpP_C2D_sp_sp{};
    Eigen::MatrixXcd SpSpP_C2D_sp_sp{};

    // Shared kernel workspaces.
    Eigen::MatrixXcd Qp1X_C2D_qp1_X{};
    Eigen::MatrixXcd XQp2Dag_C2D_X_qp2{};
    Eigen::MatrixXcd XX_C2D_X_X{};
    Eigen::MatrixXcd kernel_C2D_cfg1_cfg2{};

public:
    /**
     * @brief  Allocate projection and Pfaffian workspaces.
     * @math   φ_k = 2πk/Nφ; k = 0,…,Nφ-1.
     * @output Stored dimensions, target number, and allocated workspaces.
     */
    HFBProjectionPNP(int TargetN_I_, int Nsp_I_, const std::vector<std::vector<int>>& config1_I2D_cfg1_cqp1_, const std::vector<std::vector<int>>& config2_I2D_cfg2_cqp2_, int Nphi_I_)
    : hfb_pfaffian(Nsp_I_, config1_I2D_cfg1_cqp1_, config2_I2D_cfg2_cqp2_) {
        assert(TargetN_I_ >= 0 && TargetN_I_ <= Nsp_I_);
        assert(Nphi_I_ > 0);

        // (N,Nsp,config1,config2,Nφ) → fixed projection parameters.
        TargetN_I = TargetN_I_;
        Nsp_I = Nsp_I_;
        Ncfg1_I = static_cast<int>(hfb_pfaffian.config1_I2D_cfg1_cqp1.size());
        Ncfg2_I = static_cast<int>(hfb_pfaffian.config2_I2D_cfg2_cqp2.size());
        Nphi_I = Nphi_I_;

        // U₁,V₁,U₂,V₂ ∈ ℂ^{Nsp×Nsp}.
        U1_C2D_sp_qp1.resize(Nsp_I, Nsp_I);
        V1_C2D_sp_qp1.resize(Nsp_I, Nsp_I);
        U2_C2D_sp_qp2.resize(Nsp_I, Nsp_I);
        V2_C2D_sp_qp2.resize(Nsp_I, Nsp_I);

        D_C2D_sp_sp.resize(Nsp_I, Nsp_I);
        D_C2D_sp_sp.setIdentity();

        result_C2D_cfg1_cfg2.resize(Ncfg1_I, Ncfg2_I);
        result_C3D_cfg1_cfg2_sp.resize(Ncfg1_I, Ncfg2_I, Nsp_I);
        DU2_C2D_sp_qp2.resize(Nsp_I, Nsp_I);
        DV2_C2D_sp_qp2.resize(Nsp_I, Nsp_I);
        gDU2_C2D_sp_qp2.resize(Nsp_I, Nsp_I);
        gDV2_C2D_sp_qp2.resize(Nsp_I, Nsp_I);

        // Q contractions.
        Qp1SpDagQDag_C2D_qp1_sp.resize(Nsp_I, Nsp_I);
        Qp1SpDagQ_C2D_qp1_sp.resize(Nsp_I, Nsp_I);
        QSpDagQp2Dag_C2D_sp_qp2.resize(Nsp_I, Nsp_I);
        QDagSpDagQp2Dag_C2D_sp_qp2.resize(Nsp_I, Nsp_I);
        QSpDagSpDagQ_C2D_sp_sp.resize(Nsp_I, Nsp_I);
        QSpDagSp_C2D_sp_sp.resize(Nsp_I, Nsp_I);
        QDagSpDagSp_C2D_sp_sp.resize(Nsp_I, Nsp_I);
        QSpDagSpDag_C2D_sp_sp.resize(Nsp_I, Nsp_I);

        // P contractions.
        Qp1SpDagP_C2D_qp1_sp.resize(Nsp_I, Nsp_I);
        Qp1SpP_C2D_qp1_sp.resize(Nsp_I, Nsp_I);
        PDagSpDagQp2Dag_C2D_sp_qp2.resize(Nsp_I, Nsp_I);
        PDagSpQp2Dag_C2D_sp_qp2.resize(Nsp_I, Nsp_I);
        PDagSpDagSpDag_C2D_sp_sp.resize(Nsp_I, Nsp_I);
        PDagSpDagSp_C2D_sp_sp.resize(Nsp_I, Nsp_I);
        PDagSpDagSpP_C2D_sp_sp.resize(Nsp_I, Nsp_I);
        SpDagSpP_C2D_sp_sp.resize(Nsp_I, Nsp_I);
        SpSpP_C2D_sp_sp.resize(Nsp_I, Nsp_I);

        // Shared kernel workspaces; XX_ii = 0.
        Qp1X_C2D_qp1_X.resize(Nsp_I, 4);
        XQp2Dag_C2D_X_qp2.resize(4, Nsp_I);
        XX_C2D_X_X.resize(4, 4);
        kernel_C2D_cfg1_cfg2.resize(Ncfg1_I, Ncfg2_I);
        XX_C2D_X_X.diagonal().setZero();

        x1_C1D_sp.resize(Nsp_I);
        x2_C1D_sp.resize(Nsp_I);
        x3_C1D_sp.resize(Nsp_I);
        x4_C1D_sp.resize(Nsp_I);
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
     * @brief  Integrate creation kernels using uniform gauge quadrature.
     * @math   result_abα = (1/Nφ) Σφ eⁱᴺφ ⟨Φ₁;a|cα† e⁻ⁱφᴺ̂ D̂|Φ₂;b⟩.
     * @output Const reference to shared result_C3D_cfg1_cfg2_sp.
     */
    const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& calc_creator_pnp();

    /**
     * @brief  Integrate annihilation kernels using uniform gauge quadrature.
     * @math   result_abα = (1/Nφ) Σφ eⁱᴺφ ⟨Φ₁;a|cα e⁻ⁱφᴺ̂ D̂|Φ₂;b⟩.
     * @output Const reference to shared result_C3D_cfg1_cfg2_sp.
     */
    const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& calc_annihilator_pnp();

    /**
     * @brief  Integrate one-body kernels using uniform gauge quadrature.
     * @math   result = (1/Nφ) Σφ eⁱᴺφ Σ_ij OneBody_ij OBTD(i,j;φ,D).
     * @output Const reference to result_C2D_cfg1_cfg2.
     */
    const Eigen::MatrixXcd& calc_one_body_pnp(const Eigen::MatrixXcd& OneBody_C2D_sp_sp);

    /**
     * @brief  Integrate two-body kernels using uniform gauge quadrature.
     * @math   result = (1/Nφ) Σφ eⁱᴺφ ½Σ_ijkl TwoBody_ijkl TBTD(i,j,k,l;φ,D).
     * @output Const reference to result_C2D_cfg1_cfg2.
     * @note   TwoBody contains unsymmetrized matrix elements ⟨ij|v|kl⟩.
     */
    const Eigen::MatrixXcd& calc_two_body_pnp(const Eigen::Tensor<double, 4, Eigen::ColMajor>& TwoBody_F4D_sp_sp_sp_sp);

    /**
     * @brief  Integrate separable kernels using Pfaffian multilinearity and quadrature.
     * @math   O = Σμ[Q̂λμ†Q̂λμ − Σαδ(Qλμ†Qλμ)αδ cα†cδ].
     * @math   Q ∈ ℝ^{Nsp × Nsp × (2λ+1)}; twoLambda_I = 2λ ≥ 0, λ ∈ ℤ; μ = −λ,…,λ.
     * @output Const reference to result_C2D_cfg1_cfg2, including vacuum overlap.
     * @note   Excludes coupling strength and extra ½.
     */
    const Eigen::MatrixXcd& calc_two_body_Q_pnp(int twoLambda_I, const Eigen::Tensor<double, 3, Eigen::ColMajor>& Q_F3D_sp_sp_mu);

    /**
     * @brief  Integrate pairing kernels using Pfaffian multilinearity and quadrature.
     * @math   P̂λμ† = ½Σab Pab c_a†c_b†; O = Σμ P̂λμ†P̂λμ.
     * @math   P ∈ ℝ^{Nsp × Nsp × (2λ+1)}; Pμᵀ = −Pμ; twoLambda_I = 2λ ≥ 0, λ ∈ ℤ; μ = −λ,…,λ.
     * @output Const reference to result_C2D_cfg1_cfg2, including vacuum overlap.
     * @note   Excludes coupling strength and attraction sign.
     */
    const Eigen::MatrixXcd& calc_two_body_P_pnp(int twoLambda_I, const Eigen::Tensor<double, 3, Eigen::ColMajor>& P_F3D_sp_sp_mu);

protected:
    /**
     * @brief  Integrate kernels using uniform gauge quadrature.
     * @math   result = (1/Nφ) Σφ exp(iNφ) K(φ,D).
     * @output Filled caller-provided contiguous result view.
     * @note   Callback accumulates weighted kernels using current contractions.
     * @note   Callback preserves D,U,V, targets, mesh, and result storage.
     * @note   Callback supplies vacuum overlap; normalization is applied here.
     */
    template <typename AccumulateFunc>
    void integrate_pnp(Eigen::Ref<Eigen::VectorXcd> result_C1D_element, const AccumulateFunc& accumulate_Func);

};

template <typename AccumulateFunc>
inline void HFBProjectionPNP::integrate_pnp(Eigen::Ref<Eigen::VectorXcd> result_C1D_element, const AccumulateFunc& accumulate_Func) {
    result_C1D_element.setZero();
    assert(D_C2D_sp_sp.rows() == Nsp_I && D_C2D_sp_sp.cols() == Nsp_I && D_C2D_sp_sp.allFinite());
    // DU₂ = D U₂; DV₂ = D* V₂; independent of φ.
    DU2_C2D_sp_qp2.noalias() = D_C2D_sp_sp * U2_C2D_sp_qp2;
    DV2_C2D_sp_qp2.noalias() = D_C2D_sp_sp.conjugate() * V2_C2D_sp_qp2;
    // Contractions are prepared once per gauge node.
    for (int phi_I = 0; phi_I < Nphi_I; ++phi_I) {
        const double phi_F = 2.0 * std::numbers::pi * phi_I / Nphi_I;
        const doubleC factorPhi_C = std::exp(doubleC(0.0, TargetN_I * phi_F));
        // U₂φ = exp(-iφ) DU₂; V₂φ = exp(iφ) DV₂.
        const doubleC gauge_C = std::exp(doubleC(0.0, -phi_F));
        gDU2_C2D_sp_qp2 = gauge_C * DU2_C2D_sp_qp2;
        gDV2_C2D_sp_qp2 = std::conj(gauge_C) * DV2_C2D_sp_qp2;
        hfb_pfaffian.update_contractions(U1_C2D_sp_qp1, V1_C2D_sp_qp1, gDU2_C2D_sp_qp2, gDV2_C2D_sp_qp2, doubleC(1.0, 0.0));
        accumulate_Func(factorPhi_C);
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
    integrate_pnp(result_C1D_element, [&](doubleC factorPhi_C) {
        result_C2D_cfg1_cfg2 += factorPhi_C * hfb_pfaffian.calc_overlap();
    });
    return result_C2D_cfg1_cfg2;
}

inline const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& HFBProjectionPNP::calc_creator_pnp() {
    Eigen::Map<Eigen::VectorXcd> result_C1D_element(result_C3D_cfg1_cfg2_sp.data(), result_C3D_cfg1_cfg2_sp.size());
    const Eigen::Index Nblock_I = static_cast<Eigen::Index>(Ncfg1_I) * Ncfg2_I;
    // ColMajor: (cfg1,cfg2) → rows; sp → columns.
    Eigen::Map<Eigen::MatrixXcd> result_C2D_cfg1cfg2_sp(result_C3D_cfg1_cfg2_sp.data(), Nblock_I, Nsp_I);

    // (1/Nφ) Σφ exp(iNφ) ⟨cα† exp(-iφN̂) D̂⟩.
    integrate_pnp(result_C1D_element, [&](doubleC factorPhi_C) {
        for (int sp_I = 0; sp_I < Nsp_I; ++sp_I) {
            // x = eα selects the original orbital.
            x1_C1D_sp.setZero();
            x1_C1D_sp(sp_I) = 1.0;
            const auto& kernel_C2D_cfg1_cfg2 = hfb_pfaffian.calc_creator(x1_C1D_sp);
            const Eigen::Map<const Eigen::VectorXcd> kernel_C1D_cfg1cfg2(kernel_C2D_cfg1_cfg2.data(), Nblock_I);
            result_C2D_cfg1cfg2_sp.col(sp_I) += factorPhi_C * kernel_C1D_cfg1cfg2;
        }
    });
    return result_C3D_cfg1_cfg2_sp;
}

inline const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& HFBProjectionPNP::calc_annihilator_pnp() {
    Eigen::Map<Eigen::VectorXcd> result_C1D_element(result_C3D_cfg1_cfg2_sp.data(), result_C3D_cfg1_cfg2_sp.size());
    const Eigen::Index Nblock_I = static_cast<Eigen::Index>(Ncfg1_I) * Ncfg2_I;
    // ColMajor: (cfg1,cfg2) → rows; sp → columns.
    Eigen::Map<Eigen::MatrixXcd> result_C2D_cfg1cfg2_sp(result_C3D_cfg1_cfg2_sp.data(), Nblock_I, Nsp_I);

    // (1/Nφ) Σφ exp(iNφ) ⟨cα exp(-iφN̂) D̂⟩.
    integrate_pnp(result_C1D_element, [&](doubleC factorPhi_C) {
        for (int sp_I = 0; sp_I < Nsp_I; ++sp_I) {
            // x = eα selects the original orbital.
            x1_C1D_sp.setZero();
            x1_C1D_sp(sp_I) = 1.0;
            const auto& kernel_C2D_cfg1_cfg2 = hfb_pfaffian.calc_annihilator(x1_C1D_sp);
            const Eigen::Map<const Eigen::VectorXcd> kernel_C1D_cfg1cfg2(kernel_C2D_cfg1_cfg2.data(), Nblock_I);
            result_C2D_cfg1cfg2_sp.col(sp_I) += factorPhi_C * kernel_C1D_cfg1cfg2;
        }
    });
    return result_C3D_cfg1_cfg2_sp;
}

inline const Eigen::MatrixXcd& HFBProjectionPNP::calc_one_body_pnp(const Eigen::MatrixXcd& OneBody_C2D_sp_sp) {
    assert(OneBody_C2D_sp_sp.rows() == Nsp_I && OneBody_C2D_sp_sp.cols() == Nsp_I && OneBody_C2D_sp_sp.allFinite());

    // H = Σβ (Σα Oαβ cα†) cβ; coefficients enter linearly.
    Eigen::Map<Eigen::VectorXcd> result_C1D_element(result_C2D_cfg1_cfg2.data(), result_C2D_cfg1_cfg2.size());
    integrate_pnp(result_C1D_element, [&](doubleC factorPhi_C) {
        for (int sp2_I = 0; sp2_I < Nsp_I; ++sp2_I) {
            // x₁ = O_:β; x₂ = eβ.
            x1_C1D_sp = OneBody_C2D_sp_sp.col(sp2_I);
            x2_C1D_sp.setZero();
            x2_C1D_sp(sp2_I) = 1.0;
            result_C2D_cfg1_cfg2 += factorPhi_C * hfb_pfaffian.calc_obtd(x1_C1D_sp, x2_C1D_sp);
        }
    });
    return result_C2D_cfg1_cfg2;
}

inline const Eigen::MatrixXcd& HFBProjectionPNP::calc_two_body_pnp(const Eigen::Tensor<double, 4, Eigen::ColMajor>& TwoBody_F4D_sp_sp_sp_sp) {
    assert(TwoBody_F4D_sp_sp_sp_sp.dimension(0) == Nsp_I && TwoBody_F4D_sp_sp_sp_sp.dimension(1) == Nsp_I && TwoBody_F4D_sp_sp_sp_sp.dimension(2) == Nsp_I && TwoBody_F4D_sp_sp_sp_sp.dimension(3) == Nsp_I);
    assert(Eigen::Map<const Eigen::VectorXd>(TwoBody_F4D_sp_sp_sp_sp.data(), TwoBody_F4D_sp_sp_sp_sp.size()).allFinite());

    // x₁,…,x₄ select the original single-particle orbitals.

    // Σφ exp(iNφ) O(φ,D).
    Eigen::Map<Eigen::VectorXcd> result_C1D_element(result_C2D_cfg1_cfg2.data(), result_C2D_cfg1_cfg2.size());
    integrate_pnp(result_C1D_element, [&](doubleC factorPhi_C) {
        // H²_ab(g) = ½Σ_ijkl v_ijkl TBTD_ab(i,j,k,l;g).
        for (int sp4_I = 0; sp4_I < Nsp_I; ++sp4_I) {
            x4_C1D_sp.setZero();
            x4_C1D_sp(sp4_I) = 1.0;
            for (int sp3_I = 0; sp3_I < Nsp_I; ++sp3_I) {
                x3_C1D_sp.setZero();
                x3_C1D_sp(sp3_I) = 1.0;
                for (int sp2_I = 0; sp2_I < Nsp_I; ++sp2_I) {
                    x2_C1D_sp.setZero();
                    x2_C1D_sp(sp2_I) = 1.0;
                    for (int sp1_I = 0; sp1_I < Nsp_I; ++sp1_I) {
                        x1_C1D_sp.setZero();
                        x1_C1D_sp(sp1_I) = 1.0;
                        result_C2D_cfg1_cfg2 += factorPhi_C * 0.5 * TwoBody_F4D_sp_sp_sp_sp(sp1_I, sp2_I, sp3_I, sp4_I) * hfb_pfaffian.calc_tbtd(x1_C1D_sp, x2_C1D_sp, x3_C1D_sp, x4_C1D_sp);
                    }
                }
            }
        }
    });
    return result_C2D_cfg1_cfg2;
}

inline const Eigen::MatrixXcd& HFBProjectionPNP::calc_two_body_Q_pnp(int twoLambda_I, const Eigen::Tensor<double, 3, Eigen::ColMajor>& Q_F3D_sp_sp_mu) {
    assert(twoLambda_I >= 0 && twoLambda_I % 2 == 0);
    assert(Q_F3D_sp_sp_mu.dimension(0) == Nsp_I && Q_F3D_sp_sp_mu.dimension(1) == Nsp_I && Q_F3D_sp_sp_mu.dimension(2) == twoLambda_I + 1);
    assert(Eigen::Map<const Eigen::VectorXd>(Q_F3D_sp_sp_mu.data(), Q_F3D_sp_sp_mu.size()).allFinite());

    // Oμ = Σbd (Σa Q_ba c†a)(Σc Q_cd c†c) cd cb.
    // X = (x₁†,x₂†,cd,cb).
    Eigen::Map<Eigen::VectorXcd> result_C1D_element(result_C2D_cfg1_cfg2.data(), result_C2D_cfg1_cfg2.size());
    integrate_pnp(result_C1D_element, [&](doubleC factorPhi_C) {
        for (int mu_I = 0; mu_I < twoLambda_I + 1; ++mu_I) {
            const Eigen::Map<const Eigen::MatrixXd> Q_F2D_sp_sp(Q_F3D_sp_sp_mu.data() + static_cast<Eigen::Index>(mu_I) * Nsp_I * Nsp_I, Nsp_I, Nsp_I);

            // ⟨β₁x₁†⟩,⟨β₁x₂†⟩; ⟨x₁†β₂†⟩,⟨x₂†β₂†⟩.
            Qp1SpDagQDag_C2D_qp1_sp.noalias() = hfb_pfaffian.Qp1SpDag_C2D_qp1_sp * Q_F2D_sp_sp.transpose().cast<doubleC>();
            Qp1SpDagQ_C2D_qp1_sp.noalias() = hfb_pfaffian.Qp1SpDag_C2D_qp1_sp * Q_F2D_sp_sp.cast<doubleC>();
            QSpDagQp2Dag_C2D_sp_qp2.noalias() = Q_F2D_sp_sp.cast<doubleC>() * hfb_pfaffian.SpDagQp2Dag_C2D_sp_qp2;
            QDagSpDagQp2Dag_C2D_sp_qp2.noalias() = Q_F2D_sp_sp.transpose().cast<doubleC>() * hfb_pfaffian.SpDagQp2Dag_C2D_sp_qp2;

            // ⟨x₁†x₂†⟩ = Qκ̄Q; ⟨x₁†c⟩ = Qρ; ⟨x₂†c⟩ = Qᵀρ.
            QSpDagSpDag_C2D_sp_sp.noalias() = Q_F2D_sp_sp.cast<doubleC>() * hfb_pfaffian.SpDagSpDag_C2D_sp_sp;
            QSpDagSpDagQ_C2D_sp_sp.noalias() = QSpDagSpDag_C2D_sp_sp * Q_F2D_sp_sp.cast<doubleC>();
            QSpDagSp_C2D_sp_sp.noalias() = Q_F2D_sp_sp.cast<doubleC>() * hfb_pfaffian.SpDagSp_C2D_sp_sp;
            QDagSpDagSp_C2D_sp_sp.noalias() = Q_F2D_sp_sp.transpose().cast<doubleC>() * hfb_pfaffian.SpDagSp_C2D_sp_sp;

            // b≠d; cb cb = 0.
            for (int spD_I = 0; spD_I < Nsp_I; ++spD_I) {
                for (int spB_I = 0; spB_I < Nsp_I; ++spB_I) {
                    if (spB_I == spD_I) {continue;}
                    // Select ⟨β₁X⟩ and ⟨Xβ₂†⟩ in physical order.
                    Qp1X_C2D_qp1_X.col(0) = Qp1SpDagQDag_C2D_qp1_sp.col(spB_I);
                    Qp1X_C2D_qp1_X.col(1) = Qp1SpDagQ_C2D_qp1_sp.col(spD_I);
                    Qp1X_C2D_qp1_X.col(2) = hfb_pfaffian.Qp1Sp_C2D_qp1_sp.col(spD_I);
                    Qp1X_C2D_qp1_X.col(3) = hfb_pfaffian.Qp1Sp_C2D_qp1_sp.col(spB_I);
                    XQp2Dag_C2D_X_qp2.row(0) = QSpDagQp2Dag_C2D_sp_qp2.row(spB_I);
                    XQp2Dag_C2D_X_qp2.row(1) = QDagSpDagQp2Dag_C2D_sp_qp2.row(spD_I);
                    XQp2Dag_C2D_X_qp2.row(2) = hfb_pfaffian.SpQp2Dag_C2D_sp_qp2.row(spD_I);
                    XQp2Dag_C2D_X_qp2.row(3) = hfb_pfaffian.SpQp2Dag_C2D_sp_qp2.row(spB_I);

                    // XX_ij = ⟨X_i X_j⟩ for i<j; XX_ji = −XX_ij.
                    XX_C2D_X_X(0, 1) = QSpDagSpDagQ_C2D_sp_sp(spB_I, spD_I);
                    XX_C2D_X_X(0, 2) = QSpDagSp_C2D_sp_sp(spB_I, spD_I);
                    XX_C2D_X_X(0, 3) = QSpDagSp_C2D_sp_sp(spB_I, spB_I);
                    XX_C2D_X_X(1, 2) = QDagSpDagSp_C2D_sp_sp(spD_I, spD_I);
                    XX_C2D_X_X(1, 3) = QDagSpDagSp_C2D_sp_sp(spD_I, spB_I);
                    XX_C2D_X_X(2, 3) = hfb_pfaffian.SpSp_C2D_sp_sp(spD_I, spB_I);
                    for (int X2_I = 1; X2_I < 4; ++X2_I) {
                        for (int X1_I = 0; X1_I < X2_I; ++X1_I) {XX_C2D_X_X(X2_I, X1_I) = -XX_C2D_X_X(X1_I, X2_I);}
                    }

                    // Kernel includes vacuum overlap; no additional ½.
                    hfb_pfaffian.calc_kernel(Qp1X_C2D_qp1_X, XQp2Dag_C2D_X_qp2, XX_C2D_X_X, kernel_C2D_cfg1_cfg2);
                    result_C2D_cfg1_cfg2 += factorPhi_C * kernel_C2D_cfg1_cfg2;
                }
            }
        }
    });
    return result_C2D_cfg1_cfg2;
}

inline const Eigen::MatrixXcd& HFBProjectionPNP::calc_two_body_P_pnp(int twoLambda_I, const Eigen::Tensor<double, 3, Eigen::ColMajor>& P_F3D_sp_sp_mu) {
    assert(twoLambda_I >= 0 && twoLambda_I % 2 == 0);
    assert(P_F3D_sp_sp_mu.dimension(0) == Nsp_I && P_F3D_sp_sp_mu.dimension(1) == Nsp_I && P_F3D_sp_sp_mu.dimension(2) == twoLambda_I + 1);
    assert(Eigen::Map<const Eigen::VectorXd>(P_F3D_sp_sp_mu.data(), P_F3D_sp_sp_mu.size()).allFinite());
    for (int mu_I = 0; mu_I <= twoLambda_I; ++mu_I) {
        const Eigen::Map<const Eigen::MatrixXd> P_F2D_sp_sp(P_F3D_sp_sp_mu.data() + static_cast<Eigen::Index>(mu_I) * Nsp_I * Nsp_I, Nsp_I, Nsp_I);
        assert(P_F2D_sp_sp.isApprox(-P_F2D_sp_sp.transpose(), 1.0e-12));
    }

    // x₁†(b) = Σa Pab c_a†; x₂(d) = Σc Pcd c_c.
    // X = (x₁†,c_b†,c_d,x₂).
    Eigen::Map<Eigen::VectorXcd> result_C1D_element(result_C2D_cfg1_cfg2.data(), result_C2D_cfg1_cfg2.size());
    integrate_pnp(result_C1D_element, [&](doubleC factorPhi_C) {
        for (int mu_I = 0; mu_I <= twoLambda_I; ++mu_I) {
            const Eigen::Map<const Eigen::MatrixXd> P_F2D_sp_sp(P_F3D_sp_sp_mu.data() + static_cast<Eigen::Index>(mu_I) * Nsp_I * Nsp_I, Nsp_I, Nsp_I);

            // ⟨β₁x₁†⟩,⟨β₁x₂⟩; ⟨x₁†β₂†⟩,⟨x₂β₂†⟩.
            Qp1SpDagP_C2D_qp1_sp.noalias() = hfb_pfaffian.Qp1SpDag_C2D_qp1_sp * P_F2D_sp_sp.cast<doubleC>();
            Qp1SpP_C2D_qp1_sp.noalias() = hfb_pfaffian.Qp1Sp_C2D_qp1_sp * P_F2D_sp_sp.cast<doubleC>();
            PDagSpDagQp2Dag_C2D_sp_qp2.noalias() = P_F2D_sp_sp.transpose().cast<doubleC>() * hfb_pfaffian.SpDagQp2Dag_C2D_sp_qp2;
            PDagSpQp2Dag_C2D_sp_qp2.noalias() = P_F2D_sp_sp.transpose().cast<doubleC>() * hfb_pfaffian.SpQp2Dag_C2D_sp_qp2;

            // ⟨x₁†c†⟩ = Pᵀκ̄; ⟨x₁†c⟩ = Pᵀρ; ⟨x₁†x₂⟩ = PᵀρP.
            PDagSpDagSpDag_C2D_sp_sp.noalias() = P_F2D_sp_sp.transpose().cast<doubleC>() * hfb_pfaffian.SpDagSpDag_C2D_sp_sp;
            PDagSpDagSp_C2D_sp_sp.noalias() = P_F2D_sp_sp.transpose().cast<doubleC>() * hfb_pfaffian.SpDagSp_C2D_sp_sp;
            PDagSpDagSpP_C2D_sp_sp.noalias() = PDagSpDagSp_C2D_sp_sp * P_F2D_sp_sp.cast<doubleC>();
            SpDagSpP_C2D_sp_sp.noalias() = hfb_pfaffian.SpDagSp_C2D_sp_sp * P_F2D_sp_sp.cast<doubleC>();
            SpSpP_C2D_sp_sp.noalias() = hfb_pfaffian.SpSp_C2D_sp_sp * P_F2D_sp_sp.cast<doubleC>();

            // P̂†P̂ = ¼Σbd x₁†(b)c_b†c_d x₂(d); b=d contributes.
            for (int spD_I = 0; spD_I < Nsp_I; ++spD_I) {
                for (int spB_I = 0; spB_I < Nsp_I; ++spB_I) {
                    // Select ⟨β₁X⟩ and ⟨Xβ₂†⟩ in physical order.
                    Qp1X_C2D_qp1_X.col(0) = Qp1SpDagP_C2D_qp1_sp.col(spB_I);
                    Qp1X_C2D_qp1_X.col(1) = hfb_pfaffian.Qp1SpDag_C2D_qp1_sp.col(spB_I);
                    Qp1X_C2D_qp1_X.col(2) = hfb_pfaffian.Qp1Sp_C2D_qp1_sp.col(spD_I);
                    Qp1X_C2D_qp1_X.col(3) = Qp1SpP_C2D_qp1_sp.col(spD_I);
                    XQp2Dag_C2D_X_qp2.row(0) = PDagSpDagQp2Dag_C2D_sp_qp2.row(spB_I);
                    XQp2Dag_C2D_X_qp2.row(1) = hfb_pfaffian.SpDagQp2Dag_C2D_sp_qp2.row(spB_I);
                    XQp2Dag_C2D_X_qp2.row(2) = hfb_pfaffian.SpQp2Dag_C2D_sp_qp2.row(spD_I);
                    XQp2Dag_C2D_X_qp2.row(3) = PDagSpQp2Dag_C2D_sp_qp2.row(spD_I);

                    // XX_ij = ⟨X_i X_j⟩ for i<j; XX_ji = −XX_ij.
                    XX_C2D_X_X(0, 1) = PDagSpDagSpDag_C2D_sp_sp(spB_I, spB_I);
                    XX_C2D_X_X(0, 2) = PDagSpDagSp_C2D_sp_sp(spB_I, spD_I);
                    XX_C2D_X_X(0, 3) = PDagSpDagSpP_C2D_sp_sp(spB_I, spD_I);
                    XX_C2D_X_X(1, 2) = hfb_pfaffian.SpDagSp_C2D_sp_sp(spB_I, spD_I);
                    XX_C2D_X_X(1, 3) = SpDagSpP_C2D_sp_sp(spB_I, spD_I);
                    XX_C2D_X_X(2, 3) = SpSpP_C2D_sp_sp(spD_I, spD_I);
                    for (int X2_I = 1; X2_I < 4; ++X2_I) {
                        for (int X1_I = 0; X1_I < X2_I; ++X1_I) {XX_C2D_X_X(X2_I, X1_I) = -XX_C2D_X_X(X1_I, X2_I);}
                    }

                    // Kernel includes overlap; ¼ follows the pair definitions.
                    hfb_pfaffian.calc_kernel(Qp1X_C2D_qp1_X, XQp2Dag_C2D_X_qp2, XX_C2D_X_X, kernel_C2D_cfg1_cfg2);
                    result_C2D_cfg1_cfg2 += (0.25 * factorPhi_C) * kernel_C2D_cfg1_cfg2;
                }
            }
        }
    });
    return result_C2D_cfg1_cfg2;
}
