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

private:
    Eigen::MatrixXcd result_C2D_cfg1_cfg2{};
    Eigen::Tensor<doubleC, 3, Eigen::ColMajor> result_C3D_cfg1_cfg2_sp{};
    Eigen::Tensor<doubleC, 4, Eigen::ColMajor> result_C4D_cfg1_cfg2_sp1_sp2{};
    Eigen::MatrixXcd DU2_C2D_sp_qp2{};
    Eigen::MatrixXcd DV2_C2D_sp_qp2{};
    Eigen::MatrixXcd gDU2_C2D_sp_qp2{};
    Eigen::MatrixXcd gDV2_C2D_sp_qp2{};

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
        result_C4D_cfg1_cfg2_sp1_sp2.resize(Ncfg1_I, Ncfg2_I, Nsp_I, Nsp_I);
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
     * @brief  Prepare transformed Bogoliubov matrices using matrix products.
     * @math   DU₂ = D U₂; DV₂ = D* V₂.
     * @output Updated DU₂,DV₂ workspaces.
     */
    void prepare_DUV();

    /**
     * @brief  Update gauge-transformed contractions using Pfaffians and LU.
     * @math   U₂φ = e⁻ⁱφ D U₂; V₂φ = eⁱφ D* V₂.
     * @output Updated U₂φ,V₂φ and Pfaffian contractions.
     * @note   Call prepare_DUV after changing D,U₂,V₂.
     */
    void prepare_contractions(int phi_I);

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
     * @brief  Integrate one-body densities using uniform gauge quadrature.
     * @math   result_abij = (1/Nφ) Σφ eⁱᴺφ ⟨Φ₁;a|ci† cj e⁻ⁱφᴺ̂ D̂|Φ₂;b⟩.
     * @output Const reference to result_C4D_cfg1_cfg2_sp1_sp2.
     */
    const Eigen::Tensor<doubleC, 4, Eigen::ColMajor>& calc_obtd_pnp();

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
};

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

inline void HFBProjectionPNP::prepare_DUV() {
    assert(D_C2D_sp_sp.rows() == Nsp_I && D_C2D_sp_sp.cols() == Nsp_I && D_C2D_sp_sp.allFinite());
    // DU₂ = D U₂; DV₂ = D* V₂.
    DU2_C2D_sp_qp2.noalias() = D_C2D_sp_sp * U2_C2D_sp_qp2;
    DV2_C2D_sp_qp2.noalias() = D_C2D_sp_sp.conjugate() * V2_C2D_sp_qp2;
}

inline void HFBProjectionPNP::prepare_contractions(int phi_I) {
    assert(phi_I >= 0 && phi_I < Nphi_I);

    // U₂φ = exp(-iφ) D U₂; V₂φ = exp(iφ) D* V₂.
    const double phi_F = 2.0 * std::numbers::pi * phi_I / Nphi_I;
    const doubleC gauge_C = std::exp(doubleC(0.0, -phi_F));
    gDU2_C2D_sp_qp2 = gauge_C * DU2_C2D_sp_qp2;
    gDV2_C2D_sp_qp2 = std::conj(gauge_C) * DV2_C2D_sp_qp2;
    hfb_pfaffian.update_contractions(U1_C2D_sp_qp1, V1_C2D_sp_qp1, gDU2_C2D_sp_qp2, gDV2_C2D_sp_qp2, doubleC(1.0, 0.0));
}

inline const Eigen::MatrixXcd& HFBProjectionPNP::calc_overlap_pnp() {
    // Σφ exp(iNφ) O(φ,D).
    result_C2D_cfg1_cfg2.setZero();
    prepare_DUV();
    for (int phi_I = 0; phi_I < Nphi_I; ++phi_I) {
        // U₂φ = exp(-iφ) D U₂; V₂φ = exp(iφ) D* V₂.
        const double phi_F = 2.0 * std::numbers::pi * phi_I / Nphi_I;
        prepare_contractions(phi_I);
        const doubleC factorPhi_C = std::exp(doubleC(0.0, TargetN_I * phi_F));
        result_C2D_cfg1_cfg2 += factorPhi_C * hfb_pfaffian.calc_overlap();
    }
    result_C2D_cfg1_cfg2 /= Nphi_I;
    return result_C2D_cfg1_cfg2;
}

inline const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& HFBProjectionPNP::calc_creator_pnp() {
    result_C3D_cfg1_cfg2_sp.setZero();
    const Eigen::Index Nblock_I = static_cast<Eigen::Index>(Ncfg1_I) * Ncfg2_I;
    // ColMajor: (cfg1,cfg2) → rows; sp → columns.
    Eigen::Map<Eigen::MatrixXcd> result_C2D_cfg1cfg2_sp(result_C3D_cfg1_cfg2_sp.data(), Nblock_I, Nsp_I);

    // (1/Nφ) Σφ exp(iNφ) ⟨cα† exp(-iφN̂) D̂⟩.
    prepare_DUV();
    for (int phi_I = 0; phi_I < Nphi_I; ++phi_I) {
        const double phi_F = 2.0 * std::numbers::pi * phi_I / Nphi_I;
        const doubleC factorPhi_C = std::exp(doubleC(0.0, TargetN_I * phi_F));
        prepare_contractions(phi_I);
        for (int sp_I = 0; sp_I < Nsp_I; ++sp_I) {
            const auto& kernel_C2D_cfg1_cfg2 = hfb_pfaffian.calc_creator(sp_I);
            const Eigen::Map<const Eigen::VectorXcd> kernel_C1D_cfg1cfg2(kernel_C2D_cfg1_cfg2.data(), Nblock_I);
            result_C2D_cfg1cfg2_sp.col(sp_I) += factorPhi_C * kernel_C1D_cfg1cfg2;
        }
    }
    result_C2D_cfg1cfg2_sp /= Nphi_I;
    return result_C3D_cfg1_cfg2_sp;
}

inline const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& HFBProjectionPNP::calc_annihilator_pnp() {
    result_C3D_cfg1_cfg2_sp.setZero();
    const Eigen::Index Nblock_I = static_cast<Eigen::Index>(Ncfg1_I) * Ncfg2_I;
    // ColMajor: (cfg1,cfg2) → rows; sp → columns.
    Eigen::Map<Eigen::MatrixXcd> result_C2D_cfg1cfg2_sp(result_C3D_cfg1_cfg2_sp.data(), Nblock_I, Nsp_I);

    // (1/Nφ) Σφ exp(iNφ) ⟨cα exp(-iφN̂) D̂⟩.
    prepare_DUV();
    for (int phi_I = 0; phi_I < Nphi_I; ++phi_I) {
        const double phi_F = 2.0 * std::numbers::pi * phi_I / Nphi_I;
        const doubleC factorPhi_C = std::exp(doubleC(0.0, TargetN_I * phi_F));
        prepare_contractions(phi_I);
        for (int sp_I = 0; sp_I < Nsp_I; ++sp_I) {
            const auto& kernel_C2D_cfg1_cfg2 = hfb_pfaffian.calc_annihilator(sp_I);
            const Eigen::Map<const Eigen::VectorXcd> kernel_C1D_cfg1cfg2(kernel_C2D_cfg1_cfg2.data(), Nblock_I);
            result_C2D_cfg1cfg2_sp.col(sp_I) += factorPhi_C * kernel_C1D_cfg1cfg2;
        }
    }
    result_C2D_cfg1cfg2_sp /= Nphi_I;
    return result_C3D_cfg1_cfg2_sp;
}

inline const Eigen::Tensor<doubleC, 4, Eigen::ColMajor>& HFBProjectionPNP::calc_obtd_pnp() {
    result_C4D_cfg1_cfg2_sp1_sp2.setZero();
    const Eigen::Index Nblock_I = static_cast<Eigen::Index>(Ncfg1_I) * Ncfg2_I;
    const Eigen::Index Npair_I = static_cast<Eigen::Index>(Nsp_I) * Nsp_I;
    // ColMajor: (cfg1,cfg2) → rows; (sp1,sp2) → columns.
    Eigen::Map<Eigen::MatrixXcd> result_C2D_cfg1cfg2_sp1sp2(result_C4D_cfg1_cfg2_sp1_sp2.data(), Nblock_I, Npair_I);

    // (1/Nφ) Σφ exp(iNφ) ⟨ci† cj exp(-iφN̂) D̂⟩.
    prepare_DUV();
    for (int phi_I = 0; phi_I < Nphi_I; ++phi_I) {
        const double phi_F = 2.0 * std::numbers::pi * phi_I / Nphi_I;
        const doubleC factorPhi_C = std::exp(doubleC(0.0, TargetN_I * phi_F));
        prepare_contractions(phi_I);
        for (int sp2_I = 0; sp2_I < Nsp_I; ++sp2_I) {
            for (int sp1_I = 0; sp1_I < Nsp_I; ++sp1_I) {
                const Eigen::Index pair_I = sp1_I + static_cast<Eigen::Index>(Nsp_I) * sp2_I;
                const auto& kernel_C2D_cfg1_cfg2 = hfb_pfaffian.calc_obtd(sp1_I, sp2_I);
                const Eigen::Map<const Eigen::VectorXcd> kernel_C1D_cfg1cfg2(kernel_C2D_cfg1_cfg2.data(), Nblock_I);
                result_C2D_cfg1cfg2_sp1sp2.col(pair_I) += factorPhi_C * kernel_C1D_cfg1cfg2;
            }
        }
    }
    result_C2D_cfg1cfg2_sp1sp2 /= Nphi_I;
    return result_C4D_cfg1_cfg2_sp1_sp2;
}

inline const Eigen::MatrixXcd& HFBProjectionPNP::calc_one_body_pnp(const Eigen::MatrixXcd& OneBody_C2D_sp_sp) {
    assert(OneBody_C2D_sp_sp.rows() == Nsp_I && OneBody_C2D_sp_sp.cols() == Nsp_I && OneBody_C2D_sp_sp.allFinite());

    const auto& OBTD_C4D_cfg1_cfg2_sp1_sp2 = calc_obtd_pnp();
    const Eigen::Index Nblock_I = static_cast<Eigen::Index>(Ncfg1_I) * Ncfg2_I;
    const Eigen::Index Npair_I = static_cast<Eigen::Index>(Nsp_I) * Nsp_I;

    // ColMajor: (cfg1,cfg2) → rows; (sp1,sp2) → columns.
    const Eigen::Map<const Eigen::MatrixXcd> OBTD_C2D_cfg1cfg2_sp1sp2(OBTD_C4D_cfg1_cfg2_sp1_sp2.data(), Nblock_I, Npair_I);
    const Eigen::Map<const Eigen::VectorXcd> OneBody_C1D_sp1sp2(OneBody_C2D_sp_sp.data(), Npair_I);
    Eigen::Map<Eigen::VectorXcd> result_C1D_cfg1cfg2(result_C2D_cfg1_cfg2.data(), Nblock_I);

    // H_ab = Σ_ij O_ij OBTD_abij; no complex conjugation.
    result_C1D_cfg1cfg2.noalias() = OBTD_C2D_cfg1cfg2_sp1sp2 * OneBody_C1D_sp1sp2;
    return result_C2D_cfg1_cfg2;
}

inline const Eigen::MatrixXcd& HFBProjectionPNP::calc_two_body_pnp(const Eigen::Tensor<double, 4, Eigen::ColMajor>& TwoBody_F4D_sp_sp_sp_sp) {
    assert(TwoBody_F4D_sp_sp_sp_sp.dimension(0) == Nsp_I && TwoBody_F4D_sp_sp_sp_sp.dimension(1) == Nsp_I && TwoBody_F4D_sp_sp_sp_sp.dimension(2) == Nsp_I && TwoBody_F4D_sp_sp_sp_sp.dimension(3) == Nsp_I);
    assert(Eigen::Map<const Eigen::VectorXd>(TwoBody_F4D_sp_sp_sp_sp.data(), TwoBody_F4D_sp_sp_sp_sp.size()).allFinite());

    // Σφ exp(iNφ) O(φ,D).
    result_C2D_cfg1_cfg2.setZero();
    prepare_DUV();
    for (int phi_I = 0; phi_I < Nphi_I; ++phi_I) {
        // U₂φ = exp(-iφ) D U₂; V₂φ = exp(iφ) D* V₂.
        const double phi_F = 2.0 * std::numbers::pi * phi_I / Nphi_I;
        prepare_contractions(phi_I);
        const doubleC factorPhi_C = std::exp(doubleC(0.0, TargetN_I * phi_F));
        // H²_ab(g) = ½Σ_ijkl v_ijkl TBTD_ab(i,j,k,l;g).
        for (int sp4_I = 0; sp4_I < Nsp_I; ++sp4_I) {
            for (int sp3_I = 0; sp3_I < Nsp_I; ++sp3_I) {
                for (int sp2_I = 0; sp2_I < Nsp_I; ++sp2_I) {
                    for (int sp1_I = 0; sp1_I < Nsp_I; ++sp1_I) {
                        result_C2D_cfg1_cfg2 += factorPhi_C * 0.5 * TwoBody_F4D_sp_sp_sp_sp(sp1_I, sp2_I, sp3_I, sp4_I) * hfb_pfaffian.calc_tbtd(sp1_I, sp2_I, sp3_I, sp4_I);
                    }
                }
            }
        }
    }
    result_C2D_cfg1_cfg2 /= Nphi_I;
    return result_C2D_cfg1_cfg2;
}
