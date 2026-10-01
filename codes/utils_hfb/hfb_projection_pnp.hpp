/**
 * @file    hfb_projection_pnp.hpp
 * @author  cailea
 * @date    2026-09-18
 * @brief   Particle-number projection of HFB configurations.
 */

#pragma once

#include <algorithm>
#include <array>
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

    std::vector<std::vector<doubleC>> F0_C2D_cfg1_cfg2{};
    std::vector<std::vector<Eigen::MatrixXcd>> F2_C4D_cfg1_cfg2_chain_chain{};
    std::vector<std::vector<Eigen::Tensor<doubleC, 4, Eigen::ColMajor>>> F4_C6D_cfg1_cfg2_chain_chain_chain_chain{};

public:
    /**
     * @brief  Construct empty particle-number projection workspaces.
     * @math   Nsp = Ncfg1 = Ncfg2 = Nφ = 0.
     * @output Empty workspaces.
     * @note   Assign a configured object before calculations.
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
        result_C3D_cfg1_cfg2_sp.resize(Ncfg1_I, Ncfg2_I, Nsp_I);
        DU2_C2D_sp_qp2.resize(Nsp_I, Nsp_I);
        DV2_C2D_sp_qp2.resize(Nsp_I, Nsp_I);
        gDU2_C2D_sp_qp2.resize(Nsp_I, Nsp_I);
        gDV2_C2D_sp_qp2.resize(Nsp_I, Nsp_I);

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
     * @brief  Integrate Q/P kernels using Pfaffians and gauge quadrature.
     * @math   K = Σi gQi(Qi†Qi − one-body contraction) + Σi gPi Pi†Pi.
     * @math   P† = ½Σab Pab c†a c†b; Pᵀ = −P.
     * @output Weighted kernel in result_C2D_cfg1_cfg2.
     * @note   gQ,gP include signs; pairing ¼ is internal.
     * @note   Q/P lists may differ in length; zero coefficients skip terms.
     */
    const Eigen::MatrixXcd& calc_two_body_pnp(const Eigen::VectorXd& gQ_F1D_i, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& Q_C3D_sp_sp_i, const Eigen::VectorXd& gP_F1D_i, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& P_C3D_sp_sp_i);

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

    const auto build_F0F2F4 = [&]() {
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

        F0_C2D_cfg1_cfg2.resize(Ncfg1_I);
        F2_C4D_cfg1_cfg2_chain_chain.resize(Ncfg1_I);
        F4_C6D_cfg1_cfg2_chain_chain_chain_chain.resize(Ncfg1_I);
        for (int cfg1_I = 0; cfg1_I < Ncfg1_I; ++cfg1_I) {
            F0_C2D_cfg1_cfg2[cfg1_I].assign(Ncfg2_I, 0.0);
            F2_C4D_cfg1_cfg2_chain_chain[cfg1_I].resize(Ncfg2_I);
            F4_C6D_cfg1_cfg2_chain_chain_chain_chain[cfg1_I].resize(Ncfg2_I);
        }

        for (int cfg2_I = 0; cfg2_I < Ncfg2_I; ++cfg2_I) {
            for (int cfg1_I = 0; cfg1_I < Ncfg1_I; ++cfg1_I) {
                // z = (reversed left β, forward right β†).
                const auto& config1_I1D_cqp1 = hfb_pfaffian.config1_I2D_cfg1_cqp1[cfg1_I];
                const auto& config2_I1D_cqp2 = hfb_pfaffian.config2_I2D_cfg2_cqp2[cfg2_I];

                std::vector<int> qp_I1D_chain{};
                for (int cqp1_I = static_cast<int>(config1_I1D_cqp1.size()) - 1; cqp1_I >= 0; --cqp1_I) { qp_I1D_chain.push_back(config1_I1D_cqp1[cqp1_I]); }
                for (int cqp2_I = 0; cqp2_I < static_cast<int>(config2_I1D_cqp2.size()); ++cqp2_I) { qp_I1D_chain.push_back(Nsp_I + config2_I1D_cqp2[cqp2_I]); }
                const int Nchain_I = static_cast<int>(qp_I1D_chain.size());

                auto& F2_C2D_chain_chain = F2_C4D_cfg1_cfg2_chain_chain[cfg1_I][cfg2_I];
                auto& F4_C4D_chain_chain_chain_chain = F4_C6D_cfg1_cfg2_chain_chain_chain_chain[cfg1_I][cfg2_I];

                F2_C2D_chain_chain.resize(Nchain_I, Nchain_I);
                F2_C2D_chain_chain.setZero();
                F4_C4D_chain_chain_chain_chain.resize(Nchain_I, Nchain_I, Nchain_I, Nchain_I);
                F4_C4D_chain_chain_chain_chain.setZero();
                if (Nchain_I % 2 != 0) { continue; }

                // Sᵀ = −S.
                Eigen::MatrixXcd S_C2D_chain_chain = Eigen::MatrixXcd::Zero(Nchain_I, Nchain_I);
                for (int chain2_I = 1; chain2_I < Nchain_I; ++chain2_I) {
                    for (int chain1_I = 0; chain1_I < chain2_I; ++chain1_I) {
                        const int qp1_I = qp_I1D_chain[chain1_I];
                        const int qp2_I = qp_I1D_chain[chain2_I];
                        if (qp1_I < Nsp_I && qp2_I < Nsp_I) { S_C2D_chain_chain(chain1_I, chain2_I) = hfb_pfaffian.Qp1Qp1_C2D_qp1_qp1(qp1_I, qp2_I); }
                        if (qp1_I < Nsp_I && qp2_I >= Nsp_I) { S_C2D_chain_chain(chain1_I, chain2_I) = hfb_pfaffian.Qp1Qp2Dag_C2D_qp1_qp2(qp1_I, qp2_I - Nsp_I); }
                        if (qp1_I >= Nsp_I && qp2_I >= Nsp_I) { S_C2D_chain_chain(chain1_I, chain2_I) = hfb_pfaffian.Qp2DagQp2Dag_C2D_qp2_qp2(qp1_I - Nsp_I, qp2_I - Nsp_I); }
                        S_C2D_chain_chain(chain2_I, chain1_I) = -S_C2D_chain_chain(chain1_I, chain2_I);
                    }
                }

                Eigen::VectorXcd workspace_C1D_element(Nchain_I * Nchain_I);
                std::vector<int> remaining_chain{};
                remaining_chain.reserve(Nchain_I);

                // Selected z move forward; four insertions cross evenly.
                const auto calc_signed_pfaffian_minor = [&](int chain1_I = -1, int chain2_I = -1, int chain3_I = -1, int chain4_I = -1) {
                    const int Nremove_I = (chain1_I >= 0) + (chain2_I >= 0) + (chain3_I >= 0) + (chain4_I >= 0);
                    int phase_I = Nremove_I * (Nremove_I - 1) / 2;
                    phase_I += std::max(chain1_I, 0) + std::max(chain2_I, 0) + std::max(chain3_I, 0) + std::max(chain4_I, 0);
                    double sign_F = 1.0 - 2.0 * (phase_I % 2);

                    remaining_chain.clear();
                    for (int chain_I = 0; chain_I < Nchain_I; ++chain_I) {
                        if (chain_I != chain1_I && chain_I != chain2_I && chain_I != chain3_I && chain_I != chain4_I) { remaining_chain.push_back(chain_I); }
                    }

                    const int Nremaining_I = Nchain_I - Nremove_I;
                    Eigen::Map<Eigen::MatrixXcd> minor_C2D_remainchain_remainchain(workspace_C1D_element.data(), Nremaining_I, Nremaining_I);
                    minor_C2D_remainchain_remainchain = S_C2D_chain_chain(remaining_chain, remaining_chain);
                    return sign_F * HFBPfaffian::calc_pfaffian(minor_C2D_remainchain_remainchain);
                };

                // F⁰ = pf(S).
                F0_C2D_cfg1_cfg2[cfg1_I][cfg2_I] = calc_signed_pfaffian_minor();

                // F²ij = −F²ji.
                for (int i_I = 0; i_I < Nchain_I; ++i_I) {
                    for (int j_I = i_I + 1; j_I < Nchain_I; ++j_I) {
                        F2_C2D_chain_chain(i_I, j_I) = calc_signed_pfaffian_minor(i_I, j_I);
                        F2_C2D_chain_chain(j_I, i_I) = -F2_C2D_chain_chain(i_I, j_I);
                    }
                }

                // F⁴: i < j < k < l → all permutations.
                for (int i_I = 0; i_I < Nchain_I; ++i_I) {
                    for (int j_I = i_I + 1; j_I < Nchain_I; ++j_I) {
                        for (int k_I = j_I + 1; k_I < Nchain_I; ++k_I) {
                            for (int l_I = k_I + 1; l_I < Nchain_I; ++l_I) {
                                const doubleC F4_C = calc_signed_pfaffian_minor(i_I, j_I, k_I, l_I);
                                const std::array<int, 4> chain_I1D_cqp{i_I, j_I, k_I, l_I};
                                // F⁴[p(i,j,k,l)] = sign(p) F⁴[i,j,k,l].
                                for (const auto& [a_I, b_I, c_I, d_I, sign_I] : permutationSign_I2D_permutation_entry) {
                                    F4_C4D_chain_chain_chain_chain(chain_I1D_cqp[a_I], chain_I1D_cqp[b_I], chain_I1D_cqp[c_I], chain_I1D_cqp[d_I]) = static_cast<double>(sign_I) * F4_C;
                                }
                            }
                        }
                    }
                }
            }
        }
    };

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

    /**
     * @brief Assemble kernels using Pfaffian minors.
     * @math K = g⟨Φ₁|Φ₂⟩(F⁰T⁰ + ½F²T² + ¼F⁴Π₁Π₂).
     * @output Accumulated result_C2D_cfg1_cfg2.
     */
    const auto calc_kernel024 = [&](double g_F, const Eigen::MatrixXcd& PI1_C2D_qp_qp, const Eigen::MatrixXcd& PI2_C2D_qp_qp, doubleC factorPhi_C) {
        for (int cfg2_I = 0; cfg2_I < Ncfg2_I; ++cfg2_I) {
            for (int cfg1_I = 0; cfg1_I < Ncfg1_I; ++cfg1_I) {
                const auto& config1_I1D_cqp = hfb_pfaffian.config1_I2D_cfg1_cqp1[cfg1_I];
                const auto& config2_I1D_cqp = hfb_pfaffian.config2_I2D_cfg2_cqp2[cfg2_I];
                const int Ncqp1_I = static_cast<int>(config1_I1D_cqp.size());
                const int Nchain_I = Ncqp1_I + static_cast<int>(config2_I1D_cqp.size());
                if (Nchain_I % 2 != 0) { continue; }
                // Chain position → left β or right β† index.
                const auto qp_index = [&](int chain_I) {
                    return chain_I < Ncqp1_I ? config1_I1D_cqp[Ncqp1_I - 1 - chain_I] : Nsp_I + config2_I1D_cqp[chain_I - Ncqp1_I];
                };
                const auto& F2_C2D_chain_chain = F2_C4D_cfg1_cfg2_chain_chain[cfg1_I][cfg2_I];
                const auto& F4_C4D_chain_chain_chain_chain = F4_C6D_cfg1_cfg2_chain_chain_chain_chain[cfg1_I][cfg2_I];
                doubleC kernel_C = F0_C2D_cfg1_cfg2[cfg1_I][cfg2_I] * T0_C;
                // K = F⁰T⁰ + ½Σab F²ab T²ab + ¼Σabcd F⁴abcd Πab Πcd.
                for (int i_I = 0; i_I < Nchain_I; ++i_I) {
                    const int qpi_I = qp_index(i_I);
                    for (int j_I = 0; j_I < Nchain_I; ++j_I) {
                        const int qpj_I = qp_index(j_I);
                        kernel_C += 0.5 * F2_C2D_chain_chain(i_I, j_I) * T2_C2D_qp_qp(qpi_I, qpj_I);
                        for (int k_I = 0; k_I < Nchain_I; ++k_I) {
                            const int qpk_I = qp_index(k_I);
                            for (int l_I = 0; l_I < Nchain_I; ++l_I) {
                                const int qpl_I = qp_index(l_I);
                                kernel_C += 0.25 * F4_C4D_chain_chain_chain_chain(i_I, j_I, k_I, l_I) * PI1_C2D_qp_qp(qpi_I, qpj_I) * PI2_C2D_qp_qp(qpk_I, qpl_I);
                            }
                        }
                    }
                }
                result_C2D_cfg1_cfg2(cfg1_I, cfg2_I) += factorPhi_C * hfb_pfaffian.overlap_C * g_F * kernel_C;
            }
        }
    };

    Eigen::Map<Eigen::VectorXcd> result_C1D_element(result_C2D_cfg1_cfg2.data(), result_C2D_cfg1_cfg2.size());
    integrate_pnp(result_C1D_element, [&](doubleC factorPhi_C) {
        build_F0F2F4();

        // Q: gQ (F⁰T⁰ + ½F²T² + ¼F⁴Π¹⁴Π²³).
        for (Eigen::Index iQ_I = 0; iQ_I < gQ_F1D_i.size(); ++iQ_I) {
            if (gQ_F1D_i(iQ_I) == 0.0) { continue; }
            const Eigen::Index offset_I = iQ_I * Nsp_I * Nsp_I;
            const Eigen::Map<const Eigen::MatrixXcd> Q_C2D_sp_sp(Q_C3D_sp_sp_i.data() + offset_I, Nsp_I, Nsp_I);
            calc_constract_Q(Q_C2D_sp_sp);
            calc_kernel024(gQ_F1D_i(iQ_I), PI14_C2D_qp_qp, PI23_C2D_qp_qp, factorPhi_C);
        }

        // P: ¼gP (F⁰T⁰ + ½F²T² + ¼F⁴Π¹²Π³⁴).
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
