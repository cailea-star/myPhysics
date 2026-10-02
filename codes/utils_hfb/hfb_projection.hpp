/**
 * @file    hfb_projection.hpp
 * @author  cailea
 * @date    2026-09-11
 * @brief   Angular-momentum projection of HFB configurations.
 */

#pragma once

#include <cassert>
#include <cmath>
#include <functional>
#include <utility>
#include <vector>

#include <Eigen/Core>
#include <unsupported/Eigen/CXX11/Tensor>

#include "group_so3_su2.hpp"
#include "hfb_projection_pnp.hpp"

/**
 * @brief  Project HFB configurations onto particle number and angular momentum.
 * @math   N = TargetN; 2I = TargetTwoI; K = -I+k, k = 0,…,2I.
 * @note   Use a common ordered basis and nonzero vacuum overlaps at quadrature nodes.
 */
class HFBProjection : public HFBProjectionPNP {
public:
    int Nalpha_I = 0;
    int Nbeta_I = 0;
    int Ngamma_I = 0;

    Eigen::Tensor<doubleC, 3, Eigen::ColMajor> Rz_C3D_sp_sp_alpha{};
    Eigen::Tensor<doubleC, 3, Eigen::ColMajor> Ry_C3D_sp_sp_beta{};
    Eigen::Tensor<doubleC, 3, Eigen::ColMajor> Rz_C3D_sp_sp_gamma{};

    Eigen::VectorXd alpha_F1D_alpha{};
    Eigen::VectorXd beta_F1D_beta{};
    Eigen::VectorXd gamma_F1D_gamma{};

    Eigen::VectorXd weight_F1D_alpha{};
    Eigen::VectorXd weight_F1D_beta{};
    Eigen::VectorXd weight_F1D_gamma{};


public:
    /**
     * @brief  Construct empty angular-momentum projection workspaces.
     * @math   Nsp = Nα = Nβ = Nγ = 0.
     * @output Empty workspaces.
     * @note   Assign a configured object before calculations.
     */
    HFBProjection() = default;

    /**
     * @brief  Allocate projection and Pfaffian workspaces.
     * @math   φ_k = 2πk/Nphi; Δφ/(2π) = 1/Nphi.
     * @output Stored dimensions, target quantum numbers, and allocated workspaces.
     * @note   Initialize projection data; left/right configuration lists may differ.
     */
    HFBProjection(int Nsp_I_, const std::vector<std::vector<int>>& config1_I2D_cfg1_cqp1_, const std::vector<std::vector<int>>& config2_I2D_cfg2_cqp2_, int TargetN_I_, int Nphi_I_, int Nalpha_I_, int Nbeta_I_, int Ngamma_I_)
    : HFBProjectionPNP(Nsp_I_, config1_I2D_cfg1_cqp1_, config2_I2D_cfg2_cqp2_, TargetN_I_, Nphi_I_) {
        assert(Nalpha_I_ > 0 && Nbeta_I_ > 0 && Ngamma_I_ > 0);

        // Nα,Nβ,Nγ → Euler mesh sizes.
        Nalpha_I = Nalpha_I_;
        Nbeta_I = Nbeta_I_;
        Ngamma_I = Ngamma_I_;

        // Rz(α),Ry(β),Rz(γ) ∈ ℂ^{Nsp×Nsp}.
        Rz_C3D_sp_sp_alpha.resize(Nsp_I, Nsp_I, Nalpha_I);
        Ry_C3D_sp_sp_beta.resize(Nsp_I, Nsp_I, Nbeta_I);
        Rz_C3D_sp_sp_gamma.resize(Nsp_I, Nsp_I, Ngamma_I);

        // {α,wα}, {β,wβ}, {γ,wγ}.
        alpha_F1D_alpha.resize(Nalpha_I);
        beta_F1D_beta.resize(Nbeta_I);
        gamma_F1D_gamma.resize(Ngamma_I);
        weight_F1D_alpha.resize(Nalpha_I);
        weight_F1D_beta.resize(Nbeta_I);
        weight_F1D_gamma.resize(Ngamma_I);

    }

    /**
     * @brief  Store alpha rotations, nodes, and quadrature weights.
     * @math   Rz ∈ ℂ^{Nsp×Nsp×Nalpha}; weights represent dα.
     * @output Updated Rz(α), alpha nodes, and weights.
     */
    void update_alpha(Eigen::Tensor<doubleC, 3, Eigen::ColMajor> Rz_C3D_sp_sp_alpha_, Eigen::VectorXd alpha_F1D_alpha_, Eigen::VectorXd weight_F1D_alpha_);

    /**
     * @brief  Store beta rotations, nodes, and quadrature weights.
     * @math   Ry ∈ ℂ^{Nsp×Nsp×Nbeta}; weights represent sinβ dβ.
     * @output Updated Ry(β), beta nodes, and weights.
     */
    void update_beta(Eigen::Tensor<doubleC, 3, Eigen::ColMajor> Ry_C3D_sp_sp_beta_, Eigen::VectorXd beta_F1D_beta_, Eigen::VectorXd weight_F1D_beta_);

    /**
     * @brief  Store gamma rotations, nodes, and quadrature weights.
     * @math   Rz ∈ ℂ^{Nsp×Nsp×Ngamma}; weights represent dγ.
     * @output Updated Rz(γ), gamma nodes, and weights.
     */
    void update_gamma(Eigen::Tensor<doubleC, 3, Eigen::ColMajor> Rz_C3D_sp_sp_gamma_, Eigen::VectorXd gamma_F1D_gamma_, Eigen::VectorXd weight_F1D_gamma_);

private:
    /**
     * @brief Build kernels using uniform gauge quadrature.
     * @math cache_ab(Ω) = ⟨Φ₁;a|O Pᴺ R(Ω)|Φ₂;b⟩.
     * @output result_C5D_cfg1_cfg2_gamma_alpha_beta, returned by value.
     * @note No K selection, Euler weights, or Wigner factors.
     */
    Eigen::Tensor<doubleC, 5, Eigen::ColMajor> build_amp(const std::function<void(Eigen::Ref<Eigen::MatrixXcd>)>& kernel_Func);

public:
    /**
     * @brief  Integrate cached kernels using Euler quadrature.
     * @math   result = (2I+1)/VΩ ΣΩ wΩ Dᴵ*(Ω) cache(Ω).
     * @output Projected tensor for the requested I, returned by value.
     * @note   Cache includes particle-number projection at each Euler node.
     */
    Eigen::Tensor<doubleC, 4, Eigen::ColMajor> integrate_amp(int TargetTwoI_I, const Eigen::Tensor<doubleC, 5, Eigen::ColMajor>& cache_C5D_cfg1_cfg2_gamma_alpha_beta);

    /**
     * @brief Integrate definite-K kernels using beta quadrature.
     * @math result_ab = (2I+1)/2 Σβ wβ dᴵ*_{KaKb}(β) cache_ab(β).
     * @output Matrix indexed by cfg1,cfg2, returned by value.
     * @note Definite K; cache at α=γ=0; wβ includes sinβ.
     */
    Eigen::MatrixXcd integrate_amp_with_K(int TargetTwoI_I, const Eigen::VectorXi& TwoK1_I1D_cfg1, const Eigen::VectorXi& TwoK2_I1D_cfg2, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& cache_C3D_cfg1_cfg2_beta);

    /**
     * @brief Build overlap kernels using uniform gauge quadrature.
     * @math cache_ab(Ω) = ⟨Φ₁;a|Pᴺ R(Ω)|Φ₂;b⟩.
     * @output result_C5D_cfg1_cfg2_gamma_alpha_beta, returned by value.
     * @note Excludes Euler weights and Wigner factors.
     */
    Eigen::Tensor<doubleC, 5, Eigen::ColMajor> build_overlap();

    /**
     * @brief Build one-body kernels using uniform gauge quadrature.
     * @math cache_ab(Ω) = ⟨Φ₁;a|H¹ Pᴺ R(Ω)|Φ₂;b⟩.
     * @output result_C5D_cfg1_cfg2_gamma_alpha_beta, returned by value.
     * @note Excludes Euler weights and Wigner factors.
     */
    Eigen::Tensor<doubleC, 5, Eigen::ColMajor> build_one_body(const Eigen::MatrixXd& OneBody_F2D_sp_sp);

    /**
     * @brief Build Q/P kernels using Pfaffians and gauge quadrature.
     * @math cache(Ω) = Σi gQi KQi(Ω) + Σi gPi KPi(Ω).
     * @output result_C5D_cfg1_cfg2_gamma_alpha_beta, returned by value.
     * @note Excludes Euler weights and Wigner factors; pairing ¼ included.
     */
    Eigen::Tensor<doubleC, 5, Eigen::ColMajor> build_two_body(const Eigen::VectorXd& gQ_F1D_i, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& Q_C3D_sp_sp_i, const Eigen::VectorXd& gP_F1D_i, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& P_C3D_sp_sp_i);

};

inline void HFBProjection::update_alpha(Eigen::Tensor<doubleC, 3, Eigen::ColMajor> Rz_C3D_sp_sp_alpha_, Eigen::VectorXd alpha_F1D_alpha_, Eigen::VectorXd weight_F1D_alpha_) {
    assert(Rz_C3D_sp_sp_alpha_.dimension(0) == Nsp_I && Rz_C3D_sp_sp_alpha_.dimension(1) == Nsp_I && Rz_C3D_sp_sp_alpha_.dimension(2) == Nalpha_I);
    assert(alpha_F1D_alpha_.size() == Nalpha_I && weight_F1D_alpha_.size() == Nalpha_I);

    // {Rz(α),α,wα} → stored quadrature data.
    Rz_C3D_sp_sp_alpha = std::move(Rz_C3D_sp_sp_alpha_);
    alpha_F1D_alpha = std::move(alpha_F1D_alpha_);
    weight_F1D_alpha = std::move(weight_F1D_alpha_);
}

inline void HFBProjection::update_beta(Eigen::Tensor<doubleC, 3, Eigen::ColMajor> Ry_C3D_sp_sp_beta_, Eigen::VectorXd beta_F1D_beta_, Eigen::VectorXd weight_F1D_beta_) {
    assert(Ry_C3D_sp_sp_beta_.dimension(0) == Nsp_I && Ry_C3D_sp_sp_beta_.dimension(1) == Nsp_I && Ry_C3D_sp_sp_beta_.dimension(2) == Nbeta_I);
    assert(beta_F1D_beta_.size() == Nbeta_I && weight_F1D_beta_.size() == Nbeta_I);

    // {Ry(β),β,wβ} → stored quadrature data.
    Ry_C3D_sp_sp_beta = std::move(Ry_C3D_sp_sp_beta_);
    beta_F1D_beta = std::move(beta_F1D_beta_);
    weight_F1D_beta = std::move(weight_F1D_beta_);
}

inline void HFBProjection::update_gamma(Eigen::Tensor<doubleC, 3, Eigen::ColMajor> Rz_C3D_sp_sp_gamma_, Eigen::VectorXd gamma_F1D_gamma_, Eigen::VectorXd weight_F1D_gamma_) {
    assert(Rz_C3D_sp_sp_gamma_.dimension(0) == Nsp_I && Rz_C3D_sp_sp_gamma_.dimension(1) == Nsp_I && Rz_C3D_sp_sp_gamma_.dimension(2) == Ngamma_I);
    assert(gamma_F1D_gamma_.size() == Ngamma_I && weight_F1D_gamma_.size() == Ngamma_I);

    // {Rz(γ),γ,wγ} → stored quadrature data.
    Rz_C3D_sp_sp_gamma = std::move(Rz_C3D_sp_sp_gamma_);
    gamma_F1D_gamma = std::move(gamma_F1D_gamma_);
    weight_F1D_gamma = std::move(weight_F1D_gamma_);
}

inline Eigen::Tensor<doubleC, 5, Eigen::ColMajor> HFBProjection::build_amp(const std::function<void(Eigen::Ref<Eigen::MatrixXcd>)>& kernel_Func) {
    assert(Nalpha_I > 0 && Nbeta_I > 0 && Ngamma_I > 0);
    assert(Rz_C3D_sp_sp_alpha.dimension(0) == Nsp_I && Rz_C3D_sp_sp_alpha.dimension(1) == Nsp_I && Rz_C3D_sp_sp_alpha.dimension(2) == Nalpha_I);
    assert(Ry_C3D_sp_sp_beta.dimension(0) == Nsp_I && Ry_C3D_sp_sp_beta.dimension(1) == Nsp_I && Ry_C3D_sp_sp_beta.dimension(2) == Nbeta_I);
    assert(Rz_C3D_sp_sp_gamma.dimension(0) == Nsp_I && Rz_C3D_sp_sp_gamma.dimension(1) == Nsp_I && Rz_C3D_sp_sp_gamma.dimension(2) == Ngamma_I);

    Eigen::Tensor<doubleC, 5, Eigen::ColMajor> cache_C5D_cfg1_cfg2_gamma_alpha_beta = Eigen::Tensor<doubleC, 5, Eigen::ColMajor>(Ncfg1_I, Ncfg2_I, Ngamma_I, Nalpha_I, Nbeta_I).setZero();


    // D = Rz(α) Ry(β) Rz(γ); cache(:,:,γ,α,β) = PNP kernel.
    for (int beta_I = 0; beta_I < Nbeta_I; ++beta_I) {
        const Eigen::Map<const Eigen::MatrixXcd> RyBeta_C2D_sp_sp(Ry_C3D_sp_sp_beta.data() + static_cast<Eigen::Index>(beta_I) * Nsp_I * Nsp_I, Nsp_I, Nsp_I);
        for (int alpha_I = 0; alpha_I < Nalpha_I; ++alpha_I) {
            const Eigen::Map<const Eigen::MatrixXcd> RzAlpha_C2D_sp_sp(Rz_C3D_sp_sp_alpha.data() + static_cast<Eigen::Index>(alpha_I) * Nsp_I * Nsp_I, Nsp_I, Nsp_I);
            for (int gamma_I = 0; gamma_I < Ngamma_I; ++gamma_I) {
                const Eigen::Map<const Eigen::MatrixXcd> RzGamma_C2D_sp_sp(Rz_C3D_sp_sp_gamma.data() + static_cast<Eigen::Index>(gamma_I) * Nsp_I * Nsp_I, Nsp_I, Nsp_I);
                D_C2D_sp_sp.noalias() = RzAlpha_C2D_sp_sp * RyBeta_C2D_sp_sp * RzGamma_C2D_sp_sp;
                Eigen::Map<Eigen::MatrixXcd> cache_C2D_cfg1_cfg2(cache_C5D_cfg1_cfg2_gamma_alpha_beta.data() + (gamma_I + static_cast<Eigen::Index>(Ngamma_I) * (alpha_I + static_cast<Eigen::Index>(Nalpha_I) * beta_I)) * Ncfg1_I * Ncfg2_I, Ncfg1_I, Ncfg2_I);
                kernel_Func(cache_C2D_cfg1_cfg2);
            }
        }
    }
    return cache_C5D_cfg1_cfg2_gamma_alpha_beta;
}

inline Eigen::Tensor<doubleC, 4, Eigen::ColMajor> HFBProjection::integrate_amp(int TargetTwoI_I, const Eigen::Tensor<doubleC, 5, Eigen::ColMajor>& cache_C5D_cfg1_cfg2_gamma_alpha_beta) {
    assert(TargetTwoI_I >= 0);
    assert(cache_C5D_cfg1_cfg2_gamma_alpha_beta.dimension(0) == Ncfg1_I && cache_C5D_cfg1_cfg2_gamma_alpha_beta.dimension(1) == Ncfg2_I && cache_C5D_cfg1_cfg2_gamma_alpha_beta.dimension(2) == Ngamma_I && cache_C5D_cfg1_cfg2_gamma_alpha_beta.dimension(3) == Nalpha_I && cache_C5D_cfg1_cfg2_gamma_alpha_beta.dimension(4) == Nbeta_I);
    RepresentationSpin wignerD_X(TargetTwoI_I);
    Eigen::Tensor<doubleC, 4, Eigen::ColMajor> result_C4D_cfg1_cfg2_K_K = Eigen::Tensor<doubleC, 4, Eigen::ColMajor>(Ncfg1_I, Ncfg2_I, TargetTwoI_I + 1, TargetTwoI_I + 1).setZero();


    // (2I+1)/VΩ; beta weights contain sinβ dβ.
    const double volume_F = weight_F1D_alpha.sum() * weight_F1D_beta.sum() * weight_F1D_gamma.sum();
    assert(std::isfinite(volume_F) && volume_F > 0.0);
    const double normalization_F = (TargetTwoI_I + 1.0) / volume_F;

    const Eigen::Index Ncfg1cfg2_I = static_cast<Eigen::Index>(Ncfg1_I) * Ncfg2_I;
    Eigen::MatrixXcd weightD_C2D_K_K(TargetTwoI_I + 1, TargetTwoI_I + 1);

    // result = (2I+1)/VΩ ΣΩ wΩ cache(Ω) vec(Dᴵ*(Ω))ᵀ.
    for (int beta_I = 0; beta_I < Nbeta_I; ++beta_I) {
        for (int alpha_I = 0; alpha_I < Nalpha_I; ++alpha_I) {
            for (int gamma_I = 0; gamma_I < Ngamma_I; ++gamma_I) {
                const double weight_F = normalization_F * weight_F1D_alpha(alpha_I) * weight_F1D_beta(beta_I) * weight_F1D_gamma(gamma_I);
                weightD_C2D_K_K = weight_F * wignerD_X.calc_R(alpha_F1D_alpha(alpha_I), beta_F1D_beta(beta_I), gamma_F1D_gamma(gamma_I)).conjugate();

                // ColMajor: cfg1 → cfg2 → γ → α → β.
                const Eigen::Index offset_I = Ncfg1cfg2_I * (gamma_I + static_cast<Eigen::Index>(Ngamma_I) * alpha_I + static_cast<Eigen::Index>(Ngamma_I) * Nalpha_I * beta_I);
                Eigen::Map<Eigen::MatrixXcd> result_C2D_cfg1cfg2_KK(result_C4D_cfg1_cfg2_K_K.data(), Ncfg1cfg2_I, static_cast<Eigen::Index>(TargetTwoI_I + 1) * (TargetTwoI_I + 1));
                const Eigen::Map<const Eigen::VectorXcd> weightD_C1D_KK(weightD_C2D_K_K.data(), weightD_C2D_K_K.size());
                const Eigen::Map<const Eigen::VectorXcd> cache_C1D_cfg1cfg2(cache_C5D_cfg1_cfg2_gamma_alpha_beta.data() + offset_I, Ncfg1cfg2_I);
                result_C2D_cfg1cfg2_KK.noalias() += cache_C1D_cfg1cfg2 * weightD_C1D_KK.transpose();
            }
        }
    }
    return result_C4D_cfg1_cfg2_K_K;
}

inline Eigen::MatrixXcd HFBProjection::integrate_amp_with_K(int TargetTwoI_I, const Eigen::VectorXi& TwoK1_I1D_cfg1, const Eigen::VectorXi& TwoK2_I1D_cfg2, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& cache_C3D_cfg1_cfg2_beta) {
    assert(TargetTwoI_I >= 0);
    assert(cache_C3D_cfg1_cfg2_beta.dimension(0) == Ncfg1_I && cache_C3D_cfg1_cfg2_beta.dimension(1) == Ncfg2_I && cache_C3D_cfg1_cfg2_beta.dimension(2) == Nbeta_I);
    assert(TwoK1_I1D_cfg1.size() == Ncfg1_I && TwoK2_I1D_cfg2.size() == Ncfg2_I);
    const double normalization_F = (TargetTwoI_I + 1.0) / weight_F1D_beta.sum();
    RepresentationSpin wignerD_X(TargetTwoI_I);
    Eigen::MatrixXcd result_C2D_cfg1_cfg2 = Eigen::MatrixXcd::Zero(Ncfg1_I, Ncfg2_I);

    // result_ab = (2I+1)/2 Σβ wβ dᴵ*_{KaKb}(β) cache_ab(β).
    for (int beta_I = 0; beta_I < Nbeta_I; ++beta_I) {
        const auto& wignerDyBeta_C2D_K_K = wignerD_X.calc_Ry(beta_F1D_beta(beta_I));
        const double weight_F = normalization_F * weight_F1D_beta(beta_I);
        for (int cfg2_I = 0; cfg2_I < Ncfg2_I; ++cfg2_I) {
            const int TwoK2_I = TwoK2_I1D_cfg2(cfg2_I);
            if (std::abs(TwoK2_I) > TargetTwoI_I || (TwoK2_I + TargetTwoI_I) % 2 != 0) { continue; }
            const int K2_I = (TwoK2_I + TargetTwoI_I) / 2;
            for (int cfg1_I = 0; cfg1_I < Ncfg1_I; ++cfg1_I) {
                const int TwoK1_I = TwoK1_I1D_cfg1(cfg1_I);
                if (std::abs(TwoK1_I) > TargetTwoI_I || (TwoK1_I + TargetTwoI_I) % 2 != 0) { continue; }
                const int K1_I = (TwoK1_I + TargetTwoI_I) / 2;
                // dᴵ(β) is real in the Jz basis.
                const double weightD_F = weight_F * wignerDyBeta_C2D_K_K(K1_I, K2_I).real();
                result_C2D_cfg1_cfg2(cfg1_I, cfg2_I) += weightD_F * cache_C3D_cfg1_cfg2_beta(cfg1_I, cfg2_I, beta_I);
            }
        }
    }
    return result_C2D_cfg1_cfg2;
}

inline Eigen::Tensor<doubleC, 5, Eigen::ColMajor> HFBProjection::build_overlap() {
    return build_amp([&](Eigen::Ref<Eigen::MatrixXcd> cache_C2D_cfg1_cfg2) {
        cache_C2D_cfg1_cfg2 = calc_overlap_pnp();
    });
}

inline Eigen::Tensor<doubleC, 5, Eigen::ColMajor> HFBProjection::build_one_body(const Eigen::MatrixXd& OneBody_F2D_sp_sp) {
    assert(OneBody_F2D_sp_sp.rows() == Nsp_I && OneBody_F2D_sp_sp.cols() == Nsp_I && OneBody_F2D_sp_sp.allFinite());
    const Eigen::MatrixXcd OneBody_C2D_sp_sp = OneBody_F2D_sp_sp.cast<doubleC>();
    return build_amp([&](Eigen::Ref<Eigen::MatrixXcd> cache_C2D_cfg1_cfg2) {
        cache_C2D_cfg1_cfg2 = calc_one_body_pnp(OneBody_C2D_sp_sp);
    });
}

inline Eigen::Tensor<doubleC, 5, Eigen::ColMajor> HFBProjection::build_two_body(const Eigen::VectorXd& gQ_F1D_i, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& Q_C3D_sp_sp_i, const Eigen::VectorXd& gP_F1D_i, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& P_C3D_sp_sp_i) {
    return build_amp([&](Eigen::Ref<Eigen::MatrixXcd> cache_C2D_cfg1_cfg2) {
        cache_C2D_cfg1_cfg2 = calc_two_body_pnp(gQ_F1D_i, Q_C3D_sp_sp_i, gP_F1D_i, P_C3D_sp_sp_i);
    });
}
