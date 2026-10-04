/**
 * @file    hfb_projection_nucleus.hpp
 * @author  cailea
 * @date    2026-09-17
 * @brief   Joint neutron-proton projection interfaces.
 */

#pragma once

#include <cassert>
#include <cmath>
#include <functional>
#include <utility>

#include "hfb_projection.hpp"

/**
 * @brief  Project neutron-proton configurations using gauge and Euler quadrature.
 * @math   |Φκ⟩ = |Φν,κν⟩ ⊗ |Φπ,κπ⟩.
 * @note   Species-local indices start at zero; Euler meshes are shared.
 * @note   Particle-number targets apply to the right state.
 * @note   Requires even reference vacua and initialized projection data.
 * @note   Vacuum overlaps must remain nonzero at quadrature nodes.
 */
class HFBProjectionNucleus {
public:
    HFBProjection projection_neutron;
    HFBProjection projection_proton;

    Eigen::Tensor<doubleC, 3, Eigen::ColMajor> workspacen_C3D_cfgn1_cfgn2_operatorn{};
    Eigen::Tensor<doubleC, 3, Eigen::ColMajor> workspacep_C3D_cfgp1_cfgp2_operatorp{};

private:
    Eigen::Index Ncfgn1cfgn2_I = 0;
    Eigen::Index Ncfgp1cfgp2_I = 0;
    Eigen::Index Ncfgn1cfgn2cfgp1cfgp2_I = 0;

public:
    /**
     * @brief  Construct empty neutron-proton projection workspaces.
     * @math   Cₙ = Cₚ = ∅.
     * @output Empty workspaces.
     * @note   Assign a configured object before calculations.
     */
    HFBProjectionNucleus() = default;

    /**
     * @brief  Move species projections into owned workspaces.
     * @math   |Φκ⟩ = |Φν,κν⟩ ⊗ |Φπ,κπ⟩.
     * @output Object owning neutron and proton projections.
     */
    HFBProjectionNucleus(HFBProjection projection_neutron_, HFBProjection projection_proton_)
    : projection_neutron(std::move(projection_neutron_)), projection_proton(std::move(projection_proton_)) {
        assert(projection_neutron.Nalpha_I == projection_proton.Nalpha_I);
        assert(projection_neutron.Nbeta_I == projection_proton.Nbeta_I);
        assert(projection_neutron.Ngamma_I == projection_proton.Ngamma_I);

        Ncfgn1cfgn2_I = static_cast<Eigen::Index>(projection_neutron.Ncfg1_I) * projection_neutron.Ncfg2_I;
        Ncfgp1cfgp2_I = static_cast<Eigen::Index>(projection_proton.Ncfg1_I) * projection_proton.Ncfg2_I;
        Ncfgn1cfgn2cfgp1cfgp2_I = Ncfgn1cfgn2_I * Ncfgp1cfgp2_I;

    }

    /**
     * @brief Build nuclear kernels using gauge quadrature.
     * @math cache_ab(Ω) = ⟨Φ₁;a|O Pνᴺ Pπᶻ R(Ω)|Φ₂;b⟩.
     * @output Overwritten caller-provided result_C8D_cfgn1_cfgn2_cfgp1_cfgp2_operator_gamma_alpha_beta.
     * @note External captures are read-only; excludes K/weights/Wigner factors.
     */
    void build_amp(Eigen::Tensor<doubleC, 8, Eigen::ColMajor>& result_C8D_cfgn1_cfgn2_cfgp1_cfgp2_operator_gamma_alpha_beta, const std::function<void(HFBProjectionNucleus&, Eigen::TensorMap<Eigen::Tensor<doubleC, 5, Eigen::ColMajor>>&)>& calc_at_Omega_Func);

    /**
     * @brief Integrate cached nuclear kernels using Euler quadrature.
     * @math result = (2I+1)/VΩ ΣΩ wΩ Dᴵ*(Ω) cache(Ω).
     * @output Overwritten caller-provided result_C7D_cfgn1_cfgn2_cfgp1_cfgp2_operator_K_K.
     * @note Cache includes both species' particle-number projections.
     */
    void integrate_amp(Eigen::Tensor<doubleC, 7, Eigen::ColMajor>& result_C7D_cfgn1_cfgn2_cfgp1_cfgp2_operator_K_K, int TargetTwoI_I, const Eigen::Tensor<doubleC, 8, Eigen::ColMajor>& cache_C8D_cfgn1_cfgn2_cfgp1_cfgp2_operator_gamma_alpha_beta);

    /**
     * @brief Integrate axial kernels using beta quadrature.
     * @math M = (2I+1)/2 Σβ wβ dᴵ*_{K₁K₂}(β) cache(β).
     * @output Overwritten caller-provided result_C5D_cfgn1_cfgn2_cfgp1_cfgp2_operator.
     * @note Definite K; cache at α=γ=0; wβ includes sinβ.
     */
    void integrate_amp_with_K(Eigen::Tensor<doubleC, 5, Eigen::ColMajor>& result_C5D_cfgn1_cfgn2_cfgp1_cfgp2_operator, int TargetTwoI_I, const Eigen::VectorXi& TwoKn1_I1D_cfgn1, const Eigen::VectorXi& TwoKn2_I1D_cfgn2, const Eigen::VectorXi& TwoKp1_I1D_cfgp1, const Eigen::VectorXi& TwoKp2_I1D_cfgp2, const Eigen::Tensor<doubleC, 6, Eigen::ColMajor>& cache_C6D_cfgn1_cfgn2_cfgp1_cfgp2_operator_beta);

    /**
     * @brief Build nuclear overlaps using gauge quadrature.
     * @math result_n₁n₂p₁p₂(Ω) = Nν_n₁n₂(Ω) Nπ_p₁p₂(Ω).
     * @output Euler-grid tensor (n₁,n₂,p₁,p₂,γ,α,β), returned by value.
     */
    Eigen::Tensor<doubleC, 7, Eigen::ColMajor> build_overlap();

    /**
     * @brief Build nn kernels using Pfaffians and gauge quadrature.
     * @math result(Ω) = Nπ(Ω) Hν¹(Ω).
     * @output Euler-grid tensor (n₁,n₂,p₁,p₂,γ,α,β), returned by value.
     */
    Eigen::Tensor<doubleC, 7, Eigen::ColMajor> build_one_body_nn(const Eigen::MatrixXcd& OneBody_C2D_spn_spn);

    /**
     * @brief Build pp kernels using Pfaffians and gauge quadrature.
     * @math result(Ω) = Hπ¹(Ω) Nν(Ω).
     * @output Euler-grid tensor (n₁,n₂,p₁,p₂,γ,α,β), returned by value.
     */
    Eigen::Tensor<doubleC, 7, Eigen::ColMajor> build_one_body_pp(const Eigen::MatrixXcd& OneBody_C2D_spp_spp);

    /**
     * @brief Build pn kernels using Pfaffians and gauge quadrature.
     * @math O_pn = Σαβ O_αβ cπα†cνβ.
     * @output Euler-grid tensor (n₁,n₂,p₁,p₂,γ,α,β), returned by value.
     */
    Eigen::Tensor<doubleC, 7, Eigen::ColMajor> build_one_body_pn(const Eigen::MatrixXcd& OneBody_C2D_spp_spn);

    /**
     * @brief Build nn kernels using Pfaffians and gauge quadrature.
     * @math result(Ω) = Nπ(Ω) Hν²(Ω).
     * @output Euler-grid tensor (n₁,n₂,p₁,p₂,γ,α,β), returned by value.
     */
    Eigen::Tensor<doubleC, 7, Eigen::ColMajor> build_two_body_nn(const Eigen::VectorXd& gQ_F1D_i, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& Q_C3D_spn_spn_i, const Eigen::VectorXd& gP_F1D_i, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& P_C3D_spn_spn_i);

    /**
     * @brief Build pp kernels using Pfaffians and gauge quadrature.
     * @math result(Ω) = Hπ²(Ω) Nν(Ω).
     * @output Euler-grid tensor (n₁,n₂,p₁,p₂,γ,α,β), returned by value.
     */
    Eigen::Tensor<doubleC, 7, Eigen::ColMajor> build_two_body_pp(const Eigen::VectorXd& gQ_F1D_i, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& Q_C3D_spp_spp_i, const Eigen::VectorXd& gP_F1D_i, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& P_C3D_spp_spp_i);

    /**
     * @brief Build pn kernels using Pfaffians and gauge quadrature.
     * @math Hpn = ½ Σi gQi (Qνi† Qπi + Qνi Qπi†).
     * @output Euler-grid tensor (n₁,n₂,p₁,p₂,γ,α,β), returned by value.
     */
    Eigen::Tensor<doubleC, 7, Eigen::ColMajor> build_two_body_pn(const Eigen::VectorXd& gQ_F1D_i, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& Q_C3D_spn_spn_i, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& Q_C3D_spp_spp_i);

};

inline void HFBProjectionNucleus::build_amp(Eigen::Tensor<doubleC, 8, Eigen::ColMajor>& result_C8D_cfgn1_cfgn2_cfgp1_cfgp2_operator_gamma_alpha_beta, const std::function<void(HFBProjectionNucleus&, Eigen::TensorMap<Eigen::Tensor<doubleC, 5, Eigen::ColMajor>>&)>& calc_at_Omega_Func) {
    // Euler grids are populated after construction.
    assert(projection_neutron.alpha_F1D_alpha.isApprox(projection_proton.alpha_F1D_alpha));
    assert(projection_neutron.beta_F1D_beta.isApprox(projection_proton.beta_F1D_beta));
    assert(projection_neutron.gamma_F1D_gamma.isApprox(projection_proton.gamma_F1D_gamma));
    assert(projection_neutron.weight_F1D_alpha.isApprox(projection_proton.weight_F1D_alpha));
    assert(projection_neutron.weight_F1D_beta.isApprox(projection_proton.weight_F1D_beta));
    assert(projection_neutron.weight_F1D_gamma.isApprox(projection_proton.weight_F1D_gamma));

    const Eigen::Index Noperator_I = result_C8D_cfgn1_cfgn2_cfgp1_cfgp2_operator_gamma_alpha_beta.dimension(4);

    // cache ∈ ℂ^{cfgn1×cfgn2×cfgp1×cfgp2×Noperator×Nγ×Nα×Nβ}.
    assert((result_C8D_cfgn1_cfgn2_cfgp1_cfgp2_operator_gamma_alpha_beta.dimensions() == Eigen::array<Eigen::Index, 8>{projection_neutron.Ncfg1_I, projection_neutron.Ncfg2_I, projection_proton.Ncfg1_I, projection_proton.Ncfg2_I, Noperator_I, projection_neutron.Ngamma_I, projection_neutron.Nalpha_I, projection_neutron.Nbeta_I}));
    result_C8D_cfgn1_cfgn2_cfgp1_cfgp2_operator_gamma_alpha_beta.setZero();

    // Ω → Dν,Dπ → cache(Ω).
    const auto calc_at_beta = [&](int beta_I, HFBProjectionNucleus& nucleus_thread) {
        for (int alpha_I = 0; alpha_I < nucleus_thread.projection_neutron.Nalpha_I; ++alpha_I) {
            for (int gamma_I = 0; gamma_I < nucleus_thread.projection_neutron.Ngamma_I; ++gamma_I) {
                for (auto* projection : {&nucleus_thread.projection_neutron, &nucleus_thread.projection_proton}) {
                    // D = Rz(α) Ry(β) Rz(γ).
                    const Eigen::Map<const Eigen::MatrixXcd> RzAlpha_C2D_sp_sp(projection->Rz_C3D_sp_sp_alpha.data() + static_cast<Eigen::Index>(alpha_I) * projection->Nsp_I * projection->Nsp_I, projection->Nsp_I, projection->Nsp_I);
                    const Eigen::Map<const Eigen::MatrixXcd> RyBeta_C2D_sp_sp(projection->Ry_C3D_sp_sp_beta.data() + static_cast<Eigen::Index>(beta_I) * projection->Nsp_I * projection->Nsp_I, projection->Nsp_I, projection->Nsp_I);
                    const Eigen::Map<const Eigen::MatrixXcd> RzGamma_C2D_sp_sp(projection->Rz_C3D_sp_sp_gamma.data() + static_cast<Eigen::Index>(gamma_I) * projection->Nsp_I * projection->Nsp_I, projection->Nsp_I, projection->Nsp_I);
                    projection->D_C2D_sp_sp.noalias() = RzAlpha_C2D_sp_sp * RyBeta_C2D_sp_sp * RzGamma_C2D_sp_sp;
                }
                // ColMajor: cfgn1 → cfgn2 → cfgp1 → cfgp2 → operator → γ → α → β.
                const Eigen::Index offset_I = Ncfgn1cfgn2cfgp1cfgp2_I * Noperator_I * (gamma_I + static_cast<Eigen::Index>(nucleus_thread.projection_neutron.Ngamma_I) * alpha_I + static_cast<Eigen::Index>(nucleus_thread.projection_neutron.Ngamma_I) * nucleus_thread.projection_neutron.Nalpha_I * beta_I);
                Eigen::TensorMap<Eigen::Tensor<doubleC, 5, Eigen::ColMajor>> resultAtOmega_C5D_cfgn1_cfgn2_cfgp1_cfgp2_operator(result_C8D_cfgn1_cfgn2_cfgp1_cfgp2_operator_gamma_alpha_beta.data() + offset_I, nucleus_thread.projection_neutron.Ncfg1_I, nucleus_thread.projection_neutron.Ncfg2_I, nucleus_thread.projection_proton.Ncfg1_I, nucleus_thread.projection_proton.Ncfg2_I, Noperator_I);
                calc_at_Omega_Func(nucleus_thread, resultAtOmega_C5D_cfgn1_cfgn2_cfgp1_cfgp2_operator);
            }
        }
    };

    #pragma omp parallel
    {
        HFBProjectionNucleus nucleus_thread = *this;

        #pragma omp for schedule(static)
        for (int beta_I = 0; beta_I < projection_neutron.Nbeta_I; ++beta_I) {
            calc_at_beta(beta_I, nucleus_thread);
        }
    }
}

inline void HFBProjectionNucleus::integrate_amp(Eigen::Tensor<doubleC, 7, Eigen::ColMajor>& result_C7D_cfgn1_cfgn2_cfgp1_cfgp2_operator_K_K, int TargetTwoI_I, const Eigen::Tensor<doubleC, 8, Eigen::ColMajor>& cache_C8D_cfgn1_cfgn2_cfgp1_cfgp2_operator_gamma_alpha_beta) {
    const Eigen::Index Noperator_I = result_C7D_cfgn1_cfgn2_cfgp1_cfgp2_operator_K_K.dimension(4);

    // VΩ = Σwα Σwβ Σwγ; wβ includes sinβ.
    assert(TargetTwoI_I >= 0);
    const double volume_F = projection_neutron.weight_F1D_alpha.sum() * projection_neutron.weight_F1D_beta.sum() * projection_neutron.weight_F1D_gamma.sum();
    const double normalization_F = (TargetTwoI_I + 1.0) / volume_F;
    assert(std::isfinite(normalization_F) && normalization_F > 0.0);

    assert((cache_C8D_cfgn1_cfgn2_cfgp1_cfgp2_operator_gamma_alpha_beta.dimensions() == Eigen::array<Eigen::Index, 8>{projection_neutron.Ncfg1_I, projection_neutron.Ncfg2_I, projection_proton.Ncfg1_I, projection_proton.Ncfg2_I, Noperator_I, projection_neutron.Ngamma_I, projection_neutron.Nalpha_I, projection_neutron.Nbeta_I}));
    assert((result_C7D_cfgn1_cfgn2_cfgp1_cfgp2_operator_K_K.dimensions() == Eigen::array<Eigen::Index, 7>{projection_neutron.Ncfg1_I, projection_neutron.Ncfg2_I, projection_proton.Ncfg1_I, projection_proton.Ncfg2_I, Noperator_I, TargetTwoI_I + 1, TargetTwoI_I + 1}));
    result_C7D_cfgn1_cfgn2_cfgp1_cfgp2_operator_K_K.setZero();
    RepresentationSpin wignerD_X(TargetTwoI_I);
    // W(K₁,K₂,Ω) = (2I+1)/VΩ wΩ Dᴵ*(Ω).
    Eigen::Tensor<doubleC, 5, Eigen::ColMajor> weightD_C5D_K_K_gamma_alpha_beta(TargetTwoI_I + 1, TargetTwoI_I + 1, projection_neutron.Ngamma_I, projection_neutron.Nalpha_I, projection_neutron.Nbeta_I);
    for (int beta_I = 0; beta_I < projection_neutron.Nbeta_I; ++beta_I) {
        for (int alpha_I = 0; alpha_I < projection_neutron.Nalpha_I; ++alpha_I) {
            for (int gamma_I = 0; gamma_I < projection_neutron.Ngamma_I; ++gamma_I) {
                const double weight_F = normalization_F * projection_neutron.weight_F1D_alpha(alpha_I) * projection_neutron.weight_F1D_beta(beta_I) * projection_neutron.weight_F1D_gamma(gamma_I);
                const Eigen::Index offset_I = (gamma_I + static_cast<Eigen::Index>(projection_neutron.Ngamma_I) * (alpha_I + static_cast<Eigen::Index>(projection_neutron.Nalpha_I) * beta_I)) * (TargetTwoI_I + 1) * (TargetTwoI_I + 1);
                Eigen::Map<Eigen::MatrixXcd> weightD_C2D_K_K(weightD_C5D_K_K_gamma_alpha_beta.data() + offset_I, TargetTwoI_I + 1, TargetTwoI_I + 1);
                weightD_C2D_K_K = weight_F * wignerD_X.calc_R(projection_neutron.alpha_F1D_alpha(alpha_I), projection_neutron.beta_F1D_beta(beta_I), projection_neutron.gamma_F1D_gamma(gamma_I)).conjugate();
            }
        }
    }

    // R = C Wᵀ; Ω = γ + Nγ(α + Nαβ).
    const Eigen::Index NOmega_I = static_cast<Eigen::Index>(projection_neutron.Ngamma_I) * projection_neutron.Nalpha_I * projection_neutron.Nbeta_I;
    const Eigen::Index NKK_I = static_cast<Eigen::Index>(TargetTwoI_I + 1) * (TargetTwoI_I + 1);
    const Eigen::Map<const Eigen::MatrixXcd> cache_C2D_cfgn1cfgn2cfgp1cfgp2operator_Omega(cache_C8D_cfgn1_cfgn2_cfgp1_cfgp2_operator_gamma_alpha_beta.data(), Ncfgn1cfgn2cfgp1cfgp2_I * Noperator_I, NOmega_I);
    const Eigen::Map<const Eigen::MatrixXcd> weightD_C2D_KK_Omega(weightD_C5D_K_K_gamma_alpha_beta.data(), NKK_I, NOmega_I);

    Eigen::Map<Eigen::MatrixXcd> result_C2D_cfgn1cfgn2cfgp1cfgp2operator_KK(result_C7D_cfgn1_cfgn2_cfgp1_cfgp2_operator_K_K.data(), Ncfgn1cfgn2cfgp1cfgp2_I * Noperator_I, NKK_I);
    result_C2D_cfgn1cfgn2cfgp1cfgp2operator_KK.noalias() = cache_C2D_cfgn1cfgn2cfgp1cfgp2operator_Omega * weightD_C2D_KK_Omega.transpose();
}

inline void HFBProjectionNucleus::integrate_amp_with_K(Eigen::Tensor<doubleC, 5, Eigen::ColMajor>& result_C5D_cfgn1_cfgn2_cfgp1_cfgp2_operator, int TargetTwoI_I, const Eigen::VectorXi& TwoKn1_I1D_cfgn1, const Eigen::VectorXi& TwoKn2_I1D_cfgn2, const Eigen::VectorXi& TwoKp1_I1D_cfgp1, const Eigen::VectorXi& TwoKp2_I1D_cfgp2, const Eigen::Tensor<doubleC, 6, Eigen::ColMajor>& cache_C6D_cfgn1_cfgn2_cfgp1_cfgp2_operator_beta) {
    const Eigen::Index Noperator_I = result_C5D_cfgn1_cfgn2_cfgp1_cfgp2_operator.dimension(4);
    assert(TargetTwoI_I >= 0);
    assert(TwoKn1_I1D_cfgn1.size() == projection_neutron.Ncfg1_I && TwoKp1_I1D_cfgp1.size() == projection_proton.Ncfg1_I);
    assert(TwoKn2_I1D_cfgn2.size() == projection_neutron.Ncfg2_I && TwoKp2_I1D_cfgp2.size() == projection_proton.Ncfg2_I);
    assert((cache_C6D_cfgn1_cfgn2_cfgp1_cfgp2_operator_beta.dimensions() == Eigen::array<Eigen::Index, 6>{projection_neutron.Ncfg1_I, projection_neutron.Ncfg2_I, projection_proton.Ncfg1_I, projection_proton.Ncfg2_I, Noperator_I, projection_neutron.Nbeta_I}));
    assert((result_C5D_cfgn1_cfgn2_cfgp1_cfgp2_operator.dimensions() == Eigen::array<Eigen::Index, 5>{projection_neutron.Ncfg1_I, projection_neutron.Ncfg2_I, projection_proton.Ncfg1_I, projection_proton.Ncfg2_I, Noperator_I}));
    result_C5D_cfgn1_cfgn2_cfgp1_cfgp2_operator.setZero();
    const double normalization_F = (TargetTwoI_I + 1.0) / projection_neutron.weight_F1D_beta.sum();
    assert(std::isfinite(normalization_F) && normalization_F > 0.0);
    RepresentationSpin wignerD_X(TargetTwoI_I);

    // K₁ = Kn₁+Kp₁; K₂ = Kn₂+Kp₂.
    Eigen::Tensor<double, 3, Eigen::ColMajor> weightD_F3D_K_K_beta(TargetTwoI_I + 1, TargetTwoI_I + 1, projection_neutron.Nbeta_I);
    for (int beta_I = 0; beta_I < projection_neutron.Nbeta_I; ++beta_I) {
        const auto& wignerDyBeta_C2D_K_K = wignerD_X.calc_Ry(projection_neutron.beta_F1D_beta(beta_I));
        Eigen::Map<Eigen::MatrixXd> weightD_F2D_K_K(weightD_F3D_K_K_beta.data() + static_cast<Eigen::Index>(beta_I) * (TargetTwoI_I + 1) * (TargetTwoI_I + 1), TargetTwoI_I + 1, TargetTwoI_I + 1);
        weightD_F2D_K_K = normalization_F * projection_neutron.weight_F1D_beta(beta_I) * wignerDyBeta_C2D_K_K.real();
    }

    // Each right configuration pair owns disjoint output elements.
    #pragma omp parallel for collapse(2) schedule(static)
    for (int cfgp2_I = 0; cfgp2_I < projection_proton.Ncfg2_I; ++cfgp2_I) {
        for (int cfgn2_I = 0; cfgn2_I < projection_neutron.Ncfg2_I; ++cfgn2_I) {
            const int TwoK2_I = TwoKn2_I1D_cfgn2(cfgn2_I) + TwoKp2_I1D_cfgp2(cfgp2_I);
            if (std::abs(TwoK2_I) > TargetTwoI_I || (TwoK2_I + TargetTwoI_I) % 2 != 0) { continue; }
            const int K2_I = (TwoK2_I + TargetTwoI_I) / 2;
            for (int beta_I = 0; beta_I < projection_neutron.Nbeta_I; ++beta_I) {
                for (int cfgp1_I = 0; cfgp1_I < projection_proton.Ncfg1_I; ++cfgp1_I) {
                    for (int cfgn1_I = 0; cfgn1_I < projection_neutron.Ncfg1_I; ++cfgn1_I) {
                        const int TwoK1_I = TwoKn1_I1D_cfgn1(cfgn1_I) + TwoKp1_I1D_cfgp1(cfgp1_I);
                        if (std::abs(TwoK1_I) > TargetTwoI_I || (TwoK1_I + TargetTwoI_I) % 2 != 0) { continue; }
                        const int K1_I = (TwoK1_I + TargetTwoI_I) / 2;
                        const double weightD_F = weightD_F3D_K_K_beta(K1_I, K2_I, beta_I);
                        for (Eigen::Index operator_I = 0; operator_I < Noperator_I; ++operator_I) { result_C5D_cfgn1_cfgn2_cfgp1_cfgp2_operator(cfgn1_I, cfgn2_I, cfgp1_I, cfgp2_I, operator_I) += weightD_F * cache_C6D_cfgn1_cfgn2_cfgp1_cfgp2_operator_beta(cfgn1_I, cfgn2_I, cfgp1_I, cfgp2_I, operator_I, beta_I); }
                    }
                }
            }
        }
    }
}

inline Eigen::Tensor<doubleC, 7, Eigen::ColMajor> HFBProjectionNucleus::build_overlap() {
    Eigen::Tensor<doubleC, 8, Eigen::ColMajor> result_C8D_cfgn1_cfgn2_cfgp1_cfgp2_operator_gamma_alpha_beta(projection_neutron.Ncfg1_I, projection_neutron.Ncfg2_I, projection_proton.Ncfg1_I, projection_proton.Ncfg2_I, 1, projection_neutron.Ngamma_I, projection_neutron.Nalpha_I, projection_neutron.Nbeta_I);
    // Size PNP workspaces before thread copies.
    workspacen_C3D_cfgn1_cfgn2_operatorn.resize(projection_neutron.Ncfg1_I, projection_neutron.Ncfg2_I, 1);
    workspacep_C3D_cfgp1_cfgp2_operatorp.resize(projection_proton.Ncfg1_I, projection_proton.Ncfg2_I, 1);
    build_amp(result_C8D_cfgn1_cfgn2_cfgp1_cfgp2_operator_gamma_alpha_beta, [&](HFBProjectionNucleus& nucleus_thread, auto& resultAtOmega_C5D_cfgn1_cfgn2_cfgp1_cfgp2_operator) {

        Eigen::Map<Eigen::MatrixXcd> result_C2D_cfgn1cfgn2_cfgp1cfgp2(resultAtOmega_C5D_cfgn1_cfgn2_cfgp1_cfgp2_operator.data(), Ncfgn1cfgn2_I, Ncfgp1cfgp2_I);
        nucleus_thread.projection_neutron.integrate_pnp(Eigen::TensorMap<Eigen::Tensor<doubleC, 3, Eigen::ColMajor>>(nucleus_thread.workspacen_C3D_cfgn1_cfgn2_operatorn.data(), nucleus_thread.workspacen_C3D_cfgn1_cfgn2_operatorn.dimensions()), nucleus_thread.projection_neutron.TargetN_I, [&](auto& resultAtPhi_C3D_cfg1_cfg2_operator) {
            Eigen::Map<Eigen::MatrixXcd> result_C2D_cfg1_cfg2(resultAtPhi_C3D_cfg1_cfg2_operator.data(), nucleus_thread.projection_neutron.Ncfg1_I, nucleus_thread.projection_neutron.Ncfg2_I);
            nucleus_thread.projection_neutron.hfb_pfaffian.calc_overlap(result_C2D_cfg1_cfg2);
        });
        nucleus_thread.projection_proton.integrate_pnp(Eigen::TensorMap<Eigen::Tensor<doubleC, 3, Eigen::ColMajor>>(nucleus_thread.workspacep_C3D_cfgp1_cfgp2_operatorp.data(), nucleus_thread.workspacep_C3D_cfgp1_cfgp2_operatorp.dimensions()), nucleus_thread.projection_proton.TargetN_I, [&](auto& resultAtPhi_C3D_cfg1_cfg2_operator) {
            Eigen::Map<Eigen::MatrixXcd> result_C2D_cfg1_cfg2(resultAtPhi_C3D_cfg1_cfg2_operator.data(), nucleus_thread.projection_proton.Ncfg1_I, nucleus_thread.projection_proton.Ncfg2_I);
            nucleus_thread.projection_proton.hfb_pfaffian.calc_overlap(result_C2D_cfg1_cfg2);
        });

        // ColMajor: (n₁,n₂) → rows; (p₁,p₂) → columns.
        const Eigen::Map<const Eigen::VectorXcd> overlap_C1D_cfgn1cfgn2(nucleus_thread.workspacen_C3D_cfgn1_cfgn2_operatorn.data(), Ncfgn1cfgn2_I);
        const Eigen::Map<const Eigen::VectorXcd> overlap_C1D_cfgp1cfgp2(nucleus_thread.workspacep_C3D_cfgp1_cfgp2_operatorp.data(), Ncfgp1cfgp2_I);
        result_C2D_cfgn1cfgn2_cfgp1cfgp2.noalias() = overlap_C1D_cfgn1cfgn2 * overlap_C1D_cfgp1cfgp2.transpose();
    });
    return result_C8D_cfgn1_cfgn2_cfgp1_cfgp2_operator_gamma_alpha_beta.chip(0, 4);
}

inline Eigen::Tensor<doubleC, 7, Eigen::ColMajor> HFBProjectionNucleus::build_one_body_nn(const Eigen::MatrixXcd& OneBody_C2D_spn_spn) {
    assert(OneBody_C2D_spn_spn.rows() == projection_neutron.Nsp_I && OneBody_C2D_spn_spn.cols() == projection_neutron.Nsp_I && OneBody_C2D_spn_spn.allFinite());

    Eigen::Tensor<doubleC, 8, Eigen::ColMajor> result_C8D_cfgn1_cfgn2_cfgp1_cfgp2_operator_gamma_alpha_beta(projection_neutron.Ncfg1_I, projection_neutron.Ncfg2_I, projection_proton.Ncfg1_I, projection_proton.Ncfg2_I, 1, projection_neutron.Ngamma_I, projection_neutron.Nalpha_I, projection_neutron.Nbeta_I);
    // Size PNP workspaces before thread copies.
    workspacen_C3D_cfgn1_cfgn2_operatorn.resize(projection_neutron.Ncfg1_I, projection_neutron.Ncfg2_I, 1);
    workspacep_C3D_cfgp1_cfgp2_operatorp.resize(projection_proton.Ncfg1_I, projection_proton.Ncfg2_I, 1);
    build_amp(result_C8D_cfgn1_cfgn2_cfgp1_cfgp2_operator_gamma_alpha_beta, [&](HFBProjectionNucleus& nucleus_thread, auto& resultAtOmega_C5D_cfgn1_cfgn2_cfgp1_cfgp2_operator) {

        Eigen::Map<Eigen::MatrixXcd> result_C2D_cfgn1cfgn2_cfgp1cfgp2(resultAtOmega_C5D_cfgn1_cfgn2_cfgp1_cfgp2_operator.data(), Ncfgn1cfgn2_I, Ncfgp1cfgp2_I);
        nucleus_thread.projection_proton.integrate_pnp(Eigen::TensorMap<Eigen::Tensor<doubleC, 3, Eigen::ColMajor>>(nucleus_thread.workspacep_C3D_cfgp1_cfgp2_operatorp.data(), nucleus_thread.workspacep_C3D_cfgp1_cfgp2_operatorp.dimensions()), nucleus_thread.projection_proton.TargetN_I, [&](auto& resultAtPhi_C3D_cfg1_cfg2_operator) {
            Eigen::Map<Eigen::MatrixXcd> result_C2D_cfg1_cfg2(resultAtPhi_C3D_cfg1_cfg2_operator.data(), nucleus_thread.projection_proton.Ncfg1_I, nucleus_thread.projection_proton.Ncfg2_I);
            nucleus_thread.projection_proton.hfb_pfaffian.calc_overlap(result_C2D_cfg1_cfg2);
        });
        nucleus_thread.projection_neutron.integrate_pnp(Eigen::TensorMap<Eigen::Tensor<doubleC, 3, Eigen::ColMajor>>(nucleus_thread.workspacen_C3D_cfgn1_cfgn2_operatorn.data(), nucleus_thread.workspacen_C3D_cfgn1_cfgn2_operatorn.dimensions()), nucleus_thread.projection_neutron.TargetN_I, [&](auto& resultAtPhi_C3D_cfg1_cfg2_operator) {
            Eigen::Map<Eigen::MatrixXcd> result_C2D_cfg1_cfg2(resultAtPhi_C3D_cfg1_cfg2_operator.data(), nucleus_thread.projection_neutron.Ncfg1_I, nucleus_thread.projection_neutron.Ncfg2_I);
            nucleus_thread.projection_neutron.hfb_pfaffian.calc_one_body(result_C2D_cfg1_cfg2, OneBody_C2D_spn_spn);
        });

        // ColMajor: (n₁,n₂) → rows; (p₁,p₂) → columns.
        const Eigen::Map<const Eigen::VectorXcd> onebody_C1D_cfgn1cfgn2(nucleus_thread.workspacen_C3D_cfgn1_cfgn2_operatorn.data(), Ncfgn1cfgn2_I);
        const Eigen::Map<const Eigen::VectorXcd> overlap_C1D_cfgp1cfgp2(nucleus_thread.workspacep_C3D_cfgp1_cfgp2_operatorp.data(), Ncfgp1cfgp2_I);
        result_C2D_cfgn1cfgn2_cfgp1cfgp2.noalias() = onebody_C1D_cfgn1cfgn2 * overlap_C1D_cfgp1cfgp2.transpose();
    });
    return result_C8D_cfgn1_cfgn2_cfgp1_cfgp2_operator_gamma_alpha_beta.chip(0, 4);
}

inline Eigen::Tensor<doubleC, 7, Eigen::ColMajor> HFBProjectionNucleus::build_one_body_pp(const Eigen::MatrixXcd& OneBody_C2D_spp_spp) {
    assert(OneBody_C2D_spp_spp.rows() == projection_proton.Nsp_I && OneBody_C2D_spp_spp.cols() == projection_proton.Nsp_I && OneBody_C2D_spp_spp.allFinite());

    Eigen::Tensor<doubleC, 8, Eigen::ColMajor> result_C8D_cfgn1_cfgn2_cfgp1_cfgp2_operator_gamma_alpha_beta(projection_neutron.Ncfg1_I, projection_neutron.Ncfg2_I, projection_proton.Ncfg1_I, projection_proton.Ncfg2_I, 1, projection_neutron.Ngamma_I, projection_neutron.Nalpha_I, projection_neutron.Nbeta_I);
    // Size PNP workspaces before thread copies.
    workspacen_C3D_cfgn1_cfgn2_operatorn.resize(projection_neutron.Ncfg1_I, projection_neutron.Ncfg2_I, 1);
    workspacep_C3D_cfgp1_cfgp2_operatorp.resize(projection_proton.Ncfg1_I, projection_proton.Ncfg2_I, 1);
    build_amp(result_C8D_cfgn1_cfgn2_cfgp1_cfgp2_operator_gamma_alpha_beta, [&](HFBProjectionNucleus& nucleus_thread, auto& resultAtOmega_C5D_cfgn1_cfgn2_cfgp1_cfgp2_operator) {

        Eigen::Map<Eigen::MatrixXcd> result_C2D_cfgn1cfgn2_cfgp1cfgp2(resultAtOmega_C5D_cfgn1_cfgn2_cfgp1_cfgp2_operator.data(), Ncfgn1cfgn2_I, Ncfgp1cfgp2_I);
        nucleus_thread.projection_neutron.integrate_pnp(Eigen::TensorMap<Eigen::Tensor<doubleC, 3, Eigen::ColMajor>>(nucleus_thread.workspacen_C3D_cfgn1_cfgn2_operatorn.data(), nucleus_thread.workspacen_C3D_cfgn1_cfgn2_operatorn.dimensions()), nucleus_thread.projection_neutron.TargetN_I, [&](auto& resultAtPhi_C3D_cfg1_cfg2_operator) {
            Eigen::Map<Eigen::MatrixXcd> result_C2D_cfg1_cfg2(resultAtPhi_C3D_cfg1_cfg2_operator.data(), nucleus_thread.projection_neutron.Ncfg1_I, nucleus_thread.projection_neutron.Ncfg2_I);
            nucleus_thread.projection_neutron.hfb_pfaffian.calc_overlap(result_C2D_cfg1_cfg2);
        });
        nucleus_thread.projection_proton.integrate_pnp(Eigen::TensorMap<Eigen::Tensor<doubleC, 3, Eigen::ColMajor>>(nucleus_thread.workspacep_C3D_cfgp1_cfgp2_operatorp.data(), nucleus_thread.workspacep_C3D_cfgp1_cfgp2_operatorp.dimensions()), nucleus_thread.projection_proton.TargetN_I, [&](auto& resultAtPhi_C3D_cfg1_cfg2_operator) {
            Eigen::Map<Eigen::MatrixXcd> result_C2D_cfg1_cfg2(resultAtPhi_C3D_cfg1_cfg2_operator.data(), nucleus_thread.projection_proton.Ncfg1_I, nucleus_thread.projection_proton.Ncfg2_I);
            nucleus_thread.projection_proton.hfb_pfaffian.calc_one_body(result_C2D_cfg1_cfg2, OneBody_C2D_spp_spp);
        });

        // ColMajor: (n₁,n₂) → rows; (p₁,p₂) → columns.
        const Eigen::Map<const Eigen::VectorXcd> overlap_C1D_cfgn1cfgn2(nucleus_thread.workspacen_C3D_cfgn1_cfgn2_operatorn.data(), Ncfgn1cfgn2_I);
        const Eigen::Map<const Eigen::VectorXcd> onebody_C1D_cfgp1cfgp2(nucleus_thread.workspacep_C3D_cfgp1_cfgp2_operatorp.data(), Ncfgp1cfgp2_I);
        result_C2D_cfgn1cfgn2_cfgp1cfgp2.noalias() = overlap_C1D_cfgn1cfgn2 * onebody_C1D_cfgp1cfgp2.transpose();
    });
    return result_C8D_cfgn1_cfgn2_cfgp1_cfgp2_operator_gamma_alpha_beta.chip(0, 4);
}

inline Eigen::Tensor<doubleC, 7, Eigen::ColMajor> HFBProjectionNucleus::build_one_body_pn(const Eigen::MatrixXcd& OneBody_C2D_spp_spn) {
    assert(OneBody_C2D_spp_spn.rows() == projection_proton.Nsp_I && OneBody_C2D_spp_spn.cols() == projection_neutron.Nsp_I && OneBody_C2D_spp_spn.allFinite());

    // cπ†cν: (−1)^r₁; r₁ counts left neutron quasiparticles.
    Eigen::VectorXd sign_F1D_cfgn1cfgn2(Ncfgn1cfgn2_I);
    for (int cfgn2_I = 0; cfgn2_I < projection_neutron.Ncfg2_I; ++cfgn2_I) {
        for (int cfgn1_I = 0; cfgn1_I < projection_neutron.Ncfg1_I; ++cfgn1_I) {
            sign_F1D_cfgn1cfgn2(cfgn1_I + projection_neutron.Ncfg1_I * cfgn2_I) = 1.0 - 2.0 * (projection_neutron.hfb_pfaffian.config1_I2D_cfg1_cqp1[cfgn1_I].size() % 2);
        }
    }

    Eigen::Tensor<doubleC, 8, Eigen::ColMajor> result_C8D_cfgn1_cfgn2_cfgp1_cfgp2_operator_gamma_alpha_beta(projection_neutron.Ncfg1_I, projection_neutron.Ncfg2_I, projection_proton.Ncfg1_I, projection_proton.Ncfg2_I, 1, projection_neutron.Ngamma_I, projection_neutron.Nalpha_I, projection_neutron.Nbeta_I);
    // Size PNP workspaces before thread copies.
    workspacen_C3D_cfgn1_cfgn2_operatorn.resize(projection_neutron.Ncfg1_I, projection_neutron.Ncfg2_I, projection_neutron.Nsp_I);
    workspacep_C3D_cfgp1_cfgp2_operatorp.resize(projection_proton.Ncfg1_I, projection_proton.Ncfg2_I, projection_proton.Nsp_I);
    build_amp(result_C8D_cfgn1_cfgn2_cfgp1_cfgp2_operator_gamma_alpha_beta, [&](HFBProjectionNucleus& nucleus_thread, auto& resultAtOmega_C5D_cfgn1_cfgn2_cfgp1_cfgp2_operator) {
        // operator = sp: Aν = ⟨cν⟩; Cπ = ⟨cπ†⟩.

        Eigen::Map<Eigen::MatrixXcd> result_C2D_cfgn1cfgn2_cfgp1cfgp2(resultAtOmega_C5D_cfgn1_cfgn2_cfgp1_cfgp2_operator.data(), Ncfgn1cfgn2_I, Ncfgp1cfgp2_I);
        nucleus_thread.projection_neutron.integrate_pnp(Eigen::TensorMap<Eigen::Tensor<doubleC, 3, Eigen::ColMajor>>(nucleus_thread.workspacen_C3D_cfgn1_cfgn2_operatorn.data(), nucleus_thread.workspacen_C3D_cfgn1_cfgn2_operatorn.dimensions()), nucleus_thread.projection_neutron.TargetN_I, [&](auto& resultAtPhi_C3D_cfg1_cfg2_operator) {
            nucleus_thread.projection_neutron.hfb_pfaffian.calc_one_fermion(resultAtPhi_C3D_cfg1_cfg2_operator, false);
        });
        nucleus_thread.projection_proton.integrate_pnp(Eigen::TensorMap<Eigen::Tensor<doubleC, 3, Eigen::ColMajor>>(nucleus_thread.workspacep_C3D_cfgp1_cfgp2_operatorp.data(), nucleus_thread.workspacep_C3D_cfgp1_cfgp2_operatorp.dimensions()), nucleus_thread.projection_proton.TargetN_I, [&](auto& resultAtPhi_C3D_cfg1_cfg2_operator) {
            nucleus_thread.projection_proton.hfb_pfaffian.calc_one_fermion(resultAtPhi_C3D_cfg1_cfg2_operator, true);
        });

        // Hpn = diag((−1)^r₁) Aν Oᵀ Cπᵀ.
        const Eigen::Map<const Eigen::MatrixXcd> annihilator_C2D_cfgn1cfgn2_spn(nucleus_thread.workspacen_C3D_cfgn1_cfgn2_operatorn.data(), Ncfgn1cfgn2_I, nucleus_thread.projection_neutron.Nsp_I);
        const Eigen::Map<const Eigen::MatrixXcd> creator_C2D_cfgp1cfgp2_spp(nucleus_thread.workspacep_C3D_cfgp1_cfgp2_operatorp.data(), Ncfgp1cfgp2_I, nucleus_thread.projection_proton.Nsp_I);
        result_C2D_cfgn1cfgn2_cfgp1cfgp2.noalias() = sign_F1D_cfgn1cfgn2.asDiagonal() * annihilator_C2D_cfgn1cfgn2_spn * OneBody_C2D_spp_spn.transpose() * creator_C2D_cfgp1cfgp2_spp.transpose();
    });
    return result_C8D_cfgn1_cfgn2_cfgp1_cfgp2_operator_gamma_alpha_beta.chip(0, 4);
}

inline Eigen::Tensor<doubleC, 7, Eigen::ColMajor> HFBProjectionNucleus::build_two_body_nn(const Eigen::VectorXd& gQ_F1D_i, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& Q_C3D_spn_spn_i, const Eigen::VectorXd& gP_F1D_i, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& P_C3D_spn_spn_i) {
    // Hν²(Ω) Nπ(Ω); weighted Q/P kernel; pairing ¼ included.

    Eigen::Tensor<doubleC, 8, Eigen::ColMajor> result_C8D_cfgn1_cfgn2_cfgp1_cfgp2_operator_gamma_alpha_beta(projection_neutron.Ncfg1_I, projection_neutron.Ncfg2_I, projection_proton.Ncfg1_I, projection_proton.Ncfg2_I, 1, projection_neutron.Ngamma_I, projection_neutron.Nalpha_I, projection_neutron.Nbeta_I);
    // Size PNP workspaces before thread copies.
    workspacen_C3D_cfgn1_cfgn2_operatorn.resize(projection_neutron.Ncfg1_I, projection_neutron.Ncfg2_I, 1);
    workspacep_C3D_cfgp1_cfgp2_operatorp.resize(projection_proton.Ncfg1_I, projection_proton.Ncfg2_I, 1);
    build_amp(result_C8D_cfgn1_cfgn2_cfgp1_cfgp2_operator_gamma_alpha_beta, [&](HFBProjectionNucleus& nucleus_thread, auto& resultAtOmega_C5D_cfgn1_cfgn2_cfgp1_cfgp2_operator) {

        Eigen::Map<Eigen::MatrixXcd> result_C2D_cfgn1cfgn2_cfgp1cfgp2(resultAtOmega_C5D_cfgn1_cfgn2_cfgp1_cfgp2_operator.data(), Ncfgn1cfgn2_I, Ncfgp1cfgp2_I);
        nucleus_thread.projection_proton.integrate_pnp(Eigen::TensorMap<Eigen::Tensor<doubleC, 3, Eigen::ColMajor>>(nucleus_thread.workspacep_C3D_cfgp1_cfgp2_operatorp.data(), nucleus_thread.workspacep_C3D_cfgp1_cfgp2_operatorp.dimensions()), nucleus_thread.projection_proton.TargetN_I, [&](auto& resultAtPhi_C3D_cfg1_cfg2_operator) {
            Eigen::Map<Eigen::MatrixXcd> result_C2D_cfg1_cfg2(resultAtPhi_C3D_cfg1_cfg2_operator.data(), nucleus_thread.projection_proton.Ncfg1_I, nucleus_thread.projection_proton.Ncfg2_I);
            nucleus_thread.projection_proton.hfb_pfaffian.calc_overlap(result_C2D_cfg1_cfg2);
        });
        nucleus_thread.projection_neutron.integrate_pnp(Eigen::TensorMap<Eigen::Tensor<doubleC, 3, Eigen::ColMajor>>(nucleus_thread.workspacen_C3D_cfgn1_cfgn2_operatorn.data(), nucleus_thread.workspacen_C3D_cfgn1_cfgn2_operatorn.dimensions()), nucleus_thread.projection_neutron.TargetN_I, [&](auto& resultAtPhi_C3D_cfg1_cfg2_operator) {
            Eigen::Map<Eigen::MatrixXcd> result_C2D_cfg1_cfg2(resultAtPhi_C3D_cfg1_cfg2_operator.data(), nucleus_thread.projection_neutron.Ncfg1_I, nucleus_thread.projection_neutron.Ncfg2_I);
            nucleus_thread.projection_neutron.hfb_pfaffian.calc_two_body(result_C2D_cfg1_cfg2, gQ_F1D_i, gP_F1D_i, Q_C3D_spn_spn_i, P_C3D_spn_spn_i);
        });

        // ColMajor: (n₁,n₂) → rows; (p₁,p₂) → columns.
        const Eigen::Map<const Eigen::VectorXcd> twobody_C1D_cfgn1cfgn2(nucleus_thread.workspacen_C3D_cfgn1_cfgn2_operatorn.data(), Ncfgn1cfgn2_I);
        const Eigen::Map<const Eigen::VectorXcd> overlap_C1D_cfgp1cfgp2(nucleus_thread.workspacep_C3D_cfgp1_cfgp2_operatorp.data(), Ncfgp1cfgp2_I);
        result_C2D_cfgn1cfgn2_cfgp1cfgp2.noalias() = twobody_C1D_cfgn1cfgn2 * overlap_C1D_cfgp1cfgp2.transpose();
    });
    return result_C8D_cfgn1_cfgn2_cfgp1_cfgp2_operator_gamma_alpha_beta.chip(0, 4);
}

inline Eigen::Tensor<doubleC, 7, Eigen::ColMajor> HFBProjectionNucleus::build_two_body_pp(const Eigen::VectorXd& gQ_F1D_i, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& Q_C3D_spp_spp_i, const Eigen::VectorXd& gP_F1D_i, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& P_C3D_spp_spp_i) {
    // Nν(Ω) Hπ²(Ω); weighted Q/P kernel; pairing ¼ included.

    Eigen::Tensor<doubleC, 8, Eigen::ColMajor> result_C8D_cfgn1_cfgn2_cfgp1_cfgp2_operator_gamma_alpha_beta(projection_neutron.Ncfg1_I, projection_neutron.Ncfg2_I, projection_proton.Ncfg1_I, projection_proton.Ncfg2_I, 1, projection_neutron.Ngamma_I, projection_neutron.Nalpha_I, projection_neutron.Nbeta_I);
    // Size PNP workspaces before thread copies.
    workspacen_C3D_cfgn1_cfgn2_operatorn.resize(projection_neutron.Ncfg1_I, projection_neutron.Ncfg2_I, 1);
    workspacep_C3D_cfgp1_cfgp2_operatorp.resize(projection_proton.Ncfg1_I, projection_proton.Ncfg2_I, 1);
    build_amp(result_C8D_cfgn1_cfgn2_cfgp1_cfgp2_operator_gamma_alpha_beta, [&](HFBProjectionNucleus& nucleus_thread, auto& resultAtOmega_C5D_cfgn1_cfgn2_cfgp1_cfgp2_operator) {

        Eigen::Map<Eigen::MatrixXcd> result_C2D_cfgn1cfgn2_cfgp1cfgp2(resultAtOmega_C5D_cfgn1_cfgn2_cfgp1_cfgp2_operator.data(), Ncfgn1cfgn2_I, Ncfgp1cfgp2_I);
        nucleus_thread.projection_neutron.integrate_pnp(Eigen::TensorMap<Eigen::Tensor<doubleC, 3, Eigen::ColMajor>>(nucleus_thread.workspacen_C3D_cfgn1_cfgn2_operatorn.data(), nucleus_thread.workspacen_C3D_cfgn1_cfgn2_operatorn.dimensions()), nucleus_thread.projection_neutron.TargetN_I, [&](auto& resultAtPhi_C3D_cfg1_cfg2_operator) {
            Eigen::Map<Eigen::MatrixXcd> result_C2D_cfg1_cfg2(resultAtPhi_C3D_cfg1_cfg2_operator.data(), nucleus_thread.projection_neutron.Ncfg1_I, nucleus_thread.projection_neutron.Ncfg2_I);
            nucleus_thread.projection_neutron.hfb_pfaffian.calc_overlap(result_C2D_cfg1_cfg2);
        });
        nucleus_thread.projection_proton.integrate_pnp(Eigen::TensorMap<Eigen::Tensor<doubleC, 3, Eigen::ColMajor>>(nucleus_thread.workspacep_C3D_cfgp1_cfgp2_operatorp.data(), nucleus_thread.workspacep_C3D_cfgp1_cfgp2_operatorp.dimensions()), nucleus_thread.projection_proton.TargetN_I, [&](auto& resultAtPhi_C3D_cfg1_cfg2_operator) {
            Eigen::Map<Eigen::MatrixXcd> result_C2D_cfg1_cfg2(resultAtPhi_C3D_cfg1_cfg2_operator.data(), nucleus_thread.projection_proton.Ncfg1_I, nucleus_thread.projection_proton.Ncfg2_I);
            nucleus_thread.projection_proton.hfb_pfaffian.calc_two_body(result_C2D_cfg1_cfg2, gQ_F1D_i, gP_F1D_i, Q_C3D_spp_spp_i, P_C3D_spp_spp_i);
        });

        // ColMajor: (n₁,n₂) → rows; (p₁,p₂) → columns.
        const Eigen::Map<const Eigen::VectorXcd> overlap_C1D_cfgn1cfgn2(nucleus_thread.workspacen_C3D_cfgn1_cfgn2_operatorn.data(), Ncfgn1cfgn2_I);
        const Eigen::Map<const Eigen::VectorXcd> twobody_C1D_cfgp1cfgp2(nucleus_thread.workspacep_C3D_cfgp1_cfgp2_operatorp.data(), Ncfgp1cfgp2_I);
        result_C2D_cfgn1cfgn2_cfgp1cfgp2.noalias() = overlap_C1D_cfgn1cfgn2 * twobody_C1D_cfgp1cfgp2.transpose();
    });
    return result_C8D_cfgn1_cfgn2_cfgp1_cfgp2_operator_gamma_alpha_beta.chip(0, 4);
}

inline Eigen::Tensor<doubleC, 7, Eigen::ColMajor> HFBProjectionNucleus::build_two_body_pn(const Eigen::VectorXd& gQ_F1D_i, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& Q_C3D_spn_spn_i, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& Q_C3D_spp_spp_i) {

    // NQν = NQπ = NgQ.
    assert(gQ_F1D_i.allFinite());
    assert(Q_C3D_spn_spn_i.dimension(0) == projection_neutron.Nsp_I && Q_C3D_spn_spn_i.dimension(1) == projection_neutron.Nsp_I && Q_C3D_spn_spn_i.dimension(2) == gQ_F1D_i.size());
    assert(Q_C3D_spp_spp_i.dimension(0) == projection_proton.Nsp_I && Q_C3D_spp_spp_i.dimension(1) == projection_proton.Nsp_I && Q_C3D_spp_spp_i.dimension(2) == gQ_F1D_i.size());
    assert(Eigen::Map<const Eigen::VectorXcd>(Q_C3D_spn_spn_i.data(), Q_C3D_spn_spn_i.size()).allFinite());
    assert(Eigen::Map<const Eigen::VectorXcd>(Q_C3D_spp_spp_i.data(), Q_C3D_spp_spp_i.size()).allFinite());

    const Eigen::Index Noperator_I = 2 * gQ_F1D_i.size();

    // Q†(α,β,i) = Q*(β,α,i); shared read-only across Ω.
    const Eigen::Tensor<doubleC, 3, Eigen::ColMajor> Qdag_C3D_spn_spn_i = Q_C3D_spn_spn_i.conjugate().shuffle(Eigen::array<int, 3>{1, 0, 2});
    const Eigen::Tensor<doubleC, 3, Eigen::ColMajor> Qdag_C3D_spp_spp_i = Q_C3D_spp_spp_i.conjugate().shuffle(Eigen::array<int, 3>{1, 0, 2});

    Eigen::Tensor<doubleC, 8, Eigen::ColMajor> result_C8D_cfgn1_cfgn2_cfgp1_cfgp2_operator_gamma_alpha_beta(projection_neutron.Ncfg1_I, projection_neutron.Ncfg2_I, projection_proton.Ncfg1_I, projection_proton.Ncfg2_I, 1, projection_neutron.Ngamma_I, projection_neutron.Nalpha_I, projection_neutron.Nbeta_I);
    // Size PNP workspaces before thread copies.
    workspacen_C3D_cfgn1_cfgn2_operatorn.resize(projection_neutron.Ncfg1_I, projection_neutron.Ncfg2_I, Noperator_I);
    workspacep_C3D_cfgp1_cfgp2_operatorp.resize(projection_proton.Ncfg1_I, projection_proton.Ncfg2_I, Noperator_I);
    build_amp(result_C8D_cfgn1_cfgn2_cfgp1_cfgp2_operator_gamma_alpha_beta, [&](HFBProjectionNucleus& nucleus_thread, auto& resultAtOmega_C5D_cfgn1_cfgn2_cfgp1_cfgp2_operator) {
        // operator (2i,2i+1) = (Qᵢ,Qᵢ†).

        // All neutron operators share contractions at (Ω,φ).
        nucleus_thread.projection_neutron.integrate_pnp(Eigen::TensorMap<Eigen::Tensor<doubleC, 3, Eigen::ColMajor>>(nucleus_thread.workspacen_C3D_cfgn1_cfgn2_operatorn.data(), nucleus_thread.workspacen_C3D_cfgn1_cfgn2_operatorn.dimensions()), nucleus_thread.projection_neutron.TargetN_I, [&](auto& resultAtPhi_C3D_cfg1_cfg2_operator) {
            resultAtPhi_C3D_cfg1_cfg2_operator.setZero();
            for (Eigen::Index iQ_I = 0; iQ_I < gQ_F1D_i.size(); ++iQ_I) {
                if (gQ_F1D_i(iQ_I) == 0.0) { continue; }

                const Eigen::Index offset_I = iQ_I * nucleus_thread.projection_neutron.Nsp_I * nucleus_thread.projection_neutron.Nsp_I;
                const Eigen::Map<const Eigen::MatrixXcd> Q_C2D_spn_spn(Q_C3D_spn_spn_i.data() + offset_I, nucleus_thread.projection_neutron.Nsp_I, nucleus_thread.projection_neutron.Nsp_I);
                const Eigen::Map<const Eigen::MatrixXcd> Qdag_C2D_spn_spn(Qdag_C3D_spn_spn_i.data() + offset_I, nucleus_thread.projection_neutron.Nsp_I, nucleus_thread.projection_neutron.Nsp_I);

                Eigen::Map<Eigen::MatrixXcd> Q_C2D_cfgn1_cfgn2(resultAtPhi_C3D_cfg1_cfg2_operator.data() + 2 * iQ_I * Ncfgn1cfgn2_I, nucleus_thread.projection_neutron.Ncfg1_I, nucleus_thread.projection_neutron.Ncfg2_I);
                Eigen::Map<Eigen::MatrixXcd> Qdag_C2D_cfgn1_cfgn2(resultAtPhi_C3D_cfg1_cfg2_operator.data() + (2 * iQ_I + 1) * Ncfgn1cfgn2_I, nucleus_thread.projection_neutron.Ncfg1_I, nucleus_thread.projection_neutron.Ncfg2_I);
                nucleus_thread.projection_neutron.hfb_pfaffian.calc_one_body(Q_C2D_cfgn1_cfgn2, Q_C2D_spn_spn);
                nucleus_thread.projection_neutron.hfb_pfaffian.calc_one_body(Qdag_C2D_cfgn1_cfgn2, Qdag_C2D_spn_spn);
            }
        });

        // All proton operators share contractions at (Ω,φ).
        nucleus_thread.projection_proton.integrate_pnp(Eigen::TensorMap<Eigen::Tensor<doubleC, 3, Eigen::ColMajor>>(nucleus_thread.workspacep_C3D_cfgp1_cfgp2_operatorp.data(), nucleus_thread.workspacep_C3D_cfgp1_cfgp2_operatorp.dimensions()), nucleus_thread.projection_proton.TargetN_I, [&](auto& resultAtPhi_C3D_cfg1_cfg2_operator) {
            resultAtPhi_C3D_cfg1_cfg2_operator.setZero();
            for (Eigen::Index iQ_I = 0; iQ_I < gQ_F1D_i.size(); ++iQ_I) {
                if (gQ_F1D_i(iQ_I) == 0.0) { continue; }

                const Eigen::Index offset_I = iQ_I * nucleus_thread.projection_proton.Nsp_I * nucleus_thread.projection_proton.Nsp_I;
                const Eigen::Map<const Eigen::MatrixXcd> Q_C2D_spp_spp(Q_C3D_spp_spp_i.data() + offset_I, nucleus_thread.projection_proton.Nsp_I, nucleus_thread.projection_proton.Nsp_I);
                const Eigen::Map<const Eigen::MatrixXcd> Qdag_C2D_spp_spp(Qdag_C3D_spp_spp_i.data() + offset_I, nucleus_thread.projection_proton.Nsp_I, nucleus_thread.projection_proton.Nsp_I);

                Eigen::Map<Eigen::MatrixXcd> Q_C2D_cfgp1_cfgp2(resultAtPhi_C3D_cfg1_cfg2_operator.data() + 2 * iQ_I * Ncfgp1cfgp2_I, nucleus_thread.projection_proton.Ncfg1_I, nucleus_thread.projection_proton.Ncfg2_I);
                Eigen::Map<Eigen::MatrixXcd> Qdag_C2D_cfgp1_cfgp2(resultAtPhi_C3D_cfg1_cfg2_operator.data() + (2 * iQ_I + 1) * Ncfgp1cfgp2_I, nucleus_thread.projection_proton.Ncfg1_I, nucleus_thread.projection_proton.Ncfg2_I);
                nucleus_thread.projection_proton.hfb_pfaffian.calc_one_body(Q_C2D_cfgp1_cfgp2, Q_C2D_spp_spp);
                nucleus_thread.projection_proton.hfb_pfaffian.calc_one_body(Qdag_C2D_cfgp1_cfgp2, Qdag_C2D_spp_spp);
            }
        });

        // Hnp(Ω) = ½Σᵢ gᵢ [N(Qᵢ†)P(Qᵢ)ᵀ + N(Qᵢ)P(Qᵢ†)ᵀ].
        const Eigen::Map<const Eigen::MatrixXcd> onebody_C2D_cfgn1cfgn2_operator(nucleus_thread.workspacen_C3D_cfgn1_cfgn2_operatorn.data(), Ncfgn1cfgn2_I, Noperator_I);
        const Eigen::Map<const Eigen::MatrixXcd> onebody_C2D_cfgp1cfgp2_operator(nucleus_thread.workspacep_C3D_cfgp1_cfgp2_operatorp.data(), Ncfgp1cfgp2_I, Noperator_I);
        Eigen::Map<Eigen::MatrixXcd> result_C2D_cfgn1cfgn2_cfgp1cfgp2(resultAtOmega_C5D_cfgn1_cfgn2_cfgp1_cfgp2_operator.data(), Ncfgn1cfgn2_I, Ncfgp1cfgp2_I);
        result_C2D_cfgn1cfgn2_cfgp1cfgp2.setZero();

        for (Eigen::Index iQ_I = 0; iQ_I < gQ_F1D_i.size(); ++iQ_I) {
            result_C2D_cfgn1cfgn2_cfgp1cfgp2.noalias() += 0.5 * gQ_F1D_i(iQ_I) * onebody_C2D_cfgn1cfgn2_operator.col(2 * iQ_I + 1) * onebody_C2D_cfgp1cfgp2_operator.col(2 * iQ_I).transpose();
            result_C2D_cfgn1cfgn2_cfgp1cfgp2.noalias() += 0.5 * gQ_F1D_i(iQ_I) * onebody_C2D_cfgn1cfgn2_operator.col(2 * iQ_I) * onebody_C2D_cfgp1cfgp2_operator.col(2 * iQ_I + 1).transpose();
        }
    });
    return result_C8D_cfgn1_cfgn2_cfgp1_cfgp2_operator_gamma_alpha_beta.chip(0, 4);
}
