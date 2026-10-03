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
 * @note   Projected tensors are allocated locally and returned by value.
 * @note   Particle-number targets apply to the right state.
 * @note   Requires even reference vacua and initialized projection data.
 * @note   Vacuum overlaps must remain nonzero at quadrature nodes.
 */
class HFBProjectionNucleus {
public:
    HFBProjection projection_neutron;
    HFBProjection projection_proton;


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

private:
    /**
     * @brief Build nuclear kernels using gauge quadrature.
     * @math cache_ab(Ω) = ⟨Φ₁;a|O Pνᴺ Pπᶻ R(Ω)|Φ₂;b⟩.
     * @output Tensor (cfgn1,cfgn2,cfgp1,cfgp2,γ,α,β), returned by value.
     * @note Excludes K selection, Euler weights, and Wigner factors.
     */
    Eigen::Tensor<doubleC, 7, Eigen::ColMajor> build_amp(const std::function<void(Eigen::TensorMap<Eigen::Tensor<doubleC, 4, Eigen::ColMajor>>&)>& kernel_Func);

public:
    /**
     * @brief Integrate cached nuclear kernels using Euler quadrature.
     * @math result = (2I+1)/VΩ ΣΩ wΩ Dᴵ*(Ω) cache(Ω).
     * @output Projected tensor for one I, returned by value.
     * @note Cache includes both species' particle-number projections.
     */
    Eigen::Tensor<doubleC, 6, Eigen::ColMajor> integrate_amp(int TargetTwoI_I, const Eigen::Tensor<doubleC, 7, Eigen::ColMajor>& cache_C7D_cfgn1_cfgn2_cfgp1_cfgp2_gamma_alpha_beta);

    /**
     * @brief Integrate axial kernels using beta quadrature.
     * @math M = (2I+1)/2 Σβ wβ dᴵ*_{K₁K₂}(β) cache(β).
     * @output Tensor (cfgn1,cfgn2,cfgp1,cfgp2).
     * @note Definite K; cache at α=γ=0; wβ includes sinβ.
     */
    Eigen::Tensor<doubleC, 4, Eigen::ColMajor> integrate_amp_with_K(int TargetTwoI_I, const Eigen::VectorXi& TwoKn1_I1D_cfgn1, const Eigen::VectorXi& TwoKn2_I1D_cfgn2, const Eigen::VectorXi& TwoKp1_I1D_cfgp1, const Eigen::VectorXi& TwoKp2_I1D_cfgp2, const Eigen::Tensor<doubleC, 5, Eigen::ColMajor>& cache_C5D_cfgn1_cfgn2_cfgp1_cfgp2_beta);

    /**
     * @brief Build nuclear overlaps using gauge quadrature.
     * @math result_n₁n₂p₁p₂(Ω) = Nν_n₁n₂(Ω) Nπ_p₁p₂(Ω).
     * @output result_C7D_cfgn1_cfgn2_cfgp1_cfgp2_gamma_alpha_beta, returned by value.
     */
    Eigen::Tensor<doubleC, 7, Eigen::ColMajor> build_overlap();

    /**
     * @brief Build nn kernels using Pfaffians and gauge quadrature.
     * @math result(Ω) = Nπ(Ω) Hν¹(Ω).
     * @output result_C7D_cfgn1_cfgn2_cfgp1_cfgp2_gamma_alpha_beta, returned by value.
     */
    Eigen::Tensor<doubleC, 7, Eigen::ColMajor> build_one_body_nn(const Eigen::MatrixXcd& OneBody_C2D_spn_spn);

    /**
     * @brief Build pp kernels using Pfaffians and gauge quadrature.
     * @math result(Ω) = Hπ¹(Ω) Nν(Ω).
     * @output result_C7D_cfgn1_cfgn2_cfgp1_cfgp2_gamma_alpha_beta, returned by value.
     */
    Eigen::Tensor<doubleC, 7, Eigen::ColMajor> build_one_body_pp(const Eigen::MatrixXcd& OneBody_C2D_spp_spp);

    /**
     * @brief Build pn kernels using Pfaffians and gauge quadrature.
     * @math result(Ω) = [(-1)^r₁ Aν(Ω) Oᵀ] Cπ(Ω)ᵀ.
     * @output result_C7D_cfgn1_cfgn2_cfgp1_cfgp2_gamma_alpha_beta, returned by value.
     */
    Eigen::Tensor<doubleC, 7, Eigen::ColMajor> build_one_body_pn(const Eigen::MatrixXcd& OneBody_C2D_spp_spn);

    /**
     * @brief Build nn kernels using Pfaffians and gauge quadrature.
     * @math result(Ω) = Nπ(Ω) Hν²(Ω).
     * @output result_C7D_cfgn1_cfgn2_cfgp1_cfgp2_gamma_alpha_beta, returned by value.
     */
    Eigen::Tensor<doubleC, 7, Eigen::ColMajor> build_two_body_nn(const Eigen::VectorXd& gQ_F1D_i, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& Q_C3D_spn_spn_i, const Eigen::VectorXd& gP_F1D_i, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& P_C3D_spn_spn_i);

    /**
     * @brief Build pp kernels using Pfaffians and gauge quadrature.
     * @math result(Ω) = Hπ²(Ω) Nν(Ω).
     * @output result_C7D_cfgn1_cfgn2_cfgp1_cfgp2_gamma_alpha_beta, returned by value.
     */
    Eigen::Tensor<doubleC, 7, Eigen::ColMajor> build_two_body_pp(const Eigen::VectorXd& gQ_F1D_i, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& Q_C3D_spp_spp_i, const Eigen::VectorXd& gP_F1D_i, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& P_C3D_spp_spp_i);

    /**
     * @brief Build pn kernels using Pfaffians and gauge quadrature.
     * @math Hpn = ½ Σi gQi (Qνi† Qπi + Qνi Qπi†).
     * @output result_C7D_cfgn1_cfgn2_cfgp1_cfgp2_gamma_alpha_beta, returned by value.
     */
    Eigen::Tensor<doubleC, 7, Eigen::ColMajor> build_two_body_pn(const Eigen::VectorXd& gQ_F1D_i, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& Q_C3D_spn_spn_i, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& Q_C3D_spp_spp_i);

};

inline Eigen::Tensor<doubleC, 7, Eigen::ColMajor> HFBProjectionNucleus::build_amp(const std::function<void(Eigen::TensorMap<Eigen::Tensor<doubleC, 4, Eigen::ColMajor>>&)>& kernel_Func) {
    auto& neutron = projection_neutron;
    auto& proton = projection_proton;

    // Ων = Ωπ.
    assert(neutron.Nalpha_I > 0 && neutron.Nbeta_I > 0 && neutron.Ngamma_I > 0);
    assert(neutron.Nalpha_I == proton.Nalpha_I && neutron.Nbeta_I == proton.Nbeta_I && neutron.Ngamma_I == proton.Ngamma_I);
    assert(neutron.alpha_F1D_alpha.isApprox(proton.alpha_F1D_alpha));
    assert(neutron.beta_F1D_beta.isApprox(proton.beta_F1D_beta));
    assert(neutron.gamma_F1D_gamma.isApprox(proton.gamma_F1D_gamma));

    // cache ∈ ℂ^{cfgn1×cfgn2×cfgp1×cfgp2×Nγ×Nα×Nβ}.
    Eigen::Tensor<doubleC, 7, Eigen::ColMajor> cache_C7D_cfgn1_cfgn2_cfgp1_cfgp2_gamma_alpha_beta = Eigen::Tensor<doubleC, 7, Eigen::ColMajor>(neutron.Ncfg1_I, neutron.Ncfg2_I, proton.Ncfg1_I, proton.Ncfg2_I, neutron.Ngamma_I, neutron.Nalpha_I, neutron.Nbeta_I).setZero();

    // Ω → Dν,Dπ → cache(Ω).
    for (int beta_I = 0; beta_I < neutron.Nbeta_I; ++beta_I) {
        for (int alpha_I = 0; alpha_I < neutron.Nalpha_I; ++alpha_I) {
            for (int gamma_I = 0; gamma_I < neutron.Ngamma_I; ++gamma_I) {
                {
                    // D = Rz(α) Ry(β) Rz(γ).
                    const Eigen::Map<const Eigen::MatrixXcd> RzAlpha_C2D_sp_sp(neutron.Rz_C3D_sp_sp_alpha.data() + static_cast<Eigen::Index>(alpha_I) * neutron.Nsp_I * neutron.Nsp_I, neutron.Nsp_I, neutron.Nsp_I);
                    const Eigen::Map<const Eigen::MatrixXcd> RyBeta_C2D_sp_sp(neutron.Ry_C3D_sp_sp_beta.data() + static_cast<Eigen::Index>(beta_I) * neutron.Nsp_I * neutron.Nsp_I, neutron.Nsp_I, neutron.Nsp_I);
                    const Eigen::Map<const Eigen::MatrixXcd> RzGamma_C2D_sp_sp(neutron.Rz_C3D_sp_sp_gamma.data() + static_cast<Eigen::Index>(gamma_I) * neutron.Nsp_I * neutron.Nsp_I, neutron.Nsp_I, neutron.Nsp_I);
                    neutron.D_C2D_sp_sp.noalias() = RzAlpha_C2D_sp_sp * RyBeta_C2D_sp_sp * RzGamma_C2D_sp_sp;
                }
                {
                    // D = Rz(α) Ry(β) Rz(γ).
                    const Eigen::Map<const Eigen::MatrixXcd> RzAlpha_C2D_sp_sp(proton.Rz_C3D_sp_sp_alpha.data() + static_cast<Eigen::Index>(alpha_I) * proton.Nsp_I * proton.Nsp_I, proton.Nsp_I, proton.Nsp_I);
                    const Eigen::Map<const Eigen::MatrixXcd> RyBeta_C2D_sp_sp(proton.Ry_C3D_sp_sp_beta.data() + static_cast<Eigen::Index>(beta_I) * proton.Nsp_I * proton.Nsp_I, proton.Nsp_I, proton.Nsp_I);
                    const Eigen::Map<const Eigen::MatrixXcd> RzGamma_C2D_sp_sp(proton.Rz_C3D_sp_sp_gamma.data() + static_cast<Eigen::Index>(gamma_I) * proton.Nsp_I * proton.Nsp_I, proton.Nsp_I, proton.Nsp_I);
                    proton.D_C2D_sp_sp.noalias() = RzAlpha_C2D_sp_sp * RyBeta_C2D_sp_sp * RzGamma_C2D_sp_sp;
                }
                // ColMajor: cfgn1 → cfgn2 → cfgp1 → cfgp2 → γ → α → β.
                const Eigen::Index offset_I = Ncfgn1cfgn2cfgp1cfgp2_I * (gamma_I + static_cast<Eigen::Index>(neutron.Ngamma_I) * alpha_I + static_cast<Eigen::Index>(neutron.Ngamma_I) * neutron.Nalpha_I * beta_I);
                Eigen::TensorMap<Eigen::Tensor<doubleC, 4, Eigen::ColMajor>> cache_C4D_cfgn1_cfgn2_cfgp1_cfgp2(cache_C7D_cfgn1_cfgn2_cfgp1_cfgp2_gamma_alpha_beta.data() + offset_I, neutron.Ncfg1_I, neutron.Ncfg2_I, proton.Ncfg1_I, proton.Ncfg2_I);
                kernel_Func(cache_C4D_cfgn1_cfgn2_cfgp1_cfgp2);
            }
        }
    }
    return cache_C7D_cfgn1_cfgn2_cfgp1_cfgp2_gamma_alpha_beta;
}

inline Eigen::Tensor<doubleC, 6, Eigen::ColMajor> HFBProjectionNucleus::integrate_amp(int TargetTwoI_I, const Eigen::Tensor<doubleC, 7, Eigen::ColMajor>& cache_C7D_cfgn1_cfgn2_cfgp1_cfgp2_gamma_alpha_beta) {
    auto& neutron = projection_neutron;
    auto& proton = projection_proton;

    // Ων = Ωπ; wν = wπ.
    assert(neutron.alpha_F1D_alpha.isApprox(proton.alpha_F1D_alpha));
    assert(neutron.beta_F1D_beta.isApprox(proton.beta_F1D_beta));
    assert(neutron.gamma_F1D_gamma.isApprox(proton.gamma_F1D_gamma));
    assert(neutron.weight_F1D_alpha.isApprox(proton.weight_F1D_alpha));
    assert(neutron.weight_F1D_beta.isApprox(proton.weight_F1D_beta));
    assert(neutron.weight_F1D_gamma.isApprox(proton.weight_F1D_gamma));


    // cache ∈ ℂ^{cfgn1×cfgn2×cfgp1×cfgp2×Nγ×Nα×Nβ}.
    assert(cache_C7D_cfgn1_cfgn2_cfgp1_cfgp2_gamma_alpha_beta.dimension(0) == neutron.Ncfg1_I);
    assert(cache_C7D_cfgn1_cfgn2_cfgp1_cfgp2_gamma_alpha_beta.dimension(1) == neutron.Ncfg2_I);
    assert(cache_C7D_cfgn1_cfgn2_cfgp1_cfgp2_gamma_alpha_beta.dimension(2) == proton.Ncfg1_I);
    assert(cache_C7D_cfgn1_cfgn2_cfgp1_cfgp2_gamma_alpha_beta.dimension(3) == proton.Ncfg2_I);
    assert(cache_C7D_cfgn1_cfgn2_cfgp1_cfgp2_gamma_alpha_beta.dimension(4) == neutron.Ngamma_I);
    assert(cache_C7D_cfgn1_cfgn2_cfgp1_cfgp2_gamma_alpha_beta.dimension(5) == neutron.Nalpha_I);
    assert(cache_C7D_cfgn1_cfgn2_cfgp1_cfgp2_gamma_alpha_beta.dimension(6) == neutron.Nbeta_I);

    // VΩ = Σwα Σwβ Σwγ; wβ includes sinβ.
    assert(TargetTwoI_I >= 0);
    const double volume_F = neutron.weight_F1D_alpha.sum() * neutron.weight_F1D_beta.sum() * neutron.weight_F1D_gamma.sum();
    const double normalization_F = (TargetTwoI_I + 1.0) / volume_F;
    assert(std::isfinite(normalization_F) && normalization_F > 0.0);
    
    Eigen::Tensor<doubleC, 6, Eigen::ColMajor> result_C6D_cfgn1_cfgn2_cfgp1_cfgp2_K_K = Eigen::Tensor<doubleC, 6, Eigen::ColMajor>(neutron.Ncfg1_I, neutron.Ncfg2_I, proton.Ncfg1_I, proton.Ncfg2_I, TargetTwoI_I + 1, TargetTwoI_I + 1).setZero();
    Eigen::MatrixXcd weightD_C2D_K_K(TargetTwoI_I + 1, TargetTwoI_I + 1);

    // cfg = (cfgn1,cfgn2,cfgp1,cfgp2); K = (K1,K2).

    // result = (2I+1)/VΩ ΣΩ wΩ cache(Ω) vec(Dᴵ*(Ω))ᵀ.
    RepresentationSpin wignerD_X(TargetTwoI_I);
    for (int beta_I = 0; beta_I < neutron.Nbeta_I; ++beta_I) {
        for (int alpha_I = 0; alpha_I < neutron.Nalpha_I; ++alpha_I) {
            for (int gamma_I = 0; gamma_I < neutron.Ngamma_I; ++gamma_I) {
                const double weight_F = normalization_F * neutron.weight_F1D_alpha(alpha_I) * neutron.weight_F1D_beta(beta_I) * neutron.weight_F1D_gamma(gamma_I);

                // W_K1K2 = w exp(iK1α) dᴵ*_{K1K2}(β) exp(iK2γ).
                weightD_C2D_K_K = weight_F * wignerD_X.calc_R(neutron.alpha_F1D_alpha(alpha_I), neutron.beta_F1D_beta(beta_I), neutron.gamma_F1D_gamma(gamma_I)).conjugate();
                const Eigen::Index offset_I = Ncfgn1cfgn2cfgp1cfgp2_I * (gamma_I + static_cast<Eigen::Index>(neutron.Ngamma_I) * alpha_I + static_cast<Eigen::Index>(neutron.Ngamma_I) * neutron.Nalpha_I * beta_I);
                Eigen::Map<Eigen::MatrixXcd> result_C2D_cfgn1cfgn2cfgp1cfgp2_KK(result_C6D_cfgn1_cfgn2_cfgp1_cfgp2_K_K.data(), Ncfgn1cfgn2cfgp1cfgp2_I, static_cast<Eigen::Index>(TargetTwoI_I + 1) * (TargetTwoI_I + 1));
                const Eigen::Map<const Eigen::VectorXcd> weightD_C1D_KK(weightD_C2D_K_K.data(), weightD_C2D_K_K.size());
                const Eigen::Map<const Eigen::VectorXcd> cache_C1D_cfgn1cfgn2cfgp1cfgp2(cache_C7D_cfgn1_cfgn2_cfgp1_cfgp2_gamma_alpha_beta.data() + offset_I, Ncfgn1cfgn2cfgp1cfgp2_I);
                result_C2D_cfgn1cfgn2cfgp1cfgp2_KK.noalias() += cache_C1D_cfgn1cfgn2cfgp1cfgp2 * weightD_C1D_KK.transpose();
            }
        }
    }
    return result_C6D_cfgn1_cfgn2_cfgp1_cfgp2_K_K;
}

inline Eigen::Tensor<doubleC, 4, Eigen::ColMajor> HFBProjectionNucleus::integrate_amp_with_K(int TargetTwoI_I, const Eigen::VectorXi& TwoKn1_I1D_cfgn1, const Eigen::VectorXi& TwoKn2_I1D_cfgn2, const Eigen::VectorXi& TwoKp1_I1D_cfgp1, const Eigen::VectorXi& TwoKp2_I1D_cfgp2, const Eigen::Tensor<doubleC, 5, Eigen::ColMajor>& cache_C5D_cfgn1_cfgn2_cfgp1_cfgp2_beta) {
    auto& neutron = projection_neutron;
    auto& proton = projection_proton;
    assert(TargetTwoI_I >= 0);
    assert(neutron.beta_F1D_beta.isApprox(proton.beta_F1D_beta));
    assert(neutron.weight_F1D_beta.isApprox(proton.weight_F1D_beta));
    assert(TwoKn1_I1D_cfgn1.size() == neutron.Ncfg1_I);
    assert(TwoKn2_I1D_cfgn2.size() == neutron.Ncfg2_I);
    assert(TwoKp1_I1D_cfgp1.size() == proton.Ncfg1_I);
    assert(TwoKp2_I1D_cfgp2.size() == proton.Ncfg2_I);
    assert(cache_C5D_cfgn1_cfgn2_cfgp1_cfgp2_beta.dimension(0) == neutron.Ncfg1_I);
    assert(cache_C5D_cfgn1_cfgn2_cfgp1_cfgp2_beta.dimension(1) == neutron.Ncfg2_I);
    assert(cache_C5D_cfgn1_cfgn2_cfgp1_cfgp2_beta.dimension(2) == proton.Ncfg1_I);
    assert(cache_C5D_cfgn1_cfgn2_cfgp1_cfgp2_beta.dimension(3) == proton.Ncfg2_I);
    assert(cache_C5D_cfgn1_cfgn2_cfgp1_cfgp2_beta.dimension(4) == neutron.Nbeta_I);
    Eigen::Tensor<doubleC, 4, Eigen::ColMajor> result_C4D_cfgn1_cfgn2_cfgp1_cfgp2 = Eigen::Tensor<doubleC, 4, Eigen::ColMajor>(neutron.Ncfg1_I, neutron.Ncfg2_I, proton.Ncfg1_I, proton.Ncfg2_I).setZero();
    const double normalization_F = (TargetTwoI_I + 1.0) / neutron.weight_F1D_beta.sum();
    RepresentationSpin wignerD_X(TargetTwoI_I);

    // K₁ = Kn₁+Kp₁; K₂ = Kn₂+Kp₂.
    for (int beta_I = 0; beta_I < neutron.Nbeta_I; ++beta_I) {
        const auto& wignerDyBeta_C2D_K_K = wignerD_X.calc_Ry(neutron.beta_F1D_beta(beta_I));
        const double weight_F = normalization_F * neutron.weight_F1D_beta(beta_I);
        for (int cfgp2_I = 0; cfgp2_I < proton.Ncfg2_I; ++cfgp2_I) {
            for (int cfgn2_I = 0; cfgn2_I < neutron.Ncfg2_I; ++cfgn2_I) {
                const int TwoK2_I = TwoKn2_I1D_cfgn2(cfgn2_I) + TwoKp2_I1D_cfgp2(cfgp2_I);
                if (std::abs(TwoK2_I) > TargetTwoI_I || (TwoK2_I + TargetTwoI_I) % 2 != 0) { continue; }
                const int K2_I = (TwoK2_I + TargetTwoI_I) / 2;
                for (int cfgp1_I = 0; cfgp1_I < proton.Ncfg1_I; ++cfgp1_I) {
                    for (int cfgn1_I = 0; cfgn1_I < neutron.Ncfg1_I; ++cfgn1_I) {
                        const int TwoK1_I = TwoKn1_I1D_cfgn1(cfgn1_I) + TwoKp1_I1D_cfgp1(cfgp1_I);
                        if (std::abs(TwoK1_I) > TargetTwoI_I || (TwoK1_I + TargetTwoI_I) % 2 != 0) { continue; }
                        const int K1_I = (TwoK1_I + TargetTwoI_I) / 2;
                        // dᴵ(β) is real in the Jz basis.
                        const double weightD_F = weight_F * wignerDyBeta_C2D_K_K(K1_I, K2_I).real();
                        result_C4D_cfgn1_cfgn2_cfgp1_cfgp2(cfgn1_I, cfgn2_I, cfgp1_I, cfgp2_I) += weightD_F * cache_C5D_cfgn1_cfgn2_cfgp1_cfgp2_beta(cfgn1_I, cfgn2_I, cfgp1_I, cfgp2_I, beta_I);
                    }
                }
            }
        }
    }
    return result_C4D_cfgn1_cfgn2_cfgp1_cfgp2;
}

inline Eigen::Tensor<doubleC, 7, Eigen::ColMajor> HFBProjectionNucleus::build_overlap() {
    auto& neutron = projection_neutron;
    auto& proton = projection_proton;

    return build_amp( [&](Eigen::TensorMap<Eigen::Tensor<doubleC, 4, Eigen::ColMajor>>& cache_C4D_cfgn1_cfgn2_cfgp1_cfgp2) {
        Eigen::Map<Eigen::MatrixXcd> overlap_C2D_cfgn1cfgn2_cfgp1cfgp2(cache_C4D_cfgn1_cfgn2_cfgp1_cfgp2.data(), Ncfgn1cfgn2_I, Ncfgp1cfgp2_I);
        const auto& overlap_C2D_cfgn1_cfgn2 = neutron.calc_overlap_pnp();
        const auto& overlap_C2D_cfgp1_cfgp2 = proton.calc_overlap_pnp();

        // ColMajor: (n₁,n₂) → rows; (p₁,p₂) → columns.
        const Eigen::Map<const Eigen::VectorXcd> overlap_C1D_cfgn1cfgn2(overlap_C2D_cfgn1_cfgn2.data(), Ncfgn1cfgn2_I);
        const Eigen::Map<const Eigen::VectorXcd> overlap_C1D_cfgp1cfgp2(overlap_C2D_cfgp1_cfgp2.data(), Ncfgp1cfgp2_I);
        overlap_C2D_cfgn1cfgn2_cfgp1cfgp2.noalias() = overlap_C1D_cfgn1cfgn2 * overlap_C1D_cfgp1cfgp2.transpose();
    });
}

inline Eigen::Tensor<doubleC, 7, Eigen::ColMajor> HFBProjectionNucleus::build_one_body_nn(const Eigen::MatrixXcd& OneBody_C2D_spn_spn) {
    auto& neutron = projection_neutron;
    auto& proton = projection_proton;
    assert(OneBody_C2D_spn_spn.rows() == neutron.Nsp_I && OneBody_C2D_spn_spn.cols() == neutron.Nsp_I && OneBody_C2D_spn_spn.allFinite());

    return build_amp([&](Eigen::TensorMap<Eigen::Tensor<doubleC, 4, Eigen::ColMajor>>& cache_C4D_cfgn1_cfgn2_cfgp1_cfgp2) {
        Eigen::Map<Eigen::MatrixXcd> onebody_C2D_cfgn1cfgn2_cfgp1cfgp2(cache_C4D_cfgn1_cfgn2_cfgp1_cfgp2.data(), Ncfgn1cfgn2_I, Ncfgp1cfgp2_I);
        const auto& overlap_C2D_cfgp1_cfgp2 = proton.calc_overlap_pnp();
        const auto& onebody_C2D_cfgn1_cfgn2 = neutron.calc_one_body_pnp(OneBody_C2D_spn_spn);

        // ColMajor: (n₁,n₂) → rows; (p₁,p₂) → columns.
        const Eigen::Map<const Eigen::VectorXcd> onebody_C1D_cfgn1cfgn2(onebody_C2D_cfgn1_cfgn2.data(), Ncfgn1cfgn2_I);
        const Eigen::Map<const Eigen::VectorXcd> overlap_C1D_cfgp1cfgp2(overlap_C2D_cfgp1_cfgp2.data(), Ncfgp1cfgp2_I);
        onebody_C2D_cfgn1cfgn2_cfgp1cfgp2.noalias() = onebody_C1D_cfgn1cfgn2 * overlap_C1D_cfgp1cfgp2.transpose();
    });
}

inline Eigen::Tensor<doubleC, 7, Eigen::ColMajor> HFBProjectionNucleus::build_one_body_pp(const Eigen::MatrixXcd& OneBody_C2D_spp_spp) {
    auto& neutron = projection_neutron;
    auto& proton = projection_proton;
    assert(OneBody_C2D_spp_spp.rows() == proton.Nsp_I && OneBody_C2D_spp_spp.cols() == proton.Nsp_I && OneBody_C2D_spp_spp.allFinite());

    return build_amp([&](Eigen::TensorMap<Eigen::Tensor<doubleC, 4, Eigen::ColMajor>>& cache_C4D_cfgn1_cfgn2_cfgp1_cfgp2) {
        Eigen::Map<Eigen::MatrixXcd> onebody_C2D_cfgn1cfgn2_cfgp1cfgp2(cache_C4D_cfgn1_cfgn2_cfgp1_cfgp2.data(), Ncfgn1cfgn2_I, Ncfgp1cfgp2_I);
        const auto& overlap_C2D_cfgn1_cfgn2 = neutron.calc_overlap_pnp();
        const auto& onebody_C2D_cfgp1_cfgp2 = proton.calc_one_body_pnp(OneBody_C2D_spp_spp);

        // ColMajor: (n₁,n₂) → rows; (p₁,p₂) → columns.
        const Eigen::Map<const Eigen::VectorXcd> overlap_C1D_cfgn1cfgn2(overlap_C2D_cfgn1_cfgn2.data(), Ncfgn1cfgn2_I);
        const Eigen::Map<const Eigen::VectorXcd> onebody_C1D_cfgp1cfgp2(onebody_C2D_cfgp1_cfgp2.data(), Ncfgp1cfgp2_I);
        onebody_C2D_cfgn1cfgn2_cfgp1cfgp2.noalias() = overlap_C1D_cfgn1cfgn2 * onebody_C1D_cfgp1cfgp2.transpose();
    });
}

inline Eigen::Tensor<doubleC, 7, Eigen::ColMajor> HFBProjectionNucleus::build_one_body_pn(const Eigen::MatrixXcd& OneBody_C2D_spp_spn) {
    auto& neutron = projection_neutron;
    auto& proton = projection_proton;
    assert(OneBody_C2D_spp_spn.rows() == proton.Nsp_I && OneBody_C2D_spp_spn.cols() == neutron.Nsp_I && OneBody_C2D_spp_spn.allFinite());

    // ColMajor: (n₁,n₂) → rows; (p₁,p₂) → columns.

    // cπ†cν: (-1)^r₁; r₁ counts left neutron quasiparticles.
    Eigen::VectorXd sign_F1D_cfgn1cfgn2(Ncfgn1cfgn2_I);
    for (int cfgn2_I = 0; cfgn2_I < neutron.Ncfg2_I; ++cfgn2_I) {
        for (int cfgn1_I = 0; cfgn1_I < neutron.Ncfg1_I; ++cfgn1_I) {
            const double sign_F = 1.0 - 2.0 * (neutron.hfb_pfaffian.config1_I2D_cfg1_cqp1[cfgn1_I].size() % 2);
            const Eigen::Index row_I = cfgn1_I + static_cast<Eigen::Index>(neutron.Ncfg1_I) * cfgn2_I;
            sign_F1D_cfgn1cfgn2(row_I) = sign_F;
        }
    }

    return build_amp([&](Eigen::TensorMap<Eigen::Tensor<doubleC, 4, Eigen::ColMajor>>& cache_C4D_cfgn1_cfgn2_cfgp1_cfgp2) {
        Eigen::Map<Eigen::MatrixXcd> onebody_C2D_cfgn1cfgn2_cfgp1cfgp2(cache_C4D_cfgn1_cfgn2_cfgp1_cfgp2.data(), Ncfgn1cfgn2_I, Ncfgp1cfgp2_I);

        // Aν and Cπ include their gauge normalizations.
        const auto& annihilator_C3D_cfgn1_cfgn2_spn = neutron.calc_one_fermion_pnp(false);
        const auto& creator_C3D_cfgp1_cfgp2_spp = proton.calc_one_fermion_pnp(true);
        const Eigen::Map<const Eigen::MatrixXcd> annihilator_C2D_cfgn1cfgn2_spn(annihilator_C3D_cfgn1_cfgn2_spn.data(), Ncfgn1cfgn2_I, neutron.Nsp_I);
        const Eigen::Map<const Eigen::MatrixXcd> creator_C2D_cfgp1cfgp2_spp(creator_C3D_cfgp1_cfgp2_spp.data(), Ncfgp1cfgp2_I, proton.Nsp_I);

        // H¹ = diag(sign) Aν Oᵀ Cπᵀ; no conjugation.
        onebody_C2D_cfgn1cfgn2_cfgp1cfgp2.noalias() = sign_F1D_cfgn1cfgn2.asDiagonal() * annihilator_C2D_cfgn1cfgn2_spn * OneBody_C2D_spp_spn.transpose() * creator_C2D_cfgp1cfgp2_spp.transpose();
    });
}

inline Eigen::Tensor<doubleC, 7, Eigen::ColMajor> HFBProjectionNucleus::build_two_body_nn(const Eigen::VectorXd& gQ_F1D_i, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& Q_C3D_spn_spn_i, const Eigen::VectorXd& gP_F1D_i, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& P_C3D_spn_spn_i) {
    auto& neutron = projection_neutron;
    auto& proton = projection_proton;

    return build_amp([&](Eigen::TensorMap<Eigen::Tensor<doubleC, 4, Eigen::ColMajor>>& cache_C4D_cfgn1_cfgn2_cfgp1_cfgp2) {
        Eigen::Map<Eigen::MatrixXcd> twobody_C2D_cfgn1cfgn2_cfgp1cfgp2(cache_C4D_cfgn1_cfgn2_cfgp1_cfgp2.data(), Ncfgn1cfgn2_I, Ncfgp1cfgp2_I);
        // Hν²(Ω) Nπ(Ω); weighted Q/P kernel; pairing ¼ included.
        const auto& overlap_C2D_cfgp1_cfgp2 = proton.calc_overlap_pnp();
        const auto& twobody_C2D_cfgn1_cfgn2 = neutron.calc_two_body_pnp(gQ_F1D_i, Q_C3D_spn_spn_i, gP_F1D_i, P_C3D_spn_spn_i);

        // ColMajor: (n₁,n₂) → rows; (p₁,p₂) → columns.
        const Eigen::Map<const Eigen::VectorXcd> twobody_C1D_cfgn1cfgn2(twobody_C2D_cfgn1_cfgn2.data(), Ncfgn1cfgn2_I);
        const Eigen::Map<const Eigen::VectorXcd> overlap_C1D_cfgp1cfgp2(overlap_C2D_cfgp1_cfgp2.data(), Ncfgp1cfgp2_I);
        twobody_C2D_cfgn1cfgn2_cfgp1cfgp2.noalias() = twobody_C1D_cfgn1cfgn2 * overlap_C1D_cfgp1cfgp2.transpose();
    });
}

inline Eigen::Tensor<doubleC, 7, Eigen::ColMajor> HFBProjectionNucleus::build_two_body_pp(const Eigen::VectorXd& gQ_F1D_i, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& Q_C3D_spp_spp_i, const Eigen::VectorXd& gP_F1D_i, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& P_C3D_spp_spp_i) {
    auto& neutron = projection_neutron;
    auto& proton = projection_proton;

    return build_amp([&](Eigen::TensorMap<Eigen::Tensor<doubleC, 4, Eigen::ColMajor>>& cache_C4D_cfgn1_cfgn2_cfgp1_cfgp2) {
        Eigen::Map<Eigen::MatrixXcd> twobody_C2D_cfgn1cfgn2_cfgp1cfgp2(cache_C4D_cfgn1_cfgn2_cfgp1_cfgp2.data(), Ncfgn1cfgn2_I, Ncfgp1cfgp2_I);
        // Nν(Ω) Hπ²(Ω); weighted Q/P kernel; pairing ¼ included.
        const auto& overlap_C2D_cfgn1_cfgn2 = neutron.calc_overlap_pnp();
        const auto& twobody_C2D_cfgp1_cfgp2 = proton.calc_two_body_pnp(gQ_F1D_i, Q_C3D_spp_spp_i, gP_F1D_i, P_C3D_spp_spp_i);

        // ColMajor: (n₁,n₂) → rows; (p₁,p₂) → columns.
        const Eigen::Map<const Eigen::VectorXcd> overlap_C1D_cfgn1cfgn2(overlap_C2D_cfgn1_cfgn2.data(), Ncfgn1cfgn2_I);
        const Eigen::Map<const Eigen::VectorXcd> twobody_C1D_cfgp1cfgp2(twobody_C2D_cfgp1_cfgp2.data(), Ncfgp1cfgp2_I);
        twobody_C2D_cfgn1cfgn2_cfgp1cfgp2.noalias() = overlap_C1D_cfgn1cfgn2 * twobody_C1D_cfgp1cfgp2.transpose();
    });
}

inline Eigen::Tensor<doubleC, 7, Eigen::ColMajor> HFBProjectionNucleus::build_two_body_pn(const Eigen::VectorXd& gQ_F1D_i, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& Q_C3D_spn_spn_i, const Eigen::Tensor<doubleC, 3, Eigen::ColMajor>& Q_C3D_spp_spp_i) {
    auto& neutron = projection_neutron;
    auto& proton = projection_proton;

    // NQν = NQπ = NgQ.
    assert(gQ_F1D_i.allFinite());
    assert(Q_C3D_spn_spn_i.dimension(0) == neutron.Nsp_I && Q_C3D_spn_spn_i.dimension(1) == neutron.Nsp_I && Q_C3D_spn_spn_i.dimension(2) == gQ_F1D_i.size());
    assert(Q_C3D_spp_spp_i.dimension(0) == proton.Nsp_I && Q_C3D_spp_spp_i.dimension(1) == proton.Nsp_I && Q_C3D_spp_spp_i.dimension(2) == gQ_F1D_i.size());
    assert(Eigen::Map<const Eigen::VectorXcd>(Q_C3D_spn_spn_i.data(), Q_C3D_spn_spn_i.size()).allFinite());
    assert(Eigen::Map<const Eigen::VectorXcd>(Q_C3D_spp_spp_i.data(), Q_C3D_spp_spp_i.size()).allFinite());

    // (Qν†,Qπ) and (Qν,Qπ†) reuse species workspaces.
    Eigen::MatrixXcd Qdag_C2D_spn_spn{};
    Eigen::MatrixXcd Q_C2D_spp_spp{};
    Qdag_C2D_spn_spn.resize(neutron.Nsp_I, neutron.Nsp_I);
    Q_C2D_spp_spp.resize(proton.Nsp_I, proton.Nsp_I);

    return build_amp([&](Eigen::TensorMap<Eigen::Tensor<doubleC, 4, Eigen::ColMajor>>& cache_C4D_cfgn1_cfgn2_cfgp1_cfgp2) {
        Eigen::Map<Eigen::MatrixXcd> twobody_C2D_cfgn1cfgn2_cfgp1cfgp2(cache_C4D_cfgn1_cfgn2_cfgp1_cfgp2.data(), Ncfgn1cfgn2_I, Ncfgp1cfgp2_I);

        twobody_C2D_cfgn1cfgn2_cfgp1cfgp2.setZero();
        for (Eigen::Index iQ_I = 0; iQ_I < gQ_F1D_i.size(); ++iQ_I) {
            // Hi = ½ gQi (Qνi† Qπi + Qνi Qπi†).
            const Eigen::Map<const Eigen::MatrixXcd> Q_C2D_spn_spn(Q_C3D_spn_spn_i.data() + static_cast<Eigen::Index>(iQ_I) * neutron.Nsp_I * neutron.Nsp_I, neutron.Nsp_I, neutron.Nsp_I);
            const Eigen::Map<const Eigen::MatrixXcd> Qp_C2D_spp_spp(Q_C3D_spp_spp_i.data() + static_cast<Eigen::Index>(iQ_I) * proton.Nsp_I * proton.Nsp_I, proton.Nsp_I, proton.Nsp_I);
            // ½ gQi Qνi† Qπi.
            {
                Qdag_C2D_spn_spn = Q_C2D_spn_spn.adjoint();
                Q_C2D_spp_spp = Qp_C2D_spp_spp;
                const auto& onebody_C2D_cfgn1_cfgn2 = neutron.calc_one_body_pnp(Qdag_C2D_spn_spn);
                const auto& onebody_C2D_cfgp1_cfgp2 = proton.calc_one_body_pnp(Q_C2D_spp_spp);

                // Both kernels include overlaps and gauge normalizations.
                const Eigen::Map<const Eigen::VectorXcd> onebody_C1D_cfgn1cfgn2(onebody_C2D_cfgn1_cfgn2.data(), Ncfgn1cfgn2_I);
                const Eigen::Map<const Eigen::VectorXcd> onebody_C1D_cfgp1cfgp2(onebody_C2D_cfgp1_cfgp2.data(), Ncfgp1cfgp2_I);
                twobody_C2D_cfgn1cfgn2_cfgp1cfgp2.noalias() += 0.5 * gQ_F1D_i(iQ_I) * onebody_C1D_cfgn1cfgn2 * onebody_C1D_cfgp1cfgp2.transpose();
            }
            // ½ gQi Qνi Qπi†.
            {
                Qdag_C2D_spn_spn = Q_C2D_spn_spn;
                Q_C2D_spp_spp = Qp_C2D_spp_spp.adjoint();
                const auto& onebody_C2D_cfgn1_cfgn2 = neutron.calc_one_body_pnp(Qdag_C2D_spn_spn);
                const auto& onebody_C2D_cfgp1_cfgp2 = proton.calc_one_body_pnp(Q_C2D_spp_spp);

                // Both kernels include overlaps and gauge normalizations.
                const Eigen::Map<const Eigen::VectorXcd> onebody_C1D_cfgn1cfgn2(onebody_C2D_cfgn1_cfgn2.data(), Ncfgn1cfgn2_I);
                const Eigen::Map<const Eigen::VectorXcd> onebody_C1D_cfgp1cfgp2(onebody_C2D_cfgp1_cfgp2.data(), Ncfgp1cfgp2_I);
                twobody_C2D_cfgn1cfgn2_cfgp1cfgp2.noalias() += 0.5 * gQ_F1D_i(iQ_I) * onebody_C1D_cfgn1cfgn2 * onebody_C1D_cfgp1cfgp2.transpose();
            }
        }
    });
}
