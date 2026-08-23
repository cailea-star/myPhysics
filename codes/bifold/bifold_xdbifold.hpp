/**
 * @file    bifold_xdbifold.hpp
 * @author  cailea
 * @date    2026-08-22
 * @brief   Density-dependent XDM3Y double-folding potentials.
 * @ref     Computer Physics Communications 284 (2023) 108613
 */

#pragma once

#include <cassert>
#include <cmath>
#include <utility>

#include <Eigen/Core>
#include <unsupported/Eigen/CXX11/Tensor>

#include "bifold_bifold.hpp"
#include "bifold_xdm3y.hpp"

/**
 * @brief Precompute density-dependent XDM3Y folding kernels.
 */
class XDBifold : public Bifold {
public:
    XDM3YParameters xdm3y_parameters; // {C, α, β, γ, n, K}.

protected:
    Eigen::VectorXd rhopmexp_F1D_k = Eigen::VectorXd::Zero(Nk_I); // ℱ₀[ρ_{p,m}e⁻ᵝρᵖᵐ].
    Eigen::VectorXd rhopm2_F1D_k = Eigen::VectorXd::Zero(Nk_I); // ℱ₀[ρ_{p,m}²].
    Eigen::VectorXd rhopm3_F1D_k = Eigen::VectorXd::Zero(Nk_I); // ℱ₀[ρ_{p,m}³].
    Eigen::VectorXd rhopm4_F1D_k = Eigen::VectorXd::Zero(Nk_I); // ℱ₀[ρ_{p,m}⁴].
    Eigen::MatrixXd rhotmexp_F2D_lambda_k = Eigen::MatrixXd::Zero(Nlambda_I, Nk_I); // ℱ_λ[ρ_{t,m}e⁻ᵝρᵗᵐ].
    Eigen::MatrixXd rhotm2_F2D_lambda_k = Eigen::MatrixXd::Zero(Nlambda_I, Nk_I); // ℱ_λ[ρ_{t,m}²].
    Eigen::MatrixXd rhotm3_F2D_lambda_k = Eigen::MatrixXd::Zero(Nlambda_I, Nk_I); // ℱ_λ[ρ_{t,m}³].
    Eigen::MatrixXd rhotm4_F2D_lambda_k = Eigen::MatrixXd::Zero(Nlambda_I, Nk_I); // ℱ_λ[ρ_{t,m}⁴].
    Eigen::MatrixXd fpexp_F2D_s_k = Eigen::MatrixXd::Zero(Ns_I, Nk_I); // f̃_p[ρ_{p,m}e⁻ᵝρᵖᵐ].
    Eigen::MatrixXd fp2_F2D_s_k = Eigen::MatrixXd::Zero(Ns_I, Nk_I); // f̃_p[ρ_{p,m}²].
    Eigen::MatrixXd fp3_F2D_s_k = Eigen::MatrixXd::Zero(Ns_I, Nk_I); // f̃_p[ρ_{p,m}³].
    Eigen::MatrixXd fp4_F2D_s_k = Eigen::MatrixXd::Zero(Ns_I, Nk_I); // f̃_p[ρ_{p,m}⁴].
    Eigen::Tensor<double, 3, Eigen::ColMajor> ftexp_F3D_lambda_s_k = Eigen::Tensor<double, 3, Eigen::ColMajor>(Nlambda_I, Ns_I, Nk_I); // f̃_{tλ}[ρ_{t,m}e⁻ᵝρᵗᵐ].
    Eigen::Tensor<double, 3, Eigen::ColMajor> ft2_F3D_lambda_s_k = Eigen::Tensor<double, 3, Eigen::ColMajor>(Nlambda_I, Ns_I, Nk_I); // f̃_{tλ}[ρ_{t,m}²].
    Eigen::Tensor<double, 3, Eigen::ColMajor> ft3_F3D_lambda_s_k = Eigen::Tensor<double, 3, Eigen::ColMajor>(Nlambda_I, Ns_I, Nk_I); // f̃_{tλ}[ρ_{t,m}³].
    Eigen::Tensor<double, 3, Eigen::ColMajor> ft4_F3D_lambda_s_k = Eigen::Tensor<double, 3, Eigen::ColMajor>(Nlambda_I, Ns_I, Nk_I); // f̃_{tλ}[ρ_{t,m}⁴].

public:
    /**
     * @brief  Construct and precompute an XDM3Y folding model.
     * @math   ρ_pρ_tF(ρ_p+ρ_t) → {ρ̃, f̃}
     * @output Initialized XDM3Y folding model.
     * @note   Supports γ = 0 or n ∈ {1, 2, 3}.
     */
    XDBifold(double hmass_F_, double Ecm_F_, double Ze2_F_, XDM3YFunctions xdm3y_functions_, DensityProjec density_projec_, DensityTarget density_target_, double Cs_F_ = 1.0 / 36.0)
        : Bifold(hmass_F_, Ecm_F_, Ze2_F_, std::move(xdm3y_functions_.m3y_functions), std::move(density_projec_), std::move(density_target_), Cs_F_), xdm3y_parameters(std::move(xdm3y_functions_.xdm3y_parameters)) {
        assert(xdm3y_parameters.gamma_F == 0.0 || xdm3y_parameters.n_F == 1.0 || xdm3y_parameters.n_F == 2.0 || xdm3y_parameters.n_F == 3.0);

        // ρ_m → {ρ̃_m e⁻ᵝρᵐ, ρ̃_m², ρ̃_m³, ρ̃_m⁴, f̃}.
        fill_projec_exp(density_projec, xdm3y_parameters.beta_F);
        fill_projec_power(density_projec, 1, rhopm2_F1D_k, fp2_F2D_s_k);
        fill_projec_power(density_projec, 2, rhopm3_F1D_k, fp3_F2D_s_k);
        fill_projec_power(density_projec, 3, rhopm4_F1D_k, fp4_F2D_s_k);
        fill_target_exp(density_target, xdm3y_parameters.beta_F);
        fill_target_power(density_target, 1, rhotm2_F2D_lambda_k, ft2_F3D_lambda_s_k);
        fill_target_power(density_target, 2, rhotm3_F2D_lambda_k, ft3_F3D_lambda_s_k);
        fill_target_power(density_target, 3, rhotm4_F2D_lambda_k, ft4_F3D_lambda_s_k);
    }

    /**
     * @brief  Replace XDM3Y kernels and refresh exponential tables.
     * @math   {𝒱_M, 𝒳} → {𝒱_M′, 𝒳′}
     * @output Updated interaction and exponential tables.
     * @note   Supports γ = 0 or n ∈ {1, 2, 3}.
     */
    void update_potentials(const XDM3YFunctions& xdm3y_functions_);

    /**
     * @brief  Recompute projectile XDM3Y tables.
     * @math   ρ_p → {ρ̃_p, ρ̃_p e⁻ᵝρᵖ, ρ̃_p², ρ̃_p³, ρ̃_p⁴, f̃_p}
     * @output Updated projectile tables.
     */
    void update_density_projec(const DensityProjec& density_projec_) override;

    /**
     * @brief  Recompute target XDM3Y tables.
     * @math   ρ_t → {ρ̃_{tλ}, ρ̃_{tλ}e⁻ᵝρᵗ, ρ̃_{tλ}², ρ̃_{tλ}³, ρ̃_{tλ}⁴, f̃_{tλ}}
     * @output Updated target tables.
     */
    void update_density_target(const DensityTarget& density_target_) override;

    /**
     * @brief  Evaluate one density-dependent nuclear multipole.
     * @math   U_λ = Cℱ_λ⁻¹[ρ_pρ_t(1 + αe⁻ᵝ⁽ρᵖ⁺ρᵗ⁾ - γ(ρ_p+ρ_t)ⁿ)ṽ]
     * @output Nuclear-potential multipole, MeV.
     * @note   Uses the finite FDA binomial expansion.
     */
    double Unn_lambda(double r_F, int lambda_I, const Real2RealFunc& vnn_k_Func) const override;

protected:

    /**
     * @brief  Transform the projectile exponential FDA factor.
     * @math   ρ_{p,m}e⁻ᵝρᵖᵐ → {ρ̃_{p,m}^{exp}, f̃_p^{exp}}
     * @output Updated projectile exponential tables.
     */
    void fill_projec_exp(const DensityProjec& density_projec_, double beta_F);

    /**
     * @brief  Transform one projectile FDA power.
     * @math   ρ_{p,m}^{n+1} → {ρ̃_{p,m}^{n+1}, f̃_p^{n+1}}
     * @output Updated projectile power tables.
     * @note   Requires n ∈ {1, 2, 3}.
     */
    void fill_projec_power(const DensityProjec& density_projec_, int n_I, Eigen::VectorXd& rhopm_power_F1D_k_, Eigen::MatrixXd& fp_power_F2D_s_k_);

    /**
     * @brief  Transform the target exponential FDA factor.
     * @math   ρ_{t,m}e⁻ᵝρᵗᵐ → {ρ̃_{tλ,m}^{exp}, f̃_{tλ}^{exp}}
     * @output Updated target exponential tables.
     */
    void fill_target_exp(const DensityTarget& density_target_, double beta_F);

    /**
     * @brief  Transform one target FDA power.
     * @math   ρ_{t,m}^{n+1} → {ρ̃_{tλ,m}^{n+1}, f̃_{tλ}^{n+1}}
     * @output Updated target power tables.
     * @note   Requires n ∈ {1, 2, 3}.
     */
    void fill_target_power(const DensityTarget& density_target_, int n_I, Eigen::MatrixXd& rhotm_power_F2D_lambda_k_, Eigen::Tensor<double, 3, Eigen::ColMajor>& ft_power_F3D_lambda_s_k_);

    /**
     * @brief  Evaluate density-dependent exchange overlaps.
     * @math   G_λ = C𝒢_λ[ρ_pρ_t(1 + αe⁻ᵝ⁽ρᵖ⁺ρᵗ⁾ - γ(ρ_p+ρ_t)ⁿ)]
     * @output Exchange-overlap values over s.
     * @note   Uses the finite FDA binomial expansion.
     */
    Eigen::VectorXd G_lambda(double r_F, int lambda_I) const override;
};

inline void XDBifold::fill_projec_exp(const DensityProjec& density_projec_, double beta_F) {
    auto rhopmexp_Func = [&](double rp_F) {
        double rhopm_F = density_projec_.rho_matter(rp_F);
        return rhopm_F * std::exp(-beta_F * rhopm_F);
    };
    fill_rhop(rhopmexp_Func, rhopmexp_F1D_k);
    fill_fp(rhopmexp_Func, fpexp_F2D_s_k);
}

inline void XDBifold::fill_projec_power(const DensityProjec& density_projec_, int n_I, Eigen::VectorXd& rhopm_power_F1D_k_, Eigen::MatrixXd& fp_power_F2D_s_k_) {
    assert(n_I >= 1 && n_I <= 3);
    auto rhopm_power_Func = [&](double rp_F) {
        double rhopm_F = density_projec_.rho_matter(rp_F);
        return rhopm_F * std::pow(rhopm_F, n_I);
    };
    fill_rhop(rhopm_power_Func, rhopm_power_F1D_k_);
    fill_fp(rhopm_power_Func, fp_power_F2D_s_k_);
}

inline void XDBifold::update_density_projec(const DensityProjec& density_projec_) {
    // ρ_p → {ρ̃_p, ρ̃_p^{exp}, ρ̃_p², ρ̃_p³, ρ̃_p⁴, f̃_p}.
    Bifold::update_density_projec(density_projec_);
    fill_projec_exp(density_projec, xdm3y_parameters.beta_F);
    fill_projec_power(density_projec, 1, rhopm2_F1D_k, fp2_F2D_s_k);
    fill_projec_power(density_projec, 2, rhopm3_F1D_k, fp3_F2D_s_k);
    fill_projec_power(density_projec, 3, rhopm4_F1D_k, fp4_F2D_s_k);
}

inline void XDBifold::fill_target_exp(const DensityTarget& density_target_, double beta_F) {
    auto rhotmexp_Func = [&](double rt_F, double theta_F) {
        double rhotm_F = density_target_.rho_matter(rt_F, theta_F);
        return rhotm_F * std::exp(-beta_F * rhotm_F);
    };
    fill_rhot(rhotmexp_Func, rhotmexp_F2D_lambda_k);
    fill_ft(rhotmexp_Func, ftexp_F3D_lambda_s_k);
}

inline void XDBifold::fill_target_power(const DensityTarget& density_target_, int n_I, Eigen::MatrixXd& rhotm_power_F2D_lambda_k_, Eigen::Tensor<double, 3, Eigen::ColMajor>& ft_power_F3D_lambda_s_k_) {
    assert(n_I >= 1 && n_I <= 3);
    auto rhotm_power_Func = [&](double rt_F, double theta_F) {
        double rhotm_F = density_target_.rho_matter(rt_F, theta_F);
        return rhotm_F * std::pow(rhotm_F, n_I);
    };
    fill_rhot(rhotm_power_Func, rhotm_power_F2D_lambda_k_);
    fill_ft(rhotm_power_Func, ft_power_F3D_lambda_s_k_);
}

inline void XDBifold::update_density_target(const DensityTarget& density_target_) {
    // ρ_t → {ρ̃_{tλ}, ρ̃_{tλ}^{exp}, ρ̃_{tλ}², ρ̃_{tλ}³, ρ̃_{tλ}⁴, f̃_{tλ}}.
    Bifold::update_density_target(density_target_);
    fill_target_exp(density_target, xdm3y_parameters.beta_F);
    fill_target_power(density_target, 1, rhotm2_F2D_lambda_k, ft2_F3D_lambda_s_k);
    fill_target_power(density_target, 2, rhotm3_F2D_lambda_k, ft3_F3D_lambda_s_k);
    fill_target_power(density_target, 3, rhotm4_F2D_lambda_k, ft4_F3D_lambda_s_k);
}

inline void XDBifold::update_potentials(const XDM3YFunctions& xdm3y_functions_) {
    assert(xdm3y_functions_.xdm3y_parameters.gamma_F == 0.0 || xdm3y_functions_.xdm3y_parameters.n_F == 1.0 || xdm3y_functions_.xdm3y_parameters.n_F == 2.0 || xdm3y_functions_.xdm3y_parameters.n_F == 3.0);
    Bifold::update_potentials(xdm3y_functions_.m3y_functions);
    xdm3y_parameters = xdm3y_functions_.xdm3y_parameters;
    fill_projec_exp(density_projec, xdm3y_parameters.beta_F);
    fill_target_exp(density_target, xdm3y_parameters.beta_F);
}

inline double XDBifold::Unn_lambda(double r_F, int lambda_I, const Real2RealFunc& vnn_k_Func) const {
    double Up1t1_F = calc_Ud_lambda(r_F, lambda_I, rhopm_F1D_k, rhotm_F2D_lambda_k, vnn_k_Func); // ρ_pρ_t → U_{p1t1}.
    double Uexp_F = calc_Ud_lambda(r_F, lambda_I, rhopmexp_F1D_k, rhotmexp_F2D_lambda_k, vnn_k_Func);
    Up1t1_F += xdm3y_parameters.alpha_F * Uexp_F;

    if (xdm3y_parameters.n_F == 1.0) {
        double Up2t1_F = calc_Ud_lambda(r_F, lambda_I, rhopm2_F1D_k, rhotm_F2D_lambda_k, vnn_k_Func);
        double Up1t2_F = calc_Ud_lambda(r_F, lambda_I, rhopm_F1D_k, rhotm2_F2D_lambda_k, vnn_k_Func);
        double Upower1_F = Up2t1_F + Up1t2_F;
        Up1t1_F -= xdm3y_parameters.gamma_F * Upower1_F;
    }
    if (xdm3y_parameters.n_F == 2.0) {
        double Up3t1_F = calc_Ud_lambda(r_F, lambda_I, rhopm3_F1D_k, rhotm_F2D_lambda_k, vnn_k_Func);
        double Up2t2_F = calc_Ud_lambda(r_F, lambda_I, rhopm2_F1D_k, rhotm2_F2D_lambda_k, vnn_k_Func);
        double Up1t3_F = calc_Ud_lambda(r_F, lambda_I, rhopm_F1D_k, rhotm3_F2D_lambda_k, vnn_k_Func);
        double Upower2_F = Up3t1_F + 2.0 * Up2t2_F + Up1t3_F;
        Up1t1_F -= xdm3y_parameters.gamma_F * Upower2_F;
    }
    if (xdm3y_parameters.n_F == 3.0) {
        double Up4t1_F = calc_Ud_lambda(r_F, lambda_I, rhopm4_F1D_k, rhotm_F2D_lambda_k, vnn_k_Func);
        double Up3t2_F = calc_Ud_lambda(r_F, lambda_I, rhopm3_F1D_k, rhotm2_F2D_lambda_k, vnn_k_Func);
        double Up2t3_F = calc_Ud_lambda(r_F, lambda_I, rhopm2_F1D_k, rhotm3_F2D_lambda_k, vnn_k_Func);
        double Up1t4_F = calc_Ud_lambda(r_F, lambda_I, rhopm_F1D_k, rhotm4_F2D_lambda_k, vnn_k_Func);
        double Upower3_F = Up4t1_F + 3.0 * Up3t2_F + 3.0 * Up2t3_F + Up1t4_F;
        Up1t1_F -= xdm3y_parameters.gamma_F * Upower3_F;
    }
    return xdm3y_parameters.C_F * Up1t1_F;
}

inline Eigen::VectorXd XDBifold::G_lambda(double r_F, int lambda_I) const {
    Eigen::VectorXd Gp1t1_F1D_s = calc_G_lambda(r_F, lambda_I, fp_F2D_s_k, ft_F3D_lambda_s_k); // ρ_pρ_t → G_{p1t1}.
    Eigen::VectorXd Gexp_F1D_s = calc_G_lambda(r_F, lambda_I, fpexp_F2D_s_k, ftexp_F3D_lambda_s_k);
    Gp1t1_F1D_s += xdm3y_parameters.alpha_F * Gexp_F1D_s;

    if (xdm3y_parameters.n_F == 1.0) {
        Eigen::VectorXd Gp2t1_F1D_s = calc_G_lambda(r_F, lambda_I, fp2_F2D_s_k, ft_F3D_lambda_s_k);
        Eigen::VectorXd Gp1t2_F1D_s = calc_G_lambda(r_F, lambda_I, fp_F2D_s_k, ft2_F3D_lambda_s_k);
        Eigen::VectorXd Gpower1_F1D_s = Gp2t1_F1D_s + Gp1t2_F1D_s;
        Gp1t1_F1D_s -= xdm3y_parameters.gamma_F * Gpower1_F1D_s;
    }
    if (xdm3y_parameters.n_F == 2.0) {
        Eigen::VectorXd Gp3t1_F1D_s = calc_G_lambda(r_F, lambda_I, fp3_F2D_s_k, ft_F3D_lambda_s_k);
        Eigen::VectorXd Gp2t2_F1D_s = calc_G_lambda(r_F, lambda_I, fp2_F2D_s_k, ft2_F3D_lambda_s_k);
        Eigen::VectorXd Gp1t3_F1D_s = calc_G_lambda(r_F, lambda_I, fp_F2D_s_k, ft3_F3D_lambda_s_k);
        Eigen::VectorXd Gpower2_F1D_s = Gp3t1_F1D_s + 2.0 * Gp2t2_F1D_s + Gp1t3_F1D_s;
        Gp1t1_F1D_s -= xdm3y_parameters.gamma_F * Gpower2_F1D_s;
    }
    if (xdm3y_parameters.n_F == 3.0) {
        Eigen::VectorXd Gp4t1_F1D_s = calc_G_lambda(r_F, lambda_I, fp4_F2D_s_k, ft_F3D_lambda_s_k);
        Eigen::VectorXd Gp3t2_F1D_s = calc_G_lambda(r_F, lambda_I, fp3_F2D_s_k, ft2_F3D_lambda_s_k);
        Eigen::VectorXd Gp2t3_F1D_s = calc_G_lambda(r_F, lambda_I, fp2_F2D_s_k, ft3_F3D_lambda_s_k);
        Eigen::VectorXd Gp1t4_F1D_s = calc_G_lambda(r_F, lambda_I, fp_F2D_s_k, ft4_F3D_lambda_s_k);
        Eigen::VectorXd Gpower3_F1D_s = Gp4t1_F1D_s + 3.0 * Gp3t2_F1D_s + 3.0 * Gp2t3_F1D_s + Gp1t4_F1D_s;
        Gp1t1_F1D_s -= xdm3y_parameters.gamma_F * Gpower3_F1D_s;
    }
    return xdm3y_parameters.C_F * Gp1t1_F1D_s;
}
