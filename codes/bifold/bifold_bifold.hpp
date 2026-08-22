/**
 * @file    bifold_bifold.hpp
 * @author  cailea
 * @date    2026-08-22
 * @brief   Density-independent double-folding potentials.
 */

#pragma once

#include <cassert>
#include <cmath>
#include <numbers>
#include <utility>

#include <Eigen/Core>
#include <unsupported/Eigen/CXX11/Tensor>

#include "bifold_density.hpp"
#include "bifold_fermik.hpp"
#include "bifold_vm3y.hpp"
#include "integration_gauss.hpp"
#include "spherical_fourier.hpp"
#include "spherical_harmonics.hpp"

/**
 * @brief Precompute density-independent double-folding kernels.
 */
class Bifold {
public:
    double hmass_F = 0.0; // h_μ = (ℏc)²/(2μc²), MeV fm².
    double Ecm_F = 0.0; // Center-of-mass energy, MeV.
    double Ze2_F = 0.0; // Coulomb coupling, MeV fm.
    double Cs_F = 1.0 / 36.0; // Extended Thomas-Fermi gradient coefficient.
    VnnFunctions vnn_functions; // {v_C, v_D, v_E, v_E^ZR, g(E)}.
    DensityInfoProjec density_projec; // {A_p, Z_p, ρ_p}.
    DensityInfoTarget density_target; // {A_t, Z_t, ρ_t}.
    int Nlambda_I = 4; // λ ∈ {0, 2, 4, 6}.
    Eigen::VectorXi lambda_I1D_lambda = Eigen::VectorXi::LinSpaced(Nlambda_I, 0, 6); // Multipole values.

    /**
     * @brief  Construct an empty folding model.
     * @math   𝓑 = ∅
     * @output Empty folding model.
     */
    Bifold() {ft_F3D_lambda_s_k.setZero();}

    /**
     * @brief  Construct and precompute a folding model.
     * @math   {ρ_p, ρ_t, f_p, f_t} → k-space tables
     * @output Initialized folding model.
     * @note   Uses fixed 201-point Simpson grids.
     */
    Bifold(double hmass_F_, double Ecm_F_, double Ze2_F_, VnnFunctions vnn_functions_, DensityInfoProjec density_projec_, DensityInfoTarget density_target_, double Cs_F_ = 1.0 / 36.0)
        : hmass_F(hmass_F_), Ecm_F(Ecm_F_), Ze2_F(Ze2_F_), Cs_F(Cs_F_), vnn_functions(std::move(vnn_functions_)), density_projec(std::move(density_projec_)), density_target(std::move(density_target_)), rp_F1D_rp(Eigen::VectorXd::LinSpaced(Nr_I, 1.0e-10, density_projec.rmax_F)), rt_F1D_rt(Eigen::VectorXd::LinSpaced(Nr_I, 1.0e-10, density_target.rmax_F)) {
        assert(hmass_F > 0.0);
        assert(Ecm_F >= 0.0);
        assert(Ze2_F >= 0.0);
        assert(Cs_F >= 0.0);
        assert(density_projec.rmax_F > 0.0);
        assert(density_target.rmax_F > 0.0);
        assert(static_cast<bool>(vnn_functions.vcoul_k_Func));
        assert(static_cast<bool>(vnn_functions.vnnd_k_Func));
        assert(static_cast<bool>(vnn_functions.vnne_r_Func));
        assert(static_cast<bool>(vnn_functions.vnne_zr_Func));
        assert(static_cast<bool>(vnn_functions.gE_Func));
        ft_F3D_lambda_s_k.setZero();

        // {λ, m = 0} → {w_θ, Y_λ0(θ), ρ̃, f̃}.
        GaussLegendreMeshes spherical_meshes(lambda_I1D_lambda, Eigen::VectorXi::Zero(Nlambda_I), Ntheta_I);
        w_F1D_theta = spherical_meshes.w_F1D_x;
        Y_F2D_lambda_theta.noalias() = std::sqrt(0.5 / std::numbers::pi) * spherical_meshes.N_F1D_lm.asDiagonal() * spherical_meshes.P_F2D_lm_x;
        Real2RealFunc rhopm_Func = [&](double rp_F) { return density_projec.rho_matter(rp_F); };
        Real2RealFunc rhopc_Func = [&](double rp_F) { return density_projec.rho_charge(rp_F); };
        RealReal2RealFunc rhotm_Func = [&](double rt_F, double theta_F) { return density_target.rho_matter(rt_F, theta_F); };
        RealReal2RealFunc rhotc_Func = [&](double rt_F, double theta_F) { return density_target.rho_charge(rt_F, theta_F); };
        fill_fp(rhopm_Func, fp_F2D_s_k);
        fill_rhop(rhopm_Func, rhopm_F1D_k);
        fill_rhop(rhopc_Func, rhopc_F1D_k);
        fill_ft(rhotm_Func, ft_F3D_lambda_s_k);
        fill_rhot(rhotm_Func, rhotm_F2D_lambda_k);
        fill_rhot(rhotc_Func, rhotc_F2D_lambda_k);
    }

    /**
     * @brief  Replace the interaction kernels.
     * @math   𝒱 → 𝒱′
     * @output Updated interaction kernels.
     */
    virtual void update_potentials(const VnnFunctions& vnn_functions_);

    /**
     * @brief  Recompute all projectile tables.
     * @math   ρ_p → {ρ̃_p, ρ̃_{p,c}, f̃_p}
     * @output Updated projectile tables.
     */
    virtual void update_density_projec(const DensityInfoProjec& density_projec_);

    /**
     * @brief  Recompute all target tables.
     * @math   ρ_t → {ρ̃_{tλ}, ρ̃_{tλ,c}, f̃_{tλ}}
     * @output Updated target tables.
     */
    virtual void update_density_target(const DensityInfoTarget& density_target_);

    /**
     * @brief  Map an even multipole to storage.
     * @math   i_λ = λ/2
     * @output Multipole storage index.
     */
    int lambda_index(int lambda_I) const;

    /**
     * @brief  Evaluate a generic nuclear multipole.
     * @math   U_λ = ℱ_λ⁻¹[ρ̃_{p,m}ρ̃_{tλ,m}ṽ]
     * @output Nuclear-potential multipole, MeV.
     */
    virtual double Unn_lambda(double r_F, int lambda_I, const Real2RealFunc& vnn_k_Func) const;

    /**
     * @brief  Evaluate the direct nuclear multipole.
     * @math   U_λ^D = ℱ_λ⁻¹[ρ̃_{p,m}ρ̃_{tλ,m}ṽ_Dg(E)]
     * @output Direct-potential multipole, MeV.
     */
    double Ud_lambda(double r_F, int lambda_I) const;

    /**
     * @brief  Evaluate the zero-range exchange multipole.
     * @math   U_λ^{E,ZR} = ℱ_λ⁻¹[ρ̃_{p,m}ρ̃_{tλ,m}J₀₀g(E)]
     * @output Zero-range exchange multipole, MeV.
     */
    double Uexzr_lambda(double r_F, int lambda_I) const;

    /**
     * @brief  Evaluate the Coulomb multipole.
     * @math   U_λ^C = ℱ_λ⁻¹[ρ̃_{p,c}ρ̃_{tλ,c}ṽ_C]
     * @output Coulomb-potential multipole, MeV.
     */
    double Ucoul_lambda(double r_F, int lambda_I) const;

    /**
     * @brief  Iterate the finite-range exchange multipoles.
     * @math   U^E ↔ k_eff² = (E - U^C - U^D - U^E)/h_μ
     * @output Exchange-potential multipoles, MeV.
     * @note   Stops at 10⁻³ MeV or 100 iterations.
     */
    Eigen::VectorXd Uex(double r_F) const;

protected:
    int Nr_I = 201; // N_{r_p} = N_{r_t}.
    int Nk_I = 201; // Momentum-grid size.
    int Ns_I = 201; // Exchange-separation-grid size.
    Eigen::VectorXd rp_F1D_rp = Eigen::VectorXd::LinSpaced(Nr_I, 1.0e-10, 50.0); // r_p, fm.
    Eigen::VectorXd rt_F1D_rt = Eigen::VectorXd::LinSpaced(Nr_I, 1.0e-10, 50.0); // r_t, fm.
    Eigen::VectorXd k_F1D_k = Eigen::VectorXd::LinSpaced(Nk_I, 1.0e-10, 5.0); // k, fm⁻¹.
    Eigen::VectorXd s_F1D_s = Eigen::VectorXd::LinSpaced(Ns_I, 1.0e-10, 5.0); // s, fm.
    Eigen::VectorXd rhopm_F1D_k = Eigen::VectorXd::Zero(Nk_I); // ρ̃_{p,m}(k).
    Eigen::VectorXd rhopc_F1D_k = Eigen::VectorXd::Zero(Nk_I); // ρ̃_{p,c}(k).
    Eigen::MatrixXd rhotm_F2D_lambda_k = Eigen::MatrixXd::Zero(Nlambda_I, Nk_I); // ρ̃_{tλ,m}(k).
    Eigen::MatrixXd rhotc_F2D_lambda_k = Eigen::MatrixXd::Zero(Nlambda_I, Nk_I); // ρ̃_{tλ,c}(k).
    Eigen::MatrixXd fp_F2D_s_k = Eigen::MatrixXd::Zero(Ns_I, Nk_I); // f̃_p(s, k).
    Eigen::Tensor<double, 3, Eigen::ColMajor> ft_F3D_lambda_s_k = Eigen::Tensor<double, 3, Eigen::ColMajor>(Nlambda_I, Ns_I, Nk_I); // f̃_{tλ}(s, k).
    int Ntheta_I = 10; // Gauss-Legendre angular size.
    Eigen::VectorXd w_F1D_theta = Eigen::VectorXd::Zero(Ntheta_I); // w_θ.
    Eigen::MatrixXd Y_F2D_lambda_theta = Eigen::MatrixXd::Zero(Nlambda_I, Ntheta_I); // Y_λ0(θ).

    /**
     * @brief  Project a target density by spherical quadrature.
     * @math   ρ_{tλ}(r_t) = ∫Y*_{λ0}(Ω_t)ρ_t(r_t, Ω_t)dΩ_t
     * @output Target-density multipole.
     */
    double calc_rhot_lambda(const RealReal2RealFunc& rhot_Func, double rt_F, int lambda_I) const;

    /**
     * @brief  Project a target exchange factor.
     * @math   f_{tλ}(r_t,s) = ∫Y*_{λ0}ρ_tĵ₁(k_{F,t}s)dΩ_t
     * @output Target exchange-factor multipole.
     */
    double calc_ft_lambda(const RealReal2RealFunc& rhot_Func, double rt_F, double s_F, int lambda_I) const;

    /**
     * @brief  Evaluate the projectile exchange factor.
     * @math   f_p(r_p,s) = ρ_p(r_p)ĵ₁(k_{F,p}s)
     * @output Projectile exchange factor.
     */
    double calc_fp(const Real2RealFunc& rhop_Func, double rp_F, double s_F) const;

    /**
     * @brief  Transform projectile exchange factors.
     * @math   f_p(r_p,s) → f̃_p(s,k)
     * @output Updated projectile exchange table.
     */
    void fill_fp(const Real2RealFunc& rhop_Func, Eigen::MatrixXd& fp_F2D_s_k_);

    /**
     * @brief  Transform a projectile density.
     * @math   ρ_p(r_p) → ρ̃_p(k)
     * @output Updated projectile-density table.
     */
    void fill_rhop(const Real2RealFunc& rhop_Func, Eigen::VectorXd& rhop_F1D_k_);

    /**
     * @brief  Transform target exchange multipoles.
     * @math   f_{tλ}(r_t,s) → f̃_{tλ}(s,k)
     * @output Updated target exchange table.
     */
    void fill_ft(const RealReal2RealFunc& rhot_Func, Eigen::Tensor<double, 3, Eigen::ColMajor>& ft_F3D_lambda_s_k_);

    /**
     * @brief  Transform target-density multipoles.
     * @math   ρ_{tλ}(r_t) → ρ̃_{tλ}(k)
     * @output Updated target-density table.
     */
    void fill_rhot(const RealReal2RealFunc& rhot_Func, Eigen::MatrixXd& rhot_F2D_lambda_k_);

    /**
     * @brief  Transform one direct kernel.
     * @math   U_λ(r) = ℱ_λ⁻¹[ρ̃_pρ̃_{tλ}ṽ]
     * @output Potential multipole.
     */
    double calc_Ud_lambda(double r_F, int lambda_I, const Eigen::VectorXd& rhop_F1D_k_, const Eigen::MatrixXd& rhot_F2D_lambda_k_, const Real2RealFunc& v_k_Func) const;

    /**
     * @brief  Transform one exchange-overlap kernel.
     * @math   G_λ(r,s) = ℱ_λ⁻¹[f̃_p(s,k)f̃_{tλ}(s,k)]
     * @output Exchange-overlap values over s.
     */
    Eigen::VectorXd calc_G_lambda(double r_F, int lambda_I, const Eigen::MatrixXd& fp_F2D_s_k_, const Eigen::Tensor<double, 3, Eigen::ColMajor>& ft_F3D_lambda_s_k_) const;

    /**
     * @brief  Evaluate cached exchange-overlap values.
     * @math   G_λ(r,s) = G_λ[ρ_p,ρ_t]
     * @output Exchange-overlap values over s.
     */
    virtual Eigen::VectorXd G_lambda(double r_F, int lambda_I) const;

    /**
     * @brief  Perform one exchange fixed-point iteration.
     * @math   U_n^E → k_eff → U_{n+1}^E
     * @output L¹ multipole update, MeV.
     */
    double iterate_Uex(double r_F, const Eigen::VectorXd& Uc_F1D_lambda, const Eigen::VectorXd& Ud_F1D_lambda, Eigen::VectorXd& Ue_F1D_lambda) const;
};

inline double Bifold::calc_rhot_lambda(const RealReal2RealFunc& rhot_Func, double rt_F, int lambda_I) const {
    auto rhot_Omega_Func = [&](double theta_F, double) { return rhot_Func(rt_F, theta_F); }; // ρ_t(r_t, θ_t).
    return spherical_Omega2lm(rhot_Omega_Func, lambda_I); // ρ_{tλ}(r_t).
}

inline double Bifold::calc_ft_lambda(const RealReal2RealFunc& rhot_Func, double rt_F, double s_F, int lambda_I) const {
    // ρ_t(r_t, Ω_t) → ρ_tĵ₁(k_{F,t}s) → f_{tλ}(r_t,s).
    auto ft_Omega_Func = [&](double theta_F, double) {
        double rhot_F = rhot_Func(rt_F, theta_F);
        double kFt_F = calc_fermik(rhot_Func, rt_F, theta_F, Cs_F);
        return rhot_F * calc_hatj1(kFt_F * s_F);
    };
    return spherical_Omega2lm(ft_Omega_Func, lambda_I);
}

inline double Bifold::calc_fp(const Real2RealFunc& rhop_Func, double rp_F, double s_F) const {
    double rhop_F = rhop_Func(rp_F); // ρ_p(r_p).
    double kFp_F = calc_fermik(rhop_Func, rp_F, Cs_F); // k_{F,p}(r_p).
    return rhop_F * calc_hatj1(kFp_F * s_F); // f_p(r_p,s).
}

inline void Bifold::update_potentials(const VnnFunctions& vnn_functions_) {
    vnn_functions = vnn_functions_; // 𝒱 → 𝒱′.
}

inline void Bifold::fill_fp(const Real2RealFunc& rhop_Func, Eigen::MatrixXd& fp_F2D_s_k_) {
    // f_p(r_p,s) → f̃_p(s,k).
    Eigen::MatrixXd fp_F2D_rp_s(Nr_I, Ns_I);
    for (int rp_I = 0; rp_I < Nr_I; ++rp_I) {
        for (int s_I = 0; s_I < Ns_I; ++s_I) {
            fp_F2D_rp_s(rp_I, s_I) = calc_fp(rhop_Func, rp_F1D_rp(rp_I), s_F1D_s(s_I));
        }
    }
    for (int k_I = 0; k_I < Nk_I; ++k_I) {
        for (int s_I = 0; s_I < Ns_I; ++s_I) {
            fp_F2D_s_k_(s_I, k_I) = spherical_FT_r2k(fp_F2D_rp_s.col(s_I), rp_F1D_rp, k_F1D_k(k_I), 0);
        }
    }
}

inline void Bifold::fill_rhop(const Real2RealFunc& rhop_Func, Eigen::VectorXd& rhop_F1D_k_) {
    // ρ_p(r_p) → ρ̃_p(k).
    Eigen::VectorXd rhop_F1D_rp(Nr_I);
    for (int rp_I = 0; rp_I < Nr_I; ++rp_I) {
        rhop_F1D_rp(rp_I) = rhop_Func(rp_F1D_rp(rp_I));
    }
    for (int k_I = 0; k_I < Nk_I; ++k_I) {
        rhop_F1D_k_(k_I) = spherical_FT_r2k(rhop_F1D_rp, rp_F1D_rp, k_F1D_k(k_I), 0);
    }
}

inline void Bifold::update_density_projec(const DensityInfoProjec& density_projec_) {
    // ρ_p → {r_p, ρ̃_{p,m}, ρ̃_{p,c}, f̃_p}.
    density_projec = density_projec_;
    rp_F1D_rp = Eigen::VectorXd::LinSpaced(Nr_I, 1.0e-10, density_projec.rmax_F);
    Real2RealFunc rhopm_Func = [&](double rp_F) { return density_projec.rho_matter(rp_F); };
    Real2RealFunc rhopc_Func = [&](double rp_F) { return density_projec.rho_charge(rp_F); };
    fill_fp(rhopm_Func, fp_F2D_s_k);
    fill_rhop(rhopm_Func, rhopm_F1D_k);
    fill_rhop(rhopc_Func, rhopc_F1D_k);
}

inline void Bifold::fill_ft(const RealReal2RealFunc& rhot_Func, Eigen::Tensor<double, 3, Eigen::ColMajor>& ft_F3D_lambda_s_k_) {
    // f_{tλ}(r_t,s) → f̃_{tλ}(s,k).
    Eigen::MatrixXd ft_F2D_rt_s(Nr_I, Ns_I);
    for (int lambda_index_I = 0; lambda_index_I < Nlambda_I; ++lambda_index_I) {
        int lambda_I = lambda_I1D_lambda(lambda_index_I);
        for (int rt_I = 0; rt_I < Nr_I; ++rt_I) {
            for (int s_I = 0; s_I < Ns_I; ++s_I) {
                ft_F2D_rt_s(rt_I, s_I) = calc_ft_lambda(rhot_Func, rt_F1D_rt(rt_I), s_F1D_s(s_I), lambda_I);
            }
        }
        for (int k_I = 0; k_I < Nk_I; ++k_I) {
            for (int s_I = 0; s_I < Ns_I; ++s_I) {
                ft_F3D_lambda_s_k_(lambda_index_I, s_I, k_I) = spherical_FT_r2k(ft_F2D_rt_s.col(s_I), rt_F1D_rt, k_F1D_k(k_I), lambda_I);
            }
        }
    }
}

inline void Bifold::fill_rhot(const RealReal2RealFunc& rhot_Func, Eigen::MatrixXd& rhot_F2D_lambda_k_) {
    // ρ_{tλ}(r_t) → ρ̃_{tλ}(k).
    Eigen::VectorXd rhot_F1D_rt(Nr_I);
    for (int lambda_index_I = 0; lambda_index_I < Nlambda_I; ++lambda_index_I) {
        int lambda_I = lambda_I1D_lambda(lambda_index_I);
        for (int rt_I = 0; rt_I < Nr_I; ++rt_I) {
            rhot_F1D_rt(rt_I) = calc_rhot_lambda(rhot_Func, rt_F1D_rt(rt_I), lambda_I);
        }
        for (int k_I = 0; k_I < Nk_I; ++k_I) {
            rhot_F2D_lambda_k_(lambda_index_I, k_I) = spherical_FT_r2k(rhot_F1D_rt, rt_F1D_rt, k_F1D_k(k_I), lambda_I);
        }
    }
}

inline void Bifold::update_density_target(const DensityInfoTarget& density_target_) {
    // ρ_t → {r_t, ρ̃_{tλ,m}, ρ̃_{tλ,c}, f̃_{tλ}}.
    density_target = density_target_;
    rt_F1D_rt = Eigen::VectorXd::LinSpaced(Nr_I, 1.0e-10, density_target.rmax_F);
    RealReal2RealFunc rhotm_Func = [&](double rt_F, double theta_F) { return density_target.rho_matter(rt_F, theta_F); };
    RealReal2RealFunc rhotc_Func = [&](double rt_F, double theta_F) { return density_target.rho_charge(rt_F, theta_F); };
    fill_ft(rhotm_Func, ft_F3D_lambda_s_k);
    fill_rhot(rhotm_Func, rhotm_F2D_lambda_k);
    fill_rhot(rhotc_Func, rhotc_F2D_lambda_k);
}

inline int Bifold::lambda_index(int lambda_I) const {
    assert(lambda_I >= 0);
    assert(lambda_I % 2 == 0);
    int lambda_index_I = lambda_I / 2;
    assert(lambda_index_I < Nlambda_I);
    assert(lambda_I1D_lambda(lambda_index_I) == lambda_I);
    return lambda_index_I;
}

inline double Bifold::calc_Ud_lambda(double r_F, int lambda_I, const Eigen::VectorXd& rhop_F1D_k_, const Eigen::MatrixXd& rhot_F2D_lambda_k_, const Real2RealFunc& v_k_Func) const {
    // ρ̃_pρ̃_{tλ}ṽ → U_λ(r).
    int lambda_index_I = lambda_index(lambda_I);
    Eigen::VectorXd integrand_F1D_k(Nk_I);
    for (int k_I = 0; k_I < Nk_I; ++k_I) {
        integrand_F1D_k(k_I) = rhop_F1D_k_(k_I) * rhot_F2D_lambda_k_(lambda_index_I, k_I) * v_k_Func(k_F1D_k(k_I));
    }
    return spherical_FT_k2r(integrand_F1D_k, k_F1D_k, r_F, lambda_I);
}

inline double Bifold::Unn_lambda(double r_F, int lambda_I, const Real2RealFunc& vnn_k_Func) const {
    return calc_Ud_lambda(r_F, lambda_I, rhopm_F1D_k, rhotm_F2D_lambda_k, vnn_k_Func); // U_λ[ρ_{p,m},ρ_{t,m},v].
}

inline double Bifold::Ud_lambda(double r_F, int lambda_I) const {
    auto vnnd_k_Func = [&](double k_F) { return vnn_functions.vnnd_k_Func(k_F) * vnn_functions.gE_Func(Ecm_F, density_projec.Ap_I); }; // ṽ_D(k)g(E).
    return Unn_lambda(r_F, lambda_I, vnnd_k_Func); // U_λ^D(r).
}

inline double Bifold::Uexzr_lambda(double r_F, int lambda_I) const {
    auto vnne_zr_k_Func = [&](double) { return vnn_functions.vnne_zr_Func(Ecm_F, density_projec.Ap_I) * vnn_functions.gE_Func(Ecm_F, density_projec.Ap_I); }; // J₀₀(E)g(E).
    return Unn_lambda(r_F, lambda_I, vnne_zr_k_Func); // U_λ^{E,ZR}(r).
}

inline double Bifold::Ucoul_lambda(double r_F, int lambda_I) const {
    auto vcoul_k_Func = [&](double k_F) { return vnn_functions.vcoul_k_Func(k_F, Ze2_F); }; // ṽ_C(k; Ze²).
    return calc_Ud_lambda(r_F, lambda_I, rhopc_F1D_k, rhotc_F2D_lambda_k, vcoul_k_Func); // U_λ^C(r).
}

inline Eigen::VectorXd Bifold::calc_G_lambda(double r_F, int lambda_I, const Eigen::MatrixXd& fp_F2D_s_k_, const Eigen::Tensor<double, 3, Eigen::ColMajor>& ft_F3D_lambda_s_k_) const {
    // f̃_p(s,k)f̃_{tλ}(s,k) → G_λ(r,s).
    Eigen::VectorXd G_F1D_s = Eigen::VectorXd::Zero(Ns_I);
    Eigen::VectorXd integrand_F1D_k = Eigen::VectorXd::Zero(Nk_I);
    int lambda_index_I = lambda_index(lambda_I);
    for (int s_I = 0; s_I < Ns_I; ++s_I) {
        for (int k_I = 0; k_I < Nk_I; ++k_I) {
            integrand_F1D_k(k_I) = fp_F2D_s_k_(s_I, k_I) * ft_F3D_lambda_s_k_(lambda_index_I, s_I, k_I);
        }
        G_F1D_s(s_I) = spherical_FT_k2r(integrand_F1D_k, k_F1D_k, r_F, lambda_I);
    }
    return G_F1D_s;
}

inline Eigen::VectorXd Bifold::G_lambda(double r_F, int lambda_I) const {
    return calc_G_lambda(r_F, lambda_I, fp_F2D_s_k, ft_F3D_lambda_s_k); // G_λ[ρ_p,ρ_t].
}

inline double Bifold::iterate_Uex(double r_F, const Eigen::VectorXd& Uc_F1D_lambda, const Eigen::VectorXd& Ud_F1D_lambda, Eigen::VectorXd& Ue_F1D_lambda) const {
    assert(Uc_F1D_lambda.size() == Nlambda_I);
    assert(Ud_F1D_lambda.size() == Nlambda_I);
    assert(Ue_F1D_lambda.size() == Nlambda_I);

    // {U_λ^C,U_λ^D,U_λ^E} → k_eff²(θ).
    double Y00_F = 1.0 / std::sqrt(4.0 * std::numbers::pi);
    double Ecm_over_Y00_F = Ecm_F / Y00_F;
    Eigen::VectorXd keff2_F1D_lambda = -(Uc_F1D_lambda + Ud_F1D_lambda + Ue_F1D_lambda);
    keff2_F1D_lambda(0) += Ecm_over_Y00_F;
    keff2_F1D_lambda /= hmass_F;
    Eigen::VectorXd keff2_F1D_theta(Ntheta_I);
    keff2_F1D_theta.noalias() = Y_F2D_lambda_theta.transpose() * keff2_F1D_lambda;
    Eigen::VectorXd keff_F1D_theta = keff2_F1D_theta.cwiseMax(0.0).cwiseSqrt();

    // G_λ(r,s) → G(r,θ,s).
    Eigen::MatrixXd G_F2D_lambda_s(Nlambda_I, Ns_I);
    for (int lambda_index_I = 0; lambda_index_I < Nlambda_I; ++lambda_index_I) {
        G_F2D_lambda_s.row(lambda_index_I) = G_lambda(r_F, lambda_I1D_lambda(lambda_index_I)).transpose();
    }
    Eigen::MatrixXd G_F2D_theta_s(Ntheta_I, Ns_I);
    G_F2D_theta_s.noalias() = Y_F2D_lambda_theta.transpose() * G_F2D_lambda_s;

    // G(r,θ,s)v_E(s)g(E) → U^E(r,θ).
    auto vnne_r_Func = [&](double s_F) { return vnn_functions.vnne_r_Func(s_F) * vnn_functions.gE_Func(Ecm_F, density_projec.Ap_I); };
    for (int s_I = 0; s_I < Ns_I; ++s_I) {
        G_F2D_theta_s.col(s_I) *= vnne_r_Func(s_F1D_s(s_I));
    }
    Eigen::VectorXd Uex_F1D_theta = Eigen::VectorXd::Zero(Ntheta_I);
    for (int theta_I = 0; theta_I < Ntheta_I; ++theta_I) {
        Uex_F1D_theta(theta_I) = spherical_FT_r2k(G_F2D_theta_s.row(theta_I), s_F1D_s, keff_F1D_theta(theta_I), 0);
    }

    // U^E(r,θ) → U_λ^E(r) → ||ΔU^E||₁.
    Eigen::VectorXd integrand_F1D_theta = Uex_F1D_theta.array() * w_F1D_theta.array();
    Eigen::VectorXd Uex_F1D_lambda(Nlambda_I);
    Uex_F1D_lambda.noalias() = 2.0 * std::numbers::pi * Y_F2D_lambda_theta * integrand_F1D_theta;
    Eigen::VectorXd delta_Uex_F1D_lambda = Uex_F1D_lambda - Ue_F1D_lambda;
    double error_F = delta_Uex_F1D_lambda.cwiseAbs().sum();
    Ue_F1D_lambda = Uex_F1D_lambda;
    return error_F;
}

inline Eigen::VectorXd Bifold::Uex(double r_F) const {
    assert(r_F >= 0.0);

    // {U_λ^C,U_λ^D,U_λ^{E,ZR}} → fixed-point U_λ^E.
    double tolerance_F = 1.0e-3; // ||ΔU^E||₁, MeV.
    int Niter_I = 100; // Maximum fixed-point iterations.
    Eigen::VectorXd Uc_F1D_lambda(Nlambda_I);
    Eigen::VectorXd Ud_F1D_lambda(Nlambda_I);
    Eigen::VectorXd Ue_F1D_lambda(Nlambda_I);
    for (int lambda_index_I = 0; lambda_index_I < Nlambda_I; ++lambda_index_I) {
        int lambda_I = lambda_I1D_lambda(lambda_index_I);
        Uc_F1D_lambda(lambda_index_I) = Ucoul_lambda(r_F, lambda_I);
        Ud_F1D_lambda(lambda_index_I) = Ud_lambda(r_F, lambda_I);
        Ue_F1D_lambda(lambda_index_I) = Uexzr_lambda(r_F, lambda_I);
    }
    for (int iter_I = 0; iter_I < Niter_I; ++iter_I) {
        double error_F = iterate_Uex(r_F, Uc_F1D_lambda, Ud_F1D_lambda, Ue_F1D_lambda);
        if (error_F < tolerance_F) { break; }
    }
    return Ue_F1D_lambda;
}
