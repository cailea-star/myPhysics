/**
 * @file    spherical_basis.hpp
 * @author  cailea
 * @date    2026-09-14
 * @brief   Spherical radial basis on a generalized Gauss-Laguerre mesh.
 */

#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>
#include <Eigen/Core>
#include <gsl/gsl_sf_laguerre.h>

#include "spherical_setting.hpp"
#include "integration_gauss.hpp"

/** @brief Spherical radial basis on a generalized Gauss-Laguerre mesh. */
class SphericalLaguerreBasis {
public:
    double b_F = 0.0;                         // b = √[ℏ/(Mω)].
    std::vector<SphericalSPLabel> labels_S1D_sp{}; // α_sp = (n,l,j,m)_sp.
    Eigen::VectorXd r_F1D_r{};                // r_i = b√η_i.
    Eigen::VectorXd w_F1D_r{};                // w_r,i = b³ w_i^GL e^η_i / 2.
    Eigen::MatrixXd phi_F2D_sp_r{};           // φ_sp,i = φ_nl(r_i).
    Eigen::MatrixXd dphi_F2D_sp_r{};          // ∂_rφ_nl(r_i).
    Eigen::MatrixXd ddphi_F2D_sp_r{};         // ∂_r²φ_nl(r_i).

    /**
     * @brief  Generate radial functions using generalized Gauss-Laguerre quadrature.
     * @math   η = (r/b)², α_poly = l+1/2, α_weight = 1/2
     * @output Nodes, weights, basis functions, and first two derivatives.
     * @note   Requires b>0 and N_r>0.
     */
    SphericalLaguerreBasis(double b_F_, int Nr_I_, const std::vector<SphericalSPLabel>& labels_S1D_sp_) {
        b_F = b_F_;
        labels_S1D_sp = labels_S1D_sp_;
        assert(std::isfinite(b_F) && b_F > 0.0);
        assert(Nr_I_ > 0);

        // (N_sp,N_r) → radial arrays.
        r_F1D_r.resize(Nr_I_);
        w_F1D_r.resize(Nr_I_);
        phi_F2D_sp_r.resize(labels_S1D_sp.size(), Nr_I_);
        dphi_F2D_sp_r.resize(labels_S1D_sp.size(), Nr_I_);
        ddphi_F2D_sp_r.resize(labels_S1D_sp.size(), Nr_I_);

        // {α_sp} → {(n_sp,l_sp+1/2)}.
        Eigen::VectorXi n_I1D_sp(labels_S1D_sp.size());
        Eigen::VectorXd alpha_F1D_sp(labels_S1D_sp.size());
        for (int sp_I = 0; sp_I < static_cast<int>(labels_S1D_sp.size()); ++sp_I) {
            n_I1D_sp(sp_I) = labels_S1D_sp[sp_I].n_I;
            alpha_F1D_sp(sp_I) = labels_S1D_sp[sp_I].l_I + 0.5;
        }

        // α_poly=l+1/2, α_weight=1/2 → grid → basis.
        GaussLaguerreMeshes gl_meshes(n_I1D_sp, alpha_F1D_sp, Nr_I_, 0.5);
        fill_grid(gl_meshes);
    }

protected:
    /**
     * @brief  Evaluate radial functions and derivatives using Laguerre identities.
     * @math   φ_nl(r) = √[2n!/Γ(n+l+3/2)] L_n^(l+1/2)(η)η^(l/2)e^(-η/2)/b^(3/2)
     * @output Filled radial mesh and basis arrays.
     * @note   Requires α_weight=1/2 and α_poly,sp=l_sp+1/2.
     * @note   Physical weights represent r²dr.
     */
    void fill_grid(const GaussLaguerreMeshes& gl_meshes);
};

using SphericalBasis1D = SphericalLaguerreBasis;

inline void SphericalLaguerreBasis::fill_grid(const GaussLaguerreMeshes& gl_meshes) {
    assert(gl_meshes.weightAlpha_F == 0.5);
    assert(gl_meshes.n_I1D_na.size() == static_cast<int>(labels_S1D_sp.size()));
    assert((gl_meshes.x_F1D_x.array() > 0.0).all());

    // η=(r/b)²; r²dr=(b³/2)η^(1/2)dη.
    const double jacobian_F = 0.5 * b_F * b_F * b_F;
    const Eigen::VectorXd& eta_F1D_r = gl_meshes.x_F1D_x;
    r_F1D_r = b_F * eta_F1D_r.array().sqrt().matrix();
    w_F1D_r = jacobian_F * gl_meshes.w_F1D_x;

    // (n,l,N_n,L_n,{η_i}) → work arrays.
    const int Nr_I = static_cast<int>(gl_meshes.x_F1D_x.size());
    const int Nsp_I = static_cast<int>(gl_meshes.n_I1D_na.size());
    const double norm_F = std::sqrt(2.0) / std::pow(b_F, 1.5);
    const Eigen::VectorXd n_F1D_sp = gl_meshes.n_I1D_na.cast<double>();
    const Eigen::VectorXd l_F1D_sp = gl_meshes.alpha_F1D_na.array() - 0.5;
    const Eigen::MatrixXd L_F2D_sp_r = gl_meshes.L_F2D_na_x;
    const Eigen::VectorXd N_F1D_sp = norm_F * gl_meshes.N_F1D_na;

    // {η^(l/2),(η-l)/η,[(η-l)²-2l]/(4η²)}.
    const Eigen::VectorXd I_F1D_sp = Eigen::VectorXd::Ones(Nsp_I);
    Eigen::MatrixXd etaPl_F2D_sp_r(Nsp_I, Nr_I);
    Eigen::MatrixXd etaMl_F2D_sp_r(Nsp_I, Nr_I);
    Eigen::MatrixXd etaMl2_F2D_sp_r(Nsp_I, Nr_I);
    for (int r_I = 0; r_I < Nr_I; ++r_I) {etaPl_F2D_sp_r.col(r_I) = (I_F1D_sp * eta_F1D_r(r_I)).array().pow(0.5 * l_F1D_sp.array());}
    for (int r_I = 0; r_I < Nr_I; ++r_I) {etaMl_F2D_sp_r.col(r_I) = (eta_F1D_r(r_I) - l_F1D_sp.array()) / eta_F1D_r(r_I);}
    for (int r_I = 0; r_I < Nr_I; ++r_I) {etaMl2_F2D_sp_r.col(r_I) = 0.25 * etaMl_F2D_sp_r.col(r_I).array().square() - 0.5 * l_F1D_sp.array() / (eta_F1D_r(r_I) * eta_F1D_r(r_I));}

    // ∂_ηL_n^(l+1/2) = -L_(n-1)^(l+3/2).
    Eigen::MatrixXd dL_F2D_sp_r;
    dL_F2D_sp_r.resize(Nsp_I, Nr_I);
    dL_F2D_sp_r.setZero();
    for (int sp_I = 0; sp_I < Nsp_I; ++sp_I) {
        assert(gl_meshes.n_I1D_na(sp_I) == labels_S1D_sp[sp_I].n_I);
        assert(gl_meshes.alpha_F1D_na(sp_I) == labels_S1D_sp[sp_I].l_I + 0.5);
        const int nm1_I = std::max(gl_meshes.n_I1D_na(sp_I) - 1, 0);
        const double hasLowerOrder_F = static_cast<double>(gl_meshes.n_I1D_na(sp_I) > 0);
        const double alphap1_F = gl_meshes.alpha_F1D_na(sp_I) + 1.0;
        for (int r_I = 0; r_I < Nr_I; ++r_I) {
            const double eta_F = eta_F1D_r(r_I);
            dL_F2D_sp_r(sp_I, r_I) = -hasLowerOrder_F * gsl_sf_laguerre_n(nm1_I, alphap1_F, eta_F);
        }
    }

    // η∂²_ηL = (η-l-3/2)∂_ηL - nL.
    Eigen::MatrixXd ddL_F2D_sp_r = dL_F2D_sp_r * eta_F1D_r.asDiagonal() - (l_F1D_sp.array() + 1.5).matrix().asDiagonal() * dL_F2D_sp_r;
    ddL_F2D_sp_r.noalias() -= n_F1D_sp.asDiagonal() * L_F2D_sp_r;
    ddL_F2D_sp_r.array().rowwise() /= eta_F1D_r.transpose().array();

    // (L,L',L'',η,l) → (u'/w,u''/w), w=η^(l/2)e^(-η/2).
    const Eigen::MatrixXd dLterm_F2D_sp_r = dL_F2D_sp_r.array() - 0.5 * L_F2D_sp_r.array() * etaMl_F2D_sp_r.array();
    const Eigen::MatrixXd ddLterm_F2D_sp_r = ddL_F2D_sp_r.array() - dL_F2D_sp_r.array() * etaMl_F2D_sp_r.array() + L_F2D_sp_r.array() * etaMl2_F2D_sp_r.array();

    // η=(r/b)² → (∂_rη,∂²_rη)=(2r/b²,2/b²).
    const Eigen::VectorXd deta_F1D_r = 2.0 * r_F1D_r / (b_F * b_F);
    const double ddeta_F = 2.0 / (b_F * b_F);

    // NLη^(l/2) = e^(η/2)φ.
    phi_F2D_sp_r.noalias() = N_F1D_sp.asDiagonal() * L_F2D_sp_r;
    phi_F2D_sp_r.array() *= etaPl_F2D_sp_r.array();

    // e^(η/2)∂_ηφ → e^(η/2)∂_rφ.
    const Eigen::MatrixXd dphiTerm_F2D_sp_r = (N_F1D_sp.asDiagonal() * dLterm_F2D_sp_r).array() * etaPl_F2D_sp_r.array();
    dphi_F2D_sp_r.noalias() = dphiTerm_F2D_sp_r * deta_F1D_r.asDiagonal();

    // e^(η/2)∂²_ηφ → e^(η/2)∂²_rφ.
    const Eigen::MatrixXd ddphiTerm_F2D_sp_r = (N_F1D_sp.asDiagonal() * ddLterm_F2D_sp_r).array() * etaPl_F2D_sp_r.array();
    ddphi_F2D_sp_r.noalias() = ddphiTerm_F2D_sp_r * deta_F1D_r.array().square().matrix().asDiagonal();
    ddphi_F2D_sp_r += ddeta_F * dphiTerm_F2D_sp_r;

    // {φ,∂φ,∂²φ} × e^(-η/2); w ÷ e^-η.
    const Eigen::VectorXd expHalf_F1D_r = (-0.5 * eta_F1D_r.array()).exp();
    w_F1D_r.array() /= expHalf_F1D_r.array().square();
    phi_F2D_sp_r.array().rowwise() *= expHalf_F1D_r.transpose().array();
    dphi_F2D_sp_r.array().rowwise() *= expHalf_F1D_r.transpose().array();
    ddphi_F2D_sp_r.array().rowwise() *= expHalf_F1D_r.transpose().array();
}
