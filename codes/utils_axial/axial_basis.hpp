/**
 * @file    axial_basis.hpp
 * @author  cailea
 * @date    2026-05-02
 * @brief   Axial harmonic-oscillator basis functions, derivatives, and quadrature weights.
 */

#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>
#include <Eigen/Core>
#include <unsupported/Eigen/CXX11/Tensor>
#include <gsl/gsl_sf_hermite.h>
#include <gsl/gsl_sf_laguerre.h>

#include "axial_config.hpp"
#include "integration_gauss.hpp"

/** @brief Axial Hermite basis on a Gauss-Hermite mesh. */
class AxialHermiteBasis {
public:
    double bz_F = 0.0;                          // b_z = √[ℏ / (mω_z)].
    std::vector<AxialSPLabel> labels_S1D_sp{};  // α_sp = (n_z,n_r,Λ,Ω,Σ,π)_sp.
    Eigen::VectorXd z_F1D_z{};                  // z_i = b_zζ_i.
    Eigen::VectorXd w_F1D_z{};                  // w_z,i = b_z w_i^GH e^(ζ_i²).
    Eigen::VectorXd zeta_F1D_z{};               // ζ_i: Gauss-Hermite nodes.
    Eigen::MatrixXd phi_F2D_sp_z{};             // φ_sp,i = φ_nz,sp(z_i).
    Eigen::MatrixXd dphi_F2D_sp_z{};            // dφ_sp,i = ∂_zφ_nz,sp(z_i).
    Eigen::MatrixXd ddphi_F2D_sp_z{};           // ddφ_sp,i = ∂_z²φ_nz,sp(z_i).

    /**
     * @brief  Generate an axial Hermite basis by Gauss-Hermite quadrature.
     * @math   ζ = z / b_z
     * @output Nodes, weights, basis functions, and first two derivatives.
     * @note   useHalf_B retains z > 0 and doubles the weights.
     */
    AxialHermiteBasis(double bz_F_, int Nz_I_, const std::vector<AxialSPLabel>& labels_S1D_sp_, bool useHalf_B = false) {
        bz_F = bz_F_;
        labels_S1D_sp = labels_S1D_sp_;
        z_F1D_z.resize(Nz_I_);
        w_F1D_z.resize(Nz_I_);
        zeta_F1D_z.resize(Nz_I_);
        phi_F2D_sp_z.resize(labels_S1D_sp_.size(), Nz_I_);
        dphi_F2D_sp_z.resize(labels_S1D_sp_.size(), Nz_I_);
        ddphi_F2D_sp_z.resize(labels_S1D_sp_.size(), Nz_I_);
        assert(std::isfinite(bz_F_) && bz_F_ > 0.0);
        assert(Nz_I_ > 0);

        // {α_sp} → {n_{z,sp}}.
        Eigen::VectorXi nz_I1D_sp(labels_S1D_sp.size());
        for (int sp_I = 0; sp_I < labels_S1D_sp.size(); ++sp_I) {nz_I1D_sp(sp_I) = labels_S1D_sp[sp_I].nz_I;}

        // (-∞,∞) → (0,∞).
        if (useHalf_B) {
            GaussHermiteMeshes gh_meshes(nz_I1D_sp, 2 * Nz_I_);
            fill_grid(gh_meshes);
            z_F1D_z = z_F1D_z.tail(Nz_I_).eval();
            w_F1D_z = (2.0 * w_F1D_z.tail(Nz_I_)).eval();
            zeta_F1D_z = zeta_F1D_z.tail(Nz_I_).eval();
            phi_F2D_sp_z = phi_F2D_sp_z.rightCols(Nz_I_).eval();
            dphi_F2D_sp_z = dphi_F2D_sp_z.rightCols(Nz_I_).eval();
            ddphi_F2D_sp_z = ddphi_F2D_sp_z.rightCols(Nz_I_).eval();
        } else {
            GaussHermiteMeshes gh_meshes(nz_I1D_sp, Nz_I_);
            fill_grid(gh_meshes);
        }
    }

protected:
    /**
     * @brief  Evaluate axial Hermite basis functions and derivatives.
     * @math   φ_nz(z) = N_nz H_nz(ζ)e^(-ζ²/2) / √b_z
     * @output Filled axial mesh and basis arrays.
     */
    void fill_grid(const GaussHermiteMeshes& gh_meshes);
};

/** @brief Radial Laguerre basis on a Gauss-Laguerre mesh. */
class AxialLaguerreBasis {
public:
    double br_F = 0.0;                          // b_r = √[ℏ / (mω_r)].
    std::vector<AxialSPLabel> labels_S1D_sp{};  // α_sp = (n_z,n_r,Λ,Ω,Σ,π)_sp.
    Eigen::VectorXd r_F1D_r{};                  // r_j = b_r√η_j.
    Eigen::VectorXd w_F1D_r{};                  // w_r,j = b_r² w_j^GL e^(η_j) / 2.
    Eigen::VectorXd eta_F1D_r{};                // η_j: Gauss-Laguerre nodes.
    Eigen::MatrixXd phi_F2D_sp_r{};             // φ_sp,j = φ_nr,sp^Λ_sp(r_j).
    Eigen::MatrixXd dphi_F2D_sp_r{};            // dφ_sp,j = ∂_rφ_nr,sp^Λ_sp(r_j).
    Eigen::MatrixXd ddphi_F2D_sp_r{};           // ddφ_sp,j = ∂_r²φ_nr,sp^Λ_sp(r_j).

    /**
     * @brief  Generate a radial Laguerre basis by Gauss-Laguerre quadrature.
     * @math   η = (r / b_r)²
     * @output Nodes, weights, basis functions, and first two derivatives.
     */
    AxialLaguerreBasis(double br_F_, int Nr_I_, const std::vector<AxialSPLabel>& labels_S1D_sp_) {
        br_F = br_F_;
        labels_S1D_sp = labels_S1D_sp_;
        r_F1D_r.resize(Nr_I_);
        w_F1D_r.resize(Nr_I_);
        eta_F1D_r.resize(Nr_I_);
        phi_F2D_sp_r.resize(labels_S1D_sp_.size(), Nr_I_);
        dphi_F2D_sp_r.resize(labels_S1D_sp_.size(), Nr_I_);
        ddphi_F2D_sp_r.resize(labels_S1D_sp_.size(), Nr_I_);
        assert(std::isfinite(br_F_) && br_F_ > 0.0);
        assert(Nr_I_ > 0);

        // {α_sp} → {(n_{r,sp},Λ_sp)}.
        Eigen::VectorXi nr_I1D_sp(labels_S1D_sp.size());
        Eigen::VectorXd Lambda_F1D_sp(labels_S1D_sp.size());
        for (int sp_I = 0; sp_I < labels_S1D_sp.size(); ++sp_I) {nr_I1D_sp(sp_I) = labels_S1D_sp[sp_I].nr_I;}
        for (int sp_I = 0; sp_I < labels_S1D_sp.size(); ++sp_I) {Lambda_F1D_sp(sp_I) = labels_S1D_sp[sp_I].Lambda_I;}

        // {(n_{r,sp},Λ_sp)} → {η_j,w_j,L_{n_r}^{Λ}(η_j)}.
        GaussLaguerreMeshes gl_meshes(nr_I1D_sp, Lambda_F1D_sp, Nr_I_);
        fill_grid(gl_meshes);
    }

protected:
    /**
     * @brief  Evaluate radial Laguerre basis functions and derivatives.
     * @math   φ_nr^Λ(r) = √2 N_nr^Λ L_nr^Λ(η)η^(Λ/2)e^(-η/2) / b_r
     * @output Filled radial mesh and basis arrays.
     */
    void fill_grid(const GaussLaguerreMeshes& gl_meshes);
};

/** @brief Separable axial harmonic-oscillator basis. */
class AxialBasis {
public:
    using BasisTensor = Eigen::Tensor<double, 3, Eigen::ColMajor>;

    double br_F = 0.0;                          // b_r = √[ℏ / (mω_r)].
    double bz_F = 0.0;                          // b_z = √[ℏ / (mω_z)].
    std::vector<AxialSPLabel> labels_S1D_sp{};  // α_sp = (n_z,n_r,Λ,Ω,Σ,π)_sp.
    Eigen::VectorXi twoSigma_I1D_sp{};          // 2Σ_sp = 2Ω_sp - 2Λ_sp.
    Eigen::VectorXd z_F1D_z{};                  // z_i = b_zζ_i.
    Eigen::VectorXd r_F1D_r{};                  // r_j = b_r√η_j.
    Eigen::VectorXd zeta_F1D_z{};               // ζ_i: Gauss-Hermite nodes.
    Eigen::VectorXd eta_F1D_r{};                // η_j: Gauss-Laguerre nodes.
    Eigen::MatrixXd w_F2D_z_r{};                // w_ij = 2π w_z,i w_r,j.
    BasisTensor phi_F3D_sp_z_r{};               // φ_sp(z_i,r_j) = φ_nz(z_i)φ_nr^Λ(r_j) / √(2π).
    BasisTensor dphidr_F3D_sp_z_r{};            // dφ_sp/dr = φ_nz(z_i)∂_rφ_nr^Λ(r_j) / √(2π).
    BasisTensor dphidz_F3D_sp_z_r{};            // dφ_sp/dz = ∂_zφ_nz(z_i)φ_nr^Λ(r_j) / √(2π).
    BasisTensor ddphidr_F3D_sp_z_r{};           // d²φ_sp/dr² = φ_nz(z_i)∂_r²φ_nr^Λ(r_j) / √(2π).
    BasisTensor ddphidz_F3D_sp_z_r{};           // d²φ_sp/dz² = ∂_z²φ_nz(z_i)φ_nr^Λ(r_j) / √(2π).

    /**
     * @brief  Generate a separable axial harmonic-oscillator basis.
     * @math   φ_sp(z,r) = φ_nz(z)φ_nr^Λ(r) / √(2π)
     * @output Axial basis functions, derivatives, nodes, and weights.
     */
    AxialBasis(const AxialConfig& axialconfig, const std::vector<AxialSPLabel>& labels_S1D_sp_) {
        br_F = axialconfig.br_F;
        bz_F = axialconfig.bz_F;
        labels_S1D_sp = labels_S1D_sp_;
        twoSigma_I1D_sp.resize(labels_S1D_sp_.size());

        // config → (rBasis,zBasis,{2Σ_sp}) → basis.
        AxialLaguerreBasis rBasis(axialconfig.br_F, axialconfig.Nr_I, labels_S1D_sp_);

        AxialHermiteBasis zBasis(axialconfig.bz_F, axialconfig.Nz_I, labels_S1D_sp_, axialconfig.useReflection_B);

        for (int sp_I = 0; sp_I < labels_S1D_sp.size(); ++sp_I) {twoSigma_I1D_sp(sp_I) = labels_S1D_sp[sp_I].twoSigma_I;}

        fill_from_separable_basis(zBasis, rBasis);
    }

private:
    /**
     * @brief  Combine axial and radial basis functions.
     * @math   φ_sp(z,r) = φ_nz(z)φ_nr^Λ(r) / √(2π)
     * @output Filled basis tensors and product quadrature weights.
     */
    void fill_from_separable_basis(const AxialHermiteBasis& zBasis, const AxialLaguerreBasis& rBasis);
};

inline void AxialHermiteBasis::fill_grid(const GaussHermiteMeshes& gh_meshes) {
    // ζ=z/b_z; dz=b_z dζ.
    const double jacobian_F = bz_F;
    zeta_F1D_z = gh_meshes.x_F1D_x;
    z_F1D_z = bz_F * zeta_F1D_z;
    w_F1D_z = jacobian_F * gh_meshes.w_F1D_x;

    // (n,N_n,H_n,{ζ_i}) → work arrays.
    const int Nsp_I = gh_meshes.n_I1D_n.size();
    const int Nz_I = gh_meshes.x_F1D_x.size();
    const double norm_F = 1.0 / std::sqrt(bz_F);
    const Eigen::VectorXd nz_F1D_sp = gh_meshes.n_I1D_n.cast<double>();
    const Eigen::MatrixXd H_F2D_sp_z = gh_meshes.H_F2D_n_x;
    const Eigen::VectorXd N_F1D_sp = norm_F * gh_meshes.N_F1D_n;

    // ∂_ζ H_n = 2nH_{n-1}.
    Eigen::MatrixXd dH_F2D_sp_z = Eigen::MatrixXd::Zero(Nsp_I, Nz_I);
    for (int sp_I = 0; sp_I < Nsp_I; ++sp_I) {
        const int nz_I = gh_meshes.n_I1D_n(sp_I);
        const int nzm1_I = std::max(nz_I - 1, 0);
        const double hasLowerOrder_F = static_cast<double>(nz_I > 0);
        for (int z_I = 0; z_I < Nz_I; ++z_I) {
            const double zeta_F = zeta_F1D_z(z_I);
            dH_F2D_sp_z(sp_I, z_I) = 2.0 * nz_I * hasLowerOrder_F * gsl_sf_hermite(nzm1_I, zeta_F);
        }
    }

    // (H,H',ζ,n) → (H'',H'-ζH,H''-2ζH'+(ζ²-1)H).
    const Eigen::MatrixXd ddH_F2D_sp_z = 2.0 * dH_F2D_sp_z * zeta_F1D_z.asDiagonal() - 2.0 * nz_F1D_sp.asDiagonal() * H_F2D_sp_z;

    const Eigen::MatrixXd dHterm_F2D_sp_z = dH_F2D_sp_z - H_F2D_sp_z * zeta_F1D_z.asDiagonal();

    const Eigen::MatrixXd ddHterm_F2D_sp_z = ddH_F2D_sp_z - 2.0 * dH_F2D_sp_z * zeta_F1D_z.asDiagonal() + H_F2D_sp_z * (zeta_F1D_z.array().square() - 1.0).matrix().asDiagonal();

    // {H_n,∂_ζH_n,∂²_ζH_n} → e^{ζ²/2}{φ_n,∂_zφ_n,∂²_zφ_n}.
    phi_F2D_sp_z.noalias() = N_F1D_sp.asDiagonal() * H_F2D_sp_z;
    dphi_F2D_sp_z.noalias() = N_F1D_sp.asDiagonal() * dHterm_F2D_sp_z / bz_F;
    ddphi_F2D_sp_z.noalias() = N_F1D_sp.asDiagonal() * ddHterm_F2D_sp_z / (bz_F * bz_F);

    // {φ,∂φ,∂²φ} × e^(-ζ²/2); w ÷ e^(-ζ²).
    const Eigen::VectorXd expHalf_F1D_z = (-0.5 * zeta_F1D_z.array().square()).exp();
    w_F1D_z.array() /= expHalf_F1D_z.array().square();
    phi_F2D_sp_z.array().rowwise() *= expHalf_F1D_z.transpose().array();
    dphi_F2D_sp_z.array().rowwise() *= expHalf_F1D_z.transpose().array();
    ddphi_F2D_sp_z.array().rowwise() *= expHalf_F1D_z.transpose().array();
}

inline void AxialLaguerreBasis::fill_grid(const GaussLaguerreMeshes& gl_meshes) {
    // η=(r/b_r)²; r dr=(b_r²/2)dη.
    const double jacobian_F = 0.5 * br_F * br_F;
    eta_F1D_r = gl_meshes.x_F1D_x;
    r_F1D_r = br_F * eta_F1D_r.array().sqrt().matrix();
    w_F1D_r = jacobian_F * gl_meshes.w_F1D_x;

    // (n,N_n,L_n,{η_i}) → work arrays.
    const int Nr_I = gl_meshes.x_F1D_x.size();
    const int Nsp_I = gl_meshes.n_I1D_na.size();
    const double norm_F = std::sqrt(2.0) / br_F;
    const Eigen::VectorXd nr_F1D_sp = gl_meshes.n_I1D_na.cast<double>();
    const Eigen::VectorXd Lambda_F1D_sp = gl_meshes.alpha_F1D_na;
    const Eigen::MatrixXd L_F2D_sp_r = gl_meshes.L_F2D_na_x;
    const Eigen::VectorXd N_F1D_sp = norm_F * gl_meshes.N_F1D_na;

    // {η^(Λ/2),(η-Λ)/η,[(η-Λ)²-2Λ]/(4η²)}.
    const Eigen::VectorXd I_F1D_sp = Eigen::VectorXd::Ones(Nsp_I);
    Eigen::MatrixXd etaPLambda_F2D_sp_r(Nsp_I, Nr_I);
    Eigen::MatrixXd etaMLambda_F2D_sp_r(Nsp_I, Nr_I);
    Eigen::MatrixXd etaMLambda2_F2D_sp_r(Nsp_I, Nr_I);
    for (int r_I = 0; r_I < Nr_I; ++r_I) {etaPLambda_F2D_sp_r.col(r_I) = (I_F1D_sp * eta_F1D_r(r_I)).array().pow(0.5 * Lambda_F1D_sp.array());}
    for (int r_I = 0; r_I < Nr_I; ++r_I) {etaMLambda_F2D_sp_r.col(r_I) = (eta_F1D_r(r_I) - Lambda_F1D_sp.array()) / eta_F1D_r(r_I);}
    for (int r_I = 0; r_I < Nr_I; ++r_I) {etaMLambda2_F2D_sp_r.col(r_I) = 0.25 * etaMLambda_F2D_sp_r.col(r_I).array().square() - 0.5 * Lambda_F1D_sp.array() / (eta_F1D_r(r_I) * eta_F1D_r(r_I));}

    // ∂_ηL_nr^Λ = -L_nr-1^(Λ+1).
    Eigen::MatrixXd dL_F2D_sp_r = Eigen::MatrixXd::Zero(Nsp_I, Nr_I);
    for (int sp_I = 0; sp_I < Nsp_I; ++sp_I) {
        const int nrm1_I = std::max(gl_meshes.n_I1D_na(sp_I) - 1, 0);
        const double hasLowerOrder_F = static_cast<double>(gl_meshes.n_I1D_na(sp_I) > 0);
        const double Lambdap1_F = Lambda_F1D_sp(sp_I) + 1.0;
        for (int r_I = 0; r_I < Nr_I; ++r_I) {
            const double eta_F = eta_F1D_r(r_I);
            dL_F2D_sp_r(sp_I, r_I) = -hasLowerOrder_F * gsl_sf_laguerre_n(nrm1_I, Lambdap1_F, eta_F);
        }
    }

    // η∂²_ηL = (η-Λ-1)∂_ηL - n_rL.
    Eigen::MatrixXd ddL_F2D_sp_r = dL_F2D_sp_r * eta_F1D_r.asDiagonal() - (Lambda_F1D_sp.array() + 1.0).matrix().asDiagonal() * dL_F2D_sp_r;
    ddL_F2D_sp_r.noalias() -= nr_F1D_sp.asDiagonal() * L_F2D_sp_r;
    ddL_F2D_sp_r.array().rowwise() /= eta_F1D_r.transpose().array();

    // (L,L',L'',η,Λ) → (u'/w,u''/w), w=η^{Λ/2}e^{-η/2}.
    const Eigen::MatrixXd dLterm_F2D_sp_r = dL_F2D_sp_r.array() - 0.5 * L_F2D_sp_r.array() * etaMLambda_F2D_sp_r.array();

    const Eigen::MatrixXd ddLterm_F2D_sp_r = ddL_F2D_sp_r.array() - dL_F2D_sp_r.array() * etaMLambda_F2D_sp_r.array() + L_F2D_sp_r.array() * etaMLambda2_F2D_sp_r.array();

    // η=(r/b_r)² → (∂_rη,∂²_rη)=(2r/b_r²,2/b_r²).
    const Eigen::VectorXd deta_F1D_r = 2.0 * r_F1D_r / (br_F * br_F);
    const double ddeta_F = 2.0 / (br_F * br_F);

    // NLη^{Λ/2} = e^{η/2}φ.
    phi_F2D_sp_r.noalias() = N_F1D_sp.asDiagonal() * L_F2D_sp_r;
    phi_F2D_sp_r.array() *= etaPLambda_F2D_sp_r.array();

    // e^{η/2}∂_ηφ → e^{η/2}∂_rφ.
    const Eigen::MatrixXd dphiTerm_F2D_sp_r = (N_F1D_sp.asDiagonal() * dLterm_F2D_sp_r).array() * etaPLambda_F2D_sp_r.array();
    dphi_F2D_sp_r.noalias() = dphiTerm_F2D_sp_r * deta_F1D_r.asDiagonal();

    // e^{η/2}∂²_ηφ → e^{η/2}∂²_rφ.
    const Eigen::MatrixXd ddphiTerm_F2D_sp_r = (N_F1D_sp.asDiagonal() * ddLterm_F2D_sp_r).array() * etaPLambda_F2D_sp_r.array();
    ddphi_F2D_sp_r.noalias() = ddphiTerm_F2D_sp_r * deta_F1D_r.array().square().matrix().asDiagonal();
    ddphi_F2D_sp_r += ddeta_F * dphiTerm_F2D_sp_r;

    // {φ,∂φ,∂²φ} × e^(-η/2); w ÷ e^(-η).
    const Eigen::VectorXd expHalf_F1D_r = (-0.5 * eta_F1D_r.array()).exp();
    w_F1D_r.array() /= expHalf_F1D_r.array().square();
    phi_F2D_sp_r.array().rowwise() *= expHalf_F1D_r.transpose().array();
    dphi_F2D_sp_r.array().rowwise() *= expHalf_F1D_r.transpose().array();
    ddphi_F2D_sp_r.array().rowwise() *= expHalf_F1D_r.transpose().array();
}

inline void AxialBasis::fill_from_separable_basis(const AxialHermiteBasis& zBasis, const AxialLaguerreBasis& rBasis) {
    // (z,r,ζ,η,w_z,w_r) ← (zBasis,rBasis).
    z_F1D_z = zBasis.z_F1D_z;
    r_F1D_r = rBasis.r_F1D_r;
    zeta_F1D_z = zBasis.zeta_F1D_z;
    eta_F1D_r = rBasis.eta_F1D_r;

    const int Nsp_I = labels_S1D_sp.size();
    const int Nz_I = z_F1D_z.size();
    const int Nr_I = r_F1D_r.size();

    // (N_sp,N_z,N_r) → {w,φ,∂φ,∂²φ}.
    w_F2D_z_r.resize(Nz_I, Nr_I);
    phi_F3D_sp_z_r = BasisTensor(Nsp_I, Nz_I, Nr_I);
    dphidr_F3D_sp_z_r = BasisTensor(Nsp_I, Nz_I, Nr_I);
    dphidz_F3D_sp_z_r = BasisTensor(Nsp_I, Nz_I, Nr_I);
    ddphidr_F3D_sp_z_r = BasisTensor(Nsp_I, Nz_I, Nr_I);
    ddphidz_F3D_sp_z_r = BasisTensor(Nsp_I, Nz_I, Nr_I);

    // (φ_z,φ_r)/√(2π) → {φ,∂_rφ,∂_zφ,∂²_rφ,∂²_zφ}; w=2πw_zw_r.
    const double pi_F = std::acos(-1.0);
    const double Nphi_F = 1.0 / std::sqrt(2.0 * pi_F);
    for (int r_I = 0; r_I < Nr_I; ++r_I) {
        for (int z_I = 0; z_I < Nz_I; ++z_I) {
            w_F2D_z_r(z_I, r_I) = 2.0 * pi_F * zBasis.w_F1D_z(z_I) * rBasis.w_F1D_r(r_I);
            for (int sp_I = 0; sp_I < Nsp_I; ++sp_I) {
                phi_F3D_sp_z_r(sp_I, z_I, r_I) = Nphi_F * zBasis.phi_F2D_sp_z(sp_I, z_I) * rBasis.phi_F2D_sp_r(sp_I, r_I);
                dphidr_F3D_sp_z_r(sp_I, z_I, r_I) = Nphi_F * zBasis.phi_F2D_sp_z(sp_I, z_I) * rBasis.dphi_F2D_sp_r(sp_I, r_I);
                dphidz_F3D_sp_z_r(sp_I, z_I, r_I) = Nphi_F * zBasis.dphi_F2D_sp_z(sp_I, z_I) * rBasis.phi_F2D_sp_r(sp_I, r_I);
                ddphidr_F3D_sp_z_r(sp_I, z_I, r_I) = Nphi_F * zBasis.phi_F2D_sp_z(sp_I, z_I) * rBasis.ddphi_F2D_sp_r(sp_I, r_I);
                ddphidz_F3D_sp_z_r(sp_I, z_I, r_I) = Nphi_F * zBasis.ddphi_F2D_sp_z(sp_I, z_I) * rBasis.phi_F2D_sp_r(sp_I, r_I);
            }
        }
    }
}
