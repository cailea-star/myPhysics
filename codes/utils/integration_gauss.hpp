/**
 * @file    integration_gauss.hpp
 * @author  cailea
 * @date    2026-04-27
 * @brief   Gauss quadrature meshes and orthogonal polynomials.
 */

#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <unsupported/Eigen/CXX11/Tensor>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_hermite.h>
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_sf_legendre.h>

/** @brief Gauss-Hermite mesh and polynomials. */
class GaussHermiteMeshes {
public:
    Eigen::VectorXi n_I1D_n{};      // n ≥ 0.
    Eigen::VectorXd N_F1D_n{};      // ∫(N_nH_n)²e^-x²dx = 1.
    Eigen::VectorXd x_F1D_x{};      // {x_i}: Gauss-Hermite nodes.
    Eigen::VectorXd w_F1D_x{};      // {w_i}: e^-x² quadrature.
    Eigen::MatrixXd H_F2D_n_x{};    // H_n(x_i).

    /**
     * @brief  Create an empty mesh.
     * @math   N_x = 0
     * @output Empty mesh.
     */
    GaussHermiteMeshes() = default;

    /**
     * @brief  Generate Gauss-Hermite nodes and weights.
     * @math   W(x) = e^-x²
     * @output Nodes and weights.
     * @note   N_x>0.
     */
    GaussHermiteMeshes(int Nx_I) {
        n_I1D_n.resize(1);
        N_F1D_n.resize(1);
        x_F1D_x.resize(Nx_I);
        w_F1D_x.resize(Nx_I);
        H_F2D_n_x.resize(1, Nx_I);

        n_I1D_n.setZero();
        N_F1D_n.setZero();
        H_F2D_n_x.setZero();

        // (N_x,W) → {x_i,w_i}.
        assert(Nx_I > 0);
        double a_F = 0.0;           // x → x-a.
        double b_F = 1.0;           // x² → x²/b.
        double alpha_F = 0.0;       // W(x) = e^-x² |x|^2α.
        gsl_integration_fixed_workspace* workspace = gsl_integration_fixed_alloc(gsl_integration_fixed_hermite, Nx_I, a_F, b_F, alpha_F, 0.0);
        const double* nodes_F1D_x = gsl_integration_fixed_nodes(workspace);
        const double* weights_F1D_x = gsl_integration_fixed_weights(workspace);
        for (int i_I = 0; i_I < Nx_I; ++i_I) {x_F1D_x(i_I) = nodes_F1D_x[i_I];}
        for (int i_I = 0; i_I < Nx_I; ++i_I) {w_F1D_x(i_I) = weights_F1D_x[i_I];}
        gsl_integration_fixed_free(workspace);
    };

    /**
     * @brief  Generate nodes, weights, and normalized Hermite polynomials.
     * @math   ∫N_nH_nN_mH_m e^-x² dx = δ_nm
     * @output Nodes, weights, normalizations, and polynomials.
     * @note   n_i≥0; N_x>0.
     */
    GaussHermiteMeshes(const Eigen::VectorXi& n_I1D_n_, int Nx_I) {
        n_I1D_n = n_I1D_n_;
        N_F1D_n.resize(n_I1D_n_.size());
        x_F1D_x.resize(Nx_I);
        w_F1D_x.resize(Nx_I);
        H_F2D_n_x.resize(n_I1D_n_.size(), Nx_I);
        // (N_x,W) → {x_i,w_i}.
        assert(Nx_I > 0);
        assert((n_I1D_n_.array() >= 0).all());
        double a_F = 0.0;           // x → x-a.
        double b_F = 1.0;           // x² → x²/b.
        double alpha_F = 0.0;       // W(x) = e^-x² |x|^2α.
        gsl_integration_fixed_workspace* workspace = gsl_integration_fixed_alloc(gsl_integration_fixed_hermite, Nx_I, a_F, b_F, alpha_F, 0.0);
        const double* nodes_F1D_x = gsl_integration_fixed_nodes(workspace);
        const double* weights_F1D_x = gsl_integration_fixed_weights(workspace);
        for (int i_I = 0; i_I < Nx_I; ++i_I) {x_F1D_x(i_I) = nodes_F1D_x[i_I];}
        for (int i_I = 0; i_I < Nx_I; ++i_I) {w_F1D_x(i_I) = weights_F1D_x[i_I];}
        gsl_integration_fixed_free(workspace);
        // (n,x_i) → (N_n,H_n(x_i)).
        for (int in_I = 0; in_I < n_I1D_n.size(); ++in_I) {
            for (int ix_I = 0; ix_I < Nx_I; ++ix_I) {
                double x_F = x_F1D_x(ix_I);
                int n_I = n_I1D_n(in_I);
                N_F1D_n(in_I) = std::pow(M_PI, -0.25) * std::pow(2.0, -0.5 * n_I) * std::exp(-0.5 * gsl_sf_lnfact(n_I));
                H_F2D_n_x(in_I, ix_I) = gsl_sf_hermite(n_I, x_F);
            }
        }
    };
};

/**
 * @brief Gauss-Laguerre mesh and polynomials.
 * @note  Q_α[f]=Σ_iw_ix_i^αf(x_i).
 */
class GaussLaguerreMeshes {
public:
    Eigen::VectorXi n_I1D_na{};     // n ≥ 0.
    Eigen::VectorXd alpha_F1D_na{}; // α > -1.
    Eigen::VectorXd N_F1D_na{};     // ∫(N_nL_n^α)²x^αe^-xdx = 1.
    Eigen::VectorXd x_F1D_x{};      // {x_i}: Gauss-Laguerre nodes.
    Eigen::VectorXd w_F1D_x{};      // {w_i}: e^-x quadrature.
    Eigen::MatrixXd L_F2D_na_x{};   // L_n^α(x_i).

    /**
     * @brief  Create an empty mesh.
     * @math   N_x = 0
     * @output Empty mesh.
     */
    GaussLaguerreMeshes() = default;

    /**
     * @brief  Generate Gauss-Laguerre nodes and weights.
     * @math   W(x) = e^-x
     * @output Nodes and weights.
     * @note   N_x>0.
     */
    GaussLaguerreMeshes(int Nx_I) {
        n_I1D_na.resize(1);
        alpha_F1D_na.resize(1);
        N_F1D_na.resize(1);
        x_F1D_x.resize(Nx_I);
        w_F1D_x.resize(Nx_I);
        L_F2D_na_x.resize(1, Nx_I);

        n_I1D_na.setZero();
        alpha_F1D_na.setZero();
        N_F1D_na.setZero();
        L_F2D_na_x.setZero();

        // (N_x,W) → {x_i,w_i}.
        assert(Nx_I > 0);
        double b_F = 1.0;           // x → x/b.
        double weightAlpha_F = 0.0; // W(x) = x^α e^(-x/b).
        gsl_integration_fixed_workspace* workspace = gsl_integration_fixed_alloc(gsl_integration_fixed_laguerre, Nx_I, 0.0, b_F, weightAlpha_F, 0.0);
        const double* nodes_F1D_x = gsl_integration_fixed_nodes(workspace);
        const double* weights_F1D_x = gsl_integration_fixed_weights(workspace);
        for (int i_I = 0; i_I < Nx_I; ++i_I) {x_F1D_x(i_I) = nodes_F1D_x[i_I];}
        for (int i_I = 0; i_I < Nx_I; ++i_I) {w_F1D_x(i_I) = weights_F1D_x[i_I];}
        gsl_integration_fixed_free(workspace);
    };

    /**
     * @brief  Generate nodes, weights, and normalized Laguerre polynomials.
     * @math   ∫N_nL_n^αN_mL_m^α x^αe^-x dx = δ_nm
     * @output Nodes, weights, normalizations, and polynomials.
     * @note   n_i≥0, α_i>-1; |n|=|α|; N_x>0.
     */
    GaussLaguerreMeshes(const Eigen::VectorXi& n_I1D_na_, const Eigen::VectorXd& alpha_F1D_na_, int Nx_I) {
        n_I1D_na = n_I1D_na_;
        alpha_F1D_na = alpha_F1D_na_;
        N_F1D_na.resize(n_I1D_na_.size());
        x_F1D_x.resize(Nx_I);
        w_F1D_x.resize(Nx_I);
        L_F2D_na_x.resize(n_I1D_na_.size(), Nx_I);
        // (N_x,W) → {x_i,w_i}.
        assert(Nx_I > 0);
        assert(n_I1D_na_.size() == alpha_F1D_na_.size());
        assert((n_I1D_na_.array() >= 0).all());
        assert((alpha_F1D_na_.array() > -1.0).all());
        double b_F = 1.0;           // x → x/b.
        double weightAlpha_F = 0.0; // W(x) = x^α e^(-x/b).
        gsl_integration_fixed_workspace* workspace = gsl_integration_fixed_alloc(gsl_integration_fixed_laguerre, Nx_I, 0.0, b_F, weightAlpha_F, 0.0);
        const double* nodes_F1D_x = gsl_integration_fixed_nodes(workspace);
        const double* weights_F1D_x = gsl_integration_fixed_weights(workspace);
        for (int i_I = 0; i_I < Nx_I; ++i_I) {x_F1D_x(i_I) = nodes_F1D_x[i_I];}
        for (int i_I = 0; i_I < Nx_I; ++i_I) {w_F1D_x(i_I) = weights_F1D_x[i_I];}
        gsl_integration_fixed_free(workspace);
        // (n,α,x_i) → (N_n^α,L_n^α(x_i)).
        for (int ix_I = 0; ix_I < Nx_I; ++ix_I) {
            for (int ina_I = 0; ina_I < n_I1D_na.size(); ++ina_I) {
                int n_I = n_I1D_na(ina_I);
                double alpha_F = alpha_F1D_na(ina_I);
                double x_F = x_F1D_x(ix_I);
                N_F1D_na(ina_I) = std::exp(0.5 * (gsl_sf_lnfact(n_I) - gsl_sf_lngamma(n_I + alpha_F + 1)));
                L_F2D_na_x(ina_I, ix_I) = gsl_sf_laguerre_n(n_I, alpha_F, x_F);
            }
        }
    }
};

/** @brief Gauss-Legendre mesh and polynomials. */
class GaussLegendreMeshes {
public:
    Eigen::VectorXi l_I1D_lm{};     // l ≥ 0.
    Eigen::VectorXi m_I1D_lm{};     // 0 ≤ m ≤ l.
    Eigen::VectorXd N_F1D_lm{};     // ∫(N_lmP_l^m)²dx = 1.
    Eigen::VectorXd x_F1D_x{};      // x_i ∈ [-1,1].
    Eigen::VectorXd w_F1D_x{};      // {w_i}: unit-weight quadrature.
    Eigen::MatrixXd P_F2D_lm_x{};   // P_l^m(x_i).

    /**
     * @brief  Create an empty mesh.
     * @math   N_x = 0
     * @output Empty mesh.
     */
    GaussLegendreMeshes() = default;

    /**
     * @brief  Generate Gauss-Legendre nodes and weights.
     * @math   W(x) = 1, x ∈ [-1, 1]
     * @output Nodes and weights.
     * @note   N_x>0.
     */
    GaussLegendreMeshes(int Nx_I) {
        l_I1D_lm.resize(1);
        m_I1D_lm.resize(1);
        N_F1D_lm.resize(1);
        x_F1D_x.resize(Nx_I);
        w_F1D_x.resize(Nx_I);
        P_F2D_lm_x.resize(1, Nx_I);

        l_I1D_lm.setZero();
        m_I1D_lm.setZero();
        N_F1D_lm.setZero();
        P_F2D_lm_x.setZero();

        // (N_x,[-1,1]) → {x_i,w_i}.
        assert(Nx_I > 0);
        double xmin_F = -1.0;       // x_min = -1 = cos π.
        double xmax_F = 1.0;        // x_max = 1 = cos 0.
        gsl_integration_fixed_workspace* workspace = gsl_integration_fixed_alloc(gsl_integration_fixed_legendre, Nx_I, xmin_F, xmax_F, 0.0, 0.0);
        const double* nodes_F1D_x = gsl_integration_fixed_nodes(workspace);
        const double* weights_F1D_x = gsl_integration_fixed_weights(workspace);
        for (int i_I = 0; i_I < Nx_I; ++i_I) {x_F1D_x(i_I) = nodes_F1D_x[i_I];}
        for (int i_I = 0; i_I < Nx_I; ++i_I) {w_F1D_x(i_I) = weights_F1D_x[i_I];}
        gsl_integration_fixed_free(workspace);
    };

    /**
     * @brief  Generate nodes, weights, and normalized Legendre polynomials.
     * @math   ∫N_lmP_l^mN_kmP_k^m dx = δ_lk
     * @output Nodes, weights, normalizations, and polynomials.
     * @note   l_i≥m_i≥0; |l|=|m|; N_x>0.
     */
    GaussLegendreMeshes(const Eigen::VectorXi& l_I1D_lm_, const Eigen::VectorXi& m_I1D_lm_, int Nx_I) {
        l_I1D_lm = l_I1D_lm_;
        m_I1D_lm = m_I1D_lm_;
        N_F1D_lm.resize(l_I1D_lm_.size());
        x_F1D_x.resize(Nx_I);
        w_F1D_x.resize(Nx_I);
        P_F2D_lm_x.resize(l_I1D_lm_.size(), Nx_I);
        // (N_x,[-1,1]) → {x_i,w_i}.
        assert(Nx_I > 0);
        assert(l_I1D_lm_.size() == m_I1D_lm_.size());
        assert((l_I1D_lm_.array() >= 0).all());
        assert((m_I1D_lm_.array() >= 0).all());
        assert((m_I1D_lm_.array() <= l_I1D_lm_.array()).all());
        double xmin_F = -1.0;       // x_min = -1 = cos π.
        double xmax_F = 1.0;        // x_max = 1 = cos 0.
        gsl_integration_fixed_workspace* workspace = gsl_integration_fixed_alloc(gsl_integration_fixed_legendre, Nx_I, xmin_F, xmax_F, 0.0, 0.0);
        const double* nodes_F1D_x = gsl_integration_fixed_nodes(workspace);
        const double* weights_F1D_x = gsl_integration_fixed_weights(workspace);
        for (int i_I = 0; i_I < Nx_I; ++i_I) {x_F1D_x(i_I) = nodes_F1D_x[i_I];}
        for (int i_I = 0; i_I < Nx_I; ++i_I) {w_F1D_x(i_I) = weights_F1D_x[i_I];}
        gsl_integration_fixed_free(workspace);
        // (l,m,x_i) → (N_lm,P_l^m(x_i)).
        for (int ix_I = 0; ix_I < Nx_I; ++ix_I) {
            for (int ilm_I = 0; ilm_I < l_I1D_lm.size(); ++ilm_I) {
                int l_I = l_I1D_lm_(ilm_I);
                int m_I = m_I1D_lm_(ilm_I);
                double x_F = x_F1D_x(ix_I);
                N_F1D_lm(ilm_I) = std::sqrt((2 * l_I + 1) / 2.0) * std::exp(0.5 * (gsl_sf_lnfact(l_I - m_I) - gsl_sf_lnfact(l_I + m_I)));
                P_F2D_lm_x(ilm_I, ix_I) = gsl_sf_legendre_Plm(l_I, m_I, x_F);
            }
        }
    }
};
