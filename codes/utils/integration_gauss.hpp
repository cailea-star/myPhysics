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
    Eigen::VectorXi n_I1D_n;        // Orders.
    Eigen::VectorXd x_F1D_x;        // Nodes.
    Eigen::VectorXd N_F1D_n;        // Normalizations.
    Eigen::VectorXd w_F1D_x;        // Weights.
    Eigen::MatrixXd H_F2D_n_x;      // Polynomial values.

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
     */
    GaussHermiteMeshes(int nx_points) : x_F1D_x(nx_points), w_F1D_x(nx_points), 
    n_I1D_n(Eigen::VectorXi::Zero(1)), N_F1D_n(Eigen::VectorXd::Zero(1)), 
    H_F2D_n_x(Eigen::MatrixXd::Zero(1, nx_points)) {
        assert(nx_points > 0);
        double a_F = 0.0;           // x -> x - a.
        double b_F = 1.0;           // x² -> x² / b.
        double alpha_F = 0.0;       // W(x) = e^-x² |x|^2α.
        gsl_integration_fixed_workspace* wt = gsl_integration_fixed_alloc(gsl_integration_fixed_hermite, nx_points, a_F, b_F, alpha_F, 0.0);
        const double* nodes = gsl_integration_fixed_nodes(wt);
        const double* weights = gsl_integration_fixed_weights(wt);
        for (int i = 0; i < nx_points; ++i) {x_F1D_x(i) = nodes[i];}
        for (int i = 0; i < nx_points; ++i) {w_F1D_x(i) = weights[i];}
        gsl_integration_fixed_free(wt);
    };

    /**
     * @brief  Generate nodes, weights, and normalized Hermite polynomials.
     * @math   ∫N_nH_nN_mH_m e^-x² dx = δ_nm
     * @output Nodes, weights, normalizations, and polynomials.
     */
    GaussHermiteMeshes(const Eigen::VectorXi& n_I1D_n_, int nx_points) : x_F1D_x(nx_points), w_F1D_x(nx_points),
    n_I1D_n(n_I1D_n_), N_F1D_n(n_I1D_n_.size()), H_F2D_n_x(n_I1D_n_.size(), nx_points) {
        assert(nx_points > 0);
        assert((n_I1D_n_.array() >= 0).all());
        double a_F = 0.0;           // x -> x - a.
        double b_F = 1.0;           // x² -> x² / b.
        double alpha_F = 0.0;       // W(x) = e^-x² |x|^2α.
        gsl_integration_fixed_workspace* wt = gsl_integration_fixed_alloc(gsl_integration_fixed_hermite, nx_points, a_F, b_F, alpha_F, 0.0);
        const double* nodes = gsl_integration_fixed_nodes(wt);
        const double* weights = gsl_integration_fixed_weights(wt);
        for (int i = 0; i < nx_points; ++i) {x_F1D_x(i) = nodes[i];}
        for (int i = 0; i < nx_points; ++i) {w_F1D_x(i) = weights[i];}
        gsl_integration_fixed_free(wt);
        // Evaluate Hermite polynomials.
        for (int in = 0; in < n_I1D_n.size(); ++in) for (int ix = 0; ix < nx_points; ++ix) 
        {
            double x_F = x_F1D_x(ix);
            int n_I = n_I1D_n(in);
            N_F1D_n(in) = std::pow(M_PI, -0.25) * std::pow(2.0, -0.5 * n_I) * std::exp(-0.5 * gsl_sf_lnfact(n_I));
            H_F2D_n_x(in, ix) = gsl_sf_hermite(n_I, x_F);
        }
    };
};

/**
 * @brief Gauss-Laguerre mesh and polynomials.
 * @note  w_F1D_x weights e^-x; multiply by x_F1D_x^α.
 */
class GaussLaguerreMeshes {
public:
    Eigen::VectorXi n_I1D_na;       // Orders.
    Eigen::VectorXd a_F1D_na;       // Parameters.
    Eigen::VectorXd N_F1D_na;       // Normalizations.
    Eigen::VectorXd x_F1D_x;        // Nodes.
    Eigen::VectorXd w_F1D_x;        // Weights.
    Eigen::MatrixXd L_F2D_na_x;     // Polynomial values.

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
     */
    GaussLaguerreMeshes(int nx_points) : x_F1D_x(nx_points), w_F1D_x(nx_points),
    n_I1D_na(Eigen::VectorXi::Zero(1)), a_F1D_na(Eigen::VectorXd::Zero(1)),
    N_F1D_na(Eigen::VectorXd::Zero(1)), L_F2D_na_x(Eigen::MatrixXd::Zero(1, nx_points)) {
        assert(nx_points > 0);
        double b_F = 1.0;           // x -> x / b.
        double alpha_F = 0.0;       // W(x) = x^α e^(-x/b).
        gsl_integration_fixed_workspace* wt = gsl_integration_fixed_alloc(gsl_integration_fixed_laguerre, nx_points, 0.0, b_F, alpha_F, 0.0);
        const double* nodes = gsl_integration_fixed_nodes(wt);
        const double* weights = gsl_integration_fixed_weights(wt);
        for (int i = 0; i < nx_points; ++i) {x_F1D_x(i) = nodes[i];}
        for (int i = 0; i < nx_points; ++i) {w_F1D_x(i) = weights[i];}
        gsl_integration_fixed_free(wt);
    };

    /**
     * @brief  Generate nodes, weights, and normalized Laguerre polynomials.
     * @math   ∫N_nL_n^αN_mL_m^α x^αe^-x dx = δ_nm
     * @output Nodes, weights, normalizations, and polynomials.
     */
    GaussLaguerreMeshes(const Eigen::VectorXi& n_I1D_na_, const Eigen::VectorXd& a_F1D_na_, int nx_points) : 
    x_F1D_x(nx_points), w_F1D_x(nx_points), n_I1D_na(n_I1D_na_), a_F1D_na(a_F1D_na_), 
    N_F1D_na(n_I1D_na_.size()), L_F2D_na_x(n_I1D_na_.size(), nx_points) {
        assert(nx_points > 0);
        assert(n_I1D_na_.size() == a_F1D_na_.size());
        assert((n_I1D_na_.array() >= 0).all());
        assert((a_F1D_na_.array() > -1.0).all());
        double b_F = 1.0;           // x -> x / b.
        double alpha_F = 0.0;       // W(x) = x^α e^(-x/b).
        gsl_integration_fixed_workspace* wt = gsl_integration_fixed_alloc(gsl_integration_fixed_laguerre, nx_points, 0.0, b_F, alpha_F, 0.0);
        const double* nodes = gsl_integration_fixed_nodes(wt);
        const double* weights = gsl_integration_fixed_weights(wt);
        for (int i = 0; i < nx_points; ++i) {x_F1D_x(i) = nodes[i];}
        for (int i = 0; i < nx_points; ++i) {w_F1D_x(i) = weights[i];}
        gsl_integration_fixed_free(wt);
        // Evaluate Laguerre polynomials.
        for (int ix = 0; ix < nx_points; ++ix) for (int ina = 0; ina < n_I1D_na.size(); ++ina) 
        {
            int n_I = n_I1D_na(ina);
            double a_F = a_F1D_na(ina);
            double x_F = x_F1D_x(ix);
            N_F1D_na(ina) = std::exp(0.5 * (gsl_sf_lnfact(n_I) - gsl_sf_lngamma(n_I + a_F + 1)));
            L_F2D_na_x(ina, ix) = gsl_sf_laguerre_n(n_I, a_F, x_F);
        }
    }
};

/** @brief Gauss-Legendre mesh and polynomials. */
class GaussLegendreMeshes {
public:
    Eigen::VectorXi l_I1D_lm;       // Degrees.
    Eigen::VectorXi m_I1D_lm;       // Orders.
    Eigen::VectorXd N_F1D_lm;       // Normalizations.
    Eigen::VectorXd x_F1D_x;        // Nodes.
    Eigen::VectorXd w_F1D_x;        // Weights.
    Eigen::MatrixXd P_F2D_lm_x;     // Polynomial values.

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
     */
    GaussLegendreMeshes(int nx_points) : x_F1D_x(nx_points), w_F1D_x(nx_points),
    l_I1D_lm(Eigen::VectorXi::Zero(1)), m_I1D_lm(Eigen::VectorXi::Zero(1)),
    N_F1D_lm(Eigen::VectorXd::Zero(1)), P_F2D_lm_x(Eigen::MatrixXd::Zero(1, nx_points)) {
        assert(nx_points > 0);
        double xlo_F = -1.0;        // x = cos π.
        double xup_F = 1.0;         // x = cos 0.
        gsl_integration_fixed_workspace* wt = gsl_integration_fixed_alloc(gsl_integration_fixed_legendre, nx_points, xlo_F, xup_F, 0.0, 0.0);
        const double* nodes = gsl_integration_fixed_nodes(wt);
        const double* weights = gsl_integration_fixed_weights(wt);
        for (int i = 0; i < nx_points; ++i) {x_F1D_x(i) = nodes[i];}
        for (int i = 0; i < nx_points; ++i) {w_F1D_x(i) = weights[i];}
        gsl_integration_fixed_free(wt);
    };

    /**
     * @brief  Generate nodes, weights, and normalized Legendre polynomials.
     * @math   ∫N_lmP_l^mN_kmP_k^m dx = δ_lk
     * @output Nodes, weights, normalizations, and polynomials.
     */
    GaussLegendreMeshes(const Eigen::VectorXi& l_I1D_lm_, const Eigen::VectorXi& m_I1D_lm_, int nx_points) : 
    x_F1D_x(nx_points), w_F1D_x(nx_points), l_I1D_lm(l_I1D_lm_), m_I1D_lm(m_I1D_lm_), 
    N_F1D_lm(l_I1D_lm_.size()), P_F2D_lm_x(l_I1D_lm_.size(), nx_points) {
        assert(nx_points > 0);
        assert(l_I1D_lm_.size() == m_I1D_lm_.size());
        assert((l_I1D_lm_.array() >= 0).all());
        assert((m_I1D_lm_.array() >= 0).all());
        assert((m_I1D_lm_.array() <= l_I1D_lm_.array()).all());
        double xlo_F = -1.0;        // x = cos π.
        double xup_F = 1.0;         // x = cos 0.
        gsl_integration_fixed_workspace* wt = gsl_integration_fixed_alloc(gsl_integration_fixed_legendre, nx_points, xlo_F, xup_F, 0.0, 0.0);
        const double* nodes = gsl_integration_fixed_nodes(wt);
        const double* weights = gsl_integration_fixed_weights(wt);
        for (int i = 0; i < nx_points; ++i) {x_F1D_x(i) = nodes[i];}
        for (int i = 0; i < nx_points; ++i) {w_F1D_x(i) = weights[i];}
        gsl_integration_fixed_free(wt);
        // Evaluate Legendre polynomials.
        for (int ix = 0; ix < nx_points; ++ix) for (int ilm = 0; ilm < l_I1D_lm.size(); ++ilm) {
            int l_I = l_I1D_lm_(ilm);
            int m_I = m_I1D_lm_(ilm);
            double x_F = x_F1D_x(ix);
            N_F1D_lm(ilm) = std::sqrt((2 * l_I + 1) / 2.0) * std::exp(0.5 * (gsl_sf_lnfact(l_I - m_I) - gsl_sf_lnfact(l_I + m_I)));
            P_F2D_lm_x(ilm, ix) = gsl_sf_legendre_Plm(l_I, m_I, x_F);
        }
    }
};
