/**
 * @file    integration_gauss.hpp
 * @author  cailea
 * @date    2026-04-27
 * @brief   调用 gsl 生成积分格点的类定义
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

/**
 * Gauss-Hermite 积分格点x_i与权重w_i, 
 * 厄密多项式 H_n(x) 在格点x_i 的生成,
 * 归一化系数的生成 N_n = (2^n n! sqrt(pi))^(-0.5) 的生成
 */
class GaussHermiteMeshes {
public:
    Eigen::VectorXi n_I1D_n;        // n量子数列表
    Eigen::VectorXd x_F1D_x;        // 积分格点位置
    Eigen::VectorXd N_F1D_n;        // N_n 归一化常数列表, 与n_I1D_n对应
    Eigen::VectorXd w_F1D_x;        // 积分格点权重
    Eigen::MatrixXd H_F2D_n_x;      // Hn(x) 多项式矩阵, 行对应n量子数, 列对应x格点位置
    GaussHermiteMeshes() = default;

    GaussHermiteMeshes(int nx_points) : x_F1D_x(nx_points), w_F1D_x(nx_points), 
    n_I1D_n(Eigen::VectorXi::Zero(1)), N_F1D_n(Eigen::VectorXd::Zero(1)), 
    H_F2D_n_x(Eigen::MatrixXd::Zero(1, nx_points)) {
        assert(nx_points > 0);
        // 调用 GSL 生成 Gauss-Hermite 格点和权重
        double a_F = 0.0;           // 中心平移变换 x -> (x-a)
        double b_F = 1.0;           // 尺度缩放变换 x^2 -> x^2 / b 
        double alpha_F = 0.0;       // 权函数参数 W(x) = exp(-x^2) * |x|^(2α)
        gsl_integration_fixed_workspace* wt = gsl_integration_fixed_alloc(gsl_integration_fixed_hermite, nx_points, a_F, b_F, alpha_F, 0.0);
        const double* nodes = gsl_integration_fixed_nodes(wt);
        const double* weights = gsl_integration_fixed_weights(wt);
        for (int i = 0; i < nx_points; ++i) {x_F1D_x(i) = nodes[i];}
        for (int i = 0; i < nx_points; ++i) {w_F1D_x(i) = weights[i];}
        gsl_integration_fixed_free(wt);
    };

    GaussHermiteMeshes(const Eigen::VectorXi& n_I1D_n_, int nx_points) : x_F1D_x(nx_points), w_F1D_x(nx_points),
    n_I1D_n(n_I1D_n_), N_F1D_n(n_I1D_n_.size()), H_F2D_n_x(n_I1D_n_.size(), nx_points) {
        assert(nx_points > 0);
        assert((n_I1D_n_.array() >= 0).all());
        // 调用 GSL 生成 Gauss-Hermite 格点和权重
        double a_F = 0.0;           // 中心平移变换 x -> (x-a)
        double b_F = 1.0;           // 尺度缩放变换 x^2 -> x^2 / b 
        double alpha_F = 0.0;       // 权函数参数 W(x) = exp(-x^2) * |x|^(2α)
        gsl_integration_fixed_workspace* wt = gsl_integration_fixed_alloc(gsl_integration_fixed_hermite, nx_points, a_F, b_F, alpha_F, 0.0);
        const double* nodes = gsl_integration_fixed_nodes(wt);
        const double* weights = gsl_integration_fixed_weights(wt);
        for (int i = 0; i < nx_points; ++i) {x_F1D_x(i) = nodes[i];}
        for (int i = 0; i < nx_points; ++i) {w_F1D_x(i) = weights[i];}
        gsl_integration_fixed_free(wt);
        // 计算 Hermite 多项式矩阵
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
 * Gauss-Laguerre 积分格点x_i与权重w_i, 
 * 广义拉盖尔多项式 L_n^(a)(x) * N_n^a 在格点x_i 的生成, 
 * 归一化系数的生成 N = sqrt[n! / (n + a)!]
 * @note w_F1D_x is for e^-x; multiply by x_F1D_x^a for x^a e^-x.
 */
class GaussLaguerreMeshes {
public:
    Eigen::VectorXi n_I1D_na;       // n主量子数表
    Eigen::VectorXd a_F1D_na;       // a参数列表, 与n_I1D_na对应;
    Eigen::VectorXd N_F1D_na;       // N_n^a 归一化常数列表, 与n_I1D_na和a_F1D_na对应
    Eigen::VectorXd x_F1D_x;        // 积分格点位置
    Eigen::VectorXd w_F1D_x;        // 积分格点权重
    Eigen::MatrixXd L_F2D_na_x;     // L_n^(a)(x) 多项式矩阵, 行对应(n,a)量子数, 列对应x格点位置
    GaussLaguerreMeshes() = default;

    GaussLaguerreMeshes(int nx_points) : x_F1D_x(nx_points), w_F1D_x(nx_points),
    n_I1D_na(Eigen::VectorXi::Zero(1)), a_F1D_na(Eigen::VectorXd::Zero(1)),
    N_F1D_na(Eigen::VectorXd::Zero(1)), L_F2D_na_x(Eigen::MatrixXd::Zero(1, nx_points)) {
        assert(nx_points > 0);
        // 调用 GSL 生成 Gauss-Laguerre 格点和权重
        double b_F = 1.0;           // 尺度缩放变换 x -> x / b
        double alpha_F = 0.0;       // 权函数参数 W(x) = x^α * exp(-x / b)
        gsl_integration_fixed_workspace* wt = gsl_integration_fixed_alloc(gsl_integration_fixed_laguerre, nx_points, 0.0, b_F, alpha_F, 0.0);
        const double* nodes = gsl_integration_fixed_nodes(wt);
        const double* weights = gsl_integration_fixed_weights(wt);
        for (int i = 0; i < nx_points; ++i) {x_F1D_x(i) = nodes[i];}
        for (int i = 0; i < nx_points; ++i) {w_F1D_x(i) = weights[i];}
        gsl_integration_fixed_free(wt);
    };

    GaussLaguerreMeshes(const Eigen::VectorXi& n_I1D_na_, const Eigen::VectorXd& a_F1D_na_, int nx_points) : 
    x_F1D_x(nx_points), w_F1D_x(nx_points), n_I1D_na(n_I1D_na_), a_F1D_na(a_F1D_na_), 
    N_F1D_na(n_I1D_na_.size()), L_F2D_na_x(n_I1D_na_.size(), nx_points) {
        assert(nx_points > 0);
        assert(n_I1D_na_.size() == a_F1D_na_.size());
        assert((n_I1D_na_.array() >= 0).all());
        assert((a_F1D_na_.array() > -1.0).all());
        // 调用 GSL 生成 Gauss-Laguerre 格点和权重
        double b_F = 1.0;           // 尺度缩放变换 x -> x / b
        double alpha_F = 0.0;       // 权函数参数 W(x) = x^α * exp(-x / b)
        gsl_integration_fixed_workspace* wt = gsl_integration_fixed_alloc(gsl_integration_fixed_laguerre, nx_points, 0.0, b_F, alpha_F, 0.0);
        const double* nodes = gsl_integration_fixed_nodes(wt);
        const double* weights = gsl_integration_fixed_weights(wt);
        for (int i = 0; i < nx_points; ++i) {x_F1D_x(i) = nodes[i];}
        for (int i = 0; i < nx_points; ++i) {w_F1D_x(i) = weights[i];}
        gsl_integration_fixed_free(wt);
        // 计算广义拉盖尔多项式
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

/**
 * Gauss-Legendre 积分格点x_i与权重w_i,
 * 广义勒让德多项式 P_l^m(x) 在格点位置的生成, 
 * 归一化系数的生成 N_lm = sqrt[(2l+1)/(2) * (l-m)!/(l+m)!]
 */
class GaussLegendreMeshes {
public:
    Eigen::VectorXi l_I1D_lm;       // l量子数列表, 
    Eigen::VectorXi m_I1D_lm;       // m量子数列表, 与l_I1D_lm一一对应
    Eigen::VectorXd N_F1D_lm;       // N_l^m 归一化常数列表, 与l_I1D_lm和m_I1D_lm对应
    Eigen::VectorXd x_F1D_x;        // 积分格点位置, cos(θ)格点位置, (-1, 1)区间
    Eigen::VectorXd w_F1D_x;        // 积分格点权重, 对应cos(θ)格点位置的权重
    Eigen::MatrixXd P_F2D_lm_x;     // 广义 Legendre 多项式矩阵 P_l^m(x), 行对应(l,m)量子数, 列对应x格点位置
    GaussLegendreMeshes() = default;

    GaussLegendreMeshes(int nx_points) : x_F1D_x(nx_points), w_F1D_x(nx_points),
    l_I1D_lm(Eigen::VectorXi::Zero(1)), m_I1D_lm(Eigen::VectorXi::Zero(1)),
    N_F1D_lm(Eigen::VectorXd::Zero(1)), P_F2D_lm_x(Eigen::MatrixXd::Zero(1, nx_points)) {
        assert(nx_points > 0);
        // 调用 GSL 生成 Gauss-Legendre 格点和权重
        double xlo_F = -1.0;        // cos(θ)格点位置的下界, 对应θ=π
        double xup_F = 1.0;         // cos(θ)格点位置的上界, 对应θ=0
        gsl_integration_fixed_workspace* wt = gsl_integration_fixed_alloc(gsl_integration_fixed_legendre, nx_points, xlo_F, xup_F, 0.0, 0.0);
        const double* nodes = gsl_integration_fixed_nodes(wt);
        const double* weights = gsl_integration_fixed_weights(wt);
        for (int i = 0; i < nx_points; ++i) {x_F1D_x(i) = nodes[i];}
        for (int i = 0; i < nx_points; ++i) {w_F1D_x(i) = weights[i];}
        gsl_integration_fixed_free(wt);
    };

    GaussLegendreMeshes(const Eigen::VectorXi& l_I1D_lm_, const Eigen::VectorXi& m_I1D_lm_, int nx_points) : 
    x_F1D_x(nx_points), w_F1D_x(nx_points), l_I1D_lm(l_I1D_lm_), m_I1D_lm(m_I1D_lm_), 
    N_F1D_lm(l_I1D_lm_.size()), P_F2D_lm_x(l_I1D_lm_.size(), nx_points) {
        assert(nx_points > 0);
        assert(l_I1D_lm_.size() == m_I1D_lm_.size());
        assert((l_I1D_lm_.array() >= 0).all());
        assert((m_I1D_lm_.array() >= 0).all());
        assert((m_I1D_lm_.array() <= l_I1D_lm_.array()).all());
        // 调用 GSL 生成 Gauss-Legendre 格点和权重
        double xlo_F = -1.0;        // cos(θ)格点位置的下界, 对应θ=π
        double xup_F = 1.0;         // cos(θ)格点位置的上界, 对应θ=0
        gsl_integration_fixed_workspace* wt = gsl_integration_fixed_alloc(gsl_integration_fixed_legendre, nx_points, xlo_F, xup_F, 0.0, 0.0);
        const double* nodes = gsl_integration_fixed_nodes(wt);
        const double* weights = gsl_integration_fixed_weights(wt);
        for (int i = 0; i < nx_points; ++i) {x_F1D_x(i) = nodes[i];}
        for (int i = 0; i < nx_points; ++i) {w_F1D_x(i) = weights[i];}
        gsl_integration_fixed_free(wt);
        // 计算广义 Legendre 多项式矩阵
        for (int ix = 0; ix < nx_points; ++ix) for (int ilm = 0; ilm < l_I1D_lm.size(); ++ilm) {
            int l_I = l_I1D_lm_(ilm);
            int m_I = m_I1D_lm_(ilm);
            double x_F = x_F1D_x(ix);
            N_F1D_lm(ilm) = std::sqrt((2 * l_I + 1) / 2.0) * std::exp(0.5 * (gsl_sf_lnfact(l_I - m_I) - gsl_sf_lnfact(l_I + m_I)));
            P_F2D_lm_x(ilm, ix) = gsl_sf_legendre_Plm(l_I, m_I, x_F);
        }
    }
};
