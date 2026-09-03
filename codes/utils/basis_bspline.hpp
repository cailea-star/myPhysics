/**
 * @file    basis_bspline.hpp
 * @author  cailea
 * @date    2026-06-04
 * @brief   Cubic clamped B-spline basis on a finite interval.
 */

#pragma once
#include <algorithm>
#include <cassert>
#include <Eigen/Core>
#include <unsupported/Eigen/Splines>
#include "integration_gauss.hpp"

class BSplineBasisFunction {
public:
    using Spline1D = Eigen::Spline<double, 1>;
    using KnotVectorType = Spline1D::KnotVectorType;
    using BasisDerivativeType = Spline1D::BasisDerivativeType;

    Eigen::VectorXd x_F1D_x{};              // x_i: interval boundaries from x_min to x_max.
    const int p_I = 3;                      // p = 3, order = p + 1.
    KnotVectorType x_F1D_t{};               // t_α: clamped knot vector with endpoint multiplicity p + 1.

    /**
     * @brief  Construct a cubic clamped B-spline basis.
     * @math   N_B = N_x + p - 1
     * @output Initialized knot vector and basis definition.
     */
    BSplineBasisFunction(const Eigen::VectorXd& x_F1D_x_) {
        x_F1D_x = x_F1D_x_;
        fill_clamped_knots();
    }

    /**
     * @brief  Return the number of basis functions.
     * @math   N_B = N_t - p - 1
     * @output Basis count.
     */
    int size() const {
        return static_cast<int>(x_F1D_t.size()) - p_I - 1;
    }

    /**
     * @brief  Evaluate one basis function and its derivatives.
     * @math   B_i^(q)(x) = d^q B_i(x) / dx^q
     * @output Derivative vector from order 0 through q.
     */
    Eigen::VectorXd basis_function(double x_F, int b_I, int qmax_I = 0) const {
        // (x,b,{q}) → {B_b^{(q)}(x)}.
        assert(b_I >= 0 && b_I < size());
        assert(qmax_I >= 0);

        Eigen::VectorXd B_F1D_deriv = Eigen::VectorXd::Zero(qmax_I + 1);
        int span_I = static_cast<int>(Spline1D::Span(x_F, p_I, x_F1D_t));
        int b_local_I = b_I - (span_I - p_I);
        int bSafe_I = std::clamp(b_local_I, 0, p_I);
        double isSupported_F = static_cast<double>((b_local_I >= 0) & (b_local_I <= p_I));
        BasisDerivativeType B_F2D_deriv_b = Spline1D::BasisFunctionDerivatives(x_F, qmax_I, p_I, x_F1D_t);
        for (int deriv_I = 0; deriv_I < B_F2D_deriv_b.rows(); ++deriv_I) {
            B_F1D_deriv(deriv_I) = isSupported_F * B_F2D_deriv_b(deriv_I, bSafe_I);
        }
        return B_F1D_deriv;
    }

    /**
     * @brief  Evaluate a B-spline expansion and its derivatives.
     * @math   u^(q)(x) = Σ_i c_i d^q B_i(x) / dx^q
     * @output Derivative vector from order 0 through q.
     */
    Eigen::VectorXd wave_function(double x_F, const Eigen::VectorXd& coeff_F1D_b, int qmax_I = 0) const {
        // ({c_b},x) → {u^{(q)}(x)}.
        assert(coeff_F1D_b.size() == size());
        assert(qmax_I >= 0);

        Eigen::VectorXd u_F1D_deriv = Eigen::VectorXd::Zero(qmax_I + 1);
        for (int b_I = 0; b_I < size(); ++b_I) {
            u_F1D_deriv += coeff_F1D_b(b_I) * basis_function(x_F, b_I, qmax_I);
        }
        return u_F1D_deriv;
    }

protected:
    /**
     * @brief  Fill the clamped knot vector.
     * @math   t_0 = ... = t_p = x_min, t_{N_t-p-1} = ... = t_{N_t-1} = x_max
     * @output Updates x_F1D_t.
     */
    void fill_clamped_knots() {
        // {x_i} → {t_α}; t_edge repeated p+1.
        assert(x_F1D_x.size() >= 2);
        int Nx_I = static_cast<int>(x_F1D_x.size());
        int Nt_I = Nx_I + 2 * p_I;
        x_F1D_t.resize(Nt_I);

        for (int t_I = 0; t_I <= p_I; ++t_I) {
            x_F1D_t(t_I) = x_F1D_x(0);
        }

        for (int x_I = 1; x_I < Nx_I - 1; ++x_I) {
            x_F1D_t(p_I + x_I) = x_F1D_x(x_I);
        }

        int t_right_I = p_I + Nx_I - 1;
        for (int t_I = 0; t_I <= p_I; ++t_I) {
            x_F1D_t(t_right_I + t_I) = x_F1D_x(Nx_I - 1);
        }
    }

};

class BSplineBasis {
public:
    BSplineBasisFunction b_funcs;        // B_b(x): clamped cubic basis used for grid tabulation.
    Eigen::VectorXd x_F1D_grid{};        // x_g: Gauss-Legendre nodes mapped to each interval.
    Eigen::VectorXd w_F1D_grid{};        // w_g: mapped Gauss-Legendre weights on each interval.
    Eigen::MatrixXd B_F2D_grid_b{};      // B_gb = B_b(x_g).
    Eigen::MatrixXd dB_F2D_grid_b{};     // dB_gb = ∂_xB_b(x_g).

    /**
     * @brief  Build Gauss-Legendre grid and tabulate cubic B-splines.
     * @math   ∫_(-1)^1 f(ξ)dξ ≈ Σ_g ω_g f(ξ_g)
     * @math   x_c = (x_lo+x_up)/2, J = (x_up-x_lo)/2, x_g = x_c+Jξ_g, w_g = Jω_g
     * @math   u(x_min) = u(x_max) = 0 ⇒ b = 1,...,N_B-2
     * @output Quadrature grid, weights, basis values, and derivatives.
     * @note   isZeroBound_B selects the complete or homogeneous-Dirichlet basis space.
     */
    BSplineBasis(const BSplineBasisFunction& b_basis_func_, int nquad_I, bool isZeroBound_B)
    : b_funcs(b_basis_func_) {
        assert(nquad_I > 0);

        // {[x_i,x_{i+1}]} → {x_g,w_g}.
        const int Nseg_I = static_cast<int>(b_basis_func_.x_F1D_x.size()) - 1;
        GaussLegendreMeshes gauss_meshes(nquad_I);
        x_F1D_grid.resize(Nseg_I * nquad_I);
        w_F1D_grid.resize(Nseg_I * nquad_I);

        for (int iseg_I = 0; iseg_I < Nseg_I; ++iseg_I) {
            const double xlo_F = b_basis_func_.x_F1D_x(iseg_I);
            const double xup_F = b_basis_func_.x_F1D_x(iseg_I + 1);
            assert(xup_F > xlo_F);

            const double xc_F = 0.5 * (xlo_F + xup_F);
            const double jac_F = 0.5 * (xup_F - xlo_F);
            const int offset_I = iseg_I * nquad_I;
            x_F1D_grid.segment(offset_I, nquad_I) =
                Eigen::VectorXd::Constant(nquad_I, xc_F) + jac_F * gauss_meshes.x_F1D_x;
            w_F1D_grid.segment(offset_I, nquad_I) = jac_F * gauss_meshes.w_F1D_x;
        }

        // b ∈ [b_offset,b_offset+N_b) → {B_gb,B'_gb}.
        const int Ngrid_I = static_cast<int>(x_F1D_grid.size());
        const int Nb_all_I = b_basis_func_.size();
        assert(!isZeroBound_B || Nb_all_I > 2);
        const int isZeroBound_I = static_cast<int>(isZeroBound_B);
        const int b_offset_I = isZeroBound_I;
        const int Nb_I = Nb_all_I - 2 * isZeroBound_I;
        B_F2D_grid_b.resize(Ngrid_I, Nb_I);
        dB_F2D_grid_b.resize(Ngrid_I, Nb_I);
        for (int g_I = 0; g_I < Ngrid_I; ++g_I) {
            for (int b_I = 0; b_I < Nb_I; ++b_I) {
                const int b_all_I = b_I + b_offset_I;
                const Eigen::VectorXd B_F1D_deriv = b_basis_func_.basis_function(x_F1D_grid(g_I), b_all_I, 1);
                B_F2D_grid_b(g_I, b_I) = B_F1D_deriv(0);
                dB_F2D_grid_b(g_I, b_I) = B_F1D_deriv(1);
            }
        }
    }

    /**
     * @brief  Compute a local-operator matrix in the B-spline basis.
     * @math   M_ij = Σ_g w_g B_i(x_g) O(x_g) B_j(x_g)
     * @output Local-operator matrix M.
     */
    Eigen::MatrixXd calc_B_O_B(const Eigen::VectorXd& O_F1D_grid) const {
        // (w_g,O_g,B_gi) → M_ij.
        assert(w_F1D_grid.size() == O_F1D_grid.size());
        assert(B_F2D_grid_b.rows() == O_F1D_grid.size());

        const Eigen::VectorXd wO_F1D_grid = w_F1D_grid.array() * O_F1D_grid.array();
        return B_F2D_grid_b.transpose() * wO_F1D_grid.asDiagonal() * B_F2D_grid_b;
    }

    /**
     * @brief  Compute a derivative-form matrix in the B-spline basis.
     * @math   K_ij = Σ_g w_g B_i'(x_g) O(x_g) B_j'(x_g)
     * @output Derivative-form matrix K.
     */
    Eigen::MatrixXd calc_dB_O_dB(const Eigen::VectorXd& O_F1D_grid) const {
        // (w_g,O_g,B'_gi) → K_ij.
        assert(w_F1D_grid.size() == O_F1D_grid.size());
        assert(dB_F2D_grid_b.rows() == O_F1D_grid.size());

        const Eigen::VectorXd wO_F1D_grid = w_F1D_grid.array() * O_F1D_grid.array();
        return dB_F2D_grid_b.transpose() * wO_F1D_grid.asDiagonal() * dB_F2D_grid_b;
    }
};
