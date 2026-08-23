/**
 * @file    bvp_rmatrix.hpp
 * @author  cailea
 * @date    2026-06-26
 * @brief   B-spline R-matrix boundary-response tools.
 */

#pragma once
#include <cassert>
#include <complex>
#include <functional>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include "basis_bspline.hpp"

using Real2CMatFunc = std::function<Eigen::MatrixXcd(double)>;

class RMatrix {
private:
    const BSplineBasis& b_basis;
    Eigen::Tensor<std::complex<double>, 4, Eigen::ColMajor> Ginv_C4D_ch_b_ch_b;

public:
    /**
     * @brief  Construct an R-matrix response by B-spline Galerkin assembly.
     * @math   G⁻¹ = H_B - E
     * @output Initialized inverse Green matrix.
     * @note   Requires the complete endpoint-retaining B-spline basis.
     */
    RMatrix(const BSplineBasis& b_basis_, const Real2CMatFunc& F_Func, const Eigen::Ref<const Eigen::VectorXcd>& kxmin_C1D_ch, const Eigen::Ref<const Eigen::VectorXcd>& kxmax_C1D_ch)
    : b_basis(b_basis_) {
        update_Ginv(F_Func, kxmin_C1D_ch, kxmax_C1D_ch);
    }

    /**
     * @brief  Assemble the inverse Green matrix by B-spline Galerkin quadrature.
     * @math   G⁻¹_αi,βj = K_ij δ_αβ + ∫ B_i F_αβ B_j dx - ik_xmin,α δ_i0 δ_j0 - ik_xmax,α δ_iN δ_jN
     * @output Updates Ginv_C4D_ch_b_ch_b.
     * @note   h_μ = 1.
     * @note   Requires endpoint-retaining B-splines.
     */
    void update_Ginv(const Real2CMatFunc& F_Func, const Eigen::Ref<const Eigen::VectorXcd>& kxmin_C1D_ch, const Eigen::Ref<const Eigen::VectorXcd>& kxmax_C1D_ch);

    /**
     * @brief  Evaluate Green matrix by Galerkin contraction and full-pivoting LU.
     * @math   G_αβ(x₁,x₂) = Σ_ij B_i(x₁) G_αi,βj B_j(x₂)
     * @output Coordinate-space Green matrix.
     * @note   G equals R for h_μ = 1.
     */
    Eigen::MatrixXcd calc_green_function(double x1_F, double x2_F) const;

    /**
     * @brief  Compute S by R-matrix matching and full-pivoting LU.
     * @math   S = (O - RO′)⁻¹(I - RI′)
     * @output Channel-space scattering matrix.
     * @note   Asymptotic-wave matrices must be channel-diagonal.
     */
    Eigen::MatrixXcd calc_S_matrix(const Eigen::Ref<const Eigen::VectorXcd>& O_C1D_ch, const Eigen::Ref<const Eigen::VectorXcd>& dO_C1D_ch, const Eigen::Ref<const Eigen::VectorXcd>& I_C1D_ch, const Eigen::Ref<const Eigen::VectorXcd>& dI_C1D_ch) const;
};

inline void RMatrix::update_Ginv(const Real2CMatFunc& F_Func, const Eigen::Ref<const Eigen::VectorXcd>& kxmin_C1D_ch, const Eigen::Ref<const Eigen::VectorXcd>& kxmax_C1D_ch) {
    // (N_ch,N_b,N_grid) → validated dimensions.
    const int Nch_I = static_cast<int>(kxmin_C1D_ch.size());
    const int Nb_I = static_cast<int>(b_basis.B_F2D_grid_b.cols());
    const int Ngrid_I = static_cast<int>(b_basis.x_F1D_grid.size());
    assert(Nch_I > 0);
    assert(Nb_I > 0);
    assert(kxmax_C1D_ch.size() == Nch_I);
    assert(b_basis.w_F1D_grid.size() == Ngrid_I);
    assert(b_basis.B_F2D_grid_b.rows() == Ngrid_I);
    assert(b_basis.dB_F2D_grid_b.rows() == Ngrid_I);
    assert(Nb_I == b_basis.b_funcs.size());

    Ginv_C4D_ch_b_ch_b.resize(Nch_I, Nb_I, Nch_I, Nb_I);
    Ginv_C4D_ch_b_ch_b.setZero();

    // K_ij δ_αβ → G⁻¹_αi,βj.
    const Eigen::VectorXd one_F1D_grid = Eigen::VectorXd::Ones(Ngrid_I);
    const Eigen::MatrixXd T_F2D_b_b = b_basis.calc_dB_O_dB(one_F1D_grid);
    for (int ch_I = 0; ch_I < Nch_I; ++ch_I) {
        for (int b1_I = 0; b1_I < Nb_I; ++b1_I) {
            for (int b2_I = 0; b2_I < Nb_I; ++b2_I) {Ginv_C4D_ch_b_ch_b(ch_I, b1_I, ch_I, b2_I) += T_F2D_b_b(b1_I, b2_I);}
        }
    }

    // Σ_g w_g B_i(x_g)F_αβ(x_g)B_j(x_g) → G⁻¹_αi,βj.
    Eigen::MatrixXcd F_C2D_ch_ch(Nch_I, Nch_I);
    for (int grid_I = 0; grid_I < Ngrid_I; ++grid_I) {
        F_C2D_ch_ch = F_Func(b_basis.x_F1D_grid(grid_I));
        assert(F_C2D_ch_ch.rows() == Nch_I);
        assert(F_C2D_ch_ch.cols() == Nch_I);
        for (int b1_I = 0; b1_I < Nb_I; ++b1_I) {
            for (int b2_I = 0; b2_I < Nb_I; ++b2_I) {
                const double wBB_F = b_basis.w_F1D_grid(grid_I) * b_basis.B_F2D_grid_b(grid_I, b1_I) * b_basis.B_F2D_grid_b(grid_I, b2_I);
                for (int ch1_I = 0; ch1_I < Nch_I; ++ch1_I) {
                    for (int ch2_I = 0; ch2_I < Nch_I; ++ch2_I) {Ginv_C4D_ch_b_ch_b(ch1_I, b1_I, ch2_I, b2_I) += wBB_F * F_C2D_ch_ch(ch1_I, ch2_I);}
                }
            }
        }
    }

    // (Σ_min,Σ_max) → (-ik_min,-ik_max).
    const std::complex<double> minusI_C(0.0, -1.0);
    for (int ch_I = 0; ch_I < Nch_I; ++ch_I) {
        Ginv_C4D_ch_b_ch_b(ch_I, 0, ch_I, 0) += minusI_C * kxmin_C1D_ch(ch_I);
        Ginv_C4D_ch_b_ch_b(ch_I, Nb_I - 1, ch_I, Nb_I - 1) += minusI_C * kxmax_C1D_ch(ch_I);
    }
}

inline Eigen::MatrixXcd RMatrix::calc_green_function(double x1_F, double x2_F) const {
    // (G⁻¹,N_ch,N_b) → validated dimensions.
    const int Nch_I = static_cast<int>(Ginv_C4D_ch_b_ch_b.dimension(0));
    const int Nb_I = static_cast<int>(b_basis.B_F2D_grid_b.cols());
    assert(Nch_I > 0);
    assert(Nb_I == b_basis.b_funcs.size());
    assert(Ginv_C4D_ch_b_ch_b.dimension(1) == Nb_I);
    assert(Ginv_C4D_ch_b_ch_b.dimension(2) == Nch_I);
    assert(Ginv_C4D_ch_b_ch_b.dimension(3) == Nb_I);

    // (x_1,x_2) → ({B_i(x_1)},{B_i(x_2)}).
    Eigen::VectorXd Bx1_F1D_b(Nb_I);
    Eigen::VectorXd Bx2_F1D_b(Nb_I);
    for (int b_I = 0; b_I < Nb_I; ++b_I) {
        Bx1_F1D_b(b_I) = b_basis.b_funcs.basis_function(x1_F, b_I)(0);
        Bx2_F1D_b(b_I) = b_basis.b_funcs.basis_function(x2_F, b_I)(0);
    }

    // G⁻¹ X = B(x_2) → G(x_1,x_2).
    const int Nchb_I = Nch_I * Nb_I;
    Eigen::Map<const Eigen::MatrixXcd> Ginv_C2D_chb_chb(Ginv_C4D_ch_b_ch_b.data(), Nchb_I, Nchb_I);
    Eigen::MatrixXcd Bx2_C2D_chb_ch = Eigen::MatrixXcd::Zero(Nchb_I, Nch_I);
    for (int ch_I = 0; ch_I < Nch_I; ++ch_I) {
        for (int b_I = 0; b_I < Nb_I; ++b_I) {Bx2_C2D_chb_ch(ch_I + Nch_I * b_I, ch_I) = Bx2_F1D_b(b_I);}
    }

    const Eigen::MatrixXcd Gx2_C2D_chb_ch = Ginv_C2D_chb_chb.fullPivLu().solve(Bx2_C2D_chb_ch);
    Eigen::MatrixXcd G_C2D_ch_ch = Eigen::MatrixXcd::Zero(Nch_I, Nch_I);
    for (int ch_I = 0; ch_I < Nch_I; ++ch_I) {
        for (int b_I = 0; b_I < Nb_I; ++b_I) {G_C2D_ch_ch.row(ch_I) += Bx1_F1D_b(b_I) * Gx2_C2D_chb_ch.row(ch_I + Nch_I * b_I);}
    }
    return G_C2D_ch_ch;
}

inline Eigen::MatrixXcd RMatrix::calc_S_matrix(const Eigen::Ref<const Eigen::VectorXcd>& O_C1D_ch, const Eigen::Ref<const Eigen::VectorXcd>& dO_C1D_ch, const Eigen::Ref<const Eigen::VectorXcd>& I_C1D_ch, const Eigen::Ref<const Eigen::VectorXcd>& dI_C1D_ch) const {
    // (O,O',I,I') → validated channel vectors.
    const int Nch_I = static_cast<int>(Ginv_C4D_ch_b_ch_b.dimension(0));
    const Eigen::VectorXd& x_F1D_x = b_basis.b_funcs.x_F1D_x;
    assert(O_C1D_ch.size() == Nch_I);
    assert(dO_C1D_ch.size() == Nch_I);
    assert(I_C1D_ch.size() == Nch_I);
    assert(dI_C1D_ch.size() == Nch_I);
    assert(x_F1D_x.size() >= 2);

    // (R,O,O',I,I') → (O-RO',I-RI') → S.
    const double xmax_F = x_F1D_x(x_F1D_x.size() - 1);
    const Eigen::MatrixXcd R_C2D_ch_ch = calc_green_function(xmax_F, xmax_F);
    Eigen::MatrixXcd Zout_C2D_ch_ch = -R_C2D_ch_ch * dO_C1D_ch.asDiagonal();
    Eigen::MatrixXcd Zin_C2D_ch_ch = -R_C2D_ch_ch * dI_C1D_ch.asDiagonal();
    Zout_C2D_ch_ch.diagonal() += O_C1D_ch;
    Zin_C2D_ch_ch.diagonal() += I_C1D_ch;
    return Zout_C2D_ch_ch.fullPivLu().solve(Zin_C2D_ch_ch);
}
