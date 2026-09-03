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
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <unsupported/Eigen/CXX11/Tensor>
#include "basis_bspline.hpp"

using Real2CMatFunc = std::function<Eigen::MatrixXcd(double)>;

class RMatrix {
public:
    const BSplineBasis& b_basis;

private:
    Eigen::Tensor<std::complex<double>, 4, Eigen::ColMajor> Ginv_C4D_ch_b_ch_b{};

public:
    /**
     * @brief  Construct a fixed-size B-spline R-matrix workspace.
     * @math   G⁻¹ ∈ ℂ^(N_chN_b)×(N_chN_b)
     * @output Zeroed R-matrix workspace.
     * @note   Call update_Ginv before matrix queries.
     */
    RMatrix(const BSplineBasis& b_basis_, int Nch_I_)
    : b_basis(b_basis_) {
        assert(Nch_I_ > 0);
        const int Nb_I = static_cast<int>(b_basis.B_F2D_grid_b.cols());
        assert(Nb_I > 1);
        Ginv_C4D_ch_b_ch_b.resize(Nch_I_, Nb_I, Nch_I_, Nb_I);
        Ginv_C4D_ch_b_ch_b.setZero();
    }

    /**
     * @brief  Construct an R-matrix response by B-spline Galerkin assembly.
     * @math   G⁻¹ = H_B - E
     * @output Initialized inverse Green matrix.
     * @note   Requires the complete endpoint-retaining B-spline basis.
     */
    RMatrix(const BSplineBasis& b_basis_, const Real2CMatFunc& F_Func, const Eigen::Ref<const Eigen::VectorXcd>& lambdaIn_C1D_ch)
    : RMatrix(b_basis_, static_cast<int>(lambdaIn_C1D_ch.size())) {
        update_Ginv(F_Func, lambdaIn_C1D_ch);
    }

    /**
     * @brief  Assemble the inverse Green matrix by B-spline Galerkin quadrature.
     * @math   G⁻¹_αi,βj = K_ijδ_αβ + ∫B_iF_αβB_jdx + Λ_in,αδ_i0δ_j0
     * @output Read-only reference to Ginv_C4D_ch_b_ch_b.
     * @note   h_μ = 1.
     * @note   Requires endpoint-retaining B-splines.
     */
    const Eigen::Tensor<std::complex<double>, 4, Eigen::ColMajor>& update_Ginv(const Real2CMatFunc& F_Func, const Eigen::Ref<const Eigen::VectorXcd>& lambdaIn_C1D_ch);

    /**
     * @brief  Compute the inverse outer-boundary R matrix by Schur complement.
     * @math   R⁻¹=A_BB-A_BI A_II⁻¹ A_IB
     * @output Channel-space inverse boundary response.
     * @note   The outer endpoint basis block is the last channel block.
     */
    Eigen::MatrixXcd calc_Rinv_matrix() const;

    /**
     * @brief  Compute S from the inverse R matrix and full-pivoting LU.
     * @math   S = (R⁻¹O - O′)⁻¹(R⁻¹I - I′)
     * @output Channel-space scattering matrix.
     * @note   Asymptotic-wave matrices must be channel-diagonal.
     */
    Eigen::MatrixXcd calc_S_matrix(const Eigen::Ref<const Eigen::VectorXcd>& O_C1D_ch, const Eigen::Ref<const Eigen::VectorXcd>& dO_C1D_ch, const Eigen::Ref<const Eigen::VectorXcd>& I_C1D_ch, const Eigen::Ref<const Eigen::VectorXcd>& dI_C1D_ch) const;
};

inline const Eigen::Tensor<std::complex<double>, 4, Eigen::ColMajor>& RMatrix::update_Ginv(const Real2CMatFunc& F_Func, const Eigen::Ref<const Eigen::VectorXcd>& lambdaIn_C1D_ch) {
    // (N_ch,N_b,N_grid) → validated dimensions.
    const int Nch_I = static_cast<int>(Ginv_C4D_ch_b_ch_b.dimension(0));
    const int Nb_I = static_cast<int>(b_basis.B_F2D_grid_b.cols());
    const int Ngrid_I = static_cast<int>(b_basis.x_F1D_grid.size());
    assert(Nch_I > 0);
    assert(Nb_I > 0);
    assert(lambdaIn_C1D_ch.size() == Nch_I);
    assert(b_basis.w_F1D_grid.size() == Ngrid_I);
    assert(b_basis.B_F2D_grid_b.rows() == Ngrid_I);
    assert(b_basis.dB_F2D_grid_b.rows() == Ngrid_I);
    assert(Nb_I == b_basis.b_funcs.size());
    assert(Ginv_C4D_ch_b_ch_b.dimension(1) == Nb_I);
    assert(Ginv_C4D_ch_b_ch_b.dimension(2) == Nch_I);
    assert(Ginv_C4D_ch_b_ch_b.dimension(3) == Nb_I);

    Ginv_C4D_ch_b_ch_b.setZero();

    // Σ_g w_g[B'_iB'_jδ_αβ+B_iF_αβB_j] → G⁻¹_αi,βj.
    Eigen::MatrixXcd F_C2D_ch_ch(Nch_I, Nch_I);
    for (int grid_I = 0; grid_I < Ngrid_I; ++grid_I) {
        F_C2D_ch_ch = F_Func(b_basis.x_F1D_grid(grid_I));
        assert(F_C2D_ch_ch.rows() == Nch_I);
        assert(F_C2D_ch_ch.cols() == Nch_I);
        const int span_I = static_cast<int>(BSplineBasisFunction::Spline1D::Span(b_basis.x_F1D_grid(grid_I), b_basis.b_funcs.p_I, b_basis.b_funcs.x_F1D_t));
        const int bmin_I = span_I - b_basis.b_funcs.p_I;
        const int bmax_I = span_I;
        assert(bmin_I >= 0 && bmax_I < Nb_I);
        for (int b1_I = bmin_I; b1_I <= bmax_I; ++b1_I) {
            for (int b2_I = bmin_I; b2_I <= bmax_I; ++b2_I) {
                const double wBB_F = b_basis.w_F1D_grid(grid_I) * b_basis.B_F2D_grid_b(grid_I, b1_I) * b_basis.B_F2D_grid_b(grid_I, b2_I);
                const double dBwdB_F = b_basis.w_F1D_grid(grid_I) * b_basis.dB_F2D_grid_b(grid_I, b1_I) * b_basis.dB_F2D_grid_b(grid_I, b2_I);
                for (int ch1_I = 0; ch1_I < Nch_I; ++ch1_I) {
                    for (int ch2_I = 0; ch2_I < Nch_I; ++ch2_I) {Ginv_C4D_ch_b_ch_b(ch1_I, b1_I, ch2_I, b2_I) += wBB_F * F_C2D_ch_ch(ch1_I, ch2_I);}
                    Ginv_C4D_ch_b_ch_b(ch1_I, b1_I, ch1_I, b2_I) += dBwdB_F;
                }
            }
        }
    }

    // Σ_in → Λ_in.
    for (int ch_I = 0; ch_I < Nch_I; ++ch_I) {Ginv_C4D_ch_b_ch_b(ch_I, 0, ch_I, 0) += lambdaIn_C1D_ch(ch_I);}
    return Ginv_C4D_ch_b_ch_b;
}

inline Eigen::MatrixXcd RMatrix::calc_Rinv_matrix() const {
    // G_Σ⁻¹ → (G⁻¹_II,G⁻¹_IB,G⁻¹_BI,G⁻¹_BB).
    const int Nch_I = static_cast<int>(Ginv_C4D_ch_b_ch_b.dimension(0));
    const int Nb_I = static_cast<int>(b_basis.B_F2D_grid_b.cols());
    const int Nchb_I = Nch_I * Nb_I;
    const int Ninternal_I = Nch_I * (Nb_I - 1);
    assert(Nch_I > 0 && Nb_I > 1);
    assert(Nb_I == b_basis.b_funcs.size());
    assert(Ginv_C4D_ch_b_ch_b.dimension(1) == Nb_I);
    assert(Ginv_C4D_ch_b_ch_b.dimension(2) == Nch_I);
    assert(Ginv_C4D_ch_b_ch_b.dimension(3) == Nb_I);

    Eigen::Map<const Eigen::MatrixXcd> Ginv_C2D_chb_chb(Ginv_C4D_ch_b_ch_b.data(), Nchb_I, Nchb_I);
    const auto GinvII_C2D_chbi_chbi = Ginv_C2D_chb_chb.topLeftCorner(Ninternal_I, Ninternal_I);
    const auto GinvIB_C2D_chbi_ch = Ginv_C2D_chb_chb.topRightCorner(Ninternal_I, Nch_I);
    const auto GinvBI_C2D_ch_chbi = Ginv_C2D_chb_chb.bottomLeftCorner(Nch_I, Ninternal_I);
    const auto GinvBB_C2D_ch_ch = Ginv_C2D_chb_chb.bottomRightCorner(Nch_I, Nch_I);
    const Eigen::SparseMatrix<std::complex<double>> GinvII_C2D_chbi_chbi_sparse = GinvII_C2D_chbi_chbi.sparseView();
    Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double>>> internalSolver_SparseLU;
    internalSolver_SparseLU.compute(GinvII_C2D_chbi_chbi_sparse);
    assert(internalSolver_SparseLU.info() == Eigen::Success);
    const Eigen::MatrixXcd elimination_C2D_chbi_ch = internalSolver_SparseLU.solve(GinvIB_C2D_chbi_ch);
    assert(internalSolver_SparseLU.info() == Eigen::Success);
    Eigen::MatrixXcd Rinv_C2D_ch_ch = GinvBB_C2D_ch_ch;
    Rinv_C2D_ch_ch.noalias() -= GinvBI_C2D_ch_chbi * elimination_C2D_chbi_ch;
    assert(Rinv_C2D_ch_ch.allFinite());
    return Rinv_C2D_ch_ch;
}

inline Eigen::MatrixXcd RMatrix::calc_S_matrix(const Eigen::Ref<const Eigen::VectorXcd>& O_C1D_ch, const Eigen::Ref<const Eigen::VectorXcd>& dO_C1D_ch, const Eigen::Ref<const Eigen::VectorXcd>& I_C1D_ch, const Eigen::Ref<const Eigen::VectorXcd>& dI_C1D_ch) const {
    // (O,O',I,I') → validated channel vectors.
    const int Nch_I = static_cast<int>(Ginv_C4D_ch_b_ch_b.dimension(0));
    assert(O_C1D_ch.size() == Nch_I);
    assert(dO_C1D_ch.size() == Nch_I);
    assert(I_C1D_ch.size() == Nch_I);
    assert(dI_C1D_ch.size() == Nch_I);

    // (R⁻¹,O,O′,I,I′) → (R⁻¹O-O′,R⁻¹I-I′) → S.
    const Eigen::MatrixXcd Rinv_C2D_ch_ch = calc_Rinv_matrix();
    Eigen::MatrixXcd Zout_C2D_ch_ch = Rinv_C2D_ch_ch * O_C1D_ch.asDiagonal();
    Eigen::MatrixXcd Zin_C2D_ch_ch = Rinv_C2D_ch_ch * I_C1D_ch.asDiagonal();
    Zout_C2D_ch_ch.diagonal() -= dO_C1D_ch;
    Zin_C2D_ch_ch.diagonal() -= dI_C1D_ch;
    return Zout_C2D_ch_ch.fullPivLu().solve(Zin_C2D_ch_ch);
}
