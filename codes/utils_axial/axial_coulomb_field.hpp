/**
 * @file    axial_coulomb_field.hpp
 * @author  cailea
 * @date    2026-09-03
 * @brief   Build and apply the axial direct-Coulomb kernel.
 */

#pragma once

#include <cassert>
#include <cmath>
#include <iostream>

#include <Eigen/Core>
#include <gsl/gsl_sf_bessel.h>
#include <unsupported/Eigen/CXX11/Tensor>

#include "axial_basis.hpp"
#include "integration_gauss.hpp"

class AxialCoulombField {
public:
    bool isBuilt_B = false;

private:
    Eigen::VectorXd z_F1D_z{};
    Eigen::VectorXd r_F1D_r{};
    Eigen::MatrixXd w_F2D_z_r{};
    Eigen::Tensor<double, 4, Eigen::ColMajor> coulomb_F4D_zs_rs_zt_rt{};

public:
    /**
     * @brief  Construct an empty axial Coulomb field.
     * @math   (z_i,r_j,w_{ij}) → K_C
     * @output Stored grid and unbuilt kernel.
     */
    explicit AxialCoulombField(const AxialBasis& basis_) {
        z_F1D_z = basis_.z_F1D_z;
        r_F1D_r = basis_.r_F1D_r;
        w_F2D_z_r = basis_.w_F2D_z_r;
        coulomb_F4D_zs_rs_zt_rt.resize(z_F1D_z.size(), r_F1D_r.size(), z_F1D_z.size(), r_F1D_r.size());
    }

    /**
     * @brief  Build direct Coulomb by Gauss-Legendre quadrature.
     * @math   K_C(s,t)=2e²/√π∫_0^∞e^{-u²d_{st}²}I_0^e(2r_sr_tu²)du
     * @output Built direct-Coulomb kernel.
     */
    void build(bool useReflection_B, double e2_F);

    /**
     * @brief  Apply the direct-Coulomb kernel.
     * @math   v_C^{dir}(t)=∫d³r_sK_C(s,t)ρ(s)
     * @output Direct-Coulomb field grid.
     */
    Eigen::MatrixXd calc_direct_field(const Eigen::MatrixXd& rho_F2D_z_r) const;
};

inline void AxialCoulombField::build(bool useReflection_B, double e2_F) {
    if (isBuilt_B) {return;}

    constexpr int Nlegendre_I = 80;
    const double pi_F = std::acos(-1.0);
    const int Nz_I = static_cast<int>(z_F1D_z.size());
    const int Nr_I = static_cast<int>(r_F1D_r.size());

    // x ∈ [0,1].
    const GaussLegendreMeshes legendre_meshes(2 * Nlegendre_I);
    const Eigen::VectorXd x_F1D_x = legendre_meshes.x_F1D_x.tail(Nlegendre_I);
    const Eigen::VectorXd weight_F1D_x = legendre_meshes.w_F1D_x.tail(Nlegendre_I);

    // u=x/[b(1-x^β)^α].
    constexpr double b_F = 50.0;
    constexpr double beta_F = 2.0;
    constexpr double alpha_F = 0.5;
    Eigen::VectorXd dudx_F1D_x(Nlegendre_I);
    Eigen::VectorXd u2_F1D_x(Nlegendre_I);
    for (int x_I = 0; x_I < Nlegendre_I; ++x_I) {
        const double x_F = x_F1D_x(x_I);
        const double u_F = x_F / (b_F * std::pow(1.0 - std::pow(x_F, beta_F), alpha_F));
        u2_F1D_x(x_I) = u_F * u_F;
        dudx_F1D_x(x_I) = 1.0 / (b_F * std::pow(1.0 - std::pow(x_F, beta_F), alpha_F + 1.0));
    }

    // w_s contains the 2π azimuthal weight.
    const double factor_F = e2_F * 2.0 / std::sqrt(pi_F);
    const double reflectionWeight_F = static_cast<double>(useReflection_B);
    std::cout << "[AxialCoulombField] Building kernel (nleg=" << Nlegendre_I << ", grid=" << Nz_I << "x" << Nr_I << ")..." << std::endl;

    // (u_x,s,t) → K_C(s,t).
    coulomb_F4D_zs_rs_zt_rt.setZero();
    for (int x_I = 0; x_I < Nlegendre_I; ++x_I) {
        const double u2_F = u2_F1D_x(x_I);
        const double weight_x_F = weight_F1D_x(x_I) * dudx_F1D_x(x_I) * factor_F;
        #pragma omp parallel for collapse(4) schedule(static)
        for (int rt_I = 0; rt_I < Nr_I; ++rt_I) {
            for (int zt_I = 0; zt_I < Nz_I; ++zt_I) {
                for (int rs_I = 0; rs_I < Nr_I; ++rs_I) {
                    for (int zs_I = 0; zs_I < Nz_I; ++zs_I) {
                        const double rs_F = r_F1D_r(rs_I);
                        const double zs_F = z_F1D_z(zs_I);
                        const double rt_F = r_F1D_r(rt_I);
                        const double zt_F = z_F1D_z(zt_I);
                        const double distance2_F = (rs_F - rt_F) * (rs_F - rt_F) + (zs_F - zt_F) * (zs_F - zt_F);
                        const double distanceReflection2_F = (rs_F - rt_F) * (rs_F - rt_F) + (zs_F + zt_F) * (zs_F + zt_F);
                        const double expDirect_F = std::exp(-u2_F * distance2_F);
                        const double expReflection_F = std::exp(-u2_F * distanceReflection2_F);
                        const double expKernel_F = (1.0 - 0.5 * reflectionWeight_F) * expDirect_F + 0.5 * reflectionWeight_F * expReflection_F;
                        const double I0Scaled_F = gsl_sf_bessel_I0_scaled(2.0 * rs_F * rt_F * u2_F);
                        coulomb_F4D_zs_rs_zt_rt(zs_I, rs_I, zt_I, rt_I) += weight_x_F * expKernel_F * I0Scaled_F;
                    }
                }
            }
        }
    }
    isBuilt_B = true;
    std::cout << "[AxialCoulombField] Kernel built." << std::endl;
}

inline Eigen::MatrixXd AxialCoulombField::calc_direct_field(const Eigen::MatrixXd& rho_F2D_z_r) const {
    assert(isBuilt_B);
    assert(rho_F2D_z_r.rows() == z_F1D_z.size());
    assert(rho_F2D_z_r.cols() == r_F1D_r.size());
    const int Nz_I = static_cast<int>(z_F1D_z.size());
    const int Nr_I = static_cast<int>(r_F1D_r.size());
    const int Ngrid_I = Nz_I * Nr_I;
    const Eigen::Map<const Eigen::VectorXd> rho_F1D_zr(rho_F2D_z_r.data(), Ngrid_I);
    const Eigen::Map<const Eigen::VectorXd> weight_F1D_zr(w_F2D_z_r.data(), Ngrid_I);
    const Eigen::VectorXd rhoWeighted_F1D_zr = weight_F1D_zr.array() * rho_F1D_zr.array();
    const Eigen::Map<const Eigen::MatrixXd> coulomb_F2D_source_target(coulomb_F4D_zs_rs_zt_rt.data(), Ngrid_I, Ngrid_I);
    const Eigen::VectorXd vCoulombDirect_F1D_zr = coulomb_F2D_source_target.transpose() * rhoWeighted_F1D_zr;
    const Eigen::Map<const Eigen::MatrixXd> vCoulombDirect_F2D_z_r(vCoulombDirect_F1D_zr.data(), Nz_I, Nr_I);
    return Eigen::MatrixXd(vCoulombDirect_F2D_z_r);
}
