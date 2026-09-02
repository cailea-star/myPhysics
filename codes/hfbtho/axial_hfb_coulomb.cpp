/**
 * @file    axial_hfb_coulomb.cpp
 * @author  cailea
 * @date    2026-05-06
 * @brief   Calculate axial direct and Slater-exchange Coulomb fields.
 */

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

#include <gsl/gsl_sf_bessel.h>

#include "axial_hfb.hpp"
#include "integration_gauss.hpp"

namespace {

constexpr int Nlegendre_I = 80;
const double pi_F = std::acos(-1.0);

} // namespace

/**
 * @brief  Build direct Coulomb kernel by Gauss-Legendre quadrature.
 * @math   K_C(s,t)=Σ_x w_sw_xe^{-u_x²d_{st}²}I_0^e(2r_sr_tu_x²)
 * @output Updated Coulomb kernel and readiness flag.
 */
void AxialHFBFieldSystem::build_coulomb_kernel(const AxialBasis& basis_, const EDFParamsSkyrme& edf_skyrme_, bool useReflection_B) {
    if (isCoulombKernelBuilt_B) {return;}

    const int Nz_I = static_cast<int>(field_n.vcent_F2D_z_r.rows());
    const int Nr_I = static_cast<int>(field_n.vcent_F2D_z_r.cols());
    const Eigen::VectorXd& r_F1D_r = basis_.r_F1D_r;
    const Eigen::VectorXd& z_F1D_z = basis_.z_F1D_z;
    const Eigen::MatrixXd& w_F2D_z_r = basis_.w_F2D_z_r;

    // x ∈ [0,1].
    const GaussLegendreMeshes legendre_meshes(2 * Nlegendre_I);
    const Eigen::VectorXd x_F1D_x = legendre_meshes.x_F1D_x.tail(Nlegendre_I);
    const Eigen::VectorXd w_F1D_x = legendre_meshes.w_F1D_x.tail(Nlegendre_I);

    // u=x/[b(1-x^β)^α].
    constexpr double b_F = 50.0;
    constexpr double beta_F = 2.0;
    constexpr double alpha_F = 0.5;
    Eigen::VectorXd dudx_F1D_x(Nlegendre_I);
    Eigen::VectorXd u2_F1D_x(Nlegendre_I);
    for (int x_I = 0; x_I < Nlegendre_I; ++x_I) {
        const double x_F = x_F1D_x(x_I);
        const double u_F = x_F / (b_F * std::pow(1.0 - std::pow(x_F, beta_F), alpha_F));
        const double u2_F = u_F * u_F;
        const double dudx_F = 1.0 / (b_F * std::pow(1.0 - std::pow(x_F, beta_F), alpha_F + 1.0));
        u2_F1D_x(x_I) = u2_F;
        dudx_F1D_x(x_I) = dudx_F;
    }

    // w_s contains the 2π azimuthal weight.
    const double factor_F = edf_skyrme_.e2charg_F * 2.0 / std::sqrt(pi_F);
    const double reflectionWeight_F = static_cast<double>(useReflection_B);
    std::cout << "[AxialHFBFieldSystem] Building Coulomb kernel (nleg=" << Nlegendre_I << ", grid=" << Nz_I << "x" << Nr_I << ")..." << std::endl;

    // (u_x,s,t) → K_C(s,t).
    coulomb_F4D_zs_rs_zt_rt.setZero();
    for (int x_I = 0; x_I < Nlegendre_I; ++x_I) {
        const double u2_F = u2_F1D_x(x_I);
        const double dudx_F = dudx_F1D_x(x_I);
        const double weight_x_F = w_F1D_x(x_I) * dudx_F * factor_F;
        for (int rt_I = 0; rt_I < Nr_I; ++rt_I) {
            for (int zt_I = 0; zt_I < Nz_I; ++zt_I) {
                for (int rs_I = 0; rs_I < Nr_I; ++rs_I) {
                    for (int zs_I = 0; zs_I < Nz_I; ++zs_I) {
                        const double weight_rz1_F = w_F2D_z_r(zs_I, rs_I);
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
                        coulomb_F4D_zs_rs_zt_rt(zs_I, rs_I, zt_I, rt_I) += weight_rz1_F * weight_x_F * expKernel_F * I0Scaled_F;
                    }
                }
            }
        }
    }
    isCoulombKernelBuilt_B = true;
    std::cout << "[AxialHFBFieldSystem] Coulomb kernel built." << std::endl;
}

/**
 * @brief  Add direct and Slater-exchange Coulomb fields.
 * @math   v_C=K_Cρ_p-e²C_{ex}(3ρ_p/π)^{1/3}
 * @output Updated proton central field.
 */
void AxialHFBFieldSystem::add_coulomb_field(const AxialHFBDensitySystem& densities_, const EDFParamsSkyrme& edf_skyrme_, const HFBSettings& hfbsettings_) {
    if (!hfbsettings_.termSwitches.addLocalCoulomb_B) {return;}
    if (!isCoulombKernelBuilt_B) {
        throw std::runtime_error("Coulomb kernel is requested before AxialHFBFieldSystem::build_coulomb_kernel().");
    }

    const AxialHFBDensity& density_p_ = densities_.density_p;
    const int Nz_I = static_cast<int>(field_p.vcent_F2D_z_r.rows());
    const int Nr_I = static_cast<int>(field_p.vcent_F2D_z_r.cols());
    const int Nzr_I = Nz_I * Nr_I;

    // ρ_p → v_C^{dir}.
    const Eigen::Map<const Eigen::VectorXd> rho_p_F1D_zr(density_p_.rho_F2D_z_r.data(), Nzr_I);
    const Eigen::Map<const Eigen::MatrixXd> coulomb_F2D_zrs_zrt(coulomb_F4D_zs_rs_zt_rt.data(), Nzr_I, Nzr_I);
    const Eigen::VectorXd vCoulombDirect_F1D_zrt = coulomb_F2D_zrs_zrt.transpose() * rho_p_F1D_zr;
    const Eigen::Map<const Eigen::MatrixXd> vCoulombDirect_F2D_z_r(vCoulombDirect_F1D_zrt.data(), Nz_I, Nr_I);
    const double coex_F = -edf_skyrme_.e2charg_F * std::cbrt(3.0 / pi_F);

    // v_C=v_C^{dir}+v_C^{Slater}.
    #pragma omp parallel for collapse(2) schedule(static)
    for (int r_I = 0; r_I < Nr_I; ++r_I) {
        for (int z_I = 0; z_I < Nz_I; ++z_I) {
            const double rhoPositive_p_F = std::max(density_p_.rho_F2D_z_r(z_I, r_I), 0.0);
            const double vCoulomb_F = vCoulombDirect_F2D_z_r(z_I, r_I) + edf_skyrme_.CExPar_F * coex_F * std::cbrt(rhoPositive_p_F);
            field_p.vcent_F2D_z_r(z_I, r_I) += vCoulomb_F;
        }
    }
}
