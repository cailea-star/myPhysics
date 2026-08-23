/**
 * @file    cc0_potentials.hpp
 * @author  cailea
 * @date    2026-08-22
 * @brief   Coupled-channels nuclear, Coulomb, and centrifugal potentials.
 */

#pragma once

#include <cassert>
#include <cmath>
#include <Eigen/Core>

inline constexpr double mc2_F = 938.0;             // m_N c² [MeV].
inline constexpr double hbarc_F = 197.329;         // ℏc [MeV fm].
inline constexpr double e2_F = hbarc_F / 137.0;    // e² [MeV fm].

/**
 * @brief  Calculate the reference nuclear-potential depth by an empirical formula.
 * @math   V_{0,ref} = (51 - 33 I_t) A_p, I_t = (A_t - 2 Z_t) / A_t
 * @output Reference nuclear-potential depth V_{0,ref} in MeV.
 * @note   K. Heyde, Basic Ideas and Concepts, Sec. 9.5.
 */
inline double calc_V0ref(int At_I, int Zt_I, int Ap_I) {
    assert(At_I > 0 && Zt_I >= 0 && Zt_I <= At_I && Ap_I > 0);
    double It_F = static_cast<double>(At_I - 2 * Zt_I) / static_cast<double>(At_I);
    double V0ref_F = (51.0 - 33.0 * It_F) * static_cast<double>(Ap_I);
    return V0ref_F;
}

/**
 * @brief  Calculate the kinetic coefficient by the equal-nucleon-mass approximation.
 * @math   h_μ = (ℏc)²(A_t + A_p)/(2 A_t A_p m_N c²)
 * @output Kinetic coefficient h_μ in MeV fm².
 * @note   Equal neutron and proton masses.
 */
inline double calc_hmass(int At_I, int Ap_I) {
    assert(At_I > 0 && Ap_I > 0);
    double muc2_F = static_cast<double>(At_I) * static_cast<double>(Ap_I) / static_cast<double>(At_I + Ap_I) * mc2_F;
    double hmass_F = hbarc_F * hbarc_F / (2.0 * muc2_F);
    return hmass_F;
}

/**
 * @brief  Calculate the Woods-Saxon diffuseness by an empirical formula.
 * @math   a₀ = 0.50 + 0.33 I, I = (A - 2 Z) / A
 * @output Woods-Saxon diffuseness a₀ in fm.
 */
inline double calc_a0(int A_I, int Z_I) {
    assert(A_I > 0 && Z_I >= 0 && Z_I <= A_I);
    double I_F = static_cast<double>(A_I - 2 * Z_I) / static_cast<double>(A_I);
    double a0_F = 0.50 + 0.33 * I_F;
    return a0_F;
}

/**
 * @brief  Calculate the nuclear coupling radius by an empirical formula.
 * @math   R₀ = (1.00 + 0.39 I) A^(1/3), R₀(α) = 1.17 fm
 * @output Nuclear coupling radius R₀ in fm.
 * @note   The alpha-particle radius is an explicitly approved runtime-branch exception.
 */
inline double calc_R0(int A_I, int Z_I) {
    assert(A_I > 0 && Z_I >= 0 && Z_I <= A_I);
    if (A_I == 4 && Z_I == 2) {return 1.17;} // Alpha-particle radius in fm.
    double I_F = static_cast<double>(A_I - 2 * Z_I) / static_cast<double>(A_I);
    double r0_F = 1.00 + 0.39 * I_F;
    double R0_F = r0_F * std::pow(static_cast<double>(A_I), 1.0 / 3.0);
    return R0_F;
}

/**
 * @brief  Calculate the Coulomb radius by an empirical formula.
 * @math   R_C = 1.20 A^(1/3)
 * @output Coulomb radius R_C in fm.
 */
inline double calc_Rcoul(int A_I) {
    assert(A_I > 0);
    double Rcoul_F = 1.20 * std::pow(static_cast<double>(A_I), 1.0 / 3.0);
    return Rcoul_F;
}

class VARG {
public:
    double hmass_F;       // h_μ = (ℏc)²/(2μc²) [MeV fm²].
    double Ze2_F;         // Coulomb coupling Z_t Z_p e² [MeV fm].
    double V0_F;          // Nuclear-potential depth V₀ in MeV.
    double a0_F;          // Woods-Saxon diffuseness a₀ in fm.
    double Rp0_F;         // Projectile nuclear radius R_{p0} [fm].
    double Rt0_F;         // Target nuclear radius R_{t0} [fm].
    double Rcoul_p_F;     // Projectile Coulomb radius R_{Cp} [fm].
    double Rcoul_t_F;     // Target Coulomb radius R_{Ct} [fm].
    double beta2t_F;      // Target quadrupole deformation β_{2t}.
    double beta4t_F;      // Target hexadecapole deformation β_{4t}.
    double beta2tcoul_F;  // Target Coulomb quadrupole deformation β_{2t,C}.
    double beta4tcoul_F;  // Target Coulomb hexadecapole deformation β_{4t,C}.

    /**
     * @brief  Construct zero-valued potential parameters.
     * @math   p_i = 0
     * @output Zero-valued potential parameters.
     */
    VARG()
        : hmass_F(0.0), Ze2_F(0.0), V0_F(0.0), a0_F(0.0), Rp0_F(0.0), Rt0_F(0.0), Rcoul_p_F(0.0), Rcoul_t_F(0.0), beta2t_F(0.0), beta4t_F(0.0), beta2tcoul_F(0.0), beta4tcoul_F(0.0) {}

    /**
     * @brief  Construct parameters from target and projectile nuclei.
     * @math   p = p(A_t, Z_t, A_p, Z_p)
     * @output Spherical potential parameters with zero deformation.
     */
    VARG(int At_I, int Zt_I, int Ap_I, int Zp_I)
        : hmass_F(calc_hmass(At_I, Ap_I)), Ze2_F(e2_F * static_cast<double>(Zt_I * Zp_I)), V0_F(calc_V0ref(At_I, Zt_I, Ap_I)), a0_F(calc_a0(At_I, Zt_I)), Rp0_F(calc_R0(Ap_I, Zp_I)), Rt0_F(calc_R0(At_I, Zt_I)), Rcoul_p_F(calc_Rcoul(Ap_I)), Rcoul_t_F(calc_Rcoul(At_I)), beta2t_F(0.0), beta4t_F(0.0), beta2tcoul_F(0.0), beta4tcoul_F(0.0) {
        assert(Ap_I > 0 && Zp_I >= 0 && Zp_I <= Ap_I);
    }

    /**
     * @brief  Decode potential parameters from their canonical vector representation.
     * @math   p = (h_μ, Ze², V₀, a₀, R_{p0}, R_{t0}, R_{Cp}, R_{Ct}, β_{2t}, β_{4t}, β_{2t,C}, β_{4t,C})
     * @output Potential parameters decoded from p.
     */
    static VARG from_vector(const Eigen::VectorXd& Varg_F1D_p) {
        assert(Varg_F1D_p.size() == 12);
        VARG Varg_params;

        // Canonical parameter-vector ordering.
        Varg_params.hmass_F = Varg_F1D_p(0);
        Varg_params.Ze2_F = Varg_F1D_p(1);
        Varg_params.V0_F = Varg_F1D_p(2);
        Varg_params.a0_F = Varg_F1D_p(3);
        Varg_params.Rp0_F = Varg_F1D_p(4);
        Varg_params.Rt0_F = Varg_F1D_p(5);
        Varg_params.Rcoul_p_F = Varg_F1D_p(6);
        Varg_params.Rcoul_t_F = Varg_F1D_p(7);
        Varg_params.beta2t_F = Varg_F1D_p(8);
        Varg_params.beta4t_F = Varg_F1D_p(9);
        Varg_params.beta2tcoul_F = Varg_F1D_p(10);
        Varg_params.beta4tcoul_F = Varg_F1D_p(11);
        return Varg_params;
    }

    /**
     * @brief  Encode potential parameters in their canonical vector representation.
     * @math   p = (h_μ, Ze², V₀, a₀, R_{p0}, R_{t0}, R_{Cp}, R_{Ct}, β_{2t}, β_{4t}, β_{2t,C}, β_{4t,C})
     * @output Canonically ordered parameter vector p.
     */
    static Eigen::VectorXd to_vector(const VARG& Varg_params) {
        Eigen::VectorXd Varg_F1D_p(12);

        // Canonical parameter-vector ordering.
        Varg_F1D_p(0) = Varg_params.hmass_F;
        Varg_F1D_p(1) = Varg_params.Ze2_F;
        Varg_F1D_p(2) = Varg_params.V0_F;
        Varg_F1D_p(3) = Varg_params.a0_F;
        Varg_F1D_p(4) = Varg_params.Rp0_F;
        Varg_F1D_p(5) = Varg_params.Rt0_F;
        Varg_F1D_p(6) = Varg_params.Rcoul_p_F;
        Varg_F1D_p(7) = Varg_params.Rcoul_t_F;
        Varg_F1D_p(8) = Varg_params.beta2t_F;
        Varg_F1D_p(9) = Varg_params.beta4t_F;
        Varg_F1D_p(10) = Varg_params.beta2tcoul_F;
        Varg_F1D_p(11) = Varg_params.beta4tcoul_F;
        return Varg_F1D_p;
    }
};

/**
 * @brief  Evaluate the Woods-Saxon nuclear potential analytically.
 * @math   V_N(r) = -V₀ / {1 + exp[(r - R_{p0} - R_t) / a₀]}
 * @output Nuclear potential V_N(r) in MeV.
 */
inline double calc_Vnucl(double r_F, const VARG& Varg_params, double Rt_F) {
    assert(r_F >= 0.0 && Varg_params.a0_F > 0.0);
    double R0_F = Varg_params.Rp0_F + Rt_F;
    double x_F = (r_F - R0_F) / Varg_params.a0_F;
    return -Varg_params.V0_F / (1.0 + std::exp(x_F));
}

/**
 * @brief  Evaluate the finite-size Coulomb potential analytically.
 * @math   V_C(r) = Ze²(3 - r²/R²)/(2R), R = max(r, R_{Cp} + R_{Ct})
 * @output Coulomb potential V_C(r) in MeV.
 * @note   Fixed electromagnetic coupling.
 */
inline double calc_Vcoul(double r_F, const VARG& Varg_params) {
    double Rcoul_F = Varg_params.Rcoul_p_F + Varg_params.Rcoul_t_F;
    assert(r_F > 0.0 && Rcoul_F > 0.0);
    bool isOutside_B = r_F >= Rcoul_F;
    return static_cast<double>(isOutside_B) * Varg_params.Ze2_F / r_F + static_cast<double>(!isOutside_B) * Varg_params.Ze2_F * (3.0 - r_F * r_F / (Rcoul_F * Rcoul_F)) / (2.0 * Rcoul_F);
}

/**
 * @brief  Evaluate the centrifugal potential analytically.
 * @math   V_l(r) = h_μ l(l + 1) / r²
 * @output Centrifugal potential V_l(r) in MeV.
 * @note   Nonrelativistic orbital kinetic energy.
 */
inline double calc_Vcent(double r_F, const VARG& Varg_params, int l_I) {
    assert(r_F > 0.0 && Varg_params.hmass_F > 0.0 && l_I >= 0);
    return Varg_params.hmass_F * static_cast<double>(l_I * (l_I + 1)) / (r_F * r_F);
}
