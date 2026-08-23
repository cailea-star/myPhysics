/**
 * @file    cc0_params.hpp
 * @author  cailea
 * @date    2026-08-23
 * @brief   Coupled-channels parameters, validation, and barrier search.
 */

#pragma once

#include <Eigen/Core>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include "cc0_potentials.hpp"
#include "minimize.hpp"
#include "root.hpp"

/**
 * @brief One coupled-channel basis state.
 */
class Channel {
public:
    int twoJ_I;    // 2J.
    int twoIt_I;   // 2I_t.
    int twoKt_I;   // 2K_t.
    int twol_I;    // 2l.
    double Ech_F;  // E_ch [MeV].

    /**
     * @brief  Construct one coupled channel.
     * @math   c = (2J, 2I_t, 2K_t, 2l, E_ch)
     * @output Initialized channel c.
     */
    Channel(int twoJ_I_, int twoIt_I_, int twoKt_I_, int twol_I_, double Ech_F_)
        : twoJ_I(twoJ_I_), twoIt_I(twoIt_I_), twoKt_I(twoKt_I_), twol_I(twol_I_), Ech_F(Ech_F_)
    {
        assert(twol_I >= 0 && twoJ_I >= 0 && twoIt_I >= 0 && twoKt_I >= 0);
        assert(twol_I % 2 == 0 && std::abs(twoIt_I - twoKt_I) % 2 == 0);
        assert(std::abs(twoJ_I - twoIt_I) <= twol_I && twol_I <= twoJ_I + twoIt_I && twoKt_I <= twoIt_I);
    }

    /**
     * @brief  Format one coupled channel row.
     * @math   c ↦ (J, I_t, K_t, l, E_ch)
     * @output Fixed-width channel row.
     */
    std::string to_string() const;

    /**
     * @brief  Format coupled-channel table headers.
     * @math   h = (J, I_t, K_t, l, E_ch)
     * @output Fixed-width channel header.
     */
    static std::string table_header();
};

/**
 * @brief Coupled-channels reaction parameters.
 */
class CCParams {
public:
    int At_I;                            // Target mass number A_t.
    int Zt_I;                            // Target proton number Z_t.
    int Ap_I;                            // Projectile mass number A_p.
    int Zp_I;                            // Projectile proton number Z_p.
    VARG Varg_params;                    // Potential parameters p_V.
    int G_I;                             // G = 2n + l.
    std::vector<Channel> channel_1D_ch;  // {c_ch}.

    /**
     * @brief  Construct coupled-channels parameters.
     * @math   P = (A_t, Z_t, A_p, Z_p, G, p_V, {c_ch})
     * @output Initialized parameter set P.
     */
    CCParams(int At_I_, int Zt_I_, int Ap_I_, int Zp_I_, int G_I_, const VARG& Varg_params_, const std::vector<Channel>& channel_1D_ch_)
        // P ← (A_t, Z_t, A_p, Z_p, G, p_V, {c_ch}).
        : At_I(At_I_), Zt_I(Zt_I_), Ap_I(Ap_I_), Zp_I(Zp_I_), Varg_params(Varg_params_), G_I(G_I_), channel_1D_ch(channel_1D_ch_)
    {
    }

    /**
     * @brief  Destroy polymorphic coupled-channels parameters.
     * @math   P → ∅
     * @output Released parameter resources.
     */
    virtual ~CCParams() = default;

    /**
     * @brief  Validate coupled-channels parameters.
     * @math   A_x > 0; 0 ≤ Z_x ≤ A_x; G-l_ch ∈ 2ℤ₀
     * @output True iff all parameters are admissible.
     */
    bool is_valid() const;

    /**
     * @brief  Evaluate the coupled potential matrix.
     * @math   V_{ch,ch'}(r) = V_C(r) + V_N(r)
     * @output Coupled potential matrix in MeV.
     * @note   Implemented by the coupling module.
     */
    virtual Eigen::MatrixXcd calc_Vall_matrix(double r_F) const;

    /**
     * @brief  Evaluate uncoupled nuclear and Coulomb potential.
     * @math   V_0(r) = V_N(r) + V_C(r)
     * @output Scalar potential V_0(r) in MeV.
     * @note   Uses the spherical target radius.
     */
    virtual double calc_Vall_scalar(double r_F) const;

    /**
     * @brief  Find barriers using golden-section and bisection algorithms.
     * @math   V_eff(r_i) = E; r_w = argmin V_eff; r_b = argmax V_eff
     * @output (r_in, r_w, r_mid, r_b, r_out) in fm.
     * @note   Assumes one inner well and one outer barrier.
     * @note   Bracket validation uses approved runtime branches.
     */
    std::tuple<double, double, double, double, double> find_barrier(int l_I, double E_F, double rmin_F = 0.05, double rmax_F = 100.0) const;

    /**
     * @brief  Print parameters and barrier crossings.
     * @math   P ↦ table(P, r_min, r_max)
     * @output Parameter and channel tables.
    */
    void print_info(double rmin_F = 0.05, double rmax_F = 100.0) const;
};

inline std::string Channel::to_string() const
{
    std::ostringstream output_Oss;

    // (2J, 2I_t, 2K_t, 2l, E_ch) → (J, I_t, K_t, l, E_ch).
    output_Oss << std::fixed << std::setprecision(3) << std::left;
    output_Oss << std::setw(6) << 0.5 * static_cast<double>(twoJ_I) << ", ";
    output_Oss << std::setw(6) << 0.5 * static_cast<double>(twoIt_I) << ", ";
    output_Oss << std::setw(6) << 0.5 * static_cast<double>(twoKt_I) << ", ";
    output_Oss << std::setw(6) << 0.5 * static_cast<double>(twol_I) << ", ";
    output_Oss << std::setw(12) << Ech_F;
    return output_Oss.str();
}

inline std::string Channel::table_header()
{
    std::ostringstream output_Oss;

    // h = (J, I_t, K_t, l, E_ch).
    output_Oss << std::left;
    output_Oss << std::setw(6) << "J" << ", ";
    output_Oss << std::setw(6) << "I_t" << ", ";
    output_Oss << std::setw(6) << "K_t" << ", ";
    output_Oss << std::setw(6) << "l" << ", ";
    output_Oss << std::setw(12) << "E_ch [MeV]";
    return output_Oss.str();
}

inline bool CCParams::is_valid() const
{
    // 0 ≤ Z_t ≤ A_t.
    if (At_I <= 0 || Zt_I < 0 || Zt_I > At_I) {
        std::cerr << "[ERROR][CCParams] Invalid target nucleus.\n";
        return false;
    }

    // 0 ≤ Z_p ≤ A_p.
    if (Ap_I <= 0 || Zp_I < 0 || Zp_I > Ap_I) {
        std::cerr << "[ERROR][CCParams] Invalid projectile nucleus.\n";
        return false;
    }

    // G ∈ ℤ₀.
    if (G_I < 0) {
        std::cerr << "[ERROR][CCParams] Negative oscillator quantum number.\n";
        return false;
    }

    // ∀ch: G - l_ch ∈ 2ℤ₀.
    for (const Channel& channel_ch : channel_1D_ch) {
        int channel_l_I = channel_ch.twol_I / 2;
        if (G_I < channel_l_I || (G_I - channel_l_I) % 2 != 0) {
            std::cerr << "[ERROR][CCParams] Invalid channel selection rule.\n";
            return false;
        }
    }
    return true;
}

inline double CCParams::calc_Vall_scalar(double r_F) const
{
    assert(std::isfinite(r_F) && r_F > 0.0);
    return calc_Vnucl(r_F, Varg_params, Varg_params.Rt0_F) + calc_Vcoul(r_F, Varg_params);
}

// (l, E, r_min, r_max) → (r_in, r_w, r_mid, r_b, r_out).
inline std::tuple<double, double, double, double, double> CCParams::find_barrier(int l_I, double E_F, double rmin_F, double rmax_F) const
{
    assert(l_I >= 0 && std::isfinite(E_F));
    assert(rmin_F > 0.0 && rmax_F > rmin_F);
    assert(std::isfinite(rmin_F) && std::isfinite(rmax_F));

    double Rsum_F = Varg_params.Rp0_F + Varg_params.Rt0_F;
    assert(std::isfinite(Rsum_F) && Rsum_F > 0.0);
    double rbarrier_max_F = 2.0 * Rsum_F;
    double rwell_min_F = 0.1 * Rsum_F;

    // V_eff(r) = V_0(r) + V_l(r).
    auto Veff_Func = [this, l_I](double r_F) -> double {
        return calc_Vall_scalar(r_F) + calc_Vcent(r_F, Varg_params, l_I);
    };
    auto neg_Veff_Func = [&Veff_Func](double r_F) -> double {
        return -Veff_Func(r_F);
    };
    double rbarrier_F = minimize_golden(neg_Veff_Func, 0.8 * Rsum_F, rbarrier_max_F);
    double rwell_F = minimize_golden(Veff_Func, rwell_min_F, 0.9 * rbarrier_F);
    rwell_F = std::max(rwell_min_F, rwell_F);

    // f(r) = V_eff(r) - E.
    auto VminusE_Func = [&Veff_Func, E_F](double r_F) -> double {
        return Veff_Func(r_F) - E_F;
    };
    auto find_root_or_nan_Func = [&VminusE_Func](double rlo_F, double rup_F) -> double {
        double flo_F = VminusE_Func(rlo_F);
        double fup_F = VminusE_Func(rup_F);
        if (std::abs(flo_F) <= 1.0e-6) {
            return rlo_F;
        }
        if (std::abs(fup_F) <= 1.0e-6) {
            return rup_F;
        }
        if (flo_F * fup_F > 0.0) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        return root_bisection(VminusE_Func, rlo_F, rup_F);
    };

    // V_eff(r_i) = E, i ∈ {in, mid, out}.
    double rin_F = find_root_or_nan_Func(rmin_F, rwell_F);
    double rmid_F = find_root_or_nan_Func(rwell_F, rbarrier_F);
    double rout_F = find_root_or_nan_Func(rbarrier_F, rmax_F);
    return std::make_tuple(rin_F, rwell_F, rmid_F, rbarrier_F, rout_F);
}

inline void CCParams::print_info(double rmin_F, double rmax_F) const
{
    assert(rmin_F > 0.0 && rmax_F > rmin_F);
    assert(std::isfinite(rmin_F) && std::isfinite(rmax_F));
    bool isValid_B = is_valid();
    std::ostringstream output_Oss;
    output_Oss << std::fixed << std::setprecision(3);

    // P → parameter labels.
    output_Oss << "[CC] " << std::right;
    output_Oss << std::setw(7) << "valid";
    output_Oss << std::setw(7) << "N_ch";
    output_Oss << std::setw(7) << "A_t";
    output_Oss << std::setw(7) << "Z_t";
    output_Oss << std::setw(7) << "A_p";
    output_Oss << std::setw(7) << "Z_p";
    output_Oss << std::setw(7) << "G";
    output_Oss << std::setw(14) << "h_mu";
    output_Oss << std::setw(12) << "beta2t";
    output_Oss << std::setw(12) << "beta4t";
    output_Oss << std::setw(12) << "r_min";
    output_Oss << std::setw(12) << "r_max";
    output_Oss << std::setw(12) << "Ze2";
    output_Oss << std::setw(12) << "V0";
    output_Oss << std::setw(12) << "a0";
    output_Oss << std::setw(12) << "Rp0";
    output_Oss << std::setw(12) << "Rt0";
    output_Oss << std::setw(12) << "Rcoul_p";
    output_Oss << std::setw(12) << "Rcoul_t";
    output_Oss << '\n';

    // P → parameter values.
    output_Oss << "[CC] " << std::right;
    output_Oss << std::setw(7) << (isValid_B ? "yes" : "no");
    output_Oss << std::setw(7) << channel_1D_ch.size();
    output_Oss << std::setw(7) << At_I;
    output_Oss << std::setw(7) << Zt_I;
    output_Oss << std::setw(7) << Ap_I;
    output_Oss << std::setw(7) << Zp_I;
    output_Oss << std::setw(7) << G_I;
    output_Oss << std::setw(14) << Varg_params.hmass_F;
    output_Oss << std::setw(12) << Varg_params.beta2t_F;
    output_Oss << std::setw(12) << Varg_params.beta4t_F;
    output_Oss << std::setw(12) << rmin_F;
    output_Oss << std::setw(12) << rmax_F;
    output_Oss << std::setw(12) << Varg_params.Ze2_F;
    output_Oss << std::setw(12) << Varg_params.V0_F;
    output_Oss << std::setw(12) << Varg_params.a0_F;
    output_Oss << std::setw(12) << Varg_params.Rp0_F;
    output_Oss << std::setw(12) << Varg_params.Rt0_F;
    output_Oss << std::setw(12) << Varg_params.Rcoul_p_F;
    output_Oss << std::setw(12) << Varg_params.Rcoul_t_F;
    output_Oss << '\n';

    // ch → channel labels.
    output_Oss << "[CH] " << std::left;
    output_Oss << std::setw(6) << "index" << ", ";
    output_Oss << Channel::table_header() << ",  ";
    output_Oss << std::setw(12) << "r_in";
    output_Oss << std::setw(12) << "r_well";
    output_Oss << std::setw(12) << "r_mid";
    output_Oss << std::setw(12) << "r_barrier";
    output_Oss << std::setw(12) << "r_out";
    output_Oss << '\n';

    // ∀ch: c_ch → channel values.
    int channel_index_I = 0;
    for (const Channel& channel_ch : channel_1D_ch) {
        auto [rin_F, rwell_F, rmid_F, rbarrier_F, rout_F] = find_barrier(channel_ch.twol_I / 2, channel_ch.Ech_F, rmin_F, rmax_F);
        output_Oss << "[CH] " << std::left;
        output_Oss << std::setw(6) << channel_index_I << ", ";
        output_Oss << channel_ch.to_string() << ",  ";
        output_Oss << std::setw(12) << rin_F;
        output_Oss << std::setw(12) << rwell_F;
        output_Oss << std::setw(12) << rmid_F;
        output_Oss << std::setw(12) << rbarrier_F;
        output_Oss << std::setw(12) << rout_F;
        output_Oss << '\n';
        ++channel_index_I;
    }

    std::cout << output_Oss.str();
}
