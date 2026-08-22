/**
 * @file    bifold_m3y_xd.hpp
 * @author  cailea
 * @date    2026-08-22
 * @brief   XDM3Y density-dependence parameter sets.
 * @ref     Computer Physics Communications 284 (2023) 108613
 */

#pragma once

#include <cassert>
#include <utility>

#include "bifold_m3y_v.hpp"

/**
 * @brief Store one XDM3Y density-dependence parameterization.
 */
class XDM3YParameters {
public:
    double C_F = 1.0; // F(ρ) overall scale.
    double alpha_F = 0.0; // αe^(-βρ) coefficient.
    double beta_F = 0.0; // [β] = fm³.
    double gamma_F = 0.0; // [γ] = fm^(3n).
    double n_F = 0.0; // ρ^n exponent.
    double K_F = 0.0; // K, MeV.

    /**
     * @brief  Construct an empty XDM3Y parameter bundle.
     * @math   𝒳 = {C = 1, α = β = γ = n = K = 0}
     * @output Empty XDM3Y parameter bundle.
     */
    XDM3YParameters() = default;

    /**
     * @brief  Construct a validated XDM3Y parameter bundle.
     * @math   F(ρ) = C[1 + αe^(-βρ) - γρ^n]
     * @output Initialized XDM3Y parameter bundle.
     * @note   Requires C, K > 0 and α, β, γ, n ≥ 0.
     */
    XDM3YParameters(double C_F_, double alpha_F_, double beta_F_, double gamma_F_, double n_F_, double K_F_)
        : C_F(C_F_), alpha_F(alpha_F_), beta_F(beta_F_), gamma_F(gamma_F_), n_F(n_F_), K_F(K_F_) {
        assert(C_F > 0.0);
        assert(alpha_F >= 0.0);
        assert(beta_F >= 0.0);
        assert(gamma_F >= 0.0);
        assert(n_F >= 0.0);
        assert(K_F > 0.0);
    }
};

/**
 * @brief Couple M3Y kernels to one XDM3Y parameterization.
 */
class XDM3YFunctions {
public:
    M3YFunctions m3y_functions; // {v_C, v_D, v_E, v_E^ZR, g(E)}.
    XDM3YParameters xdm3y_parameters; // {C, α, β, γ, n, K}.

    /**
     * @brief  Construct an XDM3Y interaction bundle.
     * @math   𝒱_X = {𝒱_M, 𝒳}
     * @output Initialized XDM3Y interaction bundle.
     */
    XDM3YFunctions(M3YFunctions m3y_functions_, XDM3YParameters xdm3y_parameters_)
        : m3y_functions(std::move(m3y_functions_)), xdm3y_parameters(std::move(xdm3y_parameters_)) {}
};

// Published Reid and Paris XDM3Y parameterizations.
inline const XDM3YFunctions reid_ddm3y1(reid_functions, XDM3YParameters(0.2845, 3.6391, 2.9605, 0.0000, 0.0, 171.0));
inline const XDM3YFunctions reid_bdm3y0(reid_functions, XDM3YParameters(1.3827, 0.0000, 0.0000, 1.1135, 2.0 / 3.0, 232.0));
inline const XDM3YFunctions reid_bdm3y1(reid_functions, XDM3YParameters(1.2253, 0.0000, 0.0000, 1.5124, 1.0, 232.0));
inline const XDM3YFunctions reid_bdm3y2(reid_functions, XDM3YParameters(1.0678, 0.0000, 0.0000, 5.1069, 2.0, 354.0));
inline const XDM3YFunctions reid_bdm3y3(reid_functions, XDM3YParameters(1.0153, 0.0000, 0.0000, 21.073, 3.0, 475.0));

inline const XDM3YFunctions paris_ddm3y1(paris_functions, XDM3YParameters(0.2963, 3.7231, 3.7384, 0.0000, 0.0, 176.0));
inline const XDM3YFunctions paris_bdm3y1(paris_functions, XDM3YParameters(1.2521, 0.0000, 0.0000, 1.7452, 1.0, 270.0));
inline const XDM3YFunctions paris_bdm3y2(paris_functions, XDM3YParameters(1.0664, 0.0000, 0.0000, 6.0296, 2.0, 418.0));
inline const XDM3YFunctions paris_bdm3y3(paris_functions, XDM3YParameters(1.0045, 0.0000, 0.0000, 25.1150, 3.0, 566.0));
inline const XDM3YFunctions paris_cdm3y1(paris_functions, XDM3YParameters(0.3429, 3.0232, 3.5512, 0.5, 1.0, 188.0));
inline const XDM3YFunctions paris_cdm3y2(paris_functions, XDM3YParameters(0.3346, 3.0357, 3.0685, 1.0, 1.0, 204.0));
inline const XDM3YFunctions paris_cdm3y3(paris_functions, XDM3YParameters(0.2985, 3.4528, 2.6388, 1.5, 1.0, 217.0));
inline const XDM3YFunctions paris_cdm3y4(paris_functions, XDM3YParameters(0.3052, 3.2998, 2.3180, 2.0, 1.0, 228.0));
inline const XDM3YFunctions paris_cdm3y5(paris_functions, XDM3YParameters(0.2728, 3.7367, 1.8294, 3.0, 1.0, 241.0));
inline const XDM3YFunctions paris_cdm3y6(paris_functions, XDM3YParameters(0.2658, 3.8033, 1.4099, 4.0, 1.0, 252.0));
