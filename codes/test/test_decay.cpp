/**
 * @file    test_decay.cpp
 * @author  cailea
 * @date    2026-08-25
 * @brief   Validate coupled-channel alpha-decay observables.
 */

#include <Eigen/Core>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <tuple>
#include <vector>

#include "cc_params.hpp"
#include "cc_decay.hpp"
#include "minimize.hpp"

/**
 * @brief One alpha-decay reference case.
 */
struct DecayCase {
    const char* nuclideName_S;
    int At_I;
    int Zt_I;
    double beta2t_F;
    double beta4t_F;
    double rmax_F;
    Eigen::Vector4d Ech_F1D_ch;
    Eigen::Vector4d paperBranchPercent_F1D_ch;
    double paperHalfLife_F;
    Eigen::Vector4d expBranchPercent_F1D_ch;
    double expHalfLife_F;
};

/**
 * @brief  Validate alpha decay by Brent mismatch minimization.
 * @math   V₀=argmin e_match; T₁/₂=ℏ ln(2)/Σ_chΓ_ch
 * @output Reference, computed, and error decay observables.
 * @note   Uses Table I of Phys. Rev. C 81, 064318.
 * @note   Runtime search and tolerance branches are approved.
 */
int main()
{
    constexpr double hbar_F = hbarc_F / 2.99792458e23;
    constexpr double V0min_F = 150.0;
    constexpr double V0max_F = 160.0;
    constexpr double dV0_F = 0.25;

    // Table I: even-even actinide alpha decays.
    const std::vector<DecayCase> testCases = {
        {"Cf", 240, 96, 0.224, 0.087, 40.0, {7.329, 7.294, 7.203, 7.063}, {73.17, 25.24, 1.59, 0.00712}, 1.47e3, {75.00, 25.00, -1.0, -1.0}, 1.66e3},
        {"Cf", 242, 96, 0.224, 0.079, 43.0, {6.862, 6.819, 6.725, 6.574}, {75.91, 22.89, 1.19, 0.00358}, 1.19e5, {79.30, 20.60, 0.15, 0.016}, 1.29e5},
        {"Cf", 244, 96, 0.234, 0.073, 47.0, {6.361, 6.319, 6.218, 6.064}, {79.22, 19.98, 0.79, 0.00458}, 2.43e7, {80.00, 19.60, 0.40, -1.0}, 2.88e7},
        {"Cf", 246, 96, 0.234, 0.057, 49.0, {6.128, 6.086, 5.987, 5.831}, {80.74, 18.56, 0.70, 0.00391}, 3.29e8, {84.70, 15.00, 0.30, 0.010}, 4.13e8},
        {"Cf", 248, 96, 0.235, 0.040, 48.0, {6.217, 6.174, 6.073, 5.920}, {82.85, 16.45, 0.70, 0.00492}, 1.05e8, {84.20, 15.70, 0.24, 0.002}, 8.61e7},
        {"Cf", 250, 96, 0.225, 0.030, 51.0, {5.926, 5.884, 5.783, 5.625}, {85.62, 13.76, 0.62, 0.00302}, 3.17e9, {83.00, 17.00, -1.0, -1.0}, 1.69e9},
        {"Fm", 242, 98, 0.224, 0.079, 35.0, {8.373, 8.329, 8.223, 8.058}, {72.88, 25.44, 1.67, 0.0106}, 1.87, {80.00, 20.00, -1.0, -1.0}, 1.20},
        {"Fm", 244, 98, 0.234, 0.073, 37.0, {7.999, 7.958, 7.862, 7.712}, {73.90, 24.49, 1.59, 0.0220}, 31.7, {80.00, 20.00, -1.0, -1.0}, 38.7},
        {"Fm", 246, 98, 0.234, 0.057, 39.0, {7.557, 7.516, 7.410, 7.249}, {75.51, 23.16, 1.31, 0.0149}, 1.20e3, {83.00, 17.00, -1.0, -1.0}, 2.20e3},
        {"Fm", 248, 98, 0.235, 0.040, 41.0, {7.153, 7.111, 7.015, 6.868}, {79.22, 19.62, 1.15, 0.0137}, 4.82e4, {84.00, 15.00, 0.97, 0.023}, 9.14e4},
        {"Fm", 250, 98, 0.245, 0.026, 43.0, {7.307, 7.264, 7.163, 7.008}, {81.09, 17.81, 1.08, 0.0130}, 1.01e4, {85.00, 14.20, 0.82, 0.0066}, 1.17e4},
        {"Fm", 252, 98, 0.236, 0.015, 45.0, {7.027, 6.981, 6.875, 6.708}, {84.44, 14.66, 0.89, 0.00685}, 1.33e5, {85.00, 15.00, -1.0, -1.0}, 1.17e5},
        {"No", 248, 100, 0.235, 0.049, 35.0, {8.551, 8.507, 8.404, 8.243}, {75.43, 22.96, 1.59, 0.0276}, 2.68, {75.00, 25.00, -1.0, -1.0}, 3.91},
        {"No", 250, 100, 0.235, 0.033, 37.0, {8.223, 8.179, 8.076, 7.915}, {77.87, 20.63, 1.48, 0.0248}, 30.4, {-1.0, -1.0, -1.0, -1.0}, 56.7},
        {"No", 252, 100, 0.245, 0.018, 35.0, {8.582, 8.535, 8.427, 8.256}, {80.35, 18.32, 1.31, 0.0219}, 1.77, {87.00, 13.00, -1.0, -1.0}, 2.93},
        {"Rf", 252, 102, 0.236, 0.024, 35.0, {8.930, 8.883, 8.776, 8.609}, {78.89, 19.56, 1.52, 0.0265}, 0.989, {-1.0, -1.0, -1.0, -1.0}, 2.00},
        {"Rf", 254, 102, 0.246, 0.011, 35.0, {9.250, 9.206, 9.105, 8.946}, {80.26, 18.26, 1.45, 0.0271}, 0.0943, {-1.0, -1.0, -1.0, -1.0}, 0.0923},
        {"Sg", 256, 104, 0.247, -0.007, 33.0, {9.923, 9.872, 9.753, 9.566}, {83.66, 15.23, 1.09, 0.0165}, 0.00688, {83.00, 17.00, -1.0, -1.0}, 0.00720},
    };

    std::cout << std::scientific << std::setprecision(6);
    auto print_result_Func = [](const char* label_S, const auto& resultBranchPercent_F1D_ch, double resultHalfLife_F) {
        std::cout << label_S << std::fixed;
        for (int channel_I = 0; channel_I < 4; ++channel_I) {
            if (channel_I > 0) { std::cout << ", "; }
            if (resultBranchPercent_F1D_ch(channel_I) < 0.0) { std::cout << std::setw(10) << '-'; }
            else { std::cout << std::setw(9) << resultBranchPercent_F1D_ch(channel_I) << '%'; }
        }
        std::cout << ", " << std::scientific << resultHalfLife_F << " s\n";
    };
    for (const DecayCase& testCase : testCases) {
        VARG Varg_params(testCase.At_I, testCase.Zt_I, 4, 2);
        Varg_params.beta2t_F = testCase.beta2t_F;
        Varg_params.beta4t_F = testCase.beta4t_F;
        Varg_params.beta2tcoul_F = testCase.beta2t_F;
        Varg_params.beta4tcoul_F = testCase.beta4t_F;
        std::vector<Channel> channels = {
            Channel(0, 0, 0, 0, testCase.Ech_F1D_ch(0)),
            Channel(0, 4, 0, 4, testCase.Ech_F1D_ch(1)),
            Channel(0, 8, 0, 8, testCase.Ech_F1D_ch(2)),
            Channel(0, 12, 0, 12, testCase.Ech_F1D_ch(3)),
        };

        // r_m=max_ch r_mid,ch.
        Varg_params.V0_F = 0.5 * (V0min_F + V0max_F);
        CCParams barrierParams(testCase.At_I, testCase.Zt_I, 4, 2, 22, Varg_params, channels);
        double rmatch_F = 0.0;
        for (const Channel& channel : channels) {
            double rmid_F = std::get<2>(barrierParams.find_barrier(channel.twol_I / 2, channel.Ech_F));
            rmatch_F = std::max(rmatch_F, rmid_F);
        }
        assert(0.01 < rmatch_F && rmatch_F < testCase.rmax_F);

        // V₀ → e_match(V₀).
        auto match_Func = [&](double V0trial_F) {
            Varg_params.V0_F = V0trial_F;
            return decay_match(0.01, rmatch_F, testCase.rmax_F, 0.05, CCParams(testCase.At_I, testCase.Zt_I, 4, 2, 22, Varg_params, channels));
        };

        const Eigen::Vector3i targetNodes_I1D_ch(11, 10, 9);
        auto find_V0_Func = [&]() {
            double V0lo_F = V0min_F;
            double V0mid_F = V0lo_F + dV0_F;
            double matchErrorLo_F = std::get<0>(match_Func(V0lo_F));
            double matchErrorMid_F = std::get<0>(match_Func(V0mid_F));
            for (double V0up_F = V0mid_F + dV0_F; V0up_F <= V0max_F; V0up_F += dV0_F) {
                double matchErrorUp_F = std::get<0>(match_Func(V0up_F));
                if (matchErrorMid_F < matchErrorLo_F && matchErrorMid_F < matchErrorUp_F) {
                    double V0_F = minimize_brent([&](double V0trial_F) {return std::get<0>(match_Func(V0trial_F));}, V0lo_F, V0mid_F, V0up_F, 1.0e-5);
                    double matchError_F = std::get<0>(match_Func(V0_F));
                    Eigen::VectorXi nodes_I1D_ch = decay_nodes(0.01, rmatch_F, testCase.rmax_F, 0.05, CCParams(testCase.At_I, testCase.Zt_I, 4, 2, 22, Varg_params, channels));
                    if (matchError_F < 1.0e-5 && (nodes_I1D_ch.head<3>() - targetNodes_I1D_ch).isZero()) {
                        return std::tuple{V0_F, matchError_F, nodes_I1D_ch};
                    }
                }
                V0lo_F = V0mid_F;
                matchErrorLo_F = matchErrorMid_F;
                V0mid_F = V0up_F;
                matchErrorMid_F = matchErrorUp_F;
            }
            return std::tuple{std::nan(""), std::nan(""), Eigen::VectorXi()};
        };

        auto [V0_F, matchError_F, nodes_I1D_ch] = find_V0_Func();
        assert(std::isfinite(V0_F) && matchError_F < 1.0e-5);
        CCParams params(testCase.At_I, testCase.Zt_I, 4, 2, 22, Varg_params, channels);
        Eigen::VectorXd Gamma_F1D_ch = decay_width(0.01, rmatch_F, testCase.rmax_F, 0.05, params);
        Eigen::VectorXd branchPercent_F1D_ch = 100.0 * Gamma_F1D_ch / Gamma_F1D_ch.sum();
        double halfLife_F = hbar_F * std::log(2.0) / Gamma_F1D_ch.sum();
        Eigen::Vector4d branchTolerance_F1D_ch = testCase.Zt_I == 96 ? Eigen::Vector4d(3.2, 2.6, 2.4, 0.04) : Eigen::Vector4d(1.3, 1.0, 0.55, 0.02);
        double halfLifeTolerance_F = testCase.Zt_I == 96 ? 0.38 : 0.35;
        Eigen::Vector4d branchError_F1D_ch = (branchPercent_F1D_ch - testCase.paperBranchPercent_F1D_ch).cwiseAbs();
        double halfLifeRelativeError_F = std::abs(halfLife_F / testCase.paperHalfLife_F - 1.0);

        std::cout << std::fixed << std::setprecision(6) << "\n"
                  << testCase.nuclideName_S << testCase.At_I + 4 << ": V0 = " << V0_F << " MeV, error = "
                  << std::scientific << matchError_F << ", rmatch = " << std::fixed << rmatch_F
                  << " fm, nodes = " << nodes_I1D_ch.transpose() << '\n';

        print_result_Func("calculated: ", branchPercent_F1D_ch, halfLife_F);
        print_result_Func("reference:  ", testCase.paperBranchPercent_F1D_ch, testCase.paperHalfLife_F);
        print_result_Func("experiment: ", testCase.expBranchPercent_F1D_ch, testCase.expHalfLife_F);
        std::cout << std::flush;

        assert((nodes_I1D_ch.head<3>() - targetNodes_I1D_ch).isZero());
        assert((branchError_F1D_ch.array() < branchTolerance_F1D_ch.array()).all());
        assert(halfLifeRelativeError_F < halfLifeTolerance_F);
    }
    return 0;
}
