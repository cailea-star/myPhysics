/**
 * @file    hfb_setting.hpp
 * @author  cailea
 * @date    2026-09-14
 * @brief   Common HFB solver settings.
 */

#pragma once

class HFBSetting {
public:
    double temperature_F = 0.0; // Temperature [MeV].
    double EspCut_F = 60.0; // Equivalent single-particle energy cutoff [MeV].

    double accuracy_F = 1.0e-5;
    double mixingMin_F = 0.20;
    double mixingMax_F = 0.90;
    int NiterationsMax_I = 100;

    int NblockingCandidates_I = 6; // Maximum blocking candidates.
    double EblockingCut_F = 1.0; // Blocking window [MeV].

    bool useEspCut_B = false;
    bool useLipkinNogami_B = false;

public:
    /**
     * @brief Construct default HFB settings.
     * @math ∅ → P_HFB.
     * @output Initialized solver settings.
     */
    HFBSetting() = default;
};
