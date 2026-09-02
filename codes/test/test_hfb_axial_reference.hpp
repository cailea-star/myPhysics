/**
 * @file    test_hfb_axial_reference.hpp
 * @author  cailea
 * @date    2026-09-02
 * @brief   Store reference axial-HFB observables.
 */

#pragma once

/**
 * @brief Store one reference observable set.
 */
struct AxialHFBObservableReference {
    double Etot_F;
    double lambda_n_F;
    double lambda_p_F;
    double Epair_n_F;
    double Epair_p_F;
    double Delta_n_F;
    double Delta_p_F;
    double rms_n_F;
    double rms_p_F;
};
