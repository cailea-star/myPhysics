/**
 * @file    hfb_axial.cpp
 * @author  cailea
 * @date    2026-05-05
 * @brief   Implement axial HFB self-consistent iterations.
 */

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>

#include "hfb_axial.hpp"
#include "hfb_axial_observable.hpp"
#include "root.hpp"
#include "root_broyden.hpp"

namespace {

/**
 * @brief Calculate packed Γ and Δ storage.
 * @math N_p=2Σ_b n_b(n_b+1)/2
 * @output Packed element count.
 */
int calc_packed_size(const AxialHFBBlockList& blocklist_) {
    int Npacked_I = 0;

    // N_p ← Σ_b n_b(n_b+1)/2.
    for (const auto& block_ : blocklist_.blocks_S1D_block) {
        const int Nmatrix_I = static_cast<int>(block_.labels_S1D_bsp.size());
        Npacked_I += Nmatrix_I * (Nmatrix_I + 1) / 2;
    }
    return 2 * Npacked_I; // Γ ⊕ Δ.
}
/**
 * @brief Pack one symmetric matrix upper triangle.
 * @math A_{rc},r≤c → x_k
 * @output Updated packed vector and index.
 */
void pack_upper_triangle(const Eigen::MatrixXd& matrix_F2D_row_column_, Eigen::VectorXd& data_F1D_packed_, int& packed_I_) {
    const int Nmatrix_I = static_cast<int>(matrix_F2D_row_column_.rows());

    // x_k ← A_{rc}, r≤c.
    for (int column_I = 0; column_I < Nmatrix_I; ++column_I) {
        for (int row_I = 0; row_I <= column_I; ++row_I) {
            data_F1D_packed_(packed_I_++) = matrix_F2D_row_column_(row_I, column_I);
        }
    }
}

/**
 * @brief Unpack one symmetric matrix upper triangle.
 * @math x_k → A_{rc}=A_{cr},r≤c
 * @output Updated symmetric matrix and index.
 */
void unpack_upper_triangle(Eigen::MatrixXd& matrix_F2D_row_column_, const Eigen::VectorXd& data_F1D_packed_, int& packed_I_) {
    const int Nmatrix_I = static_cast<int>(matrix_F2D_row_column_.rows());

    // x_k → A_{rc}=A_{cr}.
    for (int column_I = 0; column_I < Nmatrix_I; ++column_I) {
        for (int row_I = 0; row_I <= column_I; ++row_I) {
            const double value_F = data_F1D_packed_(packed_I_++);
            matrix_F2D_row_column_(row_I, column_I) = value_F;
            matrix_F2D_row_column_(column_I, row_I) = value_F;
        }
    }
}

/**
 * @brief Pack neutron and proton HFB fields.
 * @math (Γ_n,Δ_n,Γ_p,Δ_p) → x
 * @output Packed field vector.
 */
void pack_Gamma_Delta(const AxialHFBBlockList& blocklist_n_, const AxialHFBBlockList& blocklist_p_, Eigen::VectorXd& data_F1D_packed_) {
    int packed_I = 0;

    // x ← Γ_n ⊕ Δ_n ⊕ Γ_p ⊕ Δ_p.
    for (const auto& block_ : blocklist_n_.blocks_S1D_block) {pack_upper_triangle(block_.Gamma_F2D_bsp_bsp, data_F1D_packed_, packed_I);}
    for (const auto& block_ : blocklist_n_.blocks_S1D_block) {pack_upper_triangle(block_.Delta_F2D_bsp_bsp, data_F1D_packed_, packed_I);}
    for (const auto& block_ : blocklist_p_.blocks_S1D_block) {pack_upper_triangle(block_.Gamma_F2D_bsp_bsp, data_F1D_packed_, packed_I);}
    for (const auto& block_ : blocklist_p_.blocks_S1D_block) {pack_upper_triangle(block_.Delta_F2D_bsp_bsp, data_F1D_packed_, packed_I);}
}

/**
 * @brief Unpack neutron and proton HFB fields.
 * @math x → (Γ_n,Δ_n,Γ_p,Δ_p)
 * @output Updated block fields.
 */
void unpack_Gamma_Delta(AxialHFBBlockList& blocklist_n_, AxialHFBBlockList& blocklist_p_, const Eigen::VectorXd& data_F1D_packed_) {
    int packed_I = 0;

    // x → Γ_n ⊕ Δ_n ⊕ Γ_p ⊕ Δ_p.
    for (auto& block_ : blocklist_n_.blocks_S1D_block) {unpack_upper_triangle(block_.Gamma_F2D_bsp_bsp, data_F1D_packed_, packed_I);}
    for (auto& block_ : blocklist_n_.blocks_S1D_block) {unpack_upper_triangle(block_.Delta_F2D_bsp_bsp, data_F1D_packed_, packed_I);}
    for (auto& block_ : blocklist_p_.blocks_S1D_block) {unpack_upper_triangle(block_.Gamma_F2D_bsp_bsp, data_F1D_packed_, packed_I);}
    for (auto& block_ : blocklist_p_.blocks_S1D_block) {unpack_upper_triangle(block_.Delta_F2D_bsp_bsp, data_F1D_packed_, packed_I);}
}

/**
 * @brief Print one HFB iteration summary.
 * @math (i,ε,α,O) → stdout
 * @output Formatted iteration row.
 */
void print_summary(const AxialHFBObservable& observable_, const HFBSettings& hfbsettings_, int iteration_I, double residual_F, double mixing_F) {
    if (iteration_I == 0) {
        std::cout << "  ----------------------------------------------------------------------------------------------------------------------------------\n";
        std::cout << "  i          si     mix   beta      Etot       A       rn      rp         En      Dn      Ep      Dp        Ln      Lp\n";
        std::cout << "  ----------------------------------------------------------------------------------------------------------------------------------" << std::endl;
    }

    // Δ_q ← |Δ_q|+λ_{2,q}.
    double Delta_n_F = std::abs(observable_.Delta_n_F);
    double Delta_p_F = std::abs(observable_.Delta_p_F);
    if (hfbsettings_.useLipkinNogami_B) {
        Delta_n_F += observable_.lambda2_n_F;
        Delta_p_F += observable_.lambda2_p_F;
    }

    // (i,ε,α,O) → stdout.
    std::cout << std::fixed;
    std::cout << std::setw(4) << iteration_I << (iteration_I > 0 ? "B" : "L") << " ";
    std::cout << std::setw(12) << std::setprecision(8) << residual_F;
    std::cout << std::setw(5) << std::setprecision(2) << mixing_F;
    std::cout << std::setw(7) << std::setprecision(3) << observable_.beta2_F;
    std::cout << std::setw(13) << std::setprecision(6) << observable_.Etot_F << " ";
    std::cout << std::setw(6) << std::setprecision(1) << observable_.A_F;
    std::cout << std::setw(8) << std::setprecision(3) << observable_.rms_n_F;
    std::cout << std::setw(8) << std::setprecision(3) << observable_.rms_p_F << " | ";
    std::cout << std::setw(8) << std::setprecision(3) << observable_.Epair_n_F;
    std::cout << std::setw(8) << std::setprecision(3) << Delta_n_F << " | ";
    std::cout << std::setw(8) << std::setprecision(3) << observable_.Epair_p_F;
    std::cout << std::setw(8) << std::setprecision(3) << Delta_p_F << " | ";
    std::cout << std::setw(8) << std::setprecision(3) << observable_.lambda_n_F;
    std::cout << std::setw(8) << std::setprecision(3) << observable_.lambda_p_F;
    std::cout << std::endl;
}

/**
 * @brief Print active quasiparticle blocking states.
 * @math {B_μ} → stdout
 * @output Formatted blocking rows.
 */
void print_blocking(const AxialConfig& axialconfig_, const std::vector<AxialHFBBlocking>& activeBlockings_) {
    for (const AxialHFBBlocking& blocking_ : activeBlockings_) {
        if (blocking_.block_I < 0 || blocking_.block_I >= static_cast<int>(axialconfig_.labels_S2D_block_bsp.size())) {continue;}
        if (blocking_.bqp_I < 0) {continue;}
        const auto& labels_S1D_bsp = axialconfig_.labels_S2D_block_bsp[blocking_.block_I];
        if (blocking_.blockedU_F1D_bsp.size() != labels_S1D_bsp.size()) {continue;}
        if (blocking_.blockedV_F1D_bsp.size() != labels_S1D_bsp.size()) {continue;}

        // max(|U|,|V|) → bsp_max.
        Eigen::Index bspMax_I = 0;
        (blocking_.blockedU_F1D_bsp.cwiseAbs().cwiseMax(blocking_.blockedV_F1D_bsp.cwiseAbs())).maxCoeff(&bspMax_I);
        const AxialSPLabel& label_ = labels_S1D_bsp[bspMax_I];

        // (q,block,bqp,overlap,α_bsp) → stdout.
        std::cout << "      blocking " << (blocking_.isNeutron_B ? "n" : "p");
        std::cout << " block=" << std::setw(3) << (blocking_.block_I + 1);
        std::cout << " qp=" << std::setw(3) << (blocking_.bqp_I + 1);
        std::cout << " ovlp=" << std::setw(8) << std::setprecision(4) << blocking_.overlap_F;
        std::cout << " label=" << label_.twoOmega_I << "/2";
        std::cout << (label_.isParityPositive_B ? "+" : "-") << "[";
        std::cout << label_.N_I << ",";
        std::cout << label_.nz_I << ",";
        std::cout << label_.Lambda_I << "]" << std::endl;
    }
}

/**
 * @brief Solve the particle-number constraint using Brent.
 * @math N_q(λ_q)-N_q^{target}=0
 * @output Updated λ, U, V, E, ρ, and κ.
 * @note Runtime branches implement root bracketing.
 */
void update_lambda_UV_E_rho_kappa(AxialHFBBlockList& blocklist_, int Ntarget_I, const HFBSettings& hfbsettings_, std::vector<AxialHFBBlocking>& activeBlockings_, bool isNeutron_B, double tolerance_F = 1.0e-6) {
    const double Ntolerance_F = std::max(1.0e-8, 1.0e-10 * std::max(1, Ntarget_I));

    // λ → N_q(λ)-N_q^{target}.
    const auto calc_particle_error_Func = [&](double lambda_F) {
        blocklist_.lambda_F = lambda_F;
        blocklist_.update_UV_E_rho_kappa(hfbsettings_);
        for (const AxialHFBBlocking& blocking_ : activeBlockings_) {
            if (blocking_.isNeutron_B != isNeutron_B) {continue;}
            AxialHFBBlocking blockingTrial_ = blocking_;
            blockingTrial_.apply_blocking(blocklist_);
        }

        double Ncalc_F = 0.0;
        const int Nblocks_I = static_cast<int>(blocklist_.blocks_S1D_block.size());
        for (int block_I = 0; block_I < Nblocks_I; ++block_I) {
            Ncalc_F += 2.0 * blocklist_.blocks_S1D_block[block_I].rho_F2D_bsp_bsp.trace();
        }
        return Ncalc_F - static_cast<double>(Ntarget_I);
    };

    // λ_0 ← λ_previous or -7 MeV.
    const double lambdaInitial_F = std::isfinite(blocklist_.lambda_F) && blocklist_.lambda_F != 0.0 ? blocklist_.lambda_F : -7.0;
    const double particleErrorInitial_F = calc_particle_error_Func(lambdaInitial_F);
    if (std::abs(particleErrorInitial_F) <= Ntolerance_F) {
        blocklist_.lambda_F = lambdaInitial_F;
        blocklist_.update_UV_E_rho_kappa(hfbsettings_);
        for (AxialHFBBlocking& blocking_ : activeBlockings_) {
            if (blocking_.isNeutron_B == isNeutron_B) {blocking_.apply_blocking(blocklist_);}
        }
        return;
    }

    // [λ_min,λ_max] ← [λ_0-2,λ_0+2].
    double lambdaMax_F = lambdaInitial_F + 2.0;
    double lambdaMin_F = lambdaInitial_F - 2.0;
    double particleErrorMin_F = calc_particle_error_Func(lambdaMin_F);
    double particleErrorMax_F = calc_particle_error_Func(lambdaMax_F);
    const double slope_F = (particleErrorMax_F - particleErrorMin_F) / (lambdaMax_F - lambdaMin_F);

    // f(λ_min)f(λ_max)>0 → expand bracket.
    if (particleErrorMin_F > 0.0 && particleErrorMax_F > 0.0 && slope_F > 0.0) {
        while (particleErrorMin_F > 0.0) {
            lambdaMin_F -= 5.0;
            particleErrorMin_F = calc_particle_error_Func(lambdaMin_F);
        }
    } else if (particleErrorMin_F > 0.0 && particleErrorMax_F > 0.0 && slope_F < 0.0) {
        while (particleErrorMax_F > 0.0) {
            lambdaMax_F += 5.0;
            particleErrorMax_F = calc_particle_error_Func(lambdaMax_F);
        }
    } else if (particleErrorMin_F < 0.0 && particleErrorMax_F < 0.0 && slope_F > 0.0) {
        while (particleErrorMax_F < 0.0) {
            lambdaMax_F += 5.0;
            particleErrorMax_F = calc_particle_error_Func(lambdaMax_F);
        }
    } else if (particleErrorMin_F < 0.0 && particleErrorMax_F < 0.0 && slope_F < 0.0) {
        while (particleErrorMin_F < 0.0) {
            lambdaMin_F -= 5.0;
            particleErrorMin_F = calc_particle_error_Func(lambdaMin_F);
        }
    }

    // f(λ_q)=0 → (U,V,E,ρ,κ)_q.
    const double lambdaRoot_F = root_brent(calc_particle_error_Func, lambdaMin_F, lambdaMax_F, tolerance_F);
    blocklist_.lambda_F = lambdaRoot_F;
    blocklist_.update_UV_E_rho_kappa(hfbsettings_);
    for (AxialHFBBlocking& blocking_ : activeBlockings_) {
        if (blocking_.isNeutron_B == isNeutron_B) {blocking_.apply_blocking(blocklist_);}
    }
}

} // namespace

void AxialHFB::iterate(int Ntarget_I, int Ztarget_I) {
    std::vector<AxialHFBBlocking> noBlockings_;
    iterate(Ntarget_I, Ztarget_I, noBlockings_);
}

void AxialHFB::iterate(int Ntarget_I, int Ztarget_I, AxialHFBBlocking& activeBlocking_) {
    std::vector<AxialHFBBlocking> activeBlockings_;
    activeBlockings_.push_back(activeBlocking_);
    iterate(Ntarget_I, Ztarget_I, activeBlockings_);
    activeBlocking_ = activeBlockings_.front();
}

void AxialHFB::iterate(int Ntarget_I, int Ztarget_I, std::vector<AxialHFBBlocking>& activeBlockings_) {
    const bool hasBlocking_B = !activeBlockings_.empty();
    const EDFParamsSkyrme edfActive_ = hfbsettings.make_active_edf(edf_skyrme, Ntarget_I + Ztarget_I);
    AxialHFBObservable observable_;

    if (hasBlocking_B) {
        std::cout << "[AxialHFB] Active blocking: N=" << Ntarget_I << ", Z=" << Ztarget_I << std::endl;
    }

    // Selected kernels → cached tables.
    if (hfbsettings.termSwitches.addFiniteRangeGogny_B) {gogny.build_tables();}
    if (hfbsettings.termSwitches.addFiniteRangeCoulomb_B) {coulomb.build_tables();}
    if (hfbsettings.termSwitches.addLocalCoulomb_B) {coulombField.build(axialconfig.useReflection_B, edfActive_.e2charg_F);}

    // (D_n,D_p) → local coordinate fields.
    const auto add_fields_Func = [&]() {
        if (hfbsettings.termSwitches.addLocalCoulomb_B) {AxialHFBField::add_coulomb_field(field_p, density_p, coulombField, edfActive_, hfbsettings);}
        if (hfbsettings.termSwitches.addLocalPair_B) {AxialHFBField::add_pairing_fields(field_p, field_n, density_p, density_n, edfActive_, hfbsettings, blocklist_n.lambda_F, blocklist_p.lambda_F);}
    };

    // (ρ,κ,V_G,V_C,λ_2) → (Γ,Δ).
    const auto add_block_fields_Func = [&]() {
        if (hfbsettings.termSwitches.addFiniteRangeGogny_B) {
            AxialHFBBlockList::add_Gamma_from_Gogny(blocklist_p, blocklist_n, gogny);
            AxialHFBBlockList::add_Delta_from_Gogny(blocklist_p, blocklist_n, gogny);
        }
        if (hfbsettings.termSwitches.addFiniteRangeCoulomb_B) {AxialHFBBlockList::add_coulomb_from_Gaussian(blocklist_p, coulomb);}
        if (hfbsettings.useLipkinNogami_B) {
            blocklist_n.add_lipkin_nogami();
            blocklist_p.add_lipkin_nogami();
        }
    };

    // x = Γ_n ⊕ Δ_n ⊕ Γ_p ⊕ Δ_p.
    const int Npacked_I = calc_packed_size(blocklist_n) + calc_packed_size(blocklist_p);
    Eigen::VectorXd x_F1D_packed(Npacked_I);
    Eigen::VectorXd Gx_F1D_packed(Npacked_I);

    // G:x → (λ,U,V,E,ρ,κ,D,F,Γ,Δ).
    const double lambdaToleranceMin_F = hfbsettings.accuracy_F * 1.0e-6;
    double lambdaTolerance_F = hfbsettings.accuracy_F;
    const auto calc_Gx_Func = [&](const Eigen::VectorXd& x_F1D_packed_, Eigen::VectorXd& Gx_F1D_packed_) {
        unpack_Gamma_Delta(blocklist_n, blocklist_p, x_F1D_packed_);
        update_lambda_UV_E_rho_kappa(blocklist_n, Ntarget_I, hfbsettings, activeBlockings_, true, lambdaTolerance_F);
        update_lambda_UV_E_rho_kappa(blocklist_p, Ztarget_I, hfbsettings, activeBlockings_, false, lambdaTolerance_F);
        density_p.update_density(global_basis, blocklist_p);
        density_n.update_density(global_basis, blocklist_n);
        field_p.set_zero();
        field_n.set_zero();
        AxialHFBField::add_nuclei_fields(field_p, field_n, density_p, density_n, edfActive_, hfbsettings);
        add_fields_Func();
        blocklist_p.set_zero_Gamma_Delta();
        blocklist_n.set_zero_Gamma_Delta();
        blocklist_p.add_Gamma_Delta_from_field(field_p, global_basis);
        blocklist_n.add_Gamma_Delta_from_field(field_n, global_basis);
        add_block_fields_Func();
        pack_Gamma_Delta(blocklist_n, blocklist_p, Gx_F1D_packed_);
    };

    // F_0 → G(x_0), x_0.
    blocklist_p.set_zero_Gamma_Delta();
    blocklist_n.set_zero_Gamma_Delta();
    blocklist_p.add_Gamma_Delta_from_field(field_p, global_basis);
    blocklist_n.add_Gamma_Delta_from_field(field_n, global_basis);
    add_block_fields_Func();
    pack_Gamma_Delta(blocklist_n, blocklist_p, Gx_F1D_packed);
    x_F1D_packed.setZero();
    if (hasBlocking_B) {
        x_F1D_packed = Gx_F1D_packed;
        calc_Gx_Func(x_F1D_packed, Gx_F1D_packed);
    }

    // (x_0,G(x_0)) → Broyden history.
    const double mixingInitial_F = hasBlocking_B ? 0.20 : (hfbsettings.termSwitches.addFiniteRangeGogny_B ? 0.70 : 0.30);
    double alpha_F = hasBlocking_B ? 0.20 : 0.70;
    BroydenIterator broyden_(7, calc_Gx_Func, mixingInitial_F, x_F1D_packed, Gx_F1D_packed);
    observable_.update_observable(*this, activeBlockings_);
    print_summary(observable_, hfbsettings, 0, 0.0, mixingInitial_F);
    print_blocking(axialconfig, activeBlockings_);

    // x_i → x_{i+1}, ||G(x_i)-x_i||_∞ → ε_i.
    const int NiterationsMax_I = 100;
    const double alphaMax_F = 0.90;
    const double alphaMin_F = 0.20;
    double errorPrevious_F = 1.0;
    for (int iteration_I = 1; iteration_I <= NiterationsMax_I; ++iteration_I) {
        const double error_F = broyden_.iterate(calc_Gx_Func, alpha_F);
        observable_.update_observable(*this, activeBlockings_);
        print_summary(observable_, hfbsettings, iteration_I, error_F, alpha_F);
        print_blocking(axialconfig, activeBlockings_);
        if (std::isfinite(error_F) && error_F <= hfbsettings.accuracy_F) {break;}
        if (std::isfinite(error_F) && error_F < errorPrevious_F) {
            alpha_F = std::min(alphaMax_F, alpha_F * 1.10);
        } else {
            alpha_F = alphaMin_F;
            const bool tightenLambdaTolerance_B = lambdaTolerance_F > lambdaToleranceMin_F * (1.0 + 1.0e-12);
            if (iteration_I > 1 && tightenLambdaTolerance_B) {
                lambdaTolerance_F = std::max(lambdaToleranceMin_F, lambdaTolerance_F * 0.1);
                std::cout << "[AxialHFB] λ tolerance → " << std::scientific << std::setprecision(3) << lambdaTolerance_F << std::fixed << std::endl;
            }
        }
        errorPrevious_F = error_F;
    }
}
