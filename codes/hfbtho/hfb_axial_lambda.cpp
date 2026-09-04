/**
 * @file    hfb_axial_lambda.cpp
 * @author  cailea
 * @date    2026-09-03
 * @brief   Update axial HFB chemical potentials.
 */

#include "hfb_axial.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <map>

#include "root.hpp"

namespace {

// b → {ℋ_b,eigensolver_b}.
struct SearchLambdaWorkspace {
    bool isNeutron_B;
    const HFBSettings& hfbsettings;
    AxialHFBBlockList& blocklist;
    Eigen::VectorXd N_F1D_block{};
    const std::vector<AxialHFBBlocking>& activeBlockings;
    std::vector<Eigen::MatrixXd> Hhfb_F3D_block_2bsp_2bsp{};
    std::vector<Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>> eigensolvers_X1D_block{};
};

/**
 * @brief  Calculate particle number without density matrices.
 * @math   N(λ)=2Σ_μ[(1-f_μ)||V_μ||²+f_μ||U_μ||²]
 * @output Trial particle number.
 */
double calc_N_at_lambda(double lambda_F, SearchLambdaWorkspace& lambdaWorkspace_) {
    AxialHFBBlockList& blocklist_ = lambdaWorkspace_.blocklist;
    const HFBSettings& hfbsettings_ = lambdaWorkspace_.hfbsettings;
    const std::vector<AxialHFBBlocking>& activeBlockings_ = lambdaWorkspace_.activeBlockings;
    const bool isNeutron_B = lambdaWorkspace_.isNeutron_B;
    Eigen::VectorXd& N_F1D_block_ = lambdaWorkspace_.N_F1D_block;
    blocklist_.lambda_F = lambda_F;
    const int Nblock_I = static_cast<int>(blocklist_.blocks_X1D_block.size());
    constexpr double EspCutTolerance_F = 1.0e-6;
    static const double EspCutTail_F = std::log(1.0 / EspCutTolerance_F - 1.0) / 100.0;

    // λ → {(E_μ,U_μ,V_μ),N_b}_b.
    #pragma omp parallel for schedule(dynamic, 1)
    for (int block_I = 0; block_I < Nblock_I; ++block_I) {
        AxialHFBBlock& block_ = blocklist_.blocks_X1D_block[block_I];
        const int Nbsp_I = static_cast<int>(block_.labels_S1D_bsp.size());
        Eigen::MatrixXd& Hhfb_F2D_2bsp_2bsp = lambdaWorkspace_.Hhfb_F3D_block_2bsp_2bsp[block_I];
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>& eigensolver_ = lambdaWorkspace_.eigensolvers_X1D_block[block_I];

        // (Γ-λI,Δ) → ℋ.
        Hhfb_F2D_2bsp_2bsp.topLeftCorner(Nbsp_I, Nbsp_I) = block_.Gamma_F2D_bsp_bsp;
        Hhfb_F2D_2bsp_2bsp.topLeftCorner(Nbsp_I, Nbsp_I).diagonal().array() -= lambda_F;
        Hhfb_F2D_2bsp_2bsp.topRightCorner(Nbsp_I, Nbsp_I) = block_.Delta_F2D_bsp_bsp;
        Hhfb_F2D_2bsp_2bsp.bottomLeftCorner(Nbsp_I, Nbsp_I) = block_.Delta_F2D_bsp_bsp;
        Hhfb_F2D_2bsp_2bsp.bottomRightCorner(Nbsp_I, Nbsp_I) = -block_.Gamma_F2D_bsp_bsp;
        Hhfb_F2D_2bsp_2bsp.bottomRightCorner(Nbsp_I, Nbsp_I).diagonal().array() += lambda_F;

        // ℋ → (E_+,U,V).
        eigensolver_.compute(Hhfb_F2D_2bsp_2bsp);
        assert(eigensolver_.info() == Eigen::Success);
        const Eigen::VectorXd& eigenvalues_F1D_state = eigensolver_.eigenvalues();
        const Eigen::MatrixXd& eigenvectors_F2D_2bsp_state = eigensolver_.eigenvectors();
        assert(eigenvalues_F1D_state(Nbsp_I - 1) <= 0.0);
        assert(eigenvalues_F1D_state(Nbsp_I) >= 0.0);
        block_.Eqp_F1D_bqp = eigenvalues_F1D_state.tail(Nbsp_I);
        block_.U_F2D_bsp_bqp = eigenvectors_F2D_2bsp_state.topRightCorner(Nbsp_I, Nbsp_I);
        block_.V_F2D_bsp_bqp = eigenvectors_F2D_2bsp_state.bottomRightCorner(Nbsp_I, Nbsp_I);

        // (E_μ,U_μ,V_μ) → N_b.
        double Nblock_F = 0.0;
        for (int bqp_I = 0; bqp_I < Nbsp_I; ++bqp_I) {
            const double Eqp_F = block_.Eqp_F1D_bqp(bqp_I);
            const double Vnorm2_F = block_.V_F2D_bsp_bqp.col(bqp_I).squaredNorm();
            const double Unorm2_F = block_.U_F2D_bsp_bqp.col(bqp_I).squaredNorm();
            const double Esp_F = Eqp_F * (1.0 - 2.0 * Vnorm2_F) + lambda_F;
            if (hfbsettings_.useEspCut_B && Esp_F > hfbsettings_.EspCut_F + EspCutTail_F) {continue;}
            const double occupation_F = hfbsettings_.temperature_F > 0.0 ? 0.5 * (1.0 - std::tanh(0.5 * Eqp_F / hfbsettings_.temperature_F)) : 0.0;
            Nblock_F += 2.0 * ((1.0 - occupation_F) * Vnorm2_F + occupation_F * Unorm2_F);
        }
        N_F1D_block_(block_I) = Nblock_F;
    }

    double Ncalc_F = N_F1D_block_.sum();

    // EFA: N → N+||U_μ||²-||V_μ||².
    for (const AxialHFBBlocking& blocking_ : activeBlockings_) {
        if (blocking_.isNeutron_B != isNeutron_B) {continue;}
        const AxialHFBBlock& block_ = blocklist_.blocks_X1D_block[blocking_.block_I];
        int bestBqp_I = -1;
        double bestOverlap_F = -1.0;
        for (int bqp_I = 0; bqp_I < block_.U_F2D_bsp_bqp.cols(); ++bqp_I) {
            const double overlapCandidate_F = (blocking_.blockedU_F1D_bsp.array() * block_.U_F2D_bsp_bqp.col(bqp_I).array()).abs().sum() + (blocking_.blockedV_F1D_bsp.array() * block_.V_F2D_bsp_bqp.col(bqp_I).array()).abs().sum();
            if (!(overlapCandidate_F > bestOverlap_F)) {continue;}
            bestOverlap_F = overlapCandidate_F;
            bestBqp_I = bqp_I;
        }
        assert(bestBqp_I >= 0);
        Ncalc_F += block_.U_F2D_bsp_bqp.col(bestBqp_I).squaredNorm();
        Ncalc_F -= block_.V_F2D_bsp_bqp.col(bestBqp_I).squaredNorm();
    }
    return Ncalc_F;
}

} // namespace

void AxialHFB::update_blocklist_lambda(AxialHFBBlockList& blocklist_, int Ntarget_I, std::vector<AxialHFBBlocking>& activeBlockings_, bool isNeutron_B, double lambdaTolerance_F) {
    const double Ntolerance_F = std::max(1.0e-8, 1.0e-10 * std::max(1, Ntarget_I));
    const int Nblock_I = static_cast<int>(blocklist_.blocks_X1D_block.size());
    SearchLambdaWorkspace lambdaWorkspace_{isNeutron_B, hfbsettings, blocklist_, {}, activeBlockings_};
    lambdaWorkspace_.Hhfb_F3D_block_2bsp_2bsp.reserve(Nblock_I);
    lambdaWorkspace_.eigensolvers_X1D_block.reserve(Nblock_I);
    for (const AxialHFBBlock& block_ : blocklist_.blocks_X1D_block) {
        const int Nbsp_I = static_cast<int>(block_.labels_S1D_bsp.size());
        lambdaWorkspace_.Hhfb_F3D_block_2bsp_2bsp.emplace_back(2 * Nbsp_I, 2 * Nbsp_I);
        lambdaWorkspace_.eigensolvers_X1D_block.emplace_back(2 * Nbsp_I);
    }
    lambdaWorkspace_.N_F1D_block.resize(Nblock_I);
    lambdaWorkspace_.N_F1D_block.setZero();
    std::map<double, double> NerrorByLambda_Map{};

    // {B_μ^q} → (ρ_q,κ_q).
    const auto apply_active_blockings_Func = [&]() {
        for (AxialHFBBlocking& blocking_ : activeBlockings_) {
            if (blocking_.isNeutron_B != isNeutron_B) {continue;}
            blocking_.apply_blocking(blocklist_);
        }
    };

    // λ → N_q(λ)-N_q^{target}.
    const auto calc_Nerror_Func = [&](double lambda_F) {
        const auto NerrorIt_ = NerrorByLambda_Map.find(lambda_F);
        if (NerrorIt_ != NerrorByLambda_Map.end()) {return NerrorIt_->second;}
        const double Ncalc_F = calc_N_at_lambda(lambda_F, lambdaWorkspace_);
        const double Nerror_F = Ncalc_F - static_cast<double>(Ntarget_I);
        NerrorByLambda_Map.emplace(lambda_F, Nerror_F);
        return Nerror_F;
    };

    // λ_0 ← λ_previous or -7 MeV.
    assert(std::isfinite(blocklist_.lambda_F));
    const double lambdaInitial_F = blocklist_.lambda_F;
    const double NerrorInitial_F = calc_Nerror_Func(lambdaInitial_F);
    if (std::abs(NerrorInitial_F) <= Ntolerance_F) {
        blocklist_.lambda_F = lambdaInitial_F;
        blocklist_.update_UV_E_rho_kappa(hfbsettings);
        apply_active_blockings_Func();
        return;
    }

    // [λ_min,λ_max] ← [λ_0-2,λ_0+2].
    double lambdaMax_F = lambdaInitial_F + 2.0;
    double lambdaMin_F = lambdaInitial_F - 2.0;
    double NerrorMin_F = calc_Nerror_Func(lambdaMin_F);
    double NerrorMax_F = calc_Nerror_Func(lambdaMax_F);
    const double NerrorSlope_F = (NerrorMax_F - NerrorMin_F) / (lambdaMax_F - lambdaMin_F);

    // f(λ_min)f(λ_max)>0 → expand bracket.
    const bool shouldShiftLambdaMin_B = NerrorMin_F * NerrorSlope_F > 0.0;
    double& lambdaBound_F = shouldShiftLambdaMin_B ? lambdaMin_F : lambdaMax_F;
    double& NerrorBound_F = shouldShiftLambdaMin_B ? NerrorMin_F : NerrorMax_F;
    const double lambdaStep_F = shouldShiftLambdaMin_B ? -5.0 : 5.0;
    constexpr int NexpandMax_I = 100;
    int Nexpand_I = 0;
    while (NerrorMin_F * NerrorMax_F > 0.0 && NerrorSlope_F != 0.0 && Nexpand_I < NexpandMax_I) {
        lambdaBound_F += lambdaStep_F;
        NerrorBound_F = calc_Nerror_Func(lambdaBound_F);
        ++Nexpand_I;
    }
    assert(NerrorMin_F * NerrorMax_F <= 0.0);

    // f(λ_q)=0 → (U,V,E,ρ,κ)_q.
    const double lambdaRoot_F = root_brent(calc_Nerror_Func, lambdaMin_F, lambdaMax_F, lambdaTolerance_F);
    blocklist_.lambda_F = lambdaRoot_F;
    blocklist_.update_UV_E_rho_kappa(hfbsettings);
    apply_active_blockings_Func();
}
