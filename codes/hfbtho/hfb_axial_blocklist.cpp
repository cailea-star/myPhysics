/**
 * @file    hfb_axial_blocklist.cpp
 * @author  cailea
 * @date    2026-05-06
 * @brief   Implement axial HFB block-list updates.
 */

#include <algorithm>
#include <cassert>
#include <cmath>

#include <Eigen/Eigenvalues>

#include "hfb_axial.hpp"

void AxialHFBBlockList::set_zero_Gamma_Delta() {
    for (AxialHFBBlock& block_ : blocks_X1D_block) {
        block_.set_zero_Gamma_Delta();
    }
}

/**
 * @brief  Add local fields into every symmetry block.
 * @math   {Γ_q,Δ_q}_{block} → {Γ_q+Γ_q^{loc},Δ_q+Δ_q^{loc}}_{block}
 * @output Accumulated block fields.
 */
void AxialHFBBlockList::add_Gamma_Delta_from_field(const AxialHFBField& field_, const AxialBasis& global_basis_) {
    const int Nblock_I = static_cast<int>(blocks_X1D_block.size());

    #pragma omp parallel for schedule(static)
    for (int block_I = 0; block_I < Nblock_I; ++block_I) {
        blocks_X1D_block[block_I].add_Gamma_Delta_from_field(field_, global_basis_);
    }
}

/**
 * @brief  Calculate the canonical-density Lipkin-Nogami correction.
 * @math   λ₂=-G_{eff}S_N/S_D; E_{LN}=-4λ₂Σu²v²
 * @output Updated λ₂, E_LN, and Γ matrices.
 */
void AxialHFBBlockList::add_lipkin_nogami() {
    double Su1v3_F = 0.0;
    double Su2v2_F = 0.0;
    double Su3v1_F = 0.0;
    double Su4v4_F = 0.0;
    double Epair_F = 0.0;
    double DeltaRho_F = 0.0;
    double Nparticle_F = 0.0;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver_;
    const int Nblock_I = static_cast<int>(blocks_X1D_block.size());

    // {ρ,Δ,κ} → {S,N,E_pair}.
    for (int block_I = 0; block_I < Nblock_I; ++block_I) {
        const AxialHFBBlock& block_ = blocks_X1D_block[block_I];
        const int Nbsp_I = static_cast<int>(block_.labels_S1D_bsp.size());
        Nparticle_F += block_.rho_F2D_bsp_bsp.trace();
        Epair_F += block_.Delta_F2D_bsp_bsp.cwiseProduct(block_.kappa_F2D_bsp_bsp).sum();
        DeltaRho_F += block_.Delta_F2D_bsp_bsp.cwiseProduct(block_.rho_F2D_bsp_bsp).sum();
        eigensolver_.compute(block_.rho_F2D_bsp_bsp, Eigen::EigenvaluesOnly);
        assert(eigensolver_.info() == Eigen::Success);

        // ρ → {u,v}.
        for (int bsp_I = 0; bsp_I < Nbsp_I; ++bsp_I) {
            const double v2_F = std::clamp(eigensolver_.eigenvalues()(bsp_I), 0.0, 1.0);
            const double u2_F = 1.0 - v2_F;
            const double u_F = std::sqrt(u2_F);
            const double v_F = std::sqrt(v2_F);
            const double u3_F = u2_F * u_F;
            const double v3_F = v2_F * v_F;
            const double u4v4_F = u2_F * u2_F * v2_F * v2_F;
            Su1v3_F += u_F * v3_F;
            Su2v2_F += u2_F * v2_F;
            Su3v1_F += u3_F * v_F;
            Su4v4_F += u4v4_F;
        }
    }

    // {S,N,E_pair} → (λ₂,E_LN).
    const double DeltaAverage_F = -DeltaRho_F / Nparticle_F;
    const double Geff_F = DeltaAverage_F * DeltaAverage_F / Epair_F;
    const double lambda2_numer_F = 8.0 * (Su3v1_F * Su1v3_F - Su4v4_F);
    const double lambda2_denom_F = 32.0 * (Su2v2_F * Su2v2_F - Su4v4_F);
    lambda2_F = -Geff_F * lambda2_numer_F / lambda2_denom_F;
    if (!std::isfinite(lambda2_F)) {lambda2_F = 0.0;}
    if (lambda2_F >= 10.0) {lambda2_F = 4.0;}
    ELipkinNogami_F = -4.0 * lambda2_F * Su2v2_F;
    if (!std::isfinite(ELipkinNogami_F)) {ELipkinNogami_F = 0.0;}

    // λ₂ → Γ.
    #pragma omp parallel for schedule(static)
    for (int block_I = 0; block_I < Nblock_I; ++block_I) {
        blocks_X1D_block[block_I].add_lipkin_nogami(lambda2_F);
    }
}

/**
 * @brief  Solve every block by self-adjoint diagonalization.
 * @math   {ℋ_{block}} → {U,V,E,ρ,κ}_{block}
 * @output Updated block solutions.
 */
void AxialHFBBlockList::update_UV_E_rho_kappa(const HFBSettings& hfbsettings_) {
    const int Nblock_I = static_cast<int>(blocks_X1D_block.size());

    #pragma omp parallel for schedule(static)
    for (int block_I = 0; block_I < Nblock_I; ++block_I) {
        blocks_X1D_block[block_I].update_UV_E_rho_kappa(hfbsettings_, lambda_F);
    }
}

/**
 * @brief  Contract Gogny particle-hole matrix elements.
 * @math   Γ^q_{13}=Σ_{24}(\bar v^{same}_{12;34}ρ^q_{42}+\bar v^{cross}_{12;34}ρ^{\bar q}_{42})
 * @output Updated neutron and proton Γ matrices.
 */
void AxialHFBBlockList::add_Gamma_from_Gogny(AxialHFBBlockList& blocklist_p_, AxialHFBBlockList& blocklist_n_, const AxialGaussianGogny& gogny_) {
    assert(blocklist_n_.blocks_X1D_block.size() == blocklist_p_.blocks_X1D_block.size());
    const int Nblock_I = static_cast<int>(blocklist_n_.blocks_X1D_block.size());

    // (block_{13},bsp_1,bsp_3) → (Γ^n_{13},Γ^p_{13}).
    const auto add_Gamma_at_oneBlock_Func = [&](int block13_I, int bsp1_I, int bsp3_I) {
        AxialHFBBlock& block13n_ = blocklist_n_.blocks_X1D_block[block13_I];
        AxialHFBBlock& block13p_ = blocklist_p_.blocks_X1D_block[block13_I];
        const int sp1_I = block13n_.indices_I1D_bsp[bsp1_I];
        const int sp3_I = block13n_.indices_I1D_bsp[bsp3_I];
        const int offDiagonal_I = static_cast<int>(bsp1_I != bsp3_I);
        double Gamma13n_F = 0.0;
        double Gamma13p_F = 0.0;

        for (int block24_I = 0; block24_I < Nblock_I; ++block24_I) {
            const AxialHFBBlock& block24n_ = blocklist_n_.blocks_X1D_block[block24_I];
            const AxialHFBBlock& block24p_ = blocklist_p_.blocks_X1D_block[block24_I];
            assert(block24n_.indices_I1D_bsp.size() == block24p_.indices_I1D_bsp.size());
            const int Nbsp24_I = static_cast<int>(block24n_.indices_I1D_bsp.size());

            for (int bsp2_I = 0; bsp2_I < Nbsp24_I; ++bsp2_I) {
                for (int bsp4_I = 0; bsp4_I < Nbsp24_I; ++bsp4_I) {
                    const int sp2_I = block24n_.indices_I1D_bsp[bsp2_I];
                    const int sp4_I = block24n_.indices_I1D_bsp[bsp4_I];
                    const AxialGaussianGogny::GognyElements gognyElements1234_ = gogny_.read_v(sp1_I, sp2_I, sp3_I, sp4_I);

                    // ρ^-_{42}=η_2η_4ρ^+_{42}.
                    const double rho42Posn_F = block24n_.rho_F2D_bsp_bsp(bsp4_I, bsp2_I);
                    const double rho42Posp_F = block24p_.rho_F2D_bsp_bsp(bsp4_I, bsp2_I);
                    const int eta24_I = block24n_.labels_S1D_bsp[bsp2_I].twoSigma_I * block24n_.labels_S1D_bsp[bsp4_I].twoSigma_I;
                    const double rho42Negn_F = eta24_I * rho42Posn_F;
                    const double rho42Negp_F = eta24_I * rho42Posp_F;

                    // (v̄^{++},v̄^{+-})·(ρ^+,ρ^-) → Γ^+.
                    Gamma13n_F += gognyElements1234_.vSamePosPosPosPos_F * rho42Posn_F;
                    Gamma13n_F += gognyElements1234_.vSamePosNegPosNeg_F * rho42Negn_F;
                    Gamma13n_F += gognyElements1234_.vCrossPosPosPosPos_F * rho42Posp_F;
                    Gamma13n_F += gognyElements1234_.vCrossPosNegPosNeg_F * rho42Negp_F;

                    Gamma13p_F += gognyElements1234_.vSamePosPosPosPos_F * rho42Posp_F;
                    Gamma13p_F += gognyElements1234_.vSamePosNegPosNeg_F * rho42Negp_F;
                    Gamma13p_F += gognyElements1234_.vCrossPosPosPosPos_F * rho42Posn_F;
                    Gamma13p_F += gognyElements1234_.vCrossPosNegPosNeg_F * rho42Negn_F;
                }
            }
        }

        // Γ_{31}=Γ_{13}.
        block13n_.Gamma_F2D_bsp_bsp(bsp1_I, bsp3_I) += Gamma13n_F;
        block13p_.Gamma_F2D_bsp_bsp(bsp1_I, bsp3_I) += Gamma13p_F;
        block13n_.Gamma_F2D_bsp_bsp(bsp3_I, bsp1_I) += offDiagonal_I * Gamma13n_F;
        block13p_.Gamma_F2D_bsp_bsp(bsp3_I, bsp1_I) += offDiagonal_I * Gamma13p_F;
    };

    #pragma omp parallel
    {
        #pragma omp single
        for (int block13_I = 0; block13_I < Nblock_I; ++block13_I) {
            const AxialHFBBlock& block13n_ = blocklist_n_.blocks_X1D_block[block13_I];
            const AxialHFBBlock& block13p_ = blocklist_p_.blocks_X1D_block[block13_I];
            assert(block13n_.indices_I1D_bsp.size() == block13p_.indices_I1D_bsp.size());
            const int Nbsp13_I = static_cast<int>(block13n_.indices_I1D_bsp.size());

            // Γ^q_{13}=Σ_{24}\bar v^q_{12;34}ρ^q_{42}.
            for (int bsp1_I = 0; bsp1_I < Nbsp13_I; ++bsp1_I) {
                for (int bsp3_I = 0; bsp3_I <= bsp1_I; ++bsp3_I) {
                    #pragma omp task firstprivate(block13_I, bsp1_I, bsp3_I)
                    {
                        add_Gamma_at_oneBlock_Func(block13_I, bsp1_I, bsp3_I);
                    }
                }
            }
        }
    }
}

/**
 * @brief  Contract Gogny pairing matrix elements.
 * @math   Δ^q_{12}=Σ_{34}\bar v^{same,+-+-}_{12;34}κ^q_{34}
 * @output Updated neutron and proton Δ matrices.
 */
void AxialHFBBlockList::add_Delta_from_Gogny(AxialHFBBlockList& blocklist_p_, AxialHFBBlockList& blocklist_n_, const AxialGaussianGogny& gogny_) {
    assert(blocklist_n_.blocks_X1D_block.size() == blocklist_p_.blocks_X1D_block.size());
    const int Nblock_I = static_cast<int>(blocklist_n_.blocks_X1D_block.size());

    // (block_{12},bsp_1,bsp_2) → (Δ^n_{12},Δ^p_{12}).
    const auto add_Delta_at_oneBlock_Func = [&](int block12_I, int bsp1_I, int bsp2_I) {
        AxialHFBBlock& block12n_ = blocklist_n_.blocks_X1D_block[block12_I];
        AxialHFBBlock& block12p_ = blocklist_p_.blocks_X1D_block[block12_I];
        const int sp1_I = block12n_.indices_I1D_bsp[bsp1_I];
        const int sp2_I = block12n_.indices_I1D_bsp[bsp2_I];
        const int offDiagonal_I = static_cast<int>(bsp1_I != bsp2_I);
        double Delta12n_F = 0.0;
        double Delta12p_F = 0.0;

        for (int block34_I = 0; block34_I < Nblock_I; ++block34_I) {
            const AxialHFBBlock& block34n_ = blocklist_n_.blocks_X1D_block[block34_I];
            const AxialHFBBlock& block34p_ = blocklist_p_.blocks_X1D_block[block34_I];
            assert(block34n_.indices_I1D_bsp.size() == block34p_.indices_I1D_bsp.size());
            const int Nbsp34_I = static_cast<int>(block34n_.indices_I1D_bsp.size());

            for (int bsp3_I = 0; bsp3_I < Nbsp34_I; ++bsp3_I) {
                for (int bsp4_I = 0; bsp4_I < Nbsp34_I; ++bsp4_I) {
                    const int sp3_I = block34n_.indices_I1D_bsp[bsp3_I];
                    const int sp4_I = block34n_.indices_I1D_bsp[bsp4_I];
                    const AxialGaussianGogny::GognyElements gognyElements1234_ = gogny_.read_v(sp1_I, sp2_I, sp3_I, sp4_I);
                    const double kappa34n_F = block34n_.kappa_F2D_bsp_bsp(bsp3_I, bsp4_I);
                    const double kappa34p_F = block34p_.kappa_F2D_bsp_bsp(bsp3_I, bsp4_I);
                    const double kappa43n_F = block34n_.kappa_F2D_bsp_bsp(bsp4_I, bsp3_I);
                    const double kappa43p_F = block34p_.kappa_F2D_bsp_bsp(bsp4_I, bsp3_I);
                    const int eta24_I = block12n_.labels_S1D_bsp[bsp2_I].twoSigma_I * block34n_.labels_S1D_bsp[bsp4_I].twoSigma_I;
                    const int eta23_I = block12n_.labels_S1D_bsp[bsp2_I].twoSigma_I * block34n_.labels_S1D_bsp[bsp3_I].twoSigma_I;
                    Delta12n_F += 0.5 * (eta24_I * gognyElements1234_.vSamePosNegPosNeg_F * kappa34n_F - eta23_I * gognyElements1234_.vSamePosNegNegPos_F * kappa43n_F);
                    Delta12p_F += 0.5 * (eta24_I * gognyElements1234_.vSamePosNegPosNeg_F * kappa34p_F - eta23_I * gognyElements1234_.vSamePosNegNegPos_F * kappa43p_F);
                }
            }
        }

        // Δ_{21}=Δ_{12}.
        block12n_.Delta_F2D_bsp_bsp(bsp1_I, bsp2_I) += Delta12n_F;
        block12p_.Delta_F2D_bsp_bsp(bsp1_I, bsp2_I) += Delta12p_F;
        block12n_.Delta_F2D_bsp_bsp(bsp2_I, bsp1_I) += offDiagonal_I * Delta12n_F;
        block12p_.Delta_F2D_bsp_bsp(bsp2_I, bsp1_I) += offDiagonal_I * Delta12p_F;
    };

    #pragma omp parallel
    {
        #pragma omp single
        for (int block12_I = 0; block12_I < Nblock_I; ++block12_I) {
            const AxialHFBBlock& block12n_ = blocklist_n_.blocks_X1D_block[block12_I];
            const AxialHFBBlock& block12p_ = blocklist_p_.blocks_X1D_block[block12_I];
            assert(block12n_.indices_I1D_bsp.size() == block12p_.indices_I1D_bsp.size());
            const int Nbsp12_I = static_cast<int>(block12n_.indices_I1D_bsp.size());

            // Δ^q_{12}=Σ_{34}\bar v^{same,+-+-}_{12;34}κ^q_{34}.
            for (int bsp1_I = 0; bsp1_I < Nbsp12_I; ++bsp1_I) {
                for (int bsp2_I = 0; bsp2_I <= bsp1_I; ++bsp2_I) {
                    #pragma omp task firstprivate(block12_I, bsp1_I, bsp2_I)
                    {
                        add_Delta_at_oneBlock_Func(block12_I, bsp1_I, bsp2_I);
                    }
                }
            }
        }
    }
}

/**
 * @brief  Contract Coulomb particle-hole matrix elements.
 * @math   Γ^p_{13}=Σ_{24}(\bar v^{C,++}_{12;34}ρ^{p,+}_{42}+\bar v^{C,+-}_{12;34}ρ^{p,-}_{42})
 * @output Updated proton Γ matrices.
 */
void AxialHFBBlockList::add_coulomb_from_Gaussian(AxialHFBBlockList& blocklist_p_, const AxialGaussianCoulomb& coulomb_) {
    const int Nblock_I = static_cast<int>(blocklist_p_.blocks_X1D_block.size());

    // (block_{13},bsp_1,bsp_3) → Γ^{C,p}_{13}.
    const auto add_coulomb_at_oneBlock_Func = [&](int block13_I, int bsp1_I, int bsp3_I) {
        AxialHFBBlock& block13p_ = blocklist_p_.blocks_X1D_block[block13_I];
        const int sp1_I = block13p_.indices_I1D_bsp[bsp1_I];
        const int sp3_I = block13p_.indices_I1D_bsp[bsp3_I];
        const int offDiagonal_I = static_cast<int>(bsp1_I != bsp3_I);
        double Gamma13p_F = 0.0;

        for (int block24_I = 0; block24_I < Nblock_I; ++block24_I) {
            const AxialHFBBlock& block24p_ = blocklist_p_.blocks_X1D_block[block24_I];
            const int Nbsp24_I = static_cast<int>(block24p_.indices_I1D_bsp.size());

            for (int bsp2_I = 0; bsp2_I < Nbsp24_I; ++bsp2_I) {
                for (int bsp4_I = 0; bsp4_I < Nbsp24_I; ++bsp4_I) {
                    const int sp2_I = block24p_.indices_I1D_bsp[bsp2_I];
                    const int sp4_I = block24p_.indices_I1D_bsp[bsp4_I];
                    const AxialGaussianCoulomb::GammaElements coulombElements1234_ = coulomb_.read_v(sp1_I, sp2_I, sp3_I, sp4_I);

                    // ρ^-_{42}=η_2η_4ρ^+_{42}.
                    const double rho42Posp_F = block24p_.rho_F2D_bsp_bsp(bsp4_I, bsp2_I);
                    const int eta24_I = block24p_.labels_S1D_bsp[bsp2_I].twoSigma_I * block24p_.labels_S1D_bsp[bsp4_I].twoSigma_I;
                    const double rho42Negp_F = eta24_I * rho42Posp_F;

                    // (v̄^{++},v̄^{+-})·(ρ^+,ρ^-) → Γ^+.
                    Gamma13p_F += coulombElements1234_.vPosPosPosPos_F * rho42Posp_F;
                    Gamma13p_F += coulombElements1234_.vPosNegPosNeg_F * rho42Negp_F;
                }
            }
        }

        // Γ_{31}=Γ_{13}.
        block13p_.Gamma_F2D_bsp_bsp(bsp1_I, bsp3_I) += Gamma13p_F;
        block13p_.Gamma_F2D_bsp_bsp(bsp3_I, bsp1_I) += offDiagonal_I * Gamma13p_F;
    };

    #pragma omp parallel for schedule(static)
    for (int block13_I = 0; block13_I < Nblock_I; ++block13_I) {
        const AxialHFBBlock& block13p_ = blocklist_p_.blocks_X1D_block[block13_I];
        const int Nbsp13_I = static_cast<int>(block13p_.indices_I1D_bsp.size());

        // Γ^p_{13}=Σ_{24}(v̄^{C,++}_{12;34}ρ^{p,+}_{42}+v̄^{C,+-}_{12;34}ρ^{p,-}_{42}).
        for (int bsp1_I = 0; bsp1_I < Nbsp13_I; ++bsp1_I) {
            for (int bsp3_I = 0; bsp3_I <= bsp1_I; ++bsp3_I) {
                add_coulomb_at_oneBlock_Func(block13_I, bsp1_I, bsp3_I);
            }
        }
    }
}
