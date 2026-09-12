/**
 * @file    hfb_cylindrical_element.cpp
 * @author  cailea
 * @date    2026-09-12
 * @brief   Assemble cylindrical HFB one-body and pairing matrices.
 */

#include "hfb_cylindrical.hpp"

void HFBKramersNucleusCylindrical::add_Gamma_Delta_from_field() {
    assert(hfb_neutron.Nblock_I == hfb_proton.Nblock_I);
    assert(hfb_neutron.Nblock_I == static_cast<int>(cylindricalsetting.labels_S2D_block_bsp.size()));

    // (block,q) → Γ_q⁺⁺, Δ_q⁺⁻; negative branches follow time reversal.
    const auto add_block_Func = [&](HFBKramersBlockField& block_, const HFBCylindricalField& field_, int block_I) {
        const auto& labels_S1D_bsp = cylindricalsetting.labels_S2D_block_bsp[block_I];
        const auto& indices_I1D_bsp = cylindricalsetting.indices_I2D_block_bsp[block_I];
        const auto& indices_I1D_bup = density_neutron.bup_I2D_block_bup[block_I];
        const auto& indices_I1D_bdn = density_neutron.bdn_I2D_block_bdn[block_I];
        const int twoOmega_I = labels_S1D_bsp.front().twoOmega_I;
        // Real time-even fields; local spin-singlet pairing.
        const int Nz_I = static_cast<int>(field_.vcent_F2D_z_r.rows());
        const int Nr_I = static_cast<int>(field_.vcent_F2D_z_r.cols());
        const int Nbsp_I = static_cast<int>(labels_S1D_bsp.size());
        const int LambdaUp_I = (twoOmega_I - 1) / 2;
        const int LambdaDown_I = (twoOmega_I + 1) / 2;
        const int spBegin_I = indices_I1D_bsp.front();
        for (int bsp_I = 0; bsp_I < Nbsp_I; ++bsp_I) {assert(indices_I1D_bsp[bsp_I] == spBegin_I + bsp_I);}
        const Eigen::VectorXd rInv_F1D_r = cylindricalbasis.r_F1D_r.array().inverse().matrix();
        Eigen::VectorXd laplacianPhi_F1D_bsp(Nbsp_I);

        struct FieldValues {
            double vcent_F = 0.0;
            double vmass_F = 0.0;
            double vpair_F = 0.0;
            double vD2_F = 0.0;
            double vDr_F = 0.0;
            double vDz_F = 0.0;
            double vdJ_F = 0.0;
            double vJphiz_F = 0.0;
            double vJzphi_F = 0.0;
            double vJphir_F = 0.0;
            double vJrphi_F = 0.0;

            FieldValues(const HFBCylindricalField& field_, int z_I, int r_I) {
                vcent_F = field_.vcent_F2D_z_r(z_I, r_I);
                vmass_F = field_.vmass_F2D_z_r(z_I, r_I);
                vpair_F = field_.vpair_F2D_z_r(z_I, r_I);
                vD2_F = field_.vD2_F2D_z_r(z_I, r_I);
                vDr_F = field_.vDr_F2D_z_r(z_I, r_I);
                vDz_F = field_.vDz_F2D_z_r(z_I, r_I);
                vdJ_F = field_.vdJ_F2D_z_r(z_I, r_I);
                vJphiz_F = field_.vJphiz_F2D_z_r(z_I, r_I);
                vJzphi_F = field_.vJzphi_F2D_z_r(z_I, r_I);
                vJphir_F = field_.vJphir_F2D_z_r(z_I, r_I);
                vJrphi_F = field_.vJrphi_F2D_z_r(z_I, r_I);
            }
        };

        // Σ_1=Σ_2 → (Γ_12,Δ_12).
        const auto add_same_spin_at_onePoint_Func = [&](const FieldValues& fieldValues_, int z_I, int r_I, const std::vector<int>& indices_I1D_bsp_, int Lambda_I, double spinSign_F) {
            const int Nspin_I = static_cast<int>(indices_I1D_bsp_.size());
            const double w_F = cylindricalbasis.w_F2D_z_r(z_I, r_I);
            const double rInv_F = rInv_F1D_r(r_I);
            const double rInv2_F = rInv_F * rInv_F;
            const double Lambda2_F = Lambda_I * Lambda_I;
            const Eigen::Map<const Eigen::VectorXd> phi_F1D_bsp(&cylindricalbasis.phi_F3D_sp_z_r(spBegin_I, z_I, r_I), Nbsp_I);
            const Eigen::Map<const Eigen::VectorXd> dphidr_F1D_bsp(&cylindricalbasis.dphidr_F3D_sp_z_r(spBegin_I, z_I, r_I), Nbsp_I);
            const Eigen::Map<const Eigen::VectorXd> dphidz_F1D_bsp(&cylindricalbasis.dphidz_F3D_sp_z_r(spBegin_I, z_I, r_I), Nbsp_I);
            const Eigen::Map<const Eigen::VectorXd> ddphidr_F1D_bsp(&cylindricalbasis.ddphidr_F3D_sp_z_r(spBegin_I, z_I, r_I), Nbsp_I);
            const Eigen::Map<const Eigen::VectorXd> ddphidz_F1D_bsp(&cylindricalbasis.ddphidz_F3D_sp_z_r(spBegin_I, z_I, r_I), Nbsp_I);

            // ∇²φ_bsp = ∂²_rφ+r⁻¹∂_rφ+∂²_zφ-r⁻²Λ²φ.
            for (int position_I = 0; position_I < Nspin_I; ++position_I) {
                const int bsp_I = indices_I1D_bsp_[position_I];
                laplacianPhi_F1D_bsp(bsp_I) = ddphidr_F1D_bsp(bsp_I) + rInv_F * dphidr_F1D_bsp(bsp_I) + ddphidz_F1D_bsp(bsp_I) - rInv2_F * Lambda2_F * phi_F1D_bsp(bsp_I);
            }

            for (int position2_I = 0; position2_I < Nspin_I; ++position2_I) {
                for (int position1_I = position2_I; position1_I < Nspin_I; ++position1_I) {
                    const int bsp1_I = indices_I1D_bsp_[position1_I];
                    const int bsp2_I = indices_I1D_bsp_[position2_I];

                    // (φ_1,φ_2,Λ) → transition densities.
                    const double phi1_F = phi_F1D_bsp(bsp1_I);
                    const double phi2_F = phi_F1D_bsp(bsp2_I);
                    const double dphidr1_F = dphidr_F1D_bsp(bsp1_I);
                    const double dphidr2_F = dphidr_F1D_bsp(bsp2_I);
                    const double dphidz1_F = dphidz_F1D_bsp(bsp1_I);
                    const double dphidz2_F = dphidz_F1D_bsp(bsp2_I);
                    const double laplacianPhi1_F = laplacianPhi_F1D_bsp(bsp1_I);
                    const double laplacianPhi2_F = laplacianPhi_F1D_bsp(bsp2_I);
                    const double rho12_F = phi1_F * phi2_F;
                    const double tau12_F = dphidr1_F * dphidr2_F + dphidz1_F * dphidz2_F + rInv2_F * Lambda2_F * rho12_F;
                    const double rhoD2_12_F = laplacianPhi1_F * phi2_F + phi1_F * laplacianPhi2_F + 2.0 * tau12_F;
                    const double rhoDr_12_F = dphidr1_F * phi2_F + phi1_F * dphidr2_F;
                    const double rhoDz_12_F = dphidz1_F * phi2_F + phi1_F * dphidz2_F;
                    const double Jphiz12_F = spinSign_F * Lambda_I * rInv_F * rho12_F;
                    const double dJ12_F = spinSign_F * Lambda_I * rInv_F * rhoDr_12_F;

                    // ({v},{ρ_12,τ_12,J_12}) → (Γ_12,Δ_12).
                    // Δ⁺⁻_{12}=η_2∫v_pair φ_1φ_2; η_1=η_2=2Σ.
                    const double Delta12_F = spinSign_F * fieldValues_.vpair_F * rho12_F;
                    const double Gamma12_F = fieldValues_.vcent_F * rho12_F + fieldValues_.vmass_F * tau12_F + fieldValues_.vD2_F * rhoD2_12_F + fieldValues_.vDr_F * rhoDr_12_F + fieldValues_.vDz_F * rhoDz_12_F + fieldValues_.vdJ_F * dJ12_F + fieldValues_.vJphiz_F * Jphiz12_F;
                    block_.GammaPosPos_F2D_bsp_bsp(bsp1_I, bsp2_I) += Gamma12_F * w_F;
                    block_.DeltaPosNeg_F2D_bsp_bsp(bsp1_I, bsp2_I) += Delta12_F * w_F;
                    const int offDiagonal_I = static_cast<int>(bsp1_I != bsp2_I);

                    block_.GammaPosPos_F2D_bsp_bsp(bsp2_I, bsp1_I) += Gamma12_F * w_F * offDiagonal_I;
                    block_.DeltaPosNeg_F2D_bsp_bsp(bsp2_I, bsp1_I) += Delta12_F * w_F * offDiagonal_I;
                }
            }
        };

        // Σ_↑≠Σ_↓ → Γ_↑↓.
        const auto add_opposite_spin_at_onePoint_Func = [&](const FieldValues& fieldValues_, int z_I, int r_I) {
            const double w_F = cylindricalbasis.w_F2D_z_r(z_I, r_I);
            const double rInv_F = rInv_F1D_r(r_I);
            const Eigen::Map<const Eigen::VectorXd> phi_F1D_bsp(&cylindricalbasis.phi_F3D_sp_z_r(spBegin_I, z_I, r_I), Nbsp_I);
            const Eigen::Map<const Eigen::VectorXd> dphidr_F1D_bsp(&cylindricalbasis.dphidr_F3D_sp_z_r(spBegin_I, z_I, r_I), Nbsp_I);
            const Eigen::Map<const Eigen::VectorXd> dphidz_F1D_bsp(&cylindricalbasis.dphidz_F3D_sp_z_r(spBegin_I, z_I, r_I), Nbsp_I);

            for (int up_I = 0; up_I < static_cast<int>(indices_I1D_bup.size()); ++up_I) {
                for (int down_I = 0; down_I < static_cast<int>(indices_I1D_bdn.size()); ++down_I) {
                    const int bspUp_I = indices_I1D_bup[up_I];
                    const int bspDown_I = indices_I1D_bdn[down_I];

                    // (φ_↑,φ_↓,Λ_↑,Λ_↓) → J_↑↓.
                    const double phiUp_F = phi_F1D_bsp(bspUp_I);
                    const double phiDown_F = phi_F1D_bsp(bspDown_I);
                    const double dphidrUp_F = dphidr_F1D_bsp(bspUp_I);
                    const double dphidrDown_F = dphidr_F1D_bsp(bspDown_I);
                    const double dphidzUp_F = dphidz_F1D_bsp(bspUp_I);
                    const double dphidzDown_F = dphidz_F1D_bsp(bspDown_I);
                    // δE=Γ_↑↓δρ_↓↑+Γ_↓↑δρ_↑↓; Γ_↑↓=Γ_↓↑.
                    const double JrphiUpDown_F = 0.5 * (dphidrUp_F * phiDown_F - phiUp_F * dphidrDown_F);
                    const double JzphiUpDown_F = 0.5 * (dphidzUp_F * phiDown_F - phiUp_F * dphidzDown_F);
                    const double JphirUpDown_F = 0.5 * twoOmega_I * rInv_F * phiUp_F * phiDown_F;
                    const double dJUpDown_F = dphidrUp_F * dphidzDown_F - dphidzUp_F * dphidrDown_F - rInv_F * (LambdaUp_I * phiUp_F * dphidzDown_F + LambdaDown_I * phiDown_F * dphidzUp_F);

                    // ({v},{J_↑↓}) → Γ_↑↓.
                    const double GammaUpDown_F = fieldValues_.vdJ_F * dJUpDown_F + fieldValues_.vJzphi_F * JzphiUpDown_F + fieldValues_.vJphir_F * JphirUpDown_F + fieldValues_.vJrphi_F * JrphiUpDown_F;
                    block_.GammaPosPos_F2D_bsp_bsp(bspUp_I, bspDown_I) += GammaUpDown_F * w_F;
                    block_.GammaPosPos_F2D_bsp_bsp(bspDown_I, bspUp_I) += GammaUpDown_F * w_F;

                }
            }
        };

        // (r,z) → (Γ,Δ).
        for (int r_I = 0; r_I < Nr_I; ++r_I) {
            for (int z_I = 0; z_I < Nz_I; ++z_I) {
                // (z,r) → {v}.
                const FieldValues fieldValues_(field_, z_I, r_I);
                add_same_spin_at_onePoint_Func(fieldValues_, z_I, r_I, indices_I1D_bup, LambdaUp_I, 1.0);
                add_same_spin_at_onePoint_Func(fieldValues_, z_I, r_I, indices_I1D_bdn, LambdaDown_I, -1.0);
                add_opposite_spin_at_onePoint_Func(fieldValues_, z_I, r_I);
            }
        }
    };

    #pragma omp parallel for schedule(static)
    for (int block_I = 0; block_I < hfb_neutron.Nblock_I; ++block_I) {
        add_block_Func(hfb_neutron.fields[block_I], field_neutron, block_I);
        add_block_Func(hfb_proton.fields[block_I], field_proton, block_I);
    }
}

void HFBKramersNucleusCylindrical::add_Gamma_from_Gogny() {
    assert(hfb_neutron.Nbsp_I1D_block == hfb_proton.Nbsp_I1D_block);

    // (block13,block24,bsp1,bsp2,bsp3,bsp4) → v̄_{12;34}.
    const auto read_element_Func = [&](int block13_I, int block24_I, int bsp1_I, int bsp2_I, int bsp3_I, int bsp4_I) {
        const auto& indices13_I1D_bsp = cylindricalsetting.indices_I2D_block_bsp[block13_I];
        const auto& indices24_I1D_bsp = cylindricalsetting.indices_I2D_block_bsp[block24_I];
        return gaussian_gogny.read_v(indices13_I1D_bsp[bsp1_I], indices24_I1D_bsp[bsp2_I], indices13_I1D_bsp[bsp3_I], indices24_I1D_bsp[bsp4_I]);
    };

    const HFBKramersNucleus::GammaElementFunc read_joint_Func = [&](int block13_I, int block24_I, int bsp1_I, int bsp2_I, int bsp3_I, int bsp4_I) {
        const auto elements = read_element_Func(block13_I, block24_I, bsp1_I, bsp2_I, bsp3_I, bsp4_I);
        return HFBKramersNucleus::Element{elements.vSamePosPosPosPos_F, elements.vSamePosNegPosNeg_F, elements.vCrossPosPosPosPos_F, elements.vCrossPosNegPosNeg_F};
    };
    add_Gamma_from_Element(read_joint_Func);
}

void HFBKramersNucleusCylindrical::add_Delta_from_Gogny() {
    assert(hfb_neutron.Nbsp_I1D_block == hfb_proton.Nbsp_I1D_block);

    // (block12,block34,bsp1,bsp2,bsp3,bsp4) → v̄_same⁺⁻⁺⁻_{12;34}.
    const HFBKramersNucleus::DeltaElementFunc read_element_Func = [&](int block12_I, int block34_I, int bsp1_I, int bsp2_I, int bsp3_I, int bsp4_I) {
        const auto& indices12_I1D_bsp = cylindricalsetting.indices_I2D_block_bsp[block12_I];
        const auto& indices34_I1D_bsp = cylindricalsetting.indices_I2D_block_bsp[block34_I];
        const auto elements = gaussian_gogny.read_v(indices12_I1D_bsp[bsp1_I], indices12_I1D_bsp[bsp2_I], indices34_I1D_bsp[bsp3_I], indices34_I1D_bsp[bsp4_I]);
        return HFBKramersNucleus::Element{elements.vSamePosPosPosPos_F, elements.vSamePosNegPosNeg_F, elements.vCrossPosPosPosPos_F, elements.vCrossPosNegPosNeg_F};
    };

    // Δ_q⁺⁻ += v̄_same⁺⁻⁺⁻ κ_q⁺⁻; q ∈ {n,p}.
    add_Delta_from_Element(read_element_Func);
}

void HFBKramersNucleusCylindrical::add_Gamma_from_Coulomb() {
    // (block13,block24,bsp1,bsp2,bsp3,bsp4) → v̄_C,12;34.
    const HFBKramers::GammaElementFunc read_element_Func = [&](int block13_I, int block24_I, int bsp1_I, int bsp2_I, int bsp3_I, int bsp4_I) {
        const auto& indices13_I1D_bsp = cylindricalsetting.indices_I2D_block_bsp[block13_I];
        const auto& indices24_I1D_bsp = cylindricalsetting.indices_I2D_block_bsp[block24_I];
        const auto elements = gaussian_coulomb.read_v(indices13_I1D_bsp[bsp1_I], indices24_I1D_bsp[bsp2_I], indices13_I1D_bsp[bsp3_I], indices24_I1D_bsp[bsp4_I]);
        return HFBKramers::Element{elements.vPosPosPosPos_F, elements.vPosNegPosNeg_F};
    };
    hfb_proton.add_Gamma_from_Element(read_element_Func);
}

void HFBKramersNucleusCylindrical::update_Gamma_Delta() {
    const int TargetA_I = hfb_neutron.TargetN_I + hfb_proton.TargetN_I;
    assert(hfb_neutron.TargetN_I >= 0 && hfb_proton.TargetN_I >= 0 && TargetA_I > 0);
    const EDFParamsSkyrme active_edf_ = hfbedfsetting.make_active_edf(edf_skyrme, TargetA_I);

    // Enabled interactions → reusable kernels.
    if (hfbedfsetting.termSwitches.addFiniteRangeGogny_B) {gaussian_gogny.build_tables();}
    if (hfbedfsetting.termSwitches.addFiniteRangeCoulomb_B) {gaussian_coulomb.build_tables();}
    if (hfbedfsetting.termSwitches.addLocalCoulomb_B) {coulomb_field.build(cylindricalsetting.useParity_B, active_edf_.e2charg_F);}

    std::vector<Eigen::MatrixXd> rhoPosPos_F3D_block_bsp_bsp{};
    std::vector<Eigen::MatrixXd> rhoNegNeg_F3D_block_bsp_bsp{};
    std::vector<Eigen::MatrixXd> kappaPosNeg_F3D_block_bsp_bsp{};
    // (ρ⁺⁺,κ⁺⁻,η) → (ρ⁺⁺,Dρ⁺⁺D,κ⁺⁻) → coordinate densities.
    const auto update_density_Func = [&](const HFBKramers& hfb_, CylindricalDensity& density_) {
        rhoPosPos_F3D_block_bsp_bsp.resize(hfb_.Nblock_I);
        rhoNegNeg_F3D_block_bsp_bsp.resize(hfb_.Nblock_I);
        kappaPosNeg_F3D_block_bsp_bsp.resize(hfb_.Nblock_I);
        for (int block_I = 0; block_I < hfb_.Nblock_I; ++block_I) {
            const auto& solution = hfb_.solutions[block_I];
            rhoPosPos_F3D_block_bsp_bsp[block_I] = solution.rhoPosPos_F2D_bsp_bsp;
            rhoNegNeg_F3D_block_bsp_bsp[block_I].noalias() = solution.eta_F1D_bsp.asDiagonal() * solution.rhoPosPos_F2D_bsp_bsp * solution.eta_F1D_bsp.asDiagonal();
            kappaPosNeg_F3D_block_bsp_bsp[block_I] = solution.kappaPosNeg_F2D_bsp_bsp;
        }
        density_.update_density(rhoPosPos_F3D_block_bsp_bsp, rhoNegNeg_F3D_block_bsp_bsp, kappaPosNeg_F3D_block_bsp_bsp);
    };
    update_density_Func(hfb_neutron, density_neutron);
    update_density_Func(hfb_proton, density_proton);

    // (D_n,D_p) → local fields; active EDF contains term switches.
    field_neutron.set_zero();
    field_proton.set_zero();
    HFBCylindricalField::add_nuclei_fields(field_proton, field_neutron, density_proton, density_neutron, active_edf_);
    if (hfbedfsetting.termSwitches.addLocalCoulomb_B) {HFBCylindricalField::add_coulomb_field(field_proton, density_proton, coulomb_field, active_edf_);}
    if (hfbedfsetting.termSwitches.addLocalPair_B) {HFBCylindricalField::add_pairing_fields(field_proton, field_neutron, density_proton, density_neutron, active_edf_, hfbedfsetting, hfb_neutron.lambda_F, hfb_proton.lambda_F);}

    // Rebuild (Γ,Δ); clear previous LN diagnostics.
    for (HFBKramers* hfb_ : {&hfb_neutron, &hfb_proton}) {
        hfb_->lambda2_F = 0.0;
        hfb_->ELipkinNogami_F = 0.0;
        for (auto& field_ : hfb_->fields) {
            field_.GammaPosPos_F2D_bsp_bsp.setZero();
            field_.DeltaPosNeg_F2D_bsp_bsp.setZero();
        }
    }
    add_Gamma_Delta_from_field();
    if (hfbedfsetting.termSwitches.addFiniteRangeGogny_B) {
        add_Gamma_from_Gogny();
        add_Delta_from_Gogny();
    }
    if (hfbedfsetting.termSwitches.addFiniteRangeCoulomb_B) {add_Gamma_from_Coulomb();}
    if (hfbedfsetting.useLipkinNogami_B) {add_Gamma_from_lipkin_nogami();}
}
