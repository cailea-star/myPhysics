/**
 * @file    test_bifold.cpp
 * @author  cailea
 * @date    2026-08-23
 * @brief   Validate BDM3Y1 double-folding potentials for alpha + 40Ca.
 */

#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>

#include <Eigen/Core>

#include "bifold_bifold_xd.hpp"

/**
 * @brief  Evaluate Reid and Paris BDM3Y1 potentials on a uniform radial grid.
 * @math   U(r) = U^C(r) + U^D(r) + U^E(r)
 * @output Labeled potential samples and bifold_out.csv.
 * @note   Reference: Phys. Lett. B 342 (1995) 6-12.
 */
int main() {
    int At_I = 40;
    int Zt_I = 20;
    int Ap_I = 4;
    int Zp_I = 2;
    double hmass_F = 197.329 * 197.329 * static_cast<double>(Ap_I * At_I) / (2.0 * 938.0 * static_cast<double>(Ap_I + At_I));
    double Ze2_F = 197.329 / 137.0;
    double Elab_F = 141.7;
    double Ecm_F = Elab_F * static_cast<double>(At_I) / static_cast<double>(At_I + Ap_I);
    double Y00_F = 1.0 / std::sqrt(4.0 * std::numbers::pi);

    DensityGaussian a_density(1.193);
    DensityFermi Ca40_density(3.60, 0.523, 0.0, 0.0);
    auto rhoa_Func = [&](double r_F) {return a_density.density(r_F);};
    auto rhoCa40_Func = [&](double r_F, double theta_F) {return Ca40_density.density(r_F, theta_F);};
    DensityProjec density_projec(Ap_I, Zp_I, rhoa_Func);
    DensityTarget density_target(At_I, Zt_I, rhoCa40_Func);
    BifoldXD bifold_reid(hmass_F, Ecm_F, Ze2_F, reid_bdm3y1, density_projec, density_target, 0.25);
    BifoldXD bifold_paris = bifold_reid;
    bifold_paris.update_potentials(paris_bdm3y1);

    int Nr_I = 201;
    Eigen::VectorXd r_F1D_r = Eigen::VectorXd::LinSpaced(Nr_I, 1.0e-4, 10.0);
    Eigen::VectorXd UcParis_F1D_r(Nr_I);
    Eigen::VectorXd UdParis_F1D_r(Nr_I);
    Eigen::VectorXd UeParis_F1D_r(Nr_I);
    Eigen::VectorXd UcReid_F1D_r(Nr_I);
    Eigen::VectorXd UdReid_F1D_r(Nr_I);
    Eigen::VectorXd UeReid_F1D_r(Nr_I);
    for (int r_I = 0; r_I < Nr_I; ++r_I) {
        double r_F = r_F1D_r(r_I);
        UcParis_F1D_r(r_I) = bifold_paris.Ucoul_lambda(r_F, 0) * Y00_F;
        UdParis_F1D_r(r_I) = bifold_paris.Ud_lambda(r_F, 0) * Y00_F;
        UeParis_F1D_r(r_I) = bifold_paris.Uex(r_F)(0) * Y00_F;
        UcReid_F1D_r(r_I) = bifold_reid.Ucoul_lambda(r_F, 0) * Y00_F;
        UdReid_F1D_r(r_I) = bifold_reid.Ud_lambda(r_F, 0) * Y00_F;
        UeReid_F1D_r(r_I) = bifold_reid.Uex(r_F)(0) * Y00_F;
    }

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "[Input] Ap = " << Ap_I << ", Zp = " << Zp_I << ", At = " << At_I << ", Zt = " << Zt_I << ", Elab = " << Elab_F << " MeV\n";
    std::cout << "[Computed] BDM3Y1 potentials:\n";
    std::cout << std::setw(10) << "r (fm)" << std::setw(10) << "Uc-p" << std::setw(10) << "Ud-p" << std::setw(10) << "Ue-p" << std::setw(10) << "Uc-r" << std::setw(10) << "Ud-r" << std::setw(10) << "Ue-r" << "\n";
    for (int r_I = 0; r_I < 5; ++r_I) {
        std::cout << std::setw(10) << r_F1D_r(r_I) << std::setw(10) << UcParis_F1D_r(r_I) << std::setw(10) << UdParis_F1D_r(r_I) << std::setw(10) << UeParis_F1D_r(r_I) << std::setw(10) << UcReid_F1D_r(r_I) << std::setw(10) << UdReid_F1D_r(r_I) << std::setw(10) << UeReid_F1D_r(r_I) << "\n";
    }
    std::cout << "       ...       ...       ...       ...       ...       ...       ...\n";
    for (int r_I = Nr_I - 5; r_I < Nr_I; ++r_I) {
        std::cout << std::setw(10) << r_F1D_r(r_I) << std::setw(10) << UcParis_F1D_r(r_I) << std::setw(10) << UdParis_F1D_r(r_I) << std::setw(10) << UeParis_F1D_r(r_I) << std::setw(10) << UcReid_F1D_r(r_I) << std::setw(10) << UdReid_F1D_r(r_I) << std::setw(10) << UeReid_F1D_r(r_I) << "\n";
    }

    std::ofstream output_stream("bifold_out.csv");
    assert(output_stream.is_open());
    output_stream << "r,Uc-p,Ud-p,Ue-p,Uc-r,Ud-r,Ue-r\n";
    for (int r_I = 0; r_I < Nr_I; ++r_I) {
        output_stream << r_F1D_r(r_I) << "," << UcParis_F1D_r(r_I) << "," << UdParis_F1D_r(r_I) << "," << UeParis_F1D_r(r_I) << "," << UcReid_F1D_r(r_I) << "," << UdReid_F1D_r(r_I) << "," << UeReid_F1D_r(r_I) << "\n";
    }
    assert(UcParis_F1D_r.allFinite() && UdParis_F1D_r.allFinite() && UeParis_F1D_r.allFinite());
    assert(UcReid_F1D_r.allFinite() && UdReid_F1D_r.allFinite() && UeReid_F1D_r.allFinite());
    return 0;
}
