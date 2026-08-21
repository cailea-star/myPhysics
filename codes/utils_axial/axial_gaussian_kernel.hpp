/**
 * @file    axial_gaussian_kernel.hpp
 * @author  cailea
 * @date    2026-05-27
 * @brief   Shared axial-HO Gaussian matrix-element cache and kernels.
 */

#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>
#include <gsl/gsl_sf_gamma.h>

#include "axial_config.hpp"
#include "tools_hashtable.hpp"

template <int Ng_I>
class AxialGaussianKernel {
public:
    using GKey = std::array<int, 4>;
    using GValues = std::array<double, Ng_I>;
    bool isBuilt_B = false;
    static_assert(Ng_I > 0, "AxialGaussianKernel requires at least one Gaussian.");

private:
    struct Metadata {
        std::uint32_t magic_I;
        std::uint32_t version_I;
        std::uint32_t Ngcache_I;
        double br_F;
        double bz_F;
        int nzMax_I;
        int nrMax_I;
        int LambdaMax_I;
        GValues mu_F1D_g;

        /**
         * @brief  Compare cache metadata fieldwise.
         * @math   M_1 = M_2.
         * @output Exact metadata equality.
         */
        bool operator==(const Metadata&) const = default;
    };
    static_assert(std::is_trivially_copyable_v<Metadata>, "Metadata must be trivially copyable.");

protected:
    std::vector<AxialSPLabel> labels_S1D_sp;
    Metadata metadata;
    MyHashTable<GValues> Gz_Table;
    MyHashTable<GValues> Gr_Table;

public:
    /**
     * @brief  Construct an axial Gaussian kernel.
     * @math   G = G^zG^r.
    * @output Empty configured tables.
     */
    AxialGaussianKernel(const AxialConfig& config_, const GValues& mu_F1D_g_)
    : isBuilt_B(false), labels_S1D_sp(config_.labels_S1D_sp), metadata{0x41474b31u, 5, static_cast<std::uint32_t>(Ng_I), config_.br_F, config_.bz_F, 0, 0, 0, mu_F1D_g_}, Gz_Table({0, 0, 0, 0}, {0, 0, 0, 0}), Gr_Table({0, 0, 0, 0}, {0, 0, 0, 0}) {
        assert((std::all_of(metadata.mu_F1D_g.begin(), metadata.mu_F1D_g.end(), [](double mu_F) {return std::isfinite(mu_F) && mu_F > 0.0;})));

        // {α_sp} -> (n_z^max,n_r^max,Λ^max,rorder^max).
        for (const AxialSPLabel& label_ : labels_S1D_sp) {
            metadata.nzMax_I = std::max(metadata.nzMax_I, label_.nz_I);
            metadata.nrMax_I = std::max(metadata.nrMax_I, label_.nr_I);
            metadata.LambdaMax_I = std::max(metadata.LambdaMax_I, std::abs(label_.Lambda_I));
        }
        const int rorderMax_I = pack_rkey(metadata.nrMax_I, metadata.LambdaMax_I);

        // Bounds -> (Gz_Table,Gr_Table).
        Gz_Table = MyHashTable<GValues>({0, 0, 0, 0}, {metadata.nzMax_I, metadata.nzMax_I, metadata.nzMax_I, metadata.nzMax_I});
        Gr_Table = MyHashTable<GValues>({0, 0, 0, 0}, {rorderMax_I, rorderMax_I, rorderMax_I, rorderMax_I});
    }

    /**
     * @brief  Build axial and transverse Gaussian tables.
     * @math   {μ_g} -> {G_g^z,G_g^r}.
     * @output Filled tables.
     */
    void build_tables();

    /**
     * @brief  Deserialize and validate Gaussian tables.
     * @math   file ⊕ config ⊕ {μ_g} -> G^z ⊕ G^r.
     * @output Validated built kernel.
     */
    static AxialGaussianKernel from_cache(const std::string& filepath_Str, const AxialConfig& config_, const GValues& mu_F1D_g_);

    /**
     * @brief  Serialize the Gaussian tables.
     * @math   metadata ⊕ G^z ⊕ G^r -> file.
     * @output Binary cache.
     */
    void to_cache(const std::string& filepath_Str);

    /**
     * @brief  Read axial Gaussian elements.
     * @math   (n_z1,n_z2,n_z3,n_z4) -> {G_g^z}.
     * @output Cached values.
     */
    const GValues& read_Gz(int nz1_I, int nz2_I, int nz3_I, int nz4_I) const;

    /**
     * @brief  Read transverse Gaussian elements.
     * @math   {(n_ri,Λ_i)} -> {G_g^r}.
     * @output Cached values.
     */
    const GValues& read_Gr(int nr1_I, int Lambda1_I, int nr2_I, int Lambda2_I, int nr3_I, int Lambda3_I, int nr4_I, int Lambda4_I) const;

private:
    /**
     * @brief  Pack (n_r,Λ) into one index.
     * @math   k = (n_r-n_r^min)(Λ^max-Λ^min+1)+Λ-Λ^min.
     * @output Packed index.
     */
    int pack_rkey(int nr_I, int Lambda_I) const;

    /**
     * @brief  Canonicalize a four-index key.
     * @math   k = min{1234,2143,3412,4321}.
     * @output Canonical key.
     */
    static GKey canonicalize_key(int key1_I, int key2_I, int key3_I, int key4_I);

    /**
     * @brief  Evaluate G^1D by the HFBTHO hypergeometric expansion.
     * @math   G^1D = μ/(sqrt(2π³)b)Σ_n T Fbar.
     * @output One-dimensional matrix element.
     */
    static double calc_G1D(int n1_I, int n2_I, int n3_I, int n4_I, double mu_F, double b_F);

    /**
     * @brief  Evaluate a polar-to-Cartesian HO coefficient.
     * @math   C = AΣ_q binom(n_x,k_x)binom(n_y,q)(-1)^(n_y-q).
     * @output Expansion coefficient.
     */
    static double calc_polar_cartesian_coeff(int nr_I, int Lambda_I, int ny_I);

    /**
     * @brief  Evaluate the transverse Gaussian factor.
     * @math   G^r = Σ P_yΠ_a C_aG_x^1DG_y^1D.
     * @output Transverse matrix element.
     */
    static double calc_Gr(int nr1_I, int Lambda1_I, int nr2_I, int Lambda2_I, int nr3_I, int Lambda3_I, int nr4_I, int Lambda4_I, double mu_F, double br_F);

};

template <int Ng_I>
void AxialGaussianKernel<Ng_I>::build_tables() {
    if (isBuilt_B) {return;}
    std::cout << "[AxialGaussianKernel]: Building Gaussian tables for " << Ng_I << " Gaussians..." << std::endl;

    // {α_sp} -> {n_z}.
    std::set<int> nz_Set;
    for (const AxialSPLabel& label_ : labels_S1D_sp) {
        nz_Set.insert(label_.nz_I);
    }
    // {n_z}^4 -> unique even-parity keys.
    std::set<GKey> Gzkeys_Set;
    std::vector<std::pair<GKey, GValues>> Gzjobs_S1D_job;
    for (int nz1_I : nz_Set) {
        for (int nz2_I : nz_Set) {
            for (int nz3_I : nz_Set) {
                for (int nz4_I : nz_Set) {
                    if ((nz1_I + nz2_I + nz3_I + nz4_I) % 2 != 0) {continue;}
                    const GKey key_I1D_q = canonicalize_key(nz1_I, nz2_I, nz3_I, nz4_I);
                    if (Gzkeys_Set.insert(key_I1D_q).second) {Gzjobs_S1D_job.push_back({key_I1D_q, GValues{}});}
                }
            }
        }
    }

    // (key,μ_g) -> G_g^z.
    #pragma omp parallel for schedule(dynamic)
    for (int job_I = 0; job_I < static_cast<int>(Gzjobs_S1D_job.size()); ++job_I) {
        const GKey& key_I1D_q = Gzjobs_S1D_job[job_I].first;
        GValues& G_F1D_g = Gzjobs_S1D_job[job_I].second;
        for (int g_I = 0; g_I < Ng_I; ++g_I) {
            G_F1D_g[g_I] = calc_G1D(key_I1D_q[0], key_I1D_q[1], key_I1D_q[2], key_I1D_q[3], metadata.mu_F1D_g[g_I], metadata.bz_F);
        }
    }

    // {(key,G_g^z)} -> Gz_Table.
    Gz_Table.reserve(Gzjobs_S1D_job.size());
    for (const auto& job_ : Gzjobs_S1D_job) {
        Gz_Table.write(job_.first, job_.second);
    }

    // {α_sp} -> {(n_r,±Λ)}.
    std::set<std::array<int, 2>> nrLambda_Set;
    for (const AxialSPLabel& label_ : labels_S1D_sp) {
        nrLambda_Set.insert({label_.nr_I, label_.Lambda_I});
        nrLambda_Set.insert({label_.nr_I, -label_.Lambda_I});
    }
    // {(n_r,Λ)}^4 -> unique conserving keys.
    std::set<GKey> Grkeys_Set;
    std::vector<std::pair<GKey, GValues>> Grjobs_S1D_job;
    for (const auto& nrLambda1_I1D_field : nrLambda_Set) {
        for (const auto& nrLambda2_I1D_field : nrLambda_Set) {
            for (const auto& nrLambda3_I1D_field : nrLambda_Set) {
                for (const auto& nrLambda4_I1D_field : nrLambda_Set) {
                    if (nrLambda1_I1D_field[1] + nrLambda2_I1D_field[1] != nrLambda3_I1D_field[1] + nrLambda4_I1D_field[1]) {continue;}
                    const int rkey1_I = pack_rkey(nrLambda1_I1D_field[0], nrLambda1_I1D_field[1]);
                    const int rkey2_I = pack_rkey(nrLambda2_I1D_field[0], nrLambda2_I1D_field[1]);
                    const int rkey3_I = pack_rkey(nrLambda3_I1D_field[0], nrLambda3_I1D_field[1]);
                    const int rkey4_I = pack_rkey(nrLambda4_I1D_field[0], nrLambda4_I1D_field[1]);
                    const GKey key_I1D_q = canonicalize_key(rkey1_I, rkey2_I, rkey3_I, rkey4_I);
                    if (Grkeys_Set.insert(key_I1D_q).second) {Grjobs_S1D_job.push_back({key_I1D_q, GValues{}});}
                }
            }
        }
    }

    // (key,μ_g) -> G_g^r.
    const int NLambda_I = 2 * metadata.LambdaMax_I + 1;
    #pragma omp parallel for schedule(dynamic)
    for (int job_I = 0; job_I < static_cast<int>(Grjobs_S1D_job.size()); ++job_I) {
        const GKey& key_I1D_q = Grjobs_S1D_job[job_I].first;
        const int nr1_I = key_I1D_q[0] / NLambda_I;
        const int nr2_I = key_I1D_q[1] / NLambda_I;
        const int nr3_I = key_I1D_q[2] / NLambda_I;
        const int nr4_I = key_I1D_q[3] / NLambda_I;
        const int Lambda1_I = key_I1D_q[0] % NLambda_I - metadata.LambdaMax_I;
        const int Lambda2_I = key_I1D_q[1] % NLambda_I - metadata.LambdaMax_I;
        const int Lambda3_I = key_I1D_q[2] % NLambda_I - metadata.LambdaMax_I;
        const int Lambda4_I = key_I1D_q[3] % NLambda_I - metadata.LambdaMax_I;
        GValues& G_F1D_g = Grjobs_S1D_job[job_I].second;
        for (int g_I = 0; g_I < Ng_I; ++g_I) {
            G_F1D_g[g_I] = calc_Gr(nr1_I, Lambda1_I, nr2_I, Lambda2_I, nr3_I, Lambda3_I, nr4_I, Lambda4_I, metadata.mu_F1D_g[g_I], metadata.br_F);
        }
    }

    // {(key,G_g^r)} -> Gr_Table.
    Gr_Table.reserve(Grjobs_S1D_job.size());
    for (const auto& job_ : Grjobs_S1D_job) {
        Gr_Table.write(job_.first, job_.second);
    }

    isBuilt_B = true;
}

template <int Ng_I>
void AxialGaussianKernel<Ng_I>::to_cache(const std::string& filepath_Str) {
    // G -> complete tables.
    build_tables();

    // metadata ⊕ (Gz_Table,Gr_Table) -> stream.
    std::ofstream output_(filepath_Str, std::ios::binary | std::ios::trunc);
    if (!output_) {throw std::runtime_error("[ERROR]: [AxialGaussianKernel::to_cache] cannot open " + filepath_Str);}
    output_.write(reinterpret_cast<const char*>(&metadata), sizeof(metadata));
    if (!output_) {throw std::runtime_error("[ERROR]: [AxialGaussianKernel::to_cache] metadata write failed");}
    Gz_Table.to_stream(output_);
    Gr_Table.to_stream(output_);
}

template <int Ng_I>
AxialGaussianKernel<Ng_I> AxialGaussianKernel<Ng_I>::from_cache(const std::string& filepath_Str, const AxialConfig& config_, const GValues& mu_F1D_g_) {
    // stream -> cache metadata.
    std::ifstream input_(filepath_Str, std::ios::binary);
    if (!input_) {throw std::runtime_error("[ERROR]: [AxialGaussianKernel::from_cache] cannot open " + filepath_Str);}
    Metadata cachedMetadata_{};
    input_.read(reinterpret_cast<char*>(&cachedMetadata_), sizeof(cachedMetadata_));
    if (!input_) {throw std::runtime_error("[ERROR]: [AxialGaussianKernel::from_cache] metadata read failed");}

    // cache metadata = requested configuration.
    AxialGaussianKernel kernel_(config_, mu_F1D_g_);
    if (cachedMetadata_ != kernel_.metadata) {throw std::runtime_error("[ERROR]: [AxialGaussianKernel::from_cache] metadata mismatch: " + filepath_Str);}

    // stream -> (Gz_Table,Gr_Table).
    kernel_.Gz_Table.from_stream(input_);
    kernel_.Gr_Table.from_stream(input_);
    kernel_.isBuilt_B = true;
    return kernel_;
}

template <int Ng_I>
const typename AxialGaussianKernel<Ng_I>::GValues& AxialGaussianKernel<Ng_I>::read_Gz(int nz1_I, int nz2_I, int nz3_I, int nz4_I) const {
    assert(isBuilt_B);
    const GKey key_I1D_q = canonicalize_key(nz1_I, nz2_I, nz3_I, nz4_I);
    return Gz_Table.read(key_I1D_q);
}

template <int Ng_I>
const typename AxialGaussianKernel<Ng_I>::GValues& AxialGaussianKernel<Ng_I>::read_Gr(int nr1_I, int Lambda1_I, int nr2_I, int Lambda2_I, int nr3_I, int Lambda3_I, int nr4_I, int Lambda4_I) const {
    assert(isBuilt_B);
    const int rkey1_I = pack_rkey(nr1_I, Lambda1_I);
    const int rkey2_I = pack_rkey(nr2_I, Lambda2_I);
    const int rkey3_I = pack_rkey(nr3_I, Lambda3_I);
    const int rkey4_I = pack_rkey(nr4_I, Lambda4_I);
    const GKey key_I1D_q = canonicalize_key(rkey1_I, rkey2_I, rkey3_I, rkey4_I);
    return Gr_Table.read(key_I1D_q);
}

template <int Ng_I>
int AxialGaussianKernel<Ng_I>::pack_rkey(int nr_I, int Lambda_I) const {return nr_I * (2 * metadata.LambdaMax_I + 1) + Lambda_I + metadata.LambdaMax_I;}

template <int Ng_I>
typename AxialGaussianKernel<Ng_I>::GKey AxialGaussianKernel<Ng_I>::canonicalize_key(int key1_I, int key2_I, int key3_I, int key4_I) {
    const GKey key1234_I1D_q{key1_I, key2_I, key3_I, key4_I};
    const GKey key2143_I1D_q{key2_I, key1_I, key4_I, key3_I};
    const GKey key3412_I1D_q{key3_I, key4_I, key1_I, key2_I};
    const GKey key4321_I1D_q{key4_I, key3_I, key2_I, key1_I};
    return std::min(std::min(key1234_I1D_q, key2143_I1D_q), std::min(key3412_I1D_q, key4321_I1D_q));
}

template <int Ng_I>
double AxialGaussianKernel<Ng_I>::calc_G1D(int n1_I, int n2_I, int n3_I, int n4_I, double mu_F, double b_F) {
    const bool useOriginalOrder_B = std::min(n2_I, n4_I) <= std::min(n1_I, n3_I);
    const int na_I = useOriginalOrder_B ? n1_I : n2_I;
    const int nb_I = useOriginalOrder_B ? n2_I : n1_I;
    const int nc_I = useOriginalOrder_B ? n3_I : n4_I;
    const int nd_I = useOriginalOrder_B ? n4_I : n3_I;
    const double chi_F = 1.0 + mu_F * mu_F / (2.0 * b_F * b_F);
    double G_F = 0.0;
    for (int n_I = std::abs(nb_I - nd_I); n_I <= nb_I + nd_I; n_I += 2) {
        if ((na_I + nc_I + n_I) % 2 != 0) {break;}

        // ln T(n_b,n_d;n).
        const double logT_F = 0.5 * (gsl_sf_lnfact(nb_I) + gsl_sf_lnfact(nd_I) + gsl_sf_lnfact(n_I)) - gsl_sf_lnfact((-nb_I + nd_I + n_I) / 2) - gsl_sf_lnfact((nb_I - nd_I + n_I) / 2) - gsl_sf_lnfact((nb_I + nd_I - n_I) / 2);
        const double xi_F = 0.5 * (na_I + nc_I + n_I + 1.0);
        const double gamma_F = gsl_sf_gamma(xi_F - na_I) * gsl_sf_gamma(xi_F - nc_I) * gsl_sf_gamma(xi_F - n_I);
        const double denominator_F = std::pow(chi_F, xi_F) * std::exp(0.5 * (gsl_sf_lnfact(na_I) + gsl_sf_lnfact(nc_I) + gsl_sf_lnfact(n_I)));

        // _2F_1(-n_a,-n_c;n+1-ξ;1-χ).
        const int pmin_I = std::min(na_I, nc_I);
        const int pmax_I = std::max(na_I, nc_I);
        double term_F = 1.0;
        double sum_F = 1.0;
        for (int i_I = 0; i_I < pmin_I; ++i_I) {
            term_F *= (pmin_I - i_I) * (pmax_I - i_I) * (1.0 - chi_F) / ((i_I + 1) * (n_I + 1.0 - xi_F + i_I));
            sum_F += term_F;
        }
        const double Fbar_F = gamma_F / denominator_F * sum_F;
        G_F += std::exp(logT_F) * Fbar_F;
    }
    const double pi_F = std::acos(-1.0);
    return mu_F * G_F / (std::sqrt(2.0 * pi_F * pi_F * pi_F) * b_F);
}

template <int Ng_I>
double AxialGaussianKernel<Ng_I>::calc_polar_cartesian_coeff(int nr_I, int Lambda_I, int ny_I) {
    if (ny_I < 0 || ny_I > 2 * nr_I + std::abs(Lambda_I)) {return 0.0;}
    const int nx_I = 2 * nr_I + std::abs(Lambda_I) - ny_I;
    const int qMax_I = std::min(ny_I, nr_I + (std::abs(Lambda_I) - Lambda_I) / 2);
    const double prefactor_F = std::pow(-1, nr_I) * std::pow(2.0, -nr_I - 0.5 * std::abs(Lambda_I)) * std::exp(0.5 * (gsl_sf_lnfact(nr_I + std::abs(Lambda_I)) + gsl_sf_lnfact(nr_I) - gsl_sf_lnfact(nx_I) - gsl_sf_lnfact(ny_I)));
    double sum_F = 0.0;
    for (int q_I = 0; q_I <= qMax_I; ++q_I) {
        const int kx_I = nr_I - q_I + (std::abs(Lambda_I) - Lambda_I) / 2;
        if (kx_I >= 0 && kx_I <= nx_I) {
            sum_F += gsl_sf_choose(static_cast<unsigned int>(nx_I), static_cast<unsigned int>(kx_I)) * gsl_sf_choose(static_cast<unsigned int>(ny_I), static_cast<unsigned int>(q_I)) * std::pow(-1, ny_I - q_I);
        }
    }
    return prefactor_F * sum_F;
}

template <int Ng_I>
double AxialGaussianKernel<Ng_I>::calc_Gr(int nr1_I, int Lambda1_I, int nr2_I, int Lambda2_I, int nr3_I, int Lambda3_I, int nr4_I, int Lambda4_I, double mu_F, double br_F) {
    if (Lambda1_I + Lambda2_I != Lambda3_I + Lambda4_I) {return 0.0;}
    double Gr_F = 0.0;
    for (int ny1_I = 0; ny1_I <= 2 * nr1_I + std::abs(Lambda1_I); ++ny1_I) {
        const int nx1_I = 2 * nr1_I + std::abs(Lambda1_I) - ny1_I;
        const double C1_F = calc_polar_cartesian_coeff(nr1_I, Lambda1_I, ny1_I);
        for (int ny2_I = 0; ny2_I <= 2 * nr2_I + std::abs(Lambda2_I); ++ny2_I) {
            const int nx2_I = 2 * nr2_I + std::abs(Lambda2_I) - ny2_I;
            const double C2_F = calc_polar_cartesian_coeff(nr2_I, Lambda2_I, ny2_I);
            for (int ny3_I = 0; ny3_I <= 2 * nr3_I + std::abs(Lambda3_I); ++ny3_I) {
                const int nx3_I = 2 * nr3_I + std::abs(Lambda3_I) - ny3_I;
                const double C3_F = calc_polar_cartesian_coeff(nr3_I, Lambda3_I, ny3_I);
                for (int ny4_I = (ny1_I + ny2_I + ny3_I) % 2; ny4_I <= 2 * nr4_I + std::abs(Lambda4_I); ny4_I += 2) {
                    const int nx4_I = 2 * nr4_I + std::abs(Lambda4_I) - ny4_I;
                    const double C4_F = calc_polar_cartesian_coeff(nr4_I, Lambda4_I, ny4_I);
                    const int phase_I = ny1_I + ny2_I + (ny1_I + ny2_I + ny3_I + ny4_I) / 2;
                    Gr_F += std::pow(-1, phase_I) * C1_F * C2_F * C3_F * C4_F * calc_G1D(nx1_I, nx2_I, nx3_I, nx4_I, mu_F, br_F) * calc_G1D(ny1_I, ny2_I, ny3_I, ny4_I, mu_F, br_F);
                }
            }
        }
    }
    return Gr_F;
}
