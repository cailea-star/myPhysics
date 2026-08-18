/**
 * @file    bayesian_mcmc.hpp
 * @author  cailea
 * @date    2026-05-27
 * @brief   Random-walk Metropolis sampler with reflecting parameter bounds.
 */

#pragma once

#include <Eigen/Core>
#include <cassert>
#include <chrono>
#include <cmath>
#include <functional>
#include <random>
#include <vector>

using Vec2RealFunc = std::function<double(const Eigen::VectorXd&)>;

// ==================== Random-walk Metropolis sampler ====================

class MCMCSampler {
public:
    std::vector<Eigen::VectorXd> samples_F2D_h_x;
    Eigen::VectorXd avg_F1D_x;
    Eigen::MatrixXd cov_F2D_x_x;

private:
    Vec2RealFunc logp_Func;
    const Eigen::VectorXd min_F1D_x;
    const Eigen::VectorXd max_F1D_x;
    Eigen::VectorXd step_F1D_x;
    Eigen::VectorXd current_F1D_x;
    double current_logp_F;
    int accepted_I = 0;
    std::mt19937 random_Generator;
    std::normal_distribution<double> normal_Distribution{0.0, 1.0};
    std::uniform_real_distribution<double> uniform_Distribution{0.0, 1.0};

public:
    /**
     * @brief  Initialize a bounded random-walk Metropolis sampler.
     * @math   x₀ = (x_min + x_max) / 2, s = x_max - x_min
     * @output Initializes the chain state, proposal scale, and sample statistics.
     */
    MCMCSampler(const Vec2RealFunc& logp_Func, const Eigen::VectorXd& min_F1D_x, const Eigen::VectorXd& max_F1D_x)
    : avg_F1D_x(Eigen::VectorXd::Zero(min_F1D_x.size())), cov_F2D_x_x(Eigen::MatrixXd::Zero(min_F1D_x.size(), min_F1D_x.size())), logp_Func(logp_Func), min_F1D_x(min_F1D_x), max_F1D_x(max_F1D_x), step_F1D_x(max_F1D_x - min_F1D_x), current_F1D_x(0.5 * (min_F1D_x + max_F1D_x)), current_logp_F(logp_Func(current_F1D_x)), random_Generator(static_cast<unsigned>(std::chrono::system_clock::now().time_since_epoch().count())) {
        assert(min_F1D_x.size() > 0);
        assert(min_F1D_x.size() == max_F1D_x.size());
        assert(min_F1D_x.allFinite() && max_F1D_x.allFinite());
    }

    /**
     * @brief  Append one sample and update its online statistics.
     * @math   avg_{n+1} = avg_n + δ_n / (n + 1), C_{n+1} = (n - 1) C_n / n + δ_n δ_n^T / (n + 1)
     * @output Appends sample_F1D_x and updates avg_F1D_x and cov_F2D_x_x.
     */
    void add_sample(const Eigen::VectorXd& sample_F1D_x) {
        assert(sample_F1D_x.size() == avg_F1D_x.size());
        const Eigen::VectorXd delta_F1D_x = sample_F1D_x - avg_F1D_x;
        const int n_I = static_cast<int>(samples_F2D_h_x.size());
        avg_F1D_x += delta_F1D_x / (n_I + 1);
        if (n_I > 0) {cov_F2D_x_x = static_cast<double>(n_I - 1) / n_I * cov_F2D_x_x + (delta_F1D_x * delta_F1D_x.transpose()) / (n_I + 1);}
        samples_F2D_h_x.push_back(sample_F1D_x);
    }

    /**
     * @brief  Clear recorded samples and their statistics.
     * @math   N = 0, avg = 0, C = 0
     * @output Clears samples_F2D_h_x and resets avg_F1D_x, cov_F2D_x_x, and accepted_I.
     */
    void reset_samples() {
        samples_F2D_h_x.clear();
        avg_F1D_x.setZero();
        cov_F2D_x_x.setZero();
        accepted_I = 0;
    }

    /**
     * @brief  Advance the random-walk Metropolis chain by one step.
     * @math   x' = R(x + s ξ), log(α) = logp(x') - logp(x)
     * @output Updates the current state and optionally records it.
     */
    void step(bool recordSample_B) {
        Eigen::VectorXd proposal_F1D_x = current_F1D_x;
        for (int x_I = 0; x_I < proposal_F1D_x.size(); ++x_I) {proposal_F1D_x(x_I) += step_F1D_x(x_I) * normal_Distribution(random_Generator);}
        use_reflect_bound(proposal_F1D_x);
        const double proposalLogp_F = logp_Func(proposal_F1D_x);
        const double logalpha_F = proposalLogp_F - current_logp_F;
        if (std::log(uniform_Distribution(random_Generator)) < logalpha_F) {
            current_F1D_x = proposal_F1D_x;
            current_logp_F = proposalLogp_F;
            ++accepted_I;
        }
        if (recordSample_B) {add_sample(current_F1D_x);}
    }

    /**
     * @brief  Run multiple Metropolis steps with periodic sample recording.
     * @math   record x_h when h mod N_interval = 0
     * @output Advances the chain by Nstep_I steps and records scheduled states.
     */
    void run(int Nstep_I, int sample_interval_I) {
        assert(Nstep_I >= 0);
        assert(sample_interval_I > 0);
        for (int h_I = 1; h_I < Nstep_I + 1; ++h_I) {step(h_I % sample_interval_I == 0);}
    }

private:
    /**
     * @brief  Reflect a parameter vector into its coordinate bounds.
     * @math   y = (x - x_min) mod 2L, R(x) = x_min + y for y ≤ L and x_max - (y - L) otherwise
     * @output Updates x_F1D_x in place within [min_F1D_x, max_F1D_x].
     */
    void use_reflect_bound(Eigen::VectorXd& x_F1D_x) const {
        assert(x_F1D_x.size() == min_F1D_x.size());
        for (int x_I = 0; x_I < x_F1D_x.size(); ++x_I) {
            const double L_F = max_F1D_x(x_I) - min_F1D_x(x_I);
            if (L_F <= 0.0) {
                x_F1D_x(x_I) = min_F1D_x(x_I);
                continue;
            }
            double y_F = std::fmod(x_F1D_x(x_I) - min_F1D_x(x_I), 2.0 * L_F);
            if (y_F < 0.0) {y_F += 2.0 * L_F;}
            x_F1D_x(x_I) = y_F <= L_F ? min_F1D_x(x_I) + y_F : max_F1D_x(x_I) - (y_F - L_F);
        }
    }
};
