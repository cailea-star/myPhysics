/**
 * @file    bayesian_mcmc.hpp
 * @author  cailea
 * @date    2026-05-27
 * @brief   Random-walk Metropolis sampler with reflecting parameter bounds.
 */

#pragma once

#include <Eigen/Core>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <functional>
#include <random>
#include <vector>

using Vec2RealFunc = std::function<double(const Eigen::VectorXd&)>;

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
    Eigen::VectorXd proposal_F1D_x;
    Eigen::VectorXd delta_F1D_x;
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
     * @note   Requires x_max > x_min coordinatewise.
     */
    MCMCSampler(const Vec2RealFunc& logp_Func_, const Eigen::VectorXd& min_F1D_x_, const Eigen::VectorXd& max_F1D_x_)
    : avg_F1D_x(Eigen::VectorXd::Zero(min_F1D_x_.size())), cov_F2D_x_x(Eigen::MatrixXd::Zero(min_F1D_x_.size(), min_F1D_x_.size())), logp_Func(logp_Func_), min_F1D_x(min_F1D_x_), max_F1D_x(max_F1D_x_), step_F1D_x(max_F1D_x_ - min_F1D_x_), current_F1D_x(0.5 * (min_F1D_x_ + max_F1D_x_)), proposal_F1D_x(min_F1D_x_.size()), delta_F1D_x(min_F1D_x_.size()), current_logp_F(logp_Func_(current_F1D_x)), random_Generator(static_cast<unsigned>(std::chrono::system_clock::now().time_since_epoch().count())) {
        // x_min < x_max; dim(x_min)=dim(x_max)>0.
        assert(min_F1D_x_.size() > 0);
        assert(min_F1D_x_.size() == max_F1D_x_.size());
        assert(min_F1D_x_.allFinite() && max_F1D_x_.allFinite());
        assert((max_F1D_x_.array() > min_F1D_x_.array()).all());
    }

    /**
     * @brief  Append one sample using Welford online covariance.
     * @math   avg_{n+1} = avg_n + δ_n / (n + 1), C_{n+1} = (n - 1) C_n / n + δ_n δ_n^T / (n + 1)
     * @output Appends sample_F1D_x and updates avg_F1D_x and cov_F2D_x_x.
     */
    void add_sample(const Eigen::VectorXd& sample_F1D_x) {
        // (x^{(n+1)},x̄_n,C_n) → (x̄_{n+1},C_{n+1}).
        assert(sample_F1D_x.size() == avg_F1D_x.size());
        delta_F1D_x = sample_F1D_x - avg_F1D_x;
        const int n_I = static_cast<int>(samples_F2D_h_x.size());
        avg_F1D_x += delta_F1D_x / (n_I + 1);
        const double hasPrevious_F = static_cast<double>(n_I > 0);
        const double nSafe_F = std::max(1, n_I);
        cov_F2D_x_x *= hasPrevious_F * static_cast<double>(n_I - 1) / nSafe_F;
        cov_F2D_x_x.noalias() += hasPrevious_F * delta_F1D_x * delta_F1D_x.transpose() / (n_I + 1);
        samples_F2D_h_x.push_back(sample_F1D_x);
    }

    /**
     * @brief  Clear recorded samples and their statistics.
     * @math   N = 0, avg = 0, C = 0
     * @output Clears samples_F2D_h_x and resets avg_F1D_x, cov_F2D_x_x, and accepted_I.
     */
    void reset_samples() {
        // ({x_h},x̄,C,N_accept) → (∅,0,0,0).
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
        // x_h → R(x_h+sξ) → x'.
        proposal_F1D_x = current_F1D_x;
        for (int x_I = 0; x_I < proposal_F1D_x.size(); ++x_I) {proposal_F1D_x(x_I) += step_F1D_x(x_I) * normal_Distribution(random_Generator);}
        use_reflect_bound(proposal_F1D_x);
        const double proposalLogp_F = logp_Func(proposal_F1D_x);
        const double logalpha_F = proposalLogp_F - current_logp_F;
        // log u < log α → x_{h+1}=x'.
        if (std::log(uniform_Distribution(random_Generator)) < logalpha_F) {
            current_F1D_x = proposal_F1D_x;
            current_logp_F = proposalLogp_F;
            ++accepted_I;
        }
        if (recordSample_B) {add_sample(current_F1D_x);}
    }

    /**
     * @brief  Run multiple Metropolis steps with periodic sample recording.
     * @math   {x_h | h mod N_interval = 0}
     * @output Advances Nstep_I steps and records scheduled states.
     */
    void run(int Nstep_I, int sample_interval_I) {
        assert(Nstep_I >= 0);
        assert(sample_interval_I > 0);
        for (int h_I = 1; h_I < Nstep_I + 1; ++h_I) {step(h_I % sample_interval_I == 0);}
    }

private:
    /**
     * @brief  Reflect a parameter vector into its coordinate bounds.
     * @math   y=(x-x_min) mod 2L, R=x_min+{y,y≤L;2L-y,y>L}
     * @output Updates x_F1D_x in place within [min_F1D_x, max_F1D_x].
     */
    void use_reflect_bound(Eigen::VectorXd& x_F1D_x) const {
        // x_i → x_i mod 2L_i → R_i(x_i).
        assert(x_F1D_x.size() == min_F1D_x.size());
        for (int x_I = 0; x_I < x_F1D_x.size(); ++x_I) {
            const double L_F = max_F1D_x(x_I) - min_F1D_x(x_I);
            double y_F = std::fmod(x_F1D_x(x_I) - min_F1D_x(x_I), 2.0 * L_F);
            y_F += 2.0 * L_F * static_cast<double>(y_F < 0.0);
            x_F1D_x(x_I) = min_F1D_x(x_I) + L_F - std::abs(y_F - L_F);
        }
    }
};
