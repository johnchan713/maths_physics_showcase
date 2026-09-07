#ifndef NS_CASCADE_LATE_GROWTH_HPP
#define NS_CASCADE_LATE_GROWTH_HPP

#include "ns_cascade/state_optimizer.hpp"

namespace ns_cascade {

struct LateGrowthOptions {
    double start_fraction = 0.5;
    int samples = 5;
    // Units are inverse time, just like d(log H1/2)/dt.
    double temperature = 0.1;
};

inline void validateLateGrowthOptions(const LateGrowthOptions& options) {
    if (!std::isfinite(options.start_fraction) ||
        options.start_fraction <= 0.0 || options.start_fraction >= 1.0 ||
        options.samples < 2 || options.samples > 64 ||
        !std::isfinite(options.temperature) || options.temperature <= 0.0) {
        throw std::invalid_argument("Late-growth window or temperature is invalid");
    }
}

inline std::vector<double> lateGrowthObservationTimes(
    const LateGrowthOptions& options, double horizon) {
    validateLateGrowthOptions(options);
    if (!std::isfinite(horizon) || horizon <= 0.0) {
        throw std::invalid_argument("Late-growth horizon must be positive");
    }
    std::vector<double> times;
    const double tolerance = 64.0 * std::numeric_limits<double>::epsilon() *
                             std::max(1.0, horizon);
    for (int i = 0; i < options.samples; ++i) {
        const double fraction = options.start_fraction +
            (1.0 - options.start_fraction) * static_cast<double>(i) /
                static_cast<double>(options.samples - 1);
        const double time = i + 1 == options.samples ? horizon : horizon * fraction;
        if (time - (times.empty() ? 0.0 : times.back()) <= tolerance) {
            throw std::invalid_argument("Late-growth observation times are not separable");
        }
        times.push_back(time);
    }
    return times;
}

// A=(-Delta)^(1/2) on the retained, mean-zero Fourier subspace.
inline OptimizationState criticalHalfLaplacian(
    const PseudospectralSystem& system, const OptimizationState& state) {
    requireOptimizationState(system, state);
    OptimizationState weighted = system.zeroState();
    const std::vector<WaveVector>& modes = system.gridModes();
    for (std::size_t i = 0; i < state.size(); ++i) {
        if (stateMaximumComponent(modes[i]) <= system.cutoff()) {
            weighted[i] = state[i] * std::sqrt(
                static_cast<double>(modes[i].normSquared()));
        }
    }
    return weighted;
}

// gamma=<Au,F(u)>/<u,Au>. The factor 1/2 from log(sqrt(Q))
// cancels the factor 2 in Q'=2<Au,F>. Viscosity is included in F.
inline double stateCriticalLogRate(
    const PseudospectralSystem& system, const OptimizationState& state) {
    const double q = optimizationStateSpectralSums(system, state).critical_h_half_squared;
    const double rate = stateRealInnerProduct(
        criticalHalfLaplacian(system, state), system.rightHandSide(state)) / q;
    if (!std::isfinite(rate)) {
        throw std::runtime_error("Critical logarithmic rate became non-finite");
    }
    return rate;
}

inline OptimizationState stateCriticalLogRateGradient(
    const PseudospectralSystem& system, const OptimizationState& state) {
    const double q = optimizationStateSpectralSums(system, state).critical_h_half_squared;
    const OptimizationState weighted = criticalHalfLaplacian(system, state);
    const OptimizationState rhs = system.rightHandSide(state);
    const double rate = stateRealInnerProduct(weighted, rhs) / q;
    if (!std::isfinite(rate)) {
        throw std::runtime_error("Critical rate gradient became non-finite");
    }
    // Quotient rule: (AF + F'^*Au)/Q - 2 gamma Au/Q.
    OptimizationState gradient = scaleOptimizationState(addOptimizationStates(
        criticalHalfLaplacian(system, rhs),
        system.adjointTangentRightHandSide(state, weighted), 1.0), 1.0 / q);
    return addOptimizationStates(gradient, weighted, -2.0 * rate / q);
}

struct LateGrowthSummary {
    std::vector<double> rates;
    std::vector<double> gradient_weights;
    double minimum = 0.0;
    double maximum = 0.0;
    double soft_minimum = 0.0;
};

inline LateGrowthSummary summarizeLateGrowthRates(
    const std::vector<double>& rates, double temperature) {
    if (rates.empty() || !std::isfinite(temperature) || temperature <= 0.0) {
        throw std::invalid_argument("Late-growth aggregation data are invalid");
    }
    for (double rate : rates) {
        if (!std::isfinite(rate)) {
            throw std::invalid_argument("Late-growth rate is non-finite");
        }
    }
    LateGrowthSummary result;
    result.rates = rates;
    result.minimum = *std::min_element(rates.begin(), rates.end());
    result.maximum = *std::max_element(rates.begin(), rates.end());
    double normalizer = 0.0;
    for (double rate : rates) {
        // Subtracting the minimum prevents exponential overflow. Negligible
        // weights may underflow to zero; at least one weight is exactly one.
        const double weight = std::exp(-(rate - result.minimum) / temperature);
        result.gradient_weights.push_back(weight);
        normalizer += weight;
    }
    result.soft_minimum = result.minimum - temperature *
        std::log(normalizer / static_cast<double>(rates.size()));
    for (double& weight : result.gradient_weights) weight /= normalizer;
    if (!std::isfinite(result.soft_minimum)) {
        throw std::runtime_error("Late-growth soft minimum became non-finite");
    }
    // This normalized soft minimum is >= the actual sampled minimum.
    // A positive score DOES NOT establish that every sampled rate is positive.
    return result;
}

inline LateGrowthSummary evaluateLateGrowth(
    const PseudospectralSystem& system,
    const std::vector<OptimizationState>& states,
    const LateGrowthOptions& options) {
    validateLateGrowthOptions(options);
    if (states.size() != static_cast<std::size_t>(options.samples)) {
        throw std::invalid_argument("Incomplete late-growth snapshots");
    }
    std::vector<double> rates;
    for (const OptimizationState& state : states) {
        rates.push_back(stateCriticalLogRate(system, state));
    }
    return summarizeLateGrowthRates(rates, options.temperature);
}

}  // namespace ns_cascade

#endif
