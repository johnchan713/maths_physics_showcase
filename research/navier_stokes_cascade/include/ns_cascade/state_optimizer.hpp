#ifndef NS_CASCADE_STATE_OPTIMIZER_HPP
#define NS_CASCADE_STATE_OPTIMIZER_HPP

#include "ns_cascade/pseudospectral.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>

namespace ns_cascade {

using OptimizationState = PseudospectralSystem::State;

struct StateObjectiveWeights {
    double characteristic_scale_weight = 0.15;
    double cutoff_penalty_weight = 0.04;
    double cutoff_fraction_threshold = 0.01;
};

struct StateObjectiveValue {
    double critical_log_growth = 0.0;
    double characteristic_log_growth = 0.0;
    double cutoff_penalty = 0.0;
    double total = 0.0;
};

struct StateSpectralSums {
    double energy = 0.0;
    double enstrophy = 0.0;
    double critical_h_half_squared = 0.0;
    double cutoff_energy = 0.0;
};

inline int stateMaximumComponent(const WaveVector& wave) {
    return std::max(std::abs(wave.x),
                    std::max(std::abs(wave.y), std::abs(wave.z)));
}

inline void requireOptimizationState(
    const PseudospectralSystem& system,
    const OptimizationState& state) {
    if (state.size() != system.gridPointCount()) {
        throw std::invalid_argument("Optimization state has the wrong grid size");
    }
}

inline double stateRealInnerProduct(const OptimizationState& left,
                                    const OptimizationState& right) {
    if (left.size() != right.size()) {
        throw std::invalid_argument("State inner product sizes do not match");
    }
    double value = 0.0;
    for (std::size_t i = 0; i < left.size(); ++i) {
        value += std::real(innerProduct(left[i], right[i]));
    }
    return value;
}

inline double optimizationStateNorm(const OptimizationState& state) {
    return std::sqrt(std::max(0.0, stateRealInnerProduct(state, state)));
}

inline OptimizationState addOptimizationStates(
    const OptimizationState& state,
    const OptimizationState& increment,
    double scale) {
    if (state.size() != increment.size() || !std::isfinite(scale)) {
        throw std::invalid_argument("Scaled state addition is invalid");
    }
    OptimizationState result = state;
    for (std::size_t i = 0; i < result.size(); ++i) {
        result[i] += increment[i] * scale;
    }
    return result;
}

inline OptimizationState scaleOptimizationState(
    const OptimizationState& state,
    double scale) {
    if (!std::isfinite(scale)) {
        throw std::invalid_argument("State scale must be finite");
    }
    OptimizationState result = state;
    for (std::size_t i = 0; i < result.size(); ++i) {
        result[i] = result[i] * scale;
    }
    return result;
}

inline OptimizationState projectBandLimitedRealSolenoidal(
    const PseudospectralSystem& system,
    const OptimizationState& state,
    int bandwidth) {
    requireOptimizationState(system, state);
    if (bandwidth < 1 || bandwidth > system.cutoff()) {
        throw std::invalid_argument(
            "Optimization bandwidth must lie inside the retained cutoff");
    }
    OptimizationState projected = system.zeroState();
    const std::vector<WaveVector>& modes = system.gridModes();
    for (std::size_t i = 0; i < state.size(); ++i) {
        const WaveVector& wave = modes[i];
        if (wave.normSquared() == 0 ||
            stateMaximumComponent(wave) > bandwidth) {
            continue;
        }
        projected[i] = lerayProject(wave, state[i]);
    }

    // Average each conjugate pair. Iterating over both members is harmless and
    // avoids relying on a storage-order convention for the canonical half.
    OptimizationState result = system.zeroState();
    for (std::size_t i = 0; i < projected.size(); ++i) {
        const WaveVector& wave = modes[i];
        if (wave.normSquared() == 0 ||
            stateMaximumComponent(wave) > bandwidth) {
            continue;
        }
        const std::size_t partner = system.indexOf(-wave);
        const ComplexVector average =
            (projected[i] + conjugate(projected[partner])) * 0.5;
        result[i] = average;
        result[partner] = conjugate(average);
    }
    return result;
}

inline OptimizationState normalizeOptimizationEnergy(
    const PseudospectralSystem& system,
    const OptimizationState& state,
    double target_energy) {
    requireOptimizationState(system, state);
    if (!std::isfinite(target_energy) || target_energy <= 0.0) {
        throw std::invalid_argument("Target optimization energy must be positive");
    }
    const double current_energy = system.energy(state);
    if (!std::isfinite(current_energy) || current_energy <= 0.0) {
        throw std::invalid_argument(
            "Cannot normalize a non-positive optimization state");
    }
    return scaleOptimizationState(
        state, std::sqrt(target_energy / current_energy));
}

inline OptimizationState makeBandLimitedOptimizationState(
    const PseudospectralSystem& system,
    const OptimizationState& state,
    int bandwidth,
    double target_energy) {
    return normalizeOptimizationEnergy(
        system,
        projectBandLimitedRealSolenoidal(system, state, bandwidth),
        target_energy);
}

inline StateSpectralSums optimizationStateSpectralSums(
    const PseudospectralSystem& system,
    const OptimizationState& state) {
    requireOptimizationState(system, state);
    StateSpectralSums sums;
    const std::vector<WaveVector>& modes = system.gridModes();
    for (std::size_t i = 0; i < state.size(); ++i) {
        const int wave_squared = modes[i].normSquared();
        if (wave_squared == 0 ||
            stateMaximumComponent(modes[i]) > system.cutoff()) {
            continue;
        }
        const double mode_energy = 0.5 * normSquared(state[i]);
        const double wave_magnitude =
            std::sqrt(static_cast<double>(wave_squared));
        sums.energy += mode_energy;
        sums.enstrophy += wave_squared * mode_energy;
        sums.critical_h_half_squared +=
            wave_magnitude * normSquared(state[i]);
        if (stateMaximumComponent(modes[i]) == system.cutoff()) {
            sums.cutoff_energy += mode_energy;
        }
    }
    if (!std::isfinite(sums.energy) || !std::isfinite(sums.enstrophy) ||
        !std::isfinite(sums.critical_h_half_squared) ||
        !std::isfinite(sums.cutoff_energy) || sums.energy <= 0.0 ||
        sums.enstrophy <= 0.0 || sums.critical_h_half_squared <= 0.0 ||
        sums.cutoff_energy < 0.0) {
        throw std::invalid_argument(
            "State objective requires finite positive spectral sums");
    }
    return sums;
}

inline void validateStateObjectiveWeights(
    const StateObjectiveWeights& weights) {
    if (!std::isfinite(weights.characteristic_scale_weight) ||
        weights.characteristic_scale_weight < 0.0 ||
        !std::isfinite(weights.cutoff_penalty_weight) ||
        weights.cutoff_penalty_weight < 0.0 ||
        !std::isfinite(weights.cutoff_fraction_threshold) ||
        weights.cutoff_fraction_threshold <= 0.0 ||
        weights.cutoff_fraction_threshold >= 1.0) {
        throw std::invalid_argument("State-objective weights are invalid");
    }
}

inline StateObjectiveValue evaluateStateObjective(
    const PseudospectralSystem& system,
    const OptimizationState& initial,
    const OptimizationState& final,
    const StateObjectiveWeights& weights) {
    validateStateObjectiveWeights(weights);
    const StateSpectralSums initial_sums =
        optimizationStateSpectralSums(system, initial);
    const StateSpectralSums final_sums =
        optimizationStateSpectralSums(system, final);
    StateObjectiveValue value;
    value.critical_log_growth = 0.5 * std::log(
        final_sums.critical_h_half_squared /
        initial_sums.critical_h_half_squared);
    value.characteristic_log_growth = 0.5 * std::log(
        (final_sums.enstrophy / final_sums.energy) /
        (initial_sums.enstrophy / initial_sums.energy));
    const double cutoff_fraction =
        final_sums.cutoff_energy / final_sums.energy;
    value.cutoff_penalty = std::log1p(
        cutoff_fraction / weights.cutoff_fraction_threshold);
    value.total = value.critical_log_growth +
                  weights.characteristic_scale_weight *
                      value.characteristic_log_growth -
                  weights.cutoff_penalty_weight * value.cutoff_penalty;
    if (!std::isfinite(value.total)) {
        throw std::runtime_error("State objective became non-finite");
    }
    return value;
}

inline OptimizationState stateLogCriticalGradient(
    const PseudospectralSystem& system,
    const OptimizationState& state) {
    const StateSpectralSums sums =
        optimizationStateSpectralSums(system, state);
    OptimizationState gradient = system.zeroState();
    const std::vector<WaveVector>& modes = system.gridModes();
    for (std::size_t i = 0; i < state.size(); ++i) {
        const int wave_squared = modes[i].normSquared();
        if (wave_squared == 0 ||
            stateMaximumComponent(modes[i]) > system.cutoff()) {
            continue;
        }
        const double wave_magnitude =
            std::sqrt(static_cast<double>(wave_squared));
        gradient[i] = state[i] *
                      (wave_magnitude / sums.critical_h_half_squared);
    }
    return gradient;
}

inline OptimizationState stateLogCharacteristicGradient(
    const PseudospectralSystem& system,
    const OptimizationState& state) {
    const StateSpectralSums sums =
        optimizationStateSpectralSums(system, state);
    OptimizationState gradient = system.zeroState();
    const std::vector<WaveVector>& modes = system.gridModes();
    for (std::size_t i = 0; i < state.size(); ++i) {
        const int wave_squared = modes[i].normSquared();
        if (wave_squared == 0 ||
            stateMaximumComponent(modes[i]) > system.cutoff()) {
            continue;
        }
        gradient[i] = state[i] *
            (0.5 * (static_cast<double>(wave_squared) / sums.enstrophy -
                    1.0 / sums.energy));
    }
    return gradient;
}

inline OptimizationState stateLogCutoffPenaltyGradient(
    const PseudospectralSystem& system,
    const OptimizationState& state,
    double cutoff_fraction_threshold) {
    if (!std::isfinite(cutoff_fraction_threshold) ||
        cutoff_fraction_threshold <= 0.0) {
        throw std::invalid_argument("Cutoff threshold must be positive");
    }
    const StateSpectralSums sums =
        optimizationStateSpectralSums(system, state);
    const double cutoff_fraction = sums.cutoff_energy / sums.energy;
    OptimizationState gradient = system.zeroState();
    const std::vector<WaveVector>& modes = system.gridModes();
    for (std::size_t i = 0; i < state.size(); ++i) {
        if (modes[i].normSquared() == 0 ||
            stateMaximumComponent(modes[i]) > system.cutoff()) {
            continue;
        }
        const double shell_indicator =
            stateMaximumComponent(modes[i]) == system.cutoff() ? 1.0 : 0.0;
        const double fraction_derivative =
            (shell_indicator - cutoff_fraction) / sums.energy;
        gradient[i] = state[i] *
            (fraction_derivative /
             (cutoff_fraction_threshold + cutoff_fraction));
    }
    return gradient;
}

inline OptimizationState terminalStateObjectiveGradient(
    const PseudospectralSystem& system,
    const OptimizationState& final,
    const StateObjectiveWeights& weights) {
    validateStateObjectiveWeights(weights);
    OptimizationState gradient = stateLogCriticalGradient(system, final);
    gradient = addOptimizationStates(
        gradient,
        stateLogCharacteristicGradient(system, final),
        weights.characteristic_scale_weight);
    gradient = addOptimizationStates(
        gradient,
        stateLogCutoffPenaltyGradient(
            system, final, weights.cutoff_fraction_threshold),
        -weights.cutoff_penalty_weight);
    return gradient;
}

inline OptimizationState initialStateObjectiveGradient(
    const PseudospectralSystem& system,
    const OptimizationState& initial,
    const StateObjectiveWeights& weights) {
    validateStateObjectiveWeights(weights);
    OptimizationState gradient = scaleOptimizationState(
        stateLogCriticalGradient(system, initial), -1.0);
    gradient = addOptimizationStates(
        gradient,
        stateLogCharacteristicGradient(system, initial),
        -weights.characteristic_scale_weight);
    return gradient;
}

inline OptimizationState projectStateObjectiveGradient(
    const PseudospectralSystem& system,
    const OptimizationState& state,
    const OptimizationState& gradient,
    int bandwidth) {
    requireOptimizationState(system, state);
    OptimizationState tangent = projectBandLimitedRealSolenoidal(
        system, gradient, bandwidth);
    const double state_norm_squared = stateRealInnerProduct(state, state);
    if (!std::isfinite(state_norm_squared) || state_norm_squared <= 0.0) {
        throw std::invalid_argument("Energy-sphere state has zero norm");
    }
    const double radial_coefficient =
        stateRealInnerProduct(state, tangent) / state_norm_squared;
    tangent = addOptimizationStates(tangent, state, -radial_coefficient);
    return projectBandLimitedRealSolenoidal(system, tangent, bandwidth);
}

inline OptimizationState stateEnergySphereDirection(
    const PseudospectralSystem& system,
    const OptimizationState& state,
    const OptimizationState& gradient,
    int bandwidth) {
    OptimizationState tangent = projectStateObjectiveGradient(
        system, state, gradient, bandwidth);
    const double tangent_norm = optimizationStateNorm(tangent);
    const double state_norm = optimizationStateNorm(state);
    if (!std::isfinite(tangent_norm) || tangent_norm <= 1e-14 ||
        !std::isfinite(state_norm) || state_norm <= 0.0) {
        throw std::runtime_error("Constrained state gradient is degenerate");
    }
    return scaleOptimizationState(tangent, state_norm / tangent_norm);
}

inline OptimizationState stateEnergySphereStep(
    const PseudospectralSystem& system,
    const OptimizationState& state,
    const OptimizationState& sphere_direction,
    double angle) {
    requireOptimizationState(system, state);
    requireOptimizationState(system, sphere_direction);
    if (!std::isfinite(angle)) {
        throw std::invalid_argument("Energy-sphere angle must be finite");
    }
    const double state_energy = system.energy(state);
    OptimizationState trial = scaleOptimizationState(state, std::cos(angle));
    trial = addOptimizationStates(
        trial, sphere_direction, std::sin(angle));
    return normalizeOptimizationEnergy(system, trial, state_energy);
}

inline OptimizationState liftOptimizationState(
    const PseudospectralSystem& source_system,
    const OptimizationState& source,
    const PseudospectralSystem& target_system) {
    requireOptimizationState(source_system, source);
    OptimizationState lifted = target_system.zeroState();
    const std::vector<WaveVector>& source_modes = source_system.gridModes();
    for (std::size_t i = 0; i < source.size(); ++i) {
        if (normSquared(source[i]) == 0.0) continue;
        const WaveVector& wave = source_modes[i];
        if (wave.normSquared() == 0 ||
            stateMaximumComponent(wave) > target_system.cutoff()) {
            throw std::invalid_argument(
                "Target grid cannot retain an optimization coefficient");
        }
        lifted[target_system.indexOf(wave)] = source[i];
    }
    return lifted;
}

inline std::size_t stateOptimizationDegreesOfFreedom(int bandwidth) {
    if (bandwidth < 1) {
        throw std::invalid_argument("Optimization bandwidth must be positive");
    }
    const std::size_t side =
        static_cast<std::size_t>(2 * bandwidth + 1);
    const std::size_t nonzero_modes = side * side * side - 1U;
    // Each conjugate pair has two complex solenoidal polarizations.
    return 2U * nonzero_modes;
}

}  // namespace ns_cascade

#endif
