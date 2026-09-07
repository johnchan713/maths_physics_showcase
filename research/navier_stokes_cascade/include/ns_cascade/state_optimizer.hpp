#ifndef NS_CASCADE_STATE_OPTIMIZER_HPP
#define NS_CASCADE_STATE_OPTIMIZER_HPP

#include "ns_cascade/pseudospectral.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

namespace ns_cascade {

using OptimizationState = PseudospectralSystem::State;

enum class SmoothSpectrumPathAggregation {
    Mean,
    SmoothMaximum
};

struct StateObjectiveWeights {
    double endpoint_critical_weight = 1.0;
    double characteristic_scale_weight = 0.15;
    double cutoff_penalty_weight = 0.04;
    double cutoff_fraction_threshold = 0.01;
    double profile_shape_penalty_weight = 0.05;
    double profile_path_penalty_weight = 0.05;
    double profile_path_temperature = 0.01;
    SmoothSpectrumPathAggregation profile_path_aggregation =
        SmoothSpectrumPathAggregation::SmoothMaximum;
    int profile_feature_count = 9;
    double profile_minimum_log_coordinate = -1.5;
    double profile_maximum_log_coordinate = 1.5;
    double profile_kernel_width = 0.35;
    double profile_log_floor = 1e-8;
};

struct StateObjectiveValue {
    double critical_log_growth = 0.0;
    double characteristic_log_growth = 0.0;
    double cutoff_penalty = 0.0;
    double profile_shape_penalty = 0.0;
    double profile_path_penalty = 0.0;
    double total = 0.0;
};

struct StateSpectralSums {
    double energy = 0.0;
    double enstrophy = 0.0;
    double critical_h_half_squared = 0.0;
    double cutoff_energy = 0.0;
};

struct SmoothSpectrumSignature {
    double energy = 0.0;
    double enstrophy = 0.0;
    double log_characteristic_wavenumber = 0.0;
    std::vector<double> features;
    std::vector<double> scale_derivative_sums;
};

struct SmoothSpectrumShapeComparison {
    SmoothSpectrumSignature initial;
    SmoothSpectrumSignature final;
    std::vector<double> log_feature_differences;
    double penalty = 0.0;
};

struct SmoothSpectrumPathComparison {
    std::vector<SmoothSpectrumShapeComparison> snapshots;
    std::vector<double> aggregation_weights;
    double average_penalty = 0.0;
    double smooth_maximum_penalty = 0.0;
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

inline void validateSmoothSpectrumParameters(
    const StateObjectiveWeights& weights) {
    if (!std::isfinite(weights.profile_shape_penalty_weight) ||
        weights.profile_shape_penalty_weight < 0.0 ||
        !std::isfinite(weights.profile_path_penalty_weight) ||
        weights.profile_path_penalty_weight < 0.0 ||
        !std::isfinite(weights.profile_path_temperature) ||
        weights.profile_path_temperature <= 0.0 ||
        (weights.profile_path_aggregation !=
             SmoothSpectrumPathAggregation::Mean &&
         weights.profile_path_aggregation !=
             SmoothSpectrumPathAggregation::SmoothMaximum) ||
        weights.profile_feature_count < 3 ||
        weights.profile_feature_count > 64 ||
        !std::isfinite(weights.profile_minimum_log_coordinate) ||
        !std::isfinite(weights.profile_maximum_log_coordinate) ||
        weights.profile_minimum_log_coordinate >=
            weights.profile_maximum_log_coordinate ||
        !std::isfinite(weights.profile_kernel_width) ||
        weights.profile_kernel_width <= 0.0 ||
        !std::isfinite(weights.profile_log_floor) ||
        weights.profile_log_floor <= 0.0 ||
        weights.profile_log_floor >= 1.0) {
        throw std::invalid_argument(
            "Smooth spectrum-shape parameters are invalid");
    }
}

inline double smoothSpectrumFeatureCenter(
    const StateObjectiveWeights& weights,
    int feature) {
    const double fraction = static_cast<double>(feature) /
                            (weights.profile_feature_count - 1);
    return weights.profile_minimum_log_coordinate +
           fraction * (weights.profile_maximum_log_coordinate -
                       weights.profile_minimum_log_coordinate);
}

inline SmoothSpectrumSignature smoothSpectrumSignature(
    const PseudospectralSystem& system,
    const OptimizationState& state,
    const StateObjectiveWeights& weights) {
    requireOptimizationState(system, state);
    validateSmoothSpectrumParameters(weights);
    const StateSpectralSums sums =
        optimizationStateSpectralSums(system, state);
    SmoothSpectrumSignature signature;
    signature.energy = sums.energy;
    signature.enstrophy = sums.enstrophy;
    signature.log_characteristic_wavenumber = 0.5 * std::log(
        sums.enstrophy / sums.energy);
    signature.features.assign(
        static_cast<std::size_t>(weights.profile_feature_count), 0.0);
    signature.scale_derivative_sums.assign(
        static_cast<std::size_t>(weights.profile_feature_count), 0.0);

    const double inverse_width_squared =
        1.0 / (weights.profile_kernel_width *
               weights.profile_kernel_width);
    const std::vector<WaveVector>& modes = system.gridModes();
    for (std::size_t i = 0; i < state.size(); ++i) {
        const int wave_squared = modes[i].normSquared();
        if (wave_squared == 0 ||
            stateMaximumComponent(modes[i]) > system.cutoff()) {
            continue;
        }
        const double mode_energy = 0.5 * normSquared(state[i]);
        if (mode_energy == 0.0) continue;
        const double log_rescaled_wave =
            0.5 * std::log(static_cast<double>(wave_squared)) -
            signature.log_characteristic_wavenumber;
        for (int feature = 0;
             feature < weights.profile_feature_count;
             ++feature) {
            const double centered = log_rescaled_wave -
                smoothSpectrumFeatureCenter(weights, feature);
            const double kernel = std::exp(
                -0.5 * centered * centered * inverse_width_squared);
            const std::size_t index = static_cast<std::size_t>(feature);
            signature.features[index] += mode_energy * kernel;
            signature.scale_derivative_sums[index] +=
                mode_energy * kernel * centered * inverse_width_squared;
        }
    }
    for (std::size_t feature = 0;
         feature < signature.features.size();
         ++feature) {
        signature.features[feature] /= signature.energy;
        signature.scale_derivative_sums[feature] /= signature.energy;
        if (!std::isfinite(signature.features[feature]) ||
            signature.features[feature] < 0.0 ||
            !std::isfinite(signature.scale_derivative_sums[feature])) {
            throw std::runtime_error(
                "Smooth spectrum signature became non-finite");
        }
    }
    return signature;
}

inline SmoothSpectrumShapeComparison compareSmoothSpectrumShapes(
    const PseudospectralSystem& system,
    const OptimizationState& initial,
    const OptimizationState& final,
    const StateObjectiveWeights& weights) {
    SmoothSpectrumShapeComparison comparison;
    comparison.initial = smoothSpectrumSignature(system, initial, weights);
    comparison.final = smoothSpectrumSignature(system, final, weights);
    comparison.log_feature_differences.assign(
        comparison.initial.features.size(), 0.0);
    const double inverse_count =
        1.0 / static_cast<double>(comparison.initial.features.size());
    for (std::size_t feature = 0;
         feature < comparison.log_feature_differences.size();
         ++feature) {
        const double difference = std::log(
            comparison.final.features[feature] +
            weights.profile_log_floor) -
            std::log(comparison.initial.features[feature] +
                     weights.profile_log_floor);
        comparison.log_feature_differences[feature] = difference;
        comparison.penalty += 0.5 * inverse_count * difference * difference;
    }
    if (!std::isfinite(comparison.penalty) || comparison.penalty < 0.0) {
        throw std::runtime_error(
            "Smooth spectrum-shape penalty became non-finite");
    }
    return comparison;
}

inline SmoothSpectrumPathComparison compareSmoothSpectrumPath(
    const PseudospectralSystem& system,
    const OptimizationState& initial,
    const std::vector<OptimizationState>& snapshots,
    const StateObjectiveWeights& weights) {
    if (snapshots.empty()) {
        throw std::invalid_argument(
            "Smooth spectrum path requires at least one snapshot");
    }
    SmoothSpectrumPathComparison path;
    path.snapshots.reserve(snapshots.size());
    double maximum_penalty = 0.0;
    for (std::size_t snapshot = 0; snapshot < snapshots.size(); ++snapshot) {
        path.snapshots.push_back(compareSmoothSpectrumShapes(
            system, initial, snapshots[snapshot], weights));
        path.average_penalty += path.snapshots.back().penalty;
        maximum_penalty = std::max(
            maximum_penalty, path.snapshots.back().penalty);
    }
    path.average_penalty /= static_cast<double>(path.snapshots.size());

    double exponential_sum = 0.0;
    path.aggregation_weights.reserve(path.snapshots.size());
    for (std::size_t snapshot = 0; snapshot < path.snapshots.size();
         ++snapshot) {
        const double exponential = std::exp(
            (path.snapshots[snapshot].penalty - maximum_penalty) /
            weights.profile_path_temperature);
        path.aggregation_weights.push_back(exponential);
        exponential_sum += exponential;
    }
    for (std::size_t snapshot = 0;
         snapshot < path.aggregation_weights.size(); ++snapshot) {
        path.aggregation_weights[snapshot] /= exponential_sum;
    }
    path.smooth_maximum_penalty = maximum_penalty +
        weights.profile_path_temperature *
            (std::log(exponential_sum) -
             std::log(static_cast<double>(path.snapshots.size())));
    if (path.smooth_maximum_penalty < 0.0 &&
        path.smooth_maximum_penalty >
            -64.0 * std::numeric_limits<double>::epsilon()) {
        path.smooth_maximum_penalty = 0.0;
    }
    if (!std::isfinite(path.average_penalty) ||
        path.average_penalty < 0.0 ||
        !std::isfinite(path.smooth_maximum_penalty) ||
        path.smooth_maximum_penalty < 0.0 ||
        !std::isfinite(exponential_sum) || exponential_sum <= 0.0) {
        throw std::runtime_error(
            "Smooth spectrum path penalty became non-finite");
    }
    return path;
}

inline double smoothSpectrumPathPenalty(
    const SmoothSpectrumPathComparison& path,
    const StateObjectiveWeights& weights) {
    return weights.profile_path_aggregation ==
                   SmoothSpectrumPathAggregation::Mean
               ? path.average_penalty
               : path.smooth_maximum_penalty;
}

inline double smoothSpectrumPathGradientWeight(
    const SmoothSpectrumPathComparison& path,
    const StateObjectiveWeights& weights,
    std::size_t snapshot) {
    if (snapshot >= path.snapshots.size() ||
        path.aggregation_weights.size() != path.snapshots.size()) {
        throw std::invalid_argument(
            "Smooth spectrum path gradient index is invalid");
    }
    return weights.profile_path_aggregation ==
                   SmoothSpectrumPathAggregation::Mean
               ? 1.0 / static_cast<double>(path.snapshots.size())
               : path.aggregation_weights[snapshot];
}

inline OptimizationState smoothSpectrumLogFeatureGradient(
    const PseudospectralSystem& system,
    const OptimizationState& state,
    const StateObjectiveWeights& weights,
    const SmoothSpectrumSignature& signature,
    const std::vector<double>& log_feature_weights) {
    requireOptimizationState(system, state);
    if (signature.features.size() !=
            static_cast<std::size_t>(weights.profile_feature_count) ||
        signature.scale_derivative_sums.size() !=
            signature.features.size() ||
        log_feature_weights.size() != signature.features.size()) {
        throw std::invalid_argument(
            "Smooth spectrum gradient has incompatible feature data");
    }
    OptimizationState gradient = system.zeroState();
    const std::vector<WaveVector>& modes = system.gridModes();
    const double inverse_width_squared =
        1.0 / (weights.profile_kernel_width *
               weights.profile_kernel_width);
    for (std::size_t i = 0; i < state.size(); ++i) {
        const int wave_squared = modes[i].normSquared();
        if (wave_squared == 0 ||
            stateMaximumComponent(modes[i]) > system.cutoff()) {
            continue;
        }
        const double log_rescaled_wave =
            0.5 * std::log(static_cast<double>(wave_squared)) -
            signature.log_characteristic_wavenumber;
        const double log_scale_gradient_coefficient = 0.5 *
            (static_cast<double>(wave_squared) / signature.enstrophy -
             1.0 / signature.energy);
        double coefficient = 0.0;
        for (int feature = 0;
             feature < weights.profile_feature_count;
             ++feature) {
            const std::size_t index = static_cast<std::size_t>(feature);
            const double centered = log_rescaled_wave -
                smoothSpectrumFeatureCenter(weights, feature);
            const double kernel = std::exp(
                -0.5 * centered * centered * inverse_width_squared);
            const double feature_gradient_coefficient =
                (kernel - signature.features[index]) / signature.energy +
                signature.scale_derivative_sums[index] *
                    log_scale_gradient_coefficient;
            coefficient += log_feature_weights[index] *
                feature_gradient_coefficient /
                (signature.features[index] + weights.profile_log_floor);
        }
        gradient[i] = state[i] * coefficient;
    }
    return gradient;
}

inline OptimizationState smoothSpectrumShapePenaltyGradient(
    const PseudospectralSystem& system,
    const OptimizationState& state,
    const StateObjectiveWeights& weights,
    const SmoothSpectrumShapeComparison& comparison,
    bool with_respect_to_final) {
    const double inverse_count =
        1.0 / static_cast<double>(comparison.log_feature_differences.size());
    std::vector<double> feature_weights =
        comparison.log_feature_differences;
    const double sign = with_respect_to_final ? 1.0 : -1.0;
    for (std::size_t feature = 0;
         feature < feature_weights.size();
         ++feature) {
        feature_weights[feature] *= sign * inverse_count;
    }
    return smoothSpectrumLogFeatureGradient(
        system,
        state,
        weights,
        with_respect_to_final ? comparison.final : comparison.initial,
        feature_weights);
}

inline void validateStateObjectiveWeights(
    const StateObjectiveWeights& weights) {
    if (!std::isfinite(weights.endpoint_critical_weight) ||
        weights.endpoint_critical_weight < 0.0 ||
        !std::isfinite(weights.characteristic_scale_weight) ||
        weights.characteristic_scale_weight < 0.0 ||
        !std::isfinite(weights.cutoff_penalty_weight) ||
        weights.cutoff_penalty_weight < 0.0 ||
        !std::isfinite(weights.cutoff_fraction_threshold) ||
        weights.cutoff_fraction_threshold <= 0.0 ||
        weights.cutoff_fraction_threshold >= 1.0) {
        throw std::invalid_argument("State-objective weights are invalid");
    }
    validateSmoothSpectrumParameters(weights);
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
    value.profile_shape_penalty = compareSmoothSpectrumShapes(
        system, initial, final, weights).penalty;
    value.profile_path_penalty = 0.0;
    value.total = weights.endpoint_critical_weight * value.critical_log_growth +
                  weights.characteristic_scale_weight *
                      value.characteristic_log_growth -
                  weights.cutoff_penalty_weight * value.cutoff_penalty -
                  weights.profile_shape_penalty_weight *
                      value.profile_shape_penalty;
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
    const OptimizationState& initial,
    const OptimizationState& final,
    const StateObjectiveWeights& weights) {
    validateStateObjectiveWeights(weights);
    OptimizationState gradient = scaleOptimizationState(
        stateLogCriticalGradient(system, final), weights.endpoint_critical_weight);
    gradient = addOptimizationStates(
        gradient,
        stateLogCharacteristicGradient(system, final),
        weights.characteristic_scale_weight);
    gradient = addOptimizationStates(
        gradient,
        stateLogCutoffPenaltyGradient(
            system, final, weights.cutoff_fraction_threshold),
        -weights.cutoff_penalty_weight);
    const SmoothSpectrumShapeComparison comparison =
        compareSmoothSpectrumShapes(system, initial, final, weights);
    gradient = addOptimizationStates(
        gradient,
        smoothSpectrumShapePenaltyGradient(
            system, final, weights, comparison, true),
        -weights.profile_shape_penalty_weight);
    return gradient;
}

inline OptimizationState initialStateObjectiveGradient(
    const PseudospectralSystem& system,
    const OptimizationState& initial,
    const OptimizationState& final,
    const StateObjectiveWeights& weights) {
    validateStateObjectiveWeights(weights);
    OptimizationState gradient = scaleOptimizationState(
        stateLogCriticalGradient(system, initial), -weights.endpoint_critical_weight);
    gradient = addOptimizationStates(
        gradient,
        stateLogCharacteristicGradient(system, initial),
        -weights.characteristic_scale_weight);
    const SmoothSpectrumShapeComparison comparison =
        compareSmoothSpectrumShapes(system, initial, final, weights);
    gradient = addOptimizationStates(
        gradient,
        smoothSpectrumShapePenaltyGradient(
            system, initial, weights, comparison, false),
        -weights.profile_shape_penalty_weight);
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

inline std::uint64_t deterministicStateMix(std::uint64_t value) {
    value += UINT64_C(0x9e3779b97f4a7c15);
    value = (value ^ (value >> 30U)) * UINT64_C(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27U)) * UINT64_C(0x94d049bb133111eb);
    return value ^ (value >> 31U);
}

inline double deterministicStateCoordinate(std::uint64_t key) {
    const std::uint64_t word = deterministicStateMix(key);
    const double unit = static_cast<double>(word >> 11U) /
                        9007199254740992.0;
    return 2.0 * unit - 1.0;
}

inline OptimizationState deterministicOptimizationDirection(
    const PseudospectralSystem& system,
    const OptimizationState& state,
    int bandwidth,
    int start_index) {
    requireOptimizationState(system, state);
    if (start_index < 1) {
        throw std::invalid_argument(
            "Deterministic optimization start index must be positive");
    }
    OptimizationState raw = system.zeroState();
    const std::vector<WaveVector>& modes = system.gridModes();
    const std::uint64_t start_key = deterministicStateMix(
        static_cast<std::uint64_t>(start_index));
    for (std::size_t i = 0; i < raw.size(); ++i) {
        const WaveVector& wave = modes[i];
        if (wave.normSquared() == 0 ||
            stateMaximumComponent(wave) > bandwidth) {
            continue;
        }
        std::uint64_t key = start_key;
        key ^= deterministicStateMix(
            static_cast<std::uint64_t>(wave.x + 4096));
        key ^= deterministicStateMix(
            static_cast<std::uint64_t>(wave.y + 8192));
        key ^= deterministicStateMix(
            static_cast<std::uint64_t>(wave.z + 16384));
        const double values[6] = {
            deterministicStateCoordinate(key + UINT64_C(0)),
            deterministicStateCoordinate(key + UINT64_C(1)),
            deterministicStateCoordinate(key + UINT64_C(2)),
            deterministicStateCoordinate(key + UINT64_C(3)),
            deterministicStateCoordinate(key + UINT64_C(4)),
            deterministicStateCoordinate(key + UINT64_C(5))};
        raw[i] = ComplexVector(
            Complex(values[0], values[1]),
            Complex(values[2], values[3]),
            Complex(values[4], values[5]));
    }
    return stateEnergySphereDirection(system, state, raw, bandwidth);
}

inline OptimizationState deterministicOptimizationStart(
    const PseudospectralSystem& system,
    const OptimizationState& base_state,
    int bandwidth,
    int start_index,
    double angle) {
    if (start_index < 0 || !std::isfinite(angle) || angle < 0.0 ||
        angle >= 1.0) {
        throw std::invalid_argument(
            "Deterministic optimization start parameters are invalid");
    }
    if (start_index == 0 || angle == 0.0) return base_state;
    const OptimizationState direction = deterministicOptimizationDirection(
        system, base_state, bandwidth, start_index);
    return stateEnergySphereStep(system, base_state, direction, angle);
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
