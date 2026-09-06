#include "ns_cascade/optimization_state_csv.hpp"
#include "ns_cascade/state_optimizer.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using State = ns_cascade::OptimizationState;

void expect(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

double relativeScalarDifference(double left, double right) {
    return std::abs(left - right) /
           std::max(1e-12, std::max(std::abs(left), std::abs(right)));
}

double outOfBandNorm(const ns_cascade::PseudospectralSystem& system,
                     const State& state,
                     int bandwidth) {
    double norm_squared = 0.0;
    const std::vector<ns_cascade::WaveVector>& modes = system.gridModes();
    for (std::size_t i = 0; i < state.size(); ++i) {
        if (modes[i].normSquared() == 0 ||
            ns_cascade::stateMaximumComponent(modes[i]) > bandwidth) {
            norm_squared += ns_cascade::normSquared(state[i]);
        }
    }
    return std::sqrt(norm_squared);
}

State makeInitial(const ns_cascade::PseudospectralSystem& system,
                  double energy = 4.0) {
    const State raw = system.interactingWavePacketState(
        ns_cascade::WavePacketParameters(1.1, 1, 0.75, 0.3), energy);
    return ns_cascade::makeBandLimitedOptimizationState(
        system, raw, 2, energy);
}

State makeSingleShellState(
    const ns_cascade::PseudospectralSystem& system,
    int wave_number,
    double energy = 4.0) {
    State state = system.zeroState();
    const ns_cascade::WaveVector positive(wave_number, 0, 0);
    const ns_cascade::WaveVector negative(-wave_number, 0, 0);
    const ns_cascade::ComplexVector coefficient(
        ns_cascade::Complex(), ns_cascade::Complex(1.0, 0.0),
        ns_cascade::Complex());
    state[system.indexOf(positive)] = coefficient;
    state[system.indexOf(negative)] = coefficient;
    return ns_cascade::normalizeOptimizationEnergy(system, state, energy);
}

void evolve(const ns_cascade::PseudospectralSystem& system,
            State& state,
            int steps,
            double time_step) {
    for (int step = 0; step < steps; ++step) {
        system.stepRungeKutta4(state, time_step);
    }
}

void testBandLimitAndEnergySphere() {
    const ns_cascade::PseudospectralSystem system(16, 0.02, 5);
    const State state = makeInitial(system);
    expect(std::abs(system.energy(state) - 4.0) < 5e-13,
           "Band-limited state did not preserve target energy");
    expect(outOfBandNorm(system, state, 2) == 0.0,
           "Band-limited state retained a forbidden coefficient");
    expect(system.divergenceDefect(state) < 2e-13,
           "Band-limited state is not divergence-free");
    expect(system.realityDefect(state) < 2e-13,
           "Band-limited state lost Fourier reality");
    expect(ns_cascade::stateOptimizationDegreesOfFreedom(2) == 248U,
           "The K=2 real-solenoidal degree count is wrong");
    expect(ns_cascade::stateOptimizationDegreesOfFreedom(3) == 684U,
           "The K=3 real-solenoidal degree count is wrong");

    const State raw_gradient = system.interactingWavePacketState(
        ns_cascade::WavePacketParameters(1.35, 2, 1.2, -0.7), 1.0);
    const State direction = ns_cascade::stateEnergySphereDirection(
        system, state, raw_gradient, 2);
    const double radius = ns_cascade::optimizationStateNorm(state);
    expect(std::abs(ns_cascade::stateRealInnerProduct(state, direction)) <
               2e-13 * radius * radius,
           "State-space direction is not tangent to fixed energy");
    expect(std::abs(ns_cascade::optimizationStateNorm(direction) - radius) <
               2e-13 * radius,
           "State-space direction does not have the sphere radius");
    expect(outOfBandNorm(system, direction, 2) == 0.0,
           "State-space direction left the seed band");

    const State trial = ns_cascade::stateEnergySphereStep(
        system, state, direction, 0.13);
    expect(std::abs(system.energy(trial) - system.energy(state)) < 5e-13,
           "Geodesic state step changed the fixed energy");
    expect(outOfBandNorm(system, trial, 2) == 0.0,
           "Geodesic state step left the seed band");
    expect(system.divergenceDefect(trial) < 2e-13,
           "Geodesic state step developed divergence");
    expect(system.realityDefect(trial) < 2e-13,
           "Geodesic state step lost Fourier reality");
}

void testDeterministicMultipleStarts() {
    const ns_cascade::PseudospectralSystem system(16, 0.02, 5);
    const State base = makeInitial(system);
    const State start_one = ns_cascade::deterministicOptimizationStart(
        system, base, 2, 1, 0.35);
    const State start_one_repeat =
        ns_cascade::deterministicOptimizationStart(
            system, base, 2, 1, 0.35);
    const State start_two = ns_cascade::deterministicOptimizationStart(
        system, base, 2, 2, 0.35);
    expect(ns_cascade::optimizationStateNorm(
               ns_cascade::addOptimizationStates(
                   start_one, start_one_repeat, -1.0)) == 0.0,
           "Repeated deterministic starts are not bitwise identical");
    expect(ns_cascade::optimizationStateNorm(
               ns_cascade::addOptimizationStates(
                   start_one, start_two, -1.0)) > 1e-3,
           "Distinct deterministic start indices produced the same state");
    expect(std::abs(system.energy(start_one) - system.energy(base)) < 5e-13,
           "Deterministic start changed fixed energy");
    expect(outOfBandNorm(system, start_one, 2) == 0.0,
           "Deterministic start left the seed band");
    expect(system.divergenceDefect(start_one) < 2e-13,
           "Deterministic start developed divergence");
    expect(system.realityDefect(start_one) < 2e-13,
           "Deterministic start lost Fourier reality");
}

void testStaticObjectiveGradient() {
    const ns_cascade::PseudospectralSystem system(16, 0.02, 5);
    const State initial = makeInitial(system);
    State final = initial;
    evolve(system, final, 10, 0.0004);
    const State initial_direction = system.interactingWavePacketState(
        ns_cascade::WavePacketParameters(1.4, 1, 1.1, -0.4), 1.0);
    const State final_direction = system.interactingWavePacketState(
        ns_cascade::WavePacketParameters(0.85, 2, 0.6, 0.8), 1.0);
    ns_cascade::StateObjectiveWeights weights;
    weights.characteristic_scale_weight = 0.2;
    weights.cutoff_penalty_weight = 0.07;
    weights.cutoff_fraction_threshold = 0.01;

    const ns_cascade::StateObjectiveValue value =
        ns_cascade::evaluateStateObjective(system, initial, final, weights);
    const double expected_critical = std::log(
        system.criticalHOneHalf(final) /
        system.criticalHOneHalf(initial));
    const double initial_scale = std::sqrt(
        system.enstrophy(initial) / system.energy(initial));
    const double final_scale = std::sqrt(
        system.enstrophy(final) / system.energy(final));
    expect(std::abs(value.critical_log_growth - expected_critical) < 2e-15,
           "State objective disagrees with the H1/2 diagnostic");
    expect(std::abs(value.characteristic_log_growth -
                    std::log(final_scale / initial_scale)) < 2e-15,
           "State objective disagrees with the characteristic scale");

    const State initial_gradient =
        ns_cascade::initialStateObjectiveGradient(
            system, initial, final, weights);
    const State final_gradient =
        ns_cascade::terminalStateObjectiveGradient(
            system, initial, final, weights);
    const double analytical =
        ns_cascade::stateRealInnerProduct(initial_gradient, initial_direction) +
        ns_cascade::stateRealInnerProduct(final_gradient, final_direction);
    const double epsilon = 1e-6;
    const ns_cascade::StateObjectiveValue plus =
        ns_cascade::evaluateStateObjective(
            system,
            ns_cascade::addOptimizationStates(
                initial, initial_direction, epsilon),
            ns_cascade::addOptimizationStates(final, final_direction, epsilon),
            weights);
    const ns_cascade::StateObjectiveValue minus =
        ns_cascade::evaluateStateObjective(
            system,
            ns_cascade::addOptimizationStates(
                initial, initial_direction, -epsilon),
            ns_cascade::addOptimizationStates(final, final_direction, -epsilon),
            weights);
    const double finite_difference =
        (plus.total - minus.total) / (2.0 * epsilon);
    expect(relativeScalarDifference(analytical, finite_difference) < 2e-8,
           "Smooth state-objective gradient failed finite differences");
}

void testSmoothSpectrumShapeInvariances() {
    const ns_cascade::PseudospectralSystem system(16, 0.02, 5);
    const ns_cascade::StateObjectiveWeights weights;
    const State state = makeInitial(system);
    const State scaled = ns_cascade::scaleOptimizationState(state, 2.75);
    const ns_cascade::SmoothSpectrumSignature signature =
        ns_cascade::smoothSpectrumSignature(system, state, weights);
    const ns_cascade::SmoothSpectrumSignature scaled_signature =
        ns_cascade::smoothSpectrumSignature(system, scaled, weights);
    for (std::size_t feature = 0;
         feature < signature.features.size();
         ++feature) {
        expect(std::abs(signature.features[feature] -
                        scaled_signature.features[feature]) < 3e-15,
               "Smooth spectrum signature depends on amplitude");
    }

    const State shell_one = makeSingleShellState(system, 1);
    const State shell_two = makeSingleShellState(system, 2);
    const ns_cascade::SmoothSpectrumSignature shell_one_signature =
        ns_cascade::smoothSpectrumSignature(system, shell_one, weights);
    const ns_cascade::SmoothSpectrumSignature shell_two_signature =
        ns_cascade::smoothSpectrumSignature(system, shell_two, weights);
    for (std::size_t feature = 0;
         feature < shell_one_signature.features.size();
         ++feature) {
        expect(std::abs(shell_one_signature.features[feature] -
                        shell_two_signature.features[feature]) < 3e-15,
               "Smooth spectrum signature depends on a pure scale shift");
    }
    expect(ns_cascade::compareSmoothSpectrumShapes(
               system, shell_one, shell_two, weights).penalty < 2e-29,
           "Pure shell rescaling has nonzero smooth shape penalty");
}

void testSmoothSpectrumShapeGradient() {
    const ns_cascade::PseudospectralSystem system(16, 0.02, 5);
    const ns_cascade::StateObjectiveWeights weights;
    const State initial = makeInitial(system);
    State final = initial;
    evolve(system, final, 12, 0.0004);
    const State initial_direction = system.interactingWavePacketState(
        ns_cascade::WavePacketParameters(1.45, 2, 0.95, -0.5), 1.0);
    const State final_direction = system.interactingWavePacketState(
        ns_cascade::WavePacketParameters(0.8, 1, 0.65, 0.9), 1.0);
    const ns_cascade::SmoothSpectrumShapeComparison comparison =
        ns_cascade::compareSmoothSpectrumShapes(
            system, initial, final, weights);
    const State initial_gradient =
        ns_cascade::smoothSpectrumShapePenaltyGradient(
            system, initial, weights, comparison, false);
    const State final_gradient =
        ns_cascade::smoothSpectrumShapePenaltyGradient(
            system, final, weights, comparison, true);
    const double analytical =
        ns_cascade::stateRealInnerProduct(
            initial_gradient, initial_direction) +
        ns_cascade::stateRealInnerProduct(final_gradient, final_direction);
    const double epsilon = 1e-6;
    const double plus = ns_cascade::compareSmoothSpectrumShapes(
        system,
        ns_cascade::addOptimizationStates(
            initial, initial_direction, epsilon),
        ns_cascade::addOptimizationStates(
            final, final_direction, epsilon),
        weights).penalty;
    const double minus = ns_cascade::compareSmoothSpectrumShapes(
        system,
        ns_cascade::addOptimizationStates(
            initial, initial_direction, -epsilon),
        ns_cascade::addOptimizationStates(
            final, final_direction, -epsilon),
        weights).penalty;
    const double finite_difference =
        (plus - minus) / (2.0 * epsilon);
    expect(relativeScalarDifference(analytical, finite_difference) < 3e-8,
           "Smooth spectrum-shape gradient failed finite differences");
}

void testSmoothSpectrumPathAggregation() {
    const ns_cascade::PseudospectralSystem system(16, 0.02, 5);
    ns_cascade::StateObjectiveWeights weights;
    weights.profile_path_temperature = 0.01;
    const State initial = makeInitial(system);
    std::vector<State> snapshots;
    State state = initial;
    for (int snapshot = 0; snapshot < 4; ++snapshot) {
        evolve(system, state, 4, 0.0004);
        snapshots.push_back(state);
    }
    const ns_cascade::SmoothSpectrumPathComparison path =
        ns_cascade::compareSmoothSpectrumPath(
            system, initial, snapshots, weights);
    double maximum_penalty = 0.0;
    double weight_sum = 0.0;
    std::size_t maximum_index = 0;
    for (std::size_t snapshot = 0; snapshot < path.snapshots.size();
         ++snapshot) {
        if (path.snapshots[snapshot].penalty > maximum_penalty) {
            maximum_penalty = path.snapshots[snapshot].penalty;
            maximum_index = snapshot;
        }
        expect(path.aggregation_weights[snapshot] > 0.0,
               "Smooth path maximum assigned a nonpositive weight");
        weight_sum += path.aggregation_weights[snapshot];
    }
    expect(std::abs(weight_sum - 1.0) < 2e-15,
           "Smooth path maximum weights do not sum to one");
    expect(path.smooth_maximum_penalty + 2e-15 >= path.average_penalty &&
               path.smooth_maximum_penalty <= maximum_penalty + 2e-15,
           "Smooth path maximum is outside its mean/maximum bounds");
    expect(ns_cascade::smoothSpectrumPathPenalty(path, weights) ==
               path.smooth_maximum_penalty,
           "Default path aggregation did not select the smooth maximum");
    for (std::size_t snapshot = 0; snapshot < path.snapshots.size();
         ++snapshot) {
        expect(path.aggregation_weights[maximum_index] + 2e-15 >=
                   path.aggregation_weights[snapshot],
               "Worst path snapshot did not receive the largest weight");
    }

    const std::vector<State> unchanged_snapshots(4, initial);
    const ns_cascade::SmoothSpectrumPathComparison unchanged =
        ns_cascade::compareSmoothSpectrumPath(
            system, initial, unchanged_snapshots, weights);
    expect(unchanged.smooth_maximum_penalty == 0.0,
           "Zero path changes produced a nonzero smooth maximum");
    for (std::size_t snapshot = 0;
         snapshot < unchanged.aggregation_weights.size(); ++snapshot) {
        expect(std::abs(unchanged.aggregation_weights[snapshot] - 0.25) <
                   2e-15,
               "Equal path penalties did not receive equal weights");
    }

    weights.profile_path_temperature = 0.0;
    bool rejected_zero_temperature = false;
    try {
        ns_cascade::compareSmoothSpectrumPath(
            system, initial, snapshots, weights);
    } catch (const std::invalid_argument&) {
        rejected_zero_temperature = true;
    }
    expect(rejected_zero_temperature,
           "Smooth path maximum accepted zero temperature");

    weights.profile_path_temperature = 0.01;
    weights.profile_path_aggregation =
        ns_cascade::SmoothSpectrumPathAggregation::Mean;
    expect(ns_cascade::smoothSpectrumPathPenalty(path, weights) ==
               path.average_penalty,
           "Historical mean path aggregation is not reproducible");
    for (std::size_t snapshot = 0; snapshot < path.snapshots.size();
         ++snapshot) {
        expect(std::abs(ns_cascade::smoothSpectrumPathGradientWeight(
                            path, weights, snapshot) - 0.25) < 2e-15,
               "Mean path aggregation has a nonuniform gradient weight");
    }
}

void testTrajectoryObjectiveGradient() {
    const ns_cascade::PseudospectralSystem system(16, 0.02, 5);
    const State initial = makeInitial(system);
    const State raw_direction = system.interactingWavePacketState(
        ns_cascade::WavePacketParameters(1.35, 2, 1.2, -0.7), 1.0);
    const State direction = ns_cascade::stateEnergySphereDirection(
        system, initial, raw_direction, 2);
    ns_cascade::StateObjectiveWeights weights;
    const int step_count = 8;
    const double time_step = 0.0004;
    std::vector<State> trajectory;
    std::vector<State> profile_path_states;
    std::vector<int> profile_path_state_indices;
    State final = initial;
    for (int step = 0; step < step_count; ++step) {
        trajectory.push_back(final);
        system.stepRungeKutta4(final, time_step);
        if ((step + 1) % 2 == 0) {
            profile_path_states.push_back(final);
            profile_path_state_indices.push_back(step + 1);
        }
    }

    State reverse_gradient =
        ns_cascade::terminalStateObjectiveGradient(
            system, initial, final, weights);
    const ns_cascade::SmoothSpectrumPathComparison path_comparison =
        ns_cascade::compareSmoothSpectrumPath(
            system, initial, profile_path_states, weights);
    for (int step = step_count; step-- > 0;) {
        for (std::size_t snapshot = 0;
             snapshot < profile_path_states.size();
             ++snapshot) {
            if (profile_path_state_indices[snapshot] != step + 1) continue;
            reverse_gradient = ns_cascade::addOptimizationStates(
                reverse_gradient,
                ns_cascade::smoothSpectrumShapePenaltyGradient(
                    system,
                    profile_path_states[snapshot],
                    weights,
                    path_comparison.snapshots[snapshot],
                    true),
                -weights.profile_path_penalty_weight *
                    ns_cascade::smoothSpectrumPathGradientWeight(
                        path_comparison, weights, snapshot));
        }
        reverse_gradient = system.adjointRungeKutta4Step(
            trajectory[static_cast<std::size_t>(step)],
            reverse_gradient,
            time_step);
    }
    reverse_gradient = ns_cascade::addOptimizationStates(
        reverse_gradient,
        ns_cascade::initialStateObjectiveGradient(
            system, initial, final, weights),
        1.0);
    for (std::size_t snapshot = 0;
         snapshot < profile_path_states.size();
         ++snapshot) {
        reverse_gradient = ns_cascade::addOptimizationStates(
            reverse_gradient,
            ns_cascade::smoothSpectrumShapePenaltyGradient(
                system,
                initial,
                weights,
                path_comparison.snapshots[snapshot],
                false),
            -weights.profile_path_penalty_weight *
                ns_cascade::smoothSpectrumPathGradientWeight(
                    path_comparison, weights, snapshot));
    }
    const double analytical =
        ns_cascade::stateRealInnerProduct(reverse_gradient, direction);

    const double epsilon = 2e-4;
    const State plus_initial = ns_cascade::stateEnergySphereStep(
        system, initial, direction, epsilon);
    const State minus_initial = ns_cascade::stateEnergySphereStep(
        system, initial, direction, -epsilon);
    State plus = plus_initial;
    State minus = minus_initial;
    std::vector<State> plus_path;
    std::vector<State> minus_path;
    for (int step = 0; step < step_count; ++step) {
        system.stepRungeKutta4(plus, time_step);
        system.stepRungeKutta4(minus, time_step);
        if ((step + 1) % 2 == 0) {
            plus_path.push_back(plus);
            minus_path.push_back(minus);
        }
    }
    const double plus_endpoint = ns_cascade::evaluateStateObjective(
        system, plus_initial, plus, weights).total;
    const double minus_endpoint = ns_cascade::evaluateStateObjective(
        system, minus_initial, minus, weights).total;
    const ns_cascade::SmoothSpectrumPathComparison plus_path_comparison =
        ns_cascade::compareSmoothSpectrumPath(
            system, plus_initial, plus_path, weights);
    const ns_cascade::SmoothSpectrumPathComparison minus_path_comparison =
        ns_cascade::compareSmoothSpectrumPath(
            system, minus_initial, minus_path, weights);
    const double plus_objective = plus_endpoint -
        weights.profile_path_penalty_weight *
            ns_cascade::smoothSpectrumPathPenalty(
                plus_path_comparison, weights);
    const double minus_objective = minus_endpoint -
        weights.profile_path_penalty_weight *
            ns_cascade::smoothSpectrumPathPenalty(
                minus_path_comparison, weights);
    const double finite_difference =
        (plus_objective - minus_objective) /
        (2.0 * epsilon);
    expect(relativeScalarDifference(analytical, finite_difference) < 2e-6,
           "Path-dependent adjoint gradient failed a geodesic difference");
}

void testLiftToFineGrid() {
    const ns_cascade::PseudospectralSystem coarse(16, 0.02, 5);
    const ns_cascade::PseudospectralSystem fine(32, 0.02, 10);
    const State state = makeInitial(coarse);
    const State lifted =
        ns_cascade::liftOptimizationState(coarse, state, fine);
    expect(std::abs(coarse.energy(state) - fine.energy(lifted)) < 2e-14,
           "Grid lift changed Fourier energy");
    expect(fine.divergenceDefect(lifted) < 2e-13,
           "Lifted state is not divergence-free");
    expect(fine.realityDefect(lifted) < 2e-13,
           "Lifted state lost Fourier reality");
    const std::vector<ns_cascade::WaveVector>& modes = coarse.gridModes();
    for (std::size_t i = 0; i < state.size(); ++i) {
        if (ns_cascade::normSquared(state[i]) == 0.0) continue;
        const ns_cascade::ComplexVector difference =
            state[i] - lifted[fine.indexOf(modes[i])];
        expect(ns_cascade::normSquared(difference) == 0.0,
               "Grid lift changed a Fourier coefficient");
    }
}

void testOptimizationStateCsvRoundTrip() {
    const ns_cascade::PseudospectralSystem system(16, 0.02, 5);
    const State state = makeInitial(system);
    const std::string path = "state_optimizer_round_trip_test.csv";
    ns_cascade::writeOptimizationStateCsv(
        path, system, state, "wave-packets", 2, 4.0);
    const ns_cascade::LoadedOptimizationState loaded =
        ns_cascade::readOptimizationStateCsv(path, system);
    const int remove_result = std::remove(path.c_str());
    expect(remove_result == 0,
           "Could not remove state-CSV round-trip fixture");
    expect(loaded.metadata.family == "wave-packets" &&
               loaded.metadata.source_grid == 16 &&
               loaded.metadata.simulation_cutoff == 5 &&
               loaded.metadata.seed_bandwidth == 2 &&
               loaded.metadata.target_energy == 4.0,
           "State-CSV round trip changed metadata");
    expect(ns_cascade::optimizationStateNorm(
               ns_cascade::addOptimizationStates(
                   state, loaded.state, -1.0)) == 0.0,
           "State-CSV round trip changed a Fourier coefficient");

    {
        std::ofstream corrupt(path.c_str());
        corrupt << ns_cascade::optimizationStateCsvHeader() << '\n';
    }
    bool rejected_corrupt_state = false;
    try {
        ns_cascade::readOptimizationStateCsv(path, system);
    } catch (const std::exception&) {
        rejected_corrupt_state = true;
    }
    expect(std::remove(path.c_str()) == 0,
           "Could not remove corrupt state-CSV fixture");
    expect(rejected_corrupt_state,
           "State-CSV reader accepted an incomplete coefficient file");

    const State broad_state = system.interactingWavePacketState(
        ns_cascade::WavePacketParameters(1.1, 1, 0.75, 0.3), 4.0);
    bool rejected_out_of_band_state = false;
    try {
        ns_cascade::writeOptimizationStateCsv(
            path, system, broad_state, "wave-packets", 2, 4.0);
    } catch (const std::invalid_argument&) {
        rejected_out_of_band_state = true;
    }
    expect(rejected_out_of_band_state,
           "State-CSV writer accepted an out-of-band coefficient");
}

}  // namespace

int main() {
    try {
        testBandLimitAndEnergySphere();
        testDeterministicMultipleStarts();
        testStaticObjectiveGradient();
        testSmoothSpectrumShapeInvariances();
        testSmoothSpectrumShapeGradient();
        testSmoothSpectrumPathAggregation();
        testTrajectoryObjectiveGradient();
        testLiftToFineGrid();
        testOptimizationStateCsvRoundTrip();
        std::cout << "All state-optimizer tests passed.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "test failure: " << error.what() << '\n';
        return 1;
    }
}
