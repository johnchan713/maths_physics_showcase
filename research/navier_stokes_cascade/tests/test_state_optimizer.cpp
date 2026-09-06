#include "ns_cascade/state_optimizer.hpp"

#include <algorithm>
#include <cmath>
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
        ns_cascade::initialStateObjectiveGradient(system, initial, weights);
    const State final_gradient =
        ns_cascade::terminalStateObjectiveGradient(system, final, weights);
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

void testTrajectoryObjectiveGradient() {
    const ns_cascade::PseudospectralSystem system(16, 0.02, 5);
    const State initial = makeInitial(system);
    const State raw_direction = system.interactingWavePacketState(
        ns_cascade::WavePacketParameters(1.35, 2, 1.2, -0.7), 1.0);
    const State direction = ns_cascade::stateEnergySphereDirection(
        system, initial, raw_direction, 2);
    ns_cascade::StateObjectiveWeights weights;
    const int step_count = 6;
    const double time_step = 0.0004;
    std::vector<State> trajectory;
    State final = initial;
    for (int step = 0; step < step_count; ++step) {
        trajectory.push_back(final);
        system.stepRungeKutta4(final, time_step);
    }

    State reverse_gradient =
        ns_cascade::terminalStateObjectiveGradient(system, final, weights);
    for (int step = step_count; step-- > 0;) {
        reverse_gradient = system.adjointRungeKutta4Step(
            trajectory[static_cast<std::size_t>(step)],
            reverse_gradient,
            time_step);
    }
    reverse_gradient = ns_cascade::addOptimizationStates(
        reverse_gradient,
        ns_cascade::initialStateObjectiveGradient(system, initial, weights),
        1.0);
    const double analytical =
        ns_cascade::stateRealInnerProduct(reverse_gradient, direction);

    const double epsilon = 2e-4;
    const State plus_initial = ns_cascade::stateEnergySphereStep(
        system, initial, direction, epsilon);
    const State minus_initial = ns_cascade::stateEnergySphereStep(
        system, initial, direction, -epsilon);
    State plus = plus_initial;
    State minus = minus_initial;
    evolve(system, plus, step_count, time_step);
    evolve(system, minus, step_count, time_step);
    const double finite_difference =
        (ns_cascade::evaluateStateObjective(
             system, plus_initial, plus, weights).total -
         ns_cascade::evaluateStateObjective(
             system, minus_initial, minus, weights).total) /
        (2.0 * epsilon);
    expect(relativeScalarDifference(analytical, finite_difference) < 2e-6,
           "Adjoint state-objective gradient failed a geodesic difference");
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

}  // namespace

int main() {
    try {
        testBandLimitAndEnergySphere();
        testStaticObjectiveGradient();
        testTrajectoryObjectiveGradient();
        testLiftToFineGrid();
        std::cout << "All state-optimizer tests passed.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "test failure: " << error.what() << '\n';
        return 1;
    }
}
