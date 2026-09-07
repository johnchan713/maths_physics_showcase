#include "ns_cascade/optimization_state_csv.hpp"

#include <iostream>

// Snapshot-level gradient of gamma(u)=d(log ||u||_Hhalf)/dt. This validates
// the proposed objective's local adjoint source; it is not a new optimizer.
namespace {
using namespace ns_cascade;

OptimizationState halfLaplacian(const PseudospectralSystem& system,
                                 const OptimizationState& state) {
    OptimizationState result = state;
    const auto& modes = system.gridModes();
    for (std::size_t i = 0; i < result.size(); ++i) {
        result[i] = result[i] * std::sqrt(static_cast<double>(modes[i].normSquared()));
    }
    return result;
}

double criticalRate(const PseudospectralSystem& system, const OptimizationState& state) {
    const auto weighted = halfLaplacian(system, state);
    const double q = stateRealInnerProduct(state, weighted);
    if (!std::isfinite(q) || q <= 0.0) throw std::invalid_argument("Invalid critical norm");
    return stateRealInnerProduct(weighted, system.rightHandSide(state)) / q;
}

OptimizationState criticalRateGradient(const PseudospectralSystem& system,
                                       const OptimizationState& state) {
    const auto weighted = halfLaplacian(system, state);
    const auto rhs = system.rightHandSide(state);
    const double q = stateRealInnerProduct(state, weighted);
    const double numerator = stateRealInnerProduct(weighted, rhs);
    const auto first = halfLaplacian(system, rhs);
    const auto second = system.adjointTangentRightHandSide(state, weighted);
    auto gradient = scaleOptimizationState(addOptimizationStates(first, second, 1.0), 1.0 / q);
    return addOptimizationStates(gradient, weighted, -2.0 * numerator / (q * q));
}
}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) throw std::invalid_argument("Pass the frozen initial coefficient CSV");
        std::cout << "grid,direction,epsilon,adjoint_slope,finite_difference_slope,relative_error\n" << std::setprecision(17);
        for (int grid : {16, 32}) {
            const PseudospectralSystem system(grid, 0.02);
            const auto state = readOptimizationStateCsv(argv[1], system).state;
            const auto gradient = criticalRateGradient(system, state);
            for (auto initial : {InitialCondition::Deterministic, InitialCondition::TaylorGreen, InitialCondition::ABC}) {
                const auto direction = system.initialState(initial, 0.5);
                const double predicted = stateRealInnerProduct(gradient, direction);
                double best_error = 1.0;
                for (double epsilon : {1e-3, 1e-4, 1e-5}) {
                    const double measured = (criticalRate(system, addOptimizationStates(state, direction, epsilon)) -
                                             criticalRate(system, addOptimizationStates(state, direction, -epsilon))) / (2 * epsilon);
                    const double error = std::abs(predicted - measured) / std::max(1e-8, std::max(std::abs(predicted), std::abs(measured)));
                    best_error = std::min(best_error, error);
                    std::cout << grid << ',' << initialConditionName(initial) << ',' << epsilon << ','
                        << predicted << ',' << measured << ',' << error << '\n';
                }
                if (best_error > 1e-6) throw std::runtime_error("Critical-rate adjoint source failed its difference check");
            }
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
