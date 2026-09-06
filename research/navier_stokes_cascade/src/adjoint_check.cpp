#include "ns_cascade/pseudospectral.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct Options {
    int grid_size = 16;
    int cutoff = 0;
    double viscosity = 0.02;
    double time_step = 0.0005;
    double final_time = 0.01;
    double initial_energy = 10.0;
    std::vector<double> epsilons = {4e-3, 2e-3, 1e-3};
    double duality_tolerance = 2e-11;
    double error_tolerance = 2e-7;
    double order_tolerance = 1.8;
    std::string output = "navier_stokes_adjoint_check.csv";
};

template <typename T>
T parseNumber(const std::string& text, const std::string& flag) {
    std::istringstream stream(text);
    T value = T();
    char trailing = '\0';
    if (!(stream >> value) || (stream >> trailing)) {
        throw std::invalid_argument("Invalid value for " + flag + ": " + text);
    }
    return value;
}

std::string requireValue(int& index, int argc, char** argv) {
    if (index + 1 >= argc) {
        throw std::invalid_argument(std::string("Missing value after ") +
                                    argv[index]);
    }
    return argv[++index];
}

std::vector<double> parseList(const std::string& text,
                              const std::string& flag) {
    std::vector<double> values;
    std::istringstream stream(text);
    std::string item;
    while (std::getline(stream, item, ',')) {
        values.push_back(parseNumber<double>(item, flag));
    }
    if (values.empty()) throw std::invalid_argument(flag + " cannot be empty");
    return values;
}

void printUsage(const char* program) {
    std::cout
        << "Usage: " << program << " [options]\n\n"
        << "Check the analytical RHS adjoint and the exact reverse pass through "
           "fixed-step RK4. The terminal objective is one half of the squared "
           "critical H1/2 norm.\n\n"
        << "  --grid N              FFT grid (default: 16)\n"
        << "  --cutoff K            Retained cutoff; zero chooses safe maximum\n"
        << "  --viscosity NU        Viscosity (default: 0.02)\n"
        << "  --dt DT               Fixed RK4 step (default: 0.0005)\n"
        << "  --final-time T        Comparison horizon (default: 0.01)\n"
        << "  --energy E            Base-state energy (default: 10)\n"
        << "  --epsilons LIST       Decreasing perturbations (0.004,0.002,0.001)\n"
        << "  --duality-tolerance X RHS/trajectory identity gate\n"
        << "  --error-tolerance X   Final objective-gradient gate\n"
        << "  --order-tolerance P   Minimum final observed order\n"
        << "  --output PATH         Convergence CSV\n"
        << "  --help                Show this message\n";
}

Options parseOptions(int argc, char** argv) {
    Options options;
    for (int i = 1; i < argc; ++i) {
        const std::string flag = argv[i];
        if (flag == "--help") {
            printUsage(argv[0]);
            std::exit(0);
        } else if (flag == "--grid") {
            options.grid_size = parseNumber<int>(requireValue(i, argc, argv), flag);
        } else if (flag == "--cutoff") {
            options.cutoff = parseNumber<int>(requireValue(i, argc, argv), flag);
        } else if (flag == "--viscosity") {
            options.viscosity =
                parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--dt") {
            options.time_step =
                parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--final-time") {
            options.final_time =
                parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--energy") {
            options.initial_energy =
                parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--epsilons") {
            options.epsilons = parseList(requireValue(i, argc, argv), flag);
        } else if (flag == "--duality-tolerance") {
            options.duality_tolerance =
                parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--error-tolerance") {
            options.error_tolerance =
                parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--order-tolerance") {
            options.order_tolerance =
                parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--output") {
            options.output = requireValue(i, argc, argv);
        } else {
            throw std::invalid_argument("Unknown option: " + flag);
        }
    }
    if (options.grid_size < 8 || options.cutoff < 0 || options.output.empty() ||
        !std::isfinite(options.viscosity) || options.viscosity < 0.0 ||
        options.epsilons.size() < 2) {
        throw std::invalid_argument(
            "Grid, cutoff, viscosity, epsilon count, or output is invalid");
    }
    const double positive[] = {
        options.time_step,
        options.final_time,
        options.initial_energy,
        options.duality_tolerance,
        options.error_tolerance,
        options.order_tolerance};
    for (std::size_t i = 0; i < sizeof(positive) / sizeof(positive[0]); ++i) {
        if (!std::isfinite(positive[i]) || positive[i] <= 0.0) {
            throw std::invalid_argument("Positive adjoint-check option is invalid");
        }
    }
    for (std::size_t i = 0; i < options.epsilons.size(); ++i) {
        if (!std::isfinite(options.epsilons[i]) || options.epsilons[i] <= 0.0 ||
            (i > 0 && options.epsilons[i] >= options.epsilons[i - 1])) {
            throw std::invalid_argument(
                "Epsilons must be finite, positive, and strictly decreasing");
        }
    }
    return options;
}

using State = ns_cascade::PseudospectralSystem::State;

State addScaled(const State& state, const State& direction, double scale) {
    State result = state;
    for (std::size_t i = 0; i < result.size(); ++i) {
        result[i] += direction[i] * scale;
    }
    return result;
}

double realInnerProduct(const State& left, const State& right) {
    double value = 0.0;
    for (std::size_t i = 0; i < left.size(); ++i) {
        value += std::real(ns_cascade::innerProduct(left[i], right[i]));
    }
    return value;
}

double stateNorm(const State& state) {
    return std::sqrt(std::max(0.0, realInnerProduct(state, state)));
}

double normalizedPairingError(double left,
                              double right,
                              double left_scale,
                              double right_scale) {
    return std::abs(left - right) /
           std::max(1e-300, std::max(left_scale, right_scale));
}

double relativeStateDifference(const State& left, const State& right) {
    State difference = left;
    for (std::size_t i = 0; i < difference.size(); ++i) {
        difference[i] = left[i] - right[i];
    }
    return stateNorm(difference) / std::max(stateNorm(right), 1e-300);
}

double objective(const ns_cascade::PseudospectralSystem& system,
                 const State& state) {
    const double critical_norm = system.criticalHOneHalf(state);
    return 0.5 * critical_norm * critical_norm;
}

State objectiveGradient(const ns_cascade::PseudospectralSystem& system,
                        const State& state) {
    State gradient = system.zeroState();
    const std::vector<ns_cascade::WaveVector>& modes = system.gridModes();
    for (std::size_t i = 0; i < gradient.size(); ++i) {
        const double wave_magnitude =
            std::sqrt(static_cast<double>(modes[i].normSquared()));
        gradient[i] = state[i] * wave_magnitude;
    }
    return gradient;
}

void evolve(const ns_cascade::PseudospectralSystem& system,
            State& state,
            double time_step,
            double final_time) {
    double time = 0.0;
    while (time < final_time) {
        const double dt = std::min(time_step, final_time - time);
        system.stepRungeKutta4(state, dt);
        time += dt;
    }
}

int run(const Options& options) {
    const ns_cascade::PseudospectralSystem system(
        options.grid_size, options.viscosity, options.cutoff);
    const ns_cascade::WavePacketParameters base_parameters(
        1.0877734104170571, 1, 0.7577474832313891, 0.299291384334349);
    State initial = system.interactingWavePacketState(
        base_parameters, options.initial_energy);
    State direction = system.interactingWavePacketState(
        ns_cascade::WavePacketParameters(1.3, 1, 1.1, -0.4), 1.0);

    // Keep perturbations tangent to the fixed-energy sphere used by the packet
    // optimizer, then normalize them in the Fourier-state inner product.
    const double radial_coefficient =
        realInnerProduct(initial, direction) /
        realInnerProduct(initial, initial);
    direction = addScaled(direction, initial, -radial_coefficient);
    const double direction_norm = stateNorm(direction);
    if (!std::isfinite(direction_norm) || direction_norm <= 1e-14) {
        throw std::runtime_error("Energy-tangent direction is degenerate");
    }
    for (std::size_t i = 0; i < direction.size(); ++i) {
        direction[i] = direction[i] * (1.0 / direction_norm);
    }

    std::vector<State> trajectory;
    std::vector<double> time_steps;
    State base = initial;
    State tangent = direction;
    State ordinary_base = initial;
    double time = 0.0;
    while (time < options.final_time) {
        const double dt = std::min(options.time_step, options.final_time - time);
        trajectory.push_back(base);
        time_steps.push_back(dt);
        system.stepTangentRungeKutta4(base, tangent, dt);
        system.stepRungeKutta4(ordinary_base, dt);
        time += dt;
    }
    const double base_difference = relativeStateDifference(base, ordinary_base);
    const State final_dual = objectiveGradient(system, base);

    const State tangent_rhs = system.tangentRightHandSide(initial, direction);
    const State adjoint_rhs =
        system.adjointTangentRightHandSide(initial, final_dual);
    const double rhs_forward_pairing =
        realInnerProduct(tangent_rhs, final_dual);
    const double rhs_reverse_pairing =
        realInnerProduct(direction, adjoint_rhs);
    const double rhs_duality_error = normalizedPairingError(
        rhs_forward_pairing,
        rhs_reverse_pairing,
        stateNorm(tangent_rhs) * stateNorm(final_dual),
        stateNorm(direction) * stateNorm(adjoint_rhs));

    State initial_dual = final_dual;
    for (std::size_t step = trajectory.size(); step-- > 0;) {
        initial_dual = system.adjointRungeKutta4Step(
            trajectory[step], initial_dual, time_steps[step]);
    }
    const double tangent_derivative = realInnerProduct(tangent, final_dual);
    const double adjoint_derivative = realInnerProduct(direction, initial_dual);
    const double trajectory_duality_error = normalizedPairingError(
        tangent_derivative,
        adjoint_derivative,
        stateNorm(tangent) * stateNorm(final_dual),
        stateNorm(direction) * stateNorm(initial_dual));

    std::ofstream csv(options.output.c_str());
    if (!csv) throw std::runtime_error("Could not open adjoint-check CSV");
    csv << std::setprecision(17)
        << "epsilon,finite_difference_derivative,adjoint_derivative,"
           "relative_gradient_error,observed_order,rhs_duality_error,"
           "trajectory_duality_error,base_state_difference,"
           "adjoint_divergence_defect,adjoint_reality_defect\n";

    double previous_epsilon = 0.0;
    double previous_error = 0.0;
    double final_error = 0.0;
    double final_order = 0.0;
    for (std::size_t sample = 0; sample < options.epsilons.size(); ++sample) {
        const double epsilon = options.epsilons[sample];
        State plus = addScaled(initial, direction, epsilon);
        State minus = addScaled(initial, direction, -epsilon);
        evolve(system, plus, options.time_step, options.final_time);
        evolve(system, minus, options.time_step, options.final_time);
        const double finite_difference_derivative =
            (objective(system, plus) - objective(system, minus)) /
            (2.0 * epsilon);
        const double error =
            std::abs(finite_difference_derivative - adjoint_derivative) /
            std::max(std::abs(adjoint_derivative), 1e-300);
        const double order = sample == 0
                                 ? 0.0
                                 : std::log(previous_error / error) /
                                       std::log(previous_epsilon / epsilon);
        csv << epsilon << ',' << finite_difference_derivative << ','
            << adjoint_derivative << ',' << error << ',';
        if (sample == 0) csv << ',';
        else csv << order << ',';
        csv << rhs_duality_error << ',' << trajectory_duality_error << ','
            << base_difference << ',' << system.divergenceDefect(initial_dual)
            << ',' << system.realityDefect(initial_dual) << '\n';
        previous_epsilon = epsilon;
        previous_error = error;
        final_error = error;
        final_order = order;
    }
    if (!csv) throw std::runtime_error("Failed while writing adjoint-check CSV");

    const double adjoint_divergence = system.divergenceDefect(initial_dual);
    const double adjoint_reality = system.realityDefect(initial_dual);
    const bool passed =
        rhs_duality_error <= options.duality_tolerance &&
        trajectory_duality_error <= options.duality_tolerance &&
        base_difference <= 1e-14 && final_error <= options.error_tolerance &&
        final_order >= options.order_tolerance &&
        adjoint_divergence <= 1e-10 && adjoint_reality <= 1e-10;
    std::cout << std::setprecision(12)
              << "Reverse discrete-adjoint check\n"
              << "  grid/cutoff/time: " << options.grid_size << '/'
              << system.cutoff() << '/' << options.final_time << '\n'
              << "  RHS/trajectory duality errors: " << rhs_duality_error
              << '/' << trajectory_duality_error << '\n'
              << "  tangent/adjoint derivative: " << tangent_derivative << '/'
              << adjoint_derivative << '\n'
              << "  final objective error/order: " << final_error << '/'
              << final_order << '\n'
              << "  base-state difference: " << base_difference << '\n'
              << "  adjoint divergence/reality: " << adjoint_divergence << '/'
              << adjoint_reality << '\n'
              << "  verdict: " << (passed ? "pass" : "fail") << '\n';
    return passed ? 0 : 2;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        return run(parseOptions(argc, argv));
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
