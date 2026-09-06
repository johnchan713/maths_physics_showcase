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
    std::vector<double> epsilons = {1e-3, 5e-4, 2.5e-4};
    double error_tolerance = 2e-5;
    double order_tolerance = 1.8;
    std::string output = "navier_stokes_tangent_check.csv";
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
        << "Check the tangent-linear RK4 trajectory against centered finite "
           "differences. Fixed timesteps deliberately exclude derivatives of "
           "adaptive timestep selection.\n\n"
        << "  --grid N              FFT grid (default: 16)\n"
        << "  --cutoff K            Retained cutoff; zero chooses safe maximum\n"
        << "  --viscosity NU        Viscosity (default: 0.02)\n"
        << "  --dt DT               Fixed RK4 step (default: 0.0005)\n"
        << "  --final-time T        Comparison horizon (default: 0.01)\n"
        << "  --energy E            Base-state energy (default: 10)\n"
        << "  --epsilons LIST       Decreasing perturbations\n"
        << "  --error-tolerance X   Final relative-error gate (default: 2e-5)\n"
        << "  --order-tolerance P   Minimum final observed order (default: 1.8)\n"
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
        !std::isfinite(options.viscosity) || options.viscosity < 0.0) {
        throw std::invalid_argument("Grid, cutoff, viscosity, or output is invalid");
    }
    const double positive[] = {
        options.time_step,
        options.final_time,
        options.initial_energy,
        options.error_tolerance,
        options.order_tolerance};
    for (std::size_t i = 0; i < sizeof(positive) / sizeof(positive[0]); ++i) {
        if (!std::isfinite(positive[i]) || positive[i] <= 0.0) {
            throw std::invalid_argument("Positive tangent-check option is invalid");
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

double relativeDifference(const State& left, const State& right) {
    double difference_squared = 0.0;
    double reference_squared = 0.0;
    for (std::size_t i = 0; i < left.size(); ++i) {
        difference_squared += ns_cascade::normSquared(left[i] - right[i]);
        reference_squared += ns_cascade::normSquared(right[i]);
    }
    return std::sqrt(difference_squared / std::max(reference_squared, 1e-300));
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

    // Remove the radial energy component: Re<u,v>=0 makes v tangent to the
    // fixed-initial-energy sphere used by the optimizer.
    const double radial_coefficient =
        realInnerProduct(initial, direction) /
        realInnerProduct(initial, initial);
    direction = addScaled(direction, initial, -radial_coefficient);
    const double direction_norm = std::sqrt(realInnerProduct(direction, direction));
    for (std::size_t i = 0; i < direction.size(); ++i) {
        direction[i] = direction[i] * (1.0 / direction_norm);
    }
    if (std::abs(realInnerProduct(initial, direction)) > 1e-12) {
        throw std::runtime_error("Energy-tangent projection failed");
    }

    State base = initial;
    State tangent = direction;
    State ordinary_base = initial;
    double time = 0.0;
    while (time < options.final_time) {
        const double dt = std::min(options.time_step, options.final_time - time);
        system.stepTangentRungeKutta4(base, tangent, dt);
        system.stepRungeKutta4(ordinary_base, dt);
        time += dt;
    }
    const double base_difference = relativeDifference(base, ordinary_base);

    std::ofstream csv(options.output.c_str());
    if (!csv) throw std::runtime_error("Could not open tangent-check CSV");
    csv << std::setprecision(17)
        << "epsilon,relative_tangent_error,observed_order,base_state_difference,"
           "tangent_divergence_defect,tangent_reality_defect\n";

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
        State difference = plus;
        for (std::size_t i = 0; i < difference.size(); ++i) {
            difference[i] = (plus[i] - minus[i]) * (0.5 / epsilon);
        }
        const double error = relativeDifference(difference, tangent);
        const double order = sample == 0
                                 ? 0.0
                                 : std::log(previous_error / error) /
                                       std::log(previous_epsilon / epsilon);
        csv << epsilon << ',' << error << ',';
        if (sample == 0) csv << ',';
        else csv << order << ',';
        csv << base_difference << ',' << system.divergenceDefect(tangent) << ','
            << system.realityDefect(tangent) << '\n';
        previous_epsilon = epsilon;
        previous_error = error;
        final_error = error;
        final_order = order;
    }
    if (!csv) throw std::runtime_error("Failed while writing tangent-check CSV");

    const bool passed =
        base_difference <= 1e-14 && final_error <= options.error_tolerance &&
        options.epsilons.size() > 1 && final_order >= options.order_tolerance &&
        system.divergenceDefect(tangent) <= 1e-10 &&
        system.realityDefect(tangent) <= 1e-10;
    std::cout << std::setprecision(12)
              << "Tangent-linear trajectory check\n"
              << "  grid/cutoff/time: " << options.grid_size << '/'
              << system.cutoff() << '/' << options.final_time << '\n'
              << "  final relative error/order: " << final_error << '/'
              << final_order << '\n'
              << "  base-state difference: " << base_difference << '\n'
              << "  tangent divergence/reality: "
              << system.divergenceDefect(tangent) << '/'
              << system.realityDefect(tangent) << '\n'
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
