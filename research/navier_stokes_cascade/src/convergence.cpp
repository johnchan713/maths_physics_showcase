#include "ns_cascade/galerkin.hpp"
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

enum class Backend {
    Direct,
    FFT
};

const char* backendName(Backend backend) {
    return backend == Backend::Direct ? "direct" : "fft";
}

struct Options {
    Backend backend = Backend::Direct;
    std::vector<int> cutoffs = std::vector<int>{2, 3};
    std::vector<int> grid_sizes = std::vector<int>{16, 32};
    double viscosity = 0.05;
    double time_step = 0.0005;
    int steps = 60;
    int diagnostic_every = 20;
    int sample_points = 0;
    double initial_energy = 1.0;
    ns_cascade::InitialCondition initial_condition =
        ns_cascade::InitialCondition::TaylorGreen;
    std::string output = "navier_stokes_convergence.csv";
};

struct Snapshot {
    int step;
    double time;
    ns_cascade::GalerkinSystem::Diagnostics diagnostics;
    std::vector<ns_cascade::GalerkinSystem::ShellDiagnostics> shells;
};

struct RunResult {
    std::vector<Snapshot> snapshots;
    double peak_cutoff_fraction;
};

template <typename T>
T parseNumber(const std::string& text, const std::string& flag) {
    std::istringstream stream(text);
    T value;
    char trailing;
    if (!(stream >> value) || (stream >> trailing)) {
        throw std::invalid_argument("Invalid value for " + flag + ": " + text);
    }
    return value;
}

std::string requireValue(int& index, int argc, char** argv) {
    if (index + 1 >= argc) {
        throw std::invalid_argument(std::string("Missing value after ") + argv[index]);
    }
    ++index;
    return argv[index];
}

ns_cascade::InitialCondition parseInitialCondition(const std::string& value) {
    if (value == "deterministic") return ns_cascade::InitialCondition::Deterministic;
    if (value == "taylor-green") return ns_cascade::InitialCondition::TaylorGreen;
    if (value == "abc") return ns_cascade::InitialCondition::ABC;
    if (value == "vortex-tubes") return ns_cascade::InitialCondition::VortexTubes;
    throw std::invalid_argument(
        "Unknown initial condition: " + value +
        " (expected deterministic, taylor-green, abc, or vortex-tubes)");
}

Backend parseBackend(const std::string& value) {
    if (value == "direct") return Backend::Direct;
    if (value == "fft") return Backend::FFT;
    throw std::invalid_argument(
        "Unknown backend: " + value + " (expected direct or fft)");
}

std::vector<int> parsePositiveList(const std::string& text,
                                   const std::string& flag) {
    std::vector<int> values;
    std::istringstream stream(text);
    std::string item;
    while (std::getline(stream, item, ',')) {
        if (item.empty()) throw std::invalid_argument("Empty value in " + flag);
        const int value = parseNumber<int>(item, flag);
        if (value < 1) {
            throw std::invalid_argument("Every value in " + flag + " must be >= 1");
        }
        values.push_back(value);
    }
    if (values.empty()) throw std::invalid_argument(flag + " cannot be empty");
    std::sort(values.begin(), values.end());
    values.erase(std::unique(values.begin(), values.end()), values.end());
    return values;
}

void printUsage(const char* program) {
    std::cout
        << "Usage: " << program << " [options]\n\n"
        << "Compare spatial resolutions and dt versus dt/2 at equal physical time.\n\n"
        << "Options:\n"
        << "  --backend B             direct or fft (default: direct)\n"
        << "  --cutoffs LIST          Comma-separated cutoffs (default: 2,3)\n"
        << "  --grids LIST            FFT grids (default: 16,32)\n"
        << "  --viscosity NU          Non-negative viscosity (default: 0.05)\n"
        << "  --dt DT                 Coarse RK4 step (default: 0.0005)\n"
        << "  --steps N               Coarse run steps (default: 60)\n"
        << "  --diagnostic-every N    Coarse sampling interval (default: 20)\n"
        << "  --sample-points N       Direct-backend samples; FFT uses native grids\n"
        << "  --energy E              Initial normalized energy (default: 1)\n"
        << "  --initial-condition C   deterministic, taylor-green, abc, or vortex-tubes\n"
        << "                            (default: taylor-green)\n"
        << "  --output PATH           Long-form comparison CSV path\n"
        << "  --help                  Show this message\n";
}

Options parseOptions(int argc, char** argv) {
    Options options;
    for (int i = 1; i < argc; ++i) {
        const std::string flag = argv[i];
        if (flag == "--help") {
            printUsage(argv[0]);
            std::exit(0);
        } else if (flag == "--backend") {
            options.backend = parseBackend(requireValue(i, argc, argv));
        } else if (flag == "--cutoffs") {
            options.cutoffs =
                parsePositiveList(requireValue(i, argc, argv), flag);
        } else if (flag == "--grids") {
            options.grid_sizes =
                parsePositiveList(requireValue(i, argc, argv), flag);
        } else if (flag == "--viscosity") {
            options.viscosity = parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--dt") {
            options.time_step = parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--steps") {
            options.steps = parseNumber<int>(requireValue(i, argc, argv), flag);
        } else if (flag == "--diagnostic-every") {
            options.diagnostic_every =
                parseNumber<int>(requireValue(i, argc, argv), flag);
        } else if (flag == "--sample-points") {
            options.sample_points = parseNumber<int>(requireValue(i, argc, argv), flag);
        } else if (flag == "--energy") {
            options.initial_energy =
                parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--initial-condition") {
            options.initial_condition =
                parseInitialCondition(requireValue(i, argc, argv));
        } else if (flag == "--output") {
            options.output = requireValue(i, argc, argv);
        } else {
            throw std::invalid_argument("Unknown option: " + flag);
        }
    }

    if (options.viscosity < 0.0) {
        throw std::invalid_argument("--viscosity must be non-negative");
    }
    if (options.time_step <= 0.0) throw std::invalid_argument("--dt must be positive");
    if (options.steps < 1) throw std::invalid_argument("--steps must be >= 1");
    if (options.diagnostic_every < 1) {
        throw std::invalid_argument("--diagnostic-every must be >= 1");
    }
    if (options.sample_points < 0) {
        throw std::invalid_argument("--sample-points cannot be negative");
    }
    if (options.backend == Backend::Direct && options.sample_points != 0 &&
        options.sample_points < 2 * options.cutoffs.back() + 1) {
        throw std::invalid_argument(
            "--sample-points must be 0 or >= 2*largest-cutoff+1");
    }
    if (options.backend == Backend::FFT && options.sample_points != 0) {
        throw std::invalid_argument(
            "FFT comparisons sample on each native grid; omit --sample-points");
    }
    if (options.initial_energy <= 0.0) {
        throw std::invalid_argument("--energy must be positive");
    }
    if (options.output.empty()) throw std::invalid_argument("--output cannot be empty");
    return options;
}

template <typename System>
RunResult runConfiguration(const System& system,
                           const Options& options,
                           double time_step,
                           int steps,
                           int diagnostic_every,
                           int diagnostic_sample_points) {
    typename System::State state =
        system.initialState(options.initial_condition, options.initial_energy);
    RunResult result;
    result.peak_cutoff_fraction = 0.0;

    for (int step = 0; step <= steps; ++step) {
        if (step % diagnostic_every == 0 || step == steps) {
            Snapshot snapshot;
            snapshot.step = step;
            snapshot.time = step * time_step;
            snapshot.diagnostics =
                system.diagnostics(state, diagnostic_sample_points);
            snapshot.shells = system.shellDiagnostics(state);
            result.peak_cutoff_fraction =
                std::max(result.peak_cutoff_fraction,
                         snapshot.diagnostics.high_shell_energy_fraction);
            result.snapshots.push_back(snapshot);
        }
        if (step != steps) system.stepRungeKutta4(state, time_step);
    }
    return result;
}

void writeHeader(std::ostream& output) {
    output
        << "initial_condition,backend,grid_size,cutoff,refinement,dt,"
        << "sample_points,step,time,shell,"
        << "lower_radius,upper_radius,total_energy,critical_l3_sample,"
        << "critical_l3_ratio,shell_energy,nonlinear_transfer,"
        << "viscous_dissipation,forward_flux,cutoff_shell_energy_fraction,"
        << "peak_cutoff_shell_fraction,cutoff_shell_ok,divergence_defect,"
        << "reality_defect\n";
}

void writeRun(std::ostream& output,
              const RunResult& result,
              ns_cascade::InitialCondition initial_condition,
              Backend backend,
              int grid_size,
              int cutoff,
              const std::string& refinement,
              double time_step,
              int sample_points) {
    const double initial_l3 =
        result.snapshots.front().diagnostics.critical_l3_sample;
    const bool cutoff_shell_ok = result.peak_cutoff_fraction <= 0.01;

    for (std::size_t snapshot_index = 0;
         snapshot_index < result.snapshots.size();
         ++snapshot_index) {
        const Snapshot& snapshot = result.snapshots[snapshot_index];
        const double l3_ratio = initial_l3 == 0.0
                                    ? 0.0
                                    : snapshot.diagnostics.critical_l3_sample /
                                          initial_l3;
        for (std::size_t shell_index = 0;
             shell_index < snapshot.shells.size();
             ++shell_index) {
            const ns_cascade::GalerkinSystem::ShellDiagnostics& shell =
                snapshot.shells[shell_index];
            output << ns_cascade::initialConditionName(initial_condition) << ','
                   << backendName(backend) << ',' << grid_size << ',' << cutoff
                   << ',' << refinement << ',' << time_step << ',' << sample_points
                   << ',' << snapshot.step << ',' << snapshot.time << ','
                   << shell.shell
                   << ',' << shell.lower_radius << ',' << shell.upper_radius << ','
                   << snapshot.diagnostics.energy << ','
                   << snapshot.diagnostics.critical_l3_sample << ',' << l3_ratio
                   << ',' << shell.energy << ',' << shell.nonlinear_transfer << ','
                   << shell.viscous_dissipation << ',' << shell.forward_flux << ','
                   << snapshot.diagnostics.high_shell_energy_fraction << ','
                   << result.peak_cutoff_fraction << ','
                   << (cutoff_shell_ok ? "true" : "false") << ','
                   << snapshot.diagnostics.divergence_defect << ','
                   << snapshot.diagnostics.reality_defect << '\n';
        }
    }
}

double relativeDifference(double left, double right) {
    const double scale = std::max(1e-15, std::max(std::abs(left), std::abs(right)));
    return std::abs(left - right) / scale;
}

double relativeCommonFluxDifference(const RunResult& left,
                                    const RunResult& right) {
    const std::vector<ns_cascade::GalerkinSystem::ShellDiagnostics>& left_shells =
        left.snapshots.back().shells;
    const std::vector<ns_cascade::GalerkinSystem::ShellDiagnostics>& right_shells =
        right.snapshots.back().shells;
    const std::size_t common_shells = std::min(left_shells.size(), right_shells.size());
    double largest_difference = 0.0;
    double largest_flux = 1e-15;
    for (std::size_t i = 0; i < common_shells; ++i) {
        largest_difference =
            std::max(largest_difference,
                     std::abs(left_shells[i].forward_flux -
                              right_shells[i].forward_flux));
        largest_flux = std::max(
            largest_flux,
            std::max(std::abs(left_shells[i].forward_flux),
                     std::abs(right_shells[i].forward_flux)));
    }
    return largest_difference / largest_flux;
}

template <typename System>
RunResult runResolution(const System& system,
                        const Options& options,
                        std::ostream& csv,
                        Backend backend,
                        int grid_size,
                        int cutoff,
                        int diagnostic_sample_points,
                        const std::string& resolution_label,
                        const RunResult* previous_refined,
                        const std::string& previous_label) {
    const RunResult coarse =
        runConfiguration(system,
                         options,
                         options.time_step,
                         options.steps,
                         options.diagnostic_every,
                         diagnostic_sample_points);
    const RunResult refined =
        runConfiguration(system,
                         options,
                         options.time_step / 2.0,
                         options.steps * 2,
                         options.diagnostic_every * 2,
                         diagnostic_sample_points);
    writeRun(csv,
             coarse,
             options.initial_condition,
             backend,
             grid_size,
             cutoff,
             "dt",
             options.time_step,
             diagnostic_sample_points);
    writeRun(csv,
             refined,
             options.initial_condition,
             backend,
             grid_size,
             cutoff,
             "dt/2",
             options.time_step / 2.0,
             diagnostic_sample_points);

    const ns_cascade::GalerkinSystem::Diagnostics& coarse_final =
        coarse.snapshots.back().diagnostics;
    const ns_cascade::GalerkinSystem::Diagnostics& refined_final =
        refined.snapshots.back().diagnostics;
    std::cout << std::setprecision(6)
              << resolution_label
              << ": relative final-energy dt error = "
              << relativeDifference(coarse_final.energy, refined_final.energy)
              << ", relative final-L3 dt error = "
              << relativeDifference(coarse_final.critical_l3_sample,
                                    refined_final.critical_l3_sample)
              << ", peak cutoff fraction = "
              << std::max(coarse.peak_cutoff_fraction,
                          refined.peak_cutoff_fraction)
              << '\n';

    if (previous_refined != NULL) {
        const ns_cascade::GalerkinSystem::Diagnostics& previous_final =
            previous_refined->snapshots.back().diagnostics;
        std::cout
            << previous_label << " -> " << resolution_label
            << ": relative final-energy resolution difference = "
            << relativeDifference(previous_final.energy, refined_final.energy)
            << ", relative final-L3 resolution difference = "
            << relativeDifference(previous_final.critical_l3_sample,
                                  refined_final.critical_l3_sample)
            << ", relative common-shell flux difference = "
            << relativeCommonFluxDifference(*previous_refined, refined) << '\n';
    }
    return refined;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parseOptions(argc, argv);
        std::ofstream csv(options.output.c_str());
        if (!csv) throw std::runtime_error("Cannot open output file: " + options.output);
        csv << std::setprecision(17);
        writeHeader(csv);

        std::cout << "Initial condition: "
                  << ns_cascade::initialConditionName(options.initial_condition)
                  << "\nBackend: " << backendName(options.backend) << '\n';
        bool have_previous_resolution = false;
        std::string previous_label;
        RunResult previous_refined;

        if (options.backend == Backend::Direct) {
            const int common_sample_points =
                options.sample_points == 0
                    ? std::max(4 * options.cutoffs.back() + 1, 9)
                    : options.sample_points;
            for (std::size_t i = 0; i < options.cutoffs.size(); ++i) {
                const int cutoff = options.cutoffs[i];
                const ns_cascade::GalerkinSystem system(cutoff, options.viscosity);
                const std::string label = "K=" + std::to_string(cutoff);
                const RunResult refined = runResolution(
                    system,
                    options,
                    csv,
                    Backend::Direct,
                    0,
                    cutoff,
                    common_sample_points,
                    label,
                    have_previous_resolution ? &previous_refined : NULL,
                    previous_label);
                previous_refined = refined;
                previous_label = label;
                have_previous_resolution = true;
            }
        } else {
            for (std::size_t i = 0; i < options.grid_sizes.size(); ++i) {
                const int grid_size = options.grid_sizes[i];
                const ns_cascade::PseudospectralSystem system(
                    grid_size, options.viscosity);
                const std::string label =
                    "N=" + std::to_string(grid_size) +
                    ", K=" + std::to_string(system.cutoff());
                const RunResult refined = runResolution(
                    system,
                    options,
                    csv,
                    Backend::FFT,
                    grid_size,
                    system.cutoff(),
                    grid_size,
                    label,
                    have_previous_resolution ? &previous_refined : NULL,
                    previous_label);
                previous_refined = refined;
                previous_label = label;
                have_previous_resolution = true;
            }
        }

        std::cout << "Comparison written to " << options.output << '\n'
                  << "A converged finite cascade is evidence, not a PDE proof.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
