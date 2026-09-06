#include "ns_cascade/pseudospectral.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct Options {
    int coarse_grid = 16;
    int fine_grid = 32;
    double viscosity = 0.02;
    double maximum_time_step = 0.005;
    double final_time = 0.02;
    double target_cfl = 0.35;
    double diffusion_safety = 2.0;
    int diagnostic_every = 5;
    int top_count = 3;
    double initial_energy = 1.0;
    double cutoff_fraction_threshold = 0.01;
    std::vector<double> core_radii = {0.55, 0.70};
    std::vector<double> separations = {1.2, 1.8};
    std::vector<double> bend_amplitudes = {0.0, 0.30};
    std::vector<int> axial_wavenumbers = {1, 2};
    std::string output = "navier_stokes_candidate_search.csv";
};

struct Candidate {
    int id;
    ns_cascade::VortexTubeParameters parameters;
};

struct RunResult {
    Candidate candidate;
    int grid_size;
    int cutoff;
    bool completed;
    bool resolution_ok;
    std::string failure;
    int steps;
    double minimum_time_step;
    double maximum_time_step;
    double maximum_cfl_bound;
    double maximum_viscous_number;
    double initial_three_dimensional_energy_fraction;
    double initial_energy;
    double final_energy;
    double initial_l3;
    double final_l3;
    double peak_l3;
    double initial_vorticity;
    double final_vorticity;
    double peak_vorticity;
    double initial_enstrophy;
    double final_enstrophy;
    double peak_enstrophy;
    double peak_palinstrophy;
    double maximum_positive_forward_flux;
    double peak_cutoff_fraction;
    double maximum_divergence_defect;
    double maximum_reality_defect;
    double maximum_energy_balance_residual;
    double score;
    double elapsed_seconds;

    RunResult(const Candidate& candidate_value, int grid, int cutoff_value)
        : candidate(candidate_value),
          grid_size(grid),
          cutoff(cutoff_value),
          completed(false),
          resolution_ok(false),
          failure(),
          steps(0),
          minimum_time_step(std::numeric_limits<double>::quiet_NaN()),
          maximum_time_step(std::numeric_limits<double>::quiet_NaN()),
          maximum_cfl_bound(std::numeric_limits<double>::quiet_NaN()),
          maximum_viscous_number(std::numeric_limits<double>::quiet_NaN()),
          initial_three_dimensional_energy_fraction(
              std::numeric_limits<double>::quiet_NaN()),
          initial_energy(std::numeric_limits<double>::quiet_NaN()),
          final_energy(std::numeric_limits<double>::quiet_NaN()),
          initial_l3(std::numeric_limits<double>::quiet_NaN()),
          final_l3(std::numeric_limits<double>::quiet_NaN()),
          peak_l3(std::numeric_limits<double>::quiet_NaN()),
          initial_vorticity(std::numeric_limits<double>::quiet_NaN()),
          final_vorticity(std::numeric_limits<double>::quiet_NaN()),
          peak_vorticity(std::numeric_limits<double>::quiet_NaN()),
          initial_enstrophy(std::numeric_limits<double>::quiet_NaN()),
          final_enstrophy(std::numeric_limits<double>::quiet_NaN()),
          peak_enstrophy(std::numeric_limits<double>::quiet_NaN()),
          peak_palinstrophy(std::numeric_limits<double>::quiet_NaN()),
          maximum_positive_forward_flux(
              std::numeric_limits<double>::quiet_NaN()),
          peak_cutoff_fraction(std::numeric_limits<double>::quiet_NaN()),
          maximum_divergence_defect(
              std::numeric_limits<double>::quiet_NaN()),
          maximum_reality_defect(std::numeric_limits<double>::quiet_NaN()),
          maximum_energy_balance_residual(
              std::numeric_limits<double>::quiet_NaN()),
          score(-std::numeric_limits<double>::infinity()),
          elapsed_seconds(std::numeric_limits<double>::quiet_NaN()) {}
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

template <typename T>
std::vector<T> parseList(const std::string& text, const std::string& flag) {
    std::vector<T> values;
    std::istringstream stream(text);
    std::string item;
    while (std::getline(stream, item, ',')) {
        if (item.empty()) {
            throw std::invalid_argument("Empty entry in " + flag);
        }
        values.push_back(parseNumber<T>(item, flag));
    }
    if (values.empty()) throw std::invalid_argument(flag + " cannot be empty");
    return values;
}

template <typename T>
void removeDuplicates(std::vector<T>& values) {
    std::vector<T> unique_values;
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (std::find(unique_values.begin(), unique_values.end(), values[i]) ==
            unique_values.end()) {
            unique_values.push_back(values[i]);
        }
    }
    values.swap(unique_values);
}

void printUsage(const char* program) {
    std::cout
        << "Usage: " << program << " [options]\n\n"
        << "Coarse-to-fine search over smooth periodic vortex-tube pairs.\n"
        << "This ranks finite numerical signals; it cannot prove PDE blow-up.\n\n"
        << "Options:\n"
        << "  --coarse-grid N         Coarse power-of-two grid (default: 16)\n"
        << "  --fine-grid N           Finalist grid, larger than coarse (default: 32)\n"
        << "  --viscosity NU          Positive viscosity (default: 0.02)\n"
        << "  --dt DT                 Maximum adaptive step (default: 0.005)\n"
        << "  --final-time T          Common physical end time (default: 0.02)\n"
        << "  --cfl C                 Conservative CFL target (default: 0.35)\n"
        << "  --diffusion-safety S    Bound on nu*|k|max^2*dt (default: 2)\n"
        << "  --diagnostic-every N    Accepted steps between samples (default: 5)\n"
        << "  --top N                 Coarse finalists rerun fine (default: 3)\n"
        << "  --energy E              Normalized initial energy (default: 1)\n"
        << "  --cores LIST            Core radii (default: 0.55,0.70)\n"
        << "  --separations LIST      Pair separations (default: 1.2,1.8)\n"
        << "  --bends LIST            Helical bend amplitudes (default: 0,0.30)\n"
        << "  --axial-modes LIST      Positive axial modes (default: 1,2)\n"
        << "  --cutoff-threshold F    Max cutoff-shell energy fraction (default: 0.01)\n"
        << "  --output PATH           Summary CSV path\n"
        << "  --help                  Show this message\n";
}

bool finiteAndPositive(double value) {
    return std::isfinite(value) && value > 0.0;
}

Options parseOptions(int argc, char** argv) {
    Options options;
    for (int i = 1; i < argc; ++i) {
        const std::string flag = argv[i];
        if (flag == "--help") {
            printUsage(argv[0]);
            std::exit(0);
        } else if (flag == "--coarse-grid") {
            options.coarse_grid =
                parseNumber<int>(requireValue(i, argc, argv), flag);
        } else if (flag == "--fine-grid") {
            options.fine_grid =
                parseNumber<int>(requireValue(i, argc, argv), flag);
        } else if (flag == "--viscosity") {
            options.viscosity =
                parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--dt") {
            options.maximum_time_step =
                parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--final-time") {
            options.final_time =
                parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--cfl") {
            options.target_cfl =
                parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--diffusion-safety") {
            options.diffusion_safety =
                parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--diagnostic-every") {
            options.diagnostic_every =
                parseNumber<int>(requireValue(i, argc, argv), flag);
        } else if (flag == "--top") {
            options.top_count =
                parseNumber<int>(requireValue(i, argc, argv), flag);
        } else if (flag == "--energy") {
            options.initial_energy =
                parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--cores") {
            options.core_radii =
                parseList<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--separations") {
            options.separations =
                parseList<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--bends") {
            options.bend_amplitudes =
                parseList<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--axial-modes") {
            options.axial_wavenumbers =
                parseList<int>(requireValue(i, argc, argv), flag);
        } else if (flag == "--cutoff-threshold") {
            options.cutoff_fraction_threshold =
                parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--output") {
            options.output = requireValue(i, argc, argv);
        } else {
            throw std::invalid_argument("Unknown option: " + flag);
        }
    }

    removeDuplicates(options.core_radii);
    removeDuplicates(options.separations);
    removeDuplicates(options.bend_amplitudes);
    removeDuplicates(options.axial_wavenumbers);
    if (options.coarse_grid >= options.fine_grid) {
        throw std::invalid_argument("--fine-grid must be larger than --coarse-grid");
    }
    if (!finiteAndPositive(options.viscosity)) {
        throw std::invalid_argument("--viscosity must be finite and positive");
    }
    if (!finiteAndPositive(options.maximum_time_step)) {
        throw std::invalid_argument("--dt must be finite and positive");
    }
    if (!finiteAndPositive(options.final_time)) {
        throw std::invalid_argument("--final-time must be finite and positive");
    }
    if (!finiteAndPositive(options.target_cfl) ||
        !finiteAndPositive(options.diffusion_safety)) {
        throw std::invalid_argument(
            "--cfl and --diffusion-safety must be finite and positive");
    }
    if (!finiteAndPositive(options.initial_energy)) {
        throw std::invalid_argument("--energy must be finite and positive");
    }
    if (!finiteAndPositive(options.cutoff_fraction_threshold) ||
        options.cutoff_fraction_threshold >= 1.0) {
        throw std::invalid_argument(
            "--cutoff-threshold must be finite and in (0,1)");
    }
    if (options.diagnostic_every < 1 || options.top_count < 1) {
        throw std::invalid_argument(
            "--diagnostic-every and --top must be positive integers");
    }
    if (options.output.empty()) throw std::invalid_argument("--output cannot be empty");

    const int coarse_safe_cutoff = (options.coarse_grid - 1) / 3;
    for (std::size_t i = 0; i < options.axial_wavenumbers.size(); ++i) {
        if (options.axial_wavenumbers[i] < 1 ||
            options.axial_wavenumbers[i] > coarse_safe_cutoff) {
            throw std::invalid_argument(
                "Every axial mode must fit inside the coarse safe cutoff");
        }
    }
    const double pi = 3.1415926535897932384626433832795;
    for (std::size_t i = 0; i < options.core_radii.size(); ++i) {
        if (!std::isfinite(options.core_radii[i]) ||
            options.core_radii[i] <= 0.0 || options.core_radii[i] > pi) {
            throw std::invalid_argument("Every core radius must be in (0,pi]");
        }
    }
    for (std::size_t i = 0; i < options.separations.size(); ++i) {
        if (!std::isfinite(options.separations[i]) ||
            options.separations[i] <= 0.0 || options.separations[i] >= 2.0 * pi) {
            throw std::invalid_argument("Every separation must be in (0,2*pi)");
        }
    }
    for (std::size_t i = 0; i < options.bend_amplitudes.size(); ++i) {
        if (!std::isfinite(options.bend_amplitudes[i]) ||
            options.bend_amplitudes[i] < 0.0 ||
            options.bend_amplitudes[i] >= pi) {
            throw std::invalid_argument("Every bend amplitude must be in [0,pi)");
        }
    }
    return options;
}

std::vector<Candidate> buildCandidates(const Options& options) {
    std::vector<Candidate> candidates;
    int id = 1;
    for (std::size_t core = 0; core < options.core_radii.size(); ++core) {
        for (std::size_t separation = 0;
             separation < options.separations.size();
             ++separation) {
            for (std::size_t bend = 0;
                 bend < options.bend_amplitudes.size();
                 ++bend) {
                for (std::size_t mode = 0;
                     mode < options.axial_wavenumbers.size();
                     ++mode) {
                    // Axial mode has no effect for a straight, z-invariant pair.
                    if (options.bend_amplitudes[bend] == 0.0 && mode != 0) continue;
                    Candidate candidate;
                    candidate.id = id++;
                    candidate.parameters = ns_cascade::VortexTubeParameters(
                        options.core_radii[core],
                        options.separations[separation],
                        options.bend_amplitudes[bend],
                        options.axial_wavenumbers[mode]);
                    candidates.push_back(candidate);
                }
            }
        }
    }
    if (candidates.empty()) {
        throw std::invalid_argument("The parameter lists produced no candidates");
    }
    return candidates;
}

double threeDimensionalEnergyFraction(
    const ns_cascade::PseudospectralSystem& system,
    const ns_cascade::PseudospectralSystem::State& state) {
    double energy = 0.0;
    for (std::size_t i = 0; i < state.size(); ++i) {
        if (system.gridModes()[i].z != 0) {
            energy += 0.5 * ns_cascade::normSquared(state[i]);
        }
    }
    return energy / system.energy(state);
}

double positiveForwardFlux(
    const std::vector<ns_cascade::PseudospectralSystem::ShellDiagnostics>& shells) {
    double maximum_flux = 0.0;
    for (std::size_t i = 0; i < shells.size(); ++i) {
        maximum_flux = std::max(maximum_flux, shells[i].forward_flux);
    }
    return maximum_flux;
}

bool finiteDiagnostics(
    const ns_cascade::PseudospectralSystem::Diagnostics& values) {
    return std::isfinite(values.energy) && std::isfinite(values.enstrophy) &&
           std::isfinite(values.palinstrophy) &&
           std::isfinite(values.critical_l3_sample) &&
           std::isfinite(values.sampled_vorticity_max) &&
           std::isfinite(values.high_shell_energy_fraction) &&
           std::isfinite(values.divergence_defect) &&
           std::isfinite(values.reality_defect) &&
           std::isfinite(values.energy_balance_residual);
}

void recordDiagnostics(
    const ns_cascade::PseudospectralSystem& system,
    const ns_cascade::PseudospectralSystem::State& state,
    const ns_cascade::PseudospectralSystem::Diagnostics& values,
    RunResult& result) {
    if (!finiteDiagnostics(values)) {
        throw std::runtime_error("A diagnostic became non-finite");
    }
    result.final_energy = values.energy;
    result.final_l3 = values.critical_l3_sample;
    result.final_vorticity = values.sampled_vorticity_max;
    result.final_enstrophy = values.enstrophy;
    result.peak_l3 = std::max(result.peak_l3, values.critical_l3_sample);
    result.peak_vorticity =
        std::max(result.peak_vorticity, values.sampled_vorticity_max);
    result.peak_enstrophy = std::max(result.peak_enstrophy, values.enstrophy);
    result.peak_palinstrophy =
        std::max(result.peak_palinstrophy, values.palinstrophy);
    result.peak_cutoff_fraction = std::max(
        result.peak_cutoff_fraction, values.high_shell_energy_fraction);
    result.maximum_divergence_defect =
        std::max(result.maximum_divergence_defect, values.divergence_defect);
    result.maximum_reality_defect =
        std::max(result.maximum_reality_defect, values.reality_defect);
    result.maximum_energy_balance_residual = std::max(
        result.maximum_energy_balance_residual,
        values.energy_balance_residual);
    result.maximum_positive_forward_flux = std::max(
        result.maximum_positive_forward_flux,
        positiveForwardFlux(system.shellDiagnostics(state)));
}

double ratio(double peak, double initial) {
    return initial == 0.0 ? 0.0 : peak / initial;
}

RunResult runCandidate(const ns_cascade::PseudospectralSystem& system,
                       const Candidate& candidate,
                       const Options& options) {
    RunResult result(candidate, system.gridSize(), system.cutoff());
    const std::chrono::steady_clock::time_point start =
        std::chrono::steady_clock::now();
    ns_cascade::PseudospectralSystem::State state =
        system.vortexTubePairState(candidate.parameters, options.initial_energy);
    const ns_cascade::PseudospectralSystem::Diagnostics initial =
        system.diagnostics(state);
    if (!finiteDiagnostics(initial)) {
        throw std::runtime_error("Initial diagnostics are non-finite");
    }

    result.initial_three_dimensional_energy_fraction =
        threeDimensionalEnergyFraction(system, state);
    result.initial_energy = initial.energy;
    result.final_energy = initial.energy;
    result.initial_l3 = initial.critical_l3_sample;
    result.final_l3 = initial.critical_l3_sample;
    result.peak_l3 = initial.critical_l3_sample;
    result.initial_vorticity = initial.sampled_vorticity_max;
    result.final_vorticity = initial.sampled_vorticity_max;
    result.peak_vorticity = initial.sampled_vorticity_max;
    result.initial_enstrophy = initial.enstrophy;
    result.final_enstrophy = initial.enstrophy;
    result.peak_enstrophy = initial.enstrophy;
    result.peak_palinstrophy = initial.palinstrophy;
    result.maximum_positive_forward_flux = 0.0;
    result.peak_cutoff_fraction = 0.0;
    result.maximum_divergence_defect = 0.0;
    result.maximum_reality_defect = 0.0;
    result.maximum_energy_balance_residual = 0.0;
    recordDiagnostics(system, state, initial, result);

    double time = 0.0;
    result.minimum_time_step = std::numeric_limits<double>::infinity();
    result.maximum_time_step = 0.0;
    result.maximum_cfl_bound = 0.0;
    result.maximum_viscous_number = 0.0;
    const double time_tolerance =
        16.0 * std::numeric_limits<double>::epsilon() *
        std::max(1.0, options.final_time);
    while (time < options.final_time) {
        if (result.steps >= 1000000) {
            throw std::runtime_error("Candidate exceeded 1000000 accepted steps");
        }
        const double remaining = options.final_time - time;
        const ns_cascade::AdaptiveStepInfo step_information =
            system.chooseAdaptiveTimeStep(
                state,
                std::min(options.maximum_time_step, remaining),
                options.target_cfl,
                options.diffusion_safety);
        if (time + step_information.time_step == time) {
            throw std::runtime_error("Adaptive step is too small to advance time");
        }
        system.stepRungeKutta4(state, step_information.time_step);
        time += step_information.time_step;
        if (options.final_time - time <= time_tolerance) time = options.final_time;
        ++result.steps;
        result.minimum_time_step =
            std::min(result.minimum_time_step, step_information.time_step);
        result.maximum_time_step =
            std::max(result.maximum_time_step, step_information.time_step);
        result.maximum_cfl_bound = std::max(
            result.maximum_cfl_bound,
            step_information.advective_cfl_upper_bound);
        result.maximum_viscous_number = std::max(
            result.maximum_viscous_number,
            step_information.viscous_stability_number);

        if (result.steps % options.diagnostic_every == 0 ||
            time >= options.final_time) {
            recordDiagnostics(system, state, system.diagnostics(state), result);
        }
    }

    result.completed = true;
    result.resolution_ok =
        result.peak_cutoff_fraction <= options.cutoff_fraction_threshold &&
        result.maximum_divergence_defect <= 1e-9 &&
        result.maximum_reality_defect <= 1e-9 &&
        result.maximum_energy_balance_residual <= 1e-9;
    const double l3_growth = ratio(result.peak_l3, result.initial_l3);
    const double vorticity_growth =
        ratio(result.peak_vorticity, result.initial_vorticity);
    const double final_l3_ratio = ratio(result.final_l3, result.initial_l3);
    result.score = std::log(std::max(l3_growth, 1e-300)) +
                   0.05 * std::log(std::max(vorticity_growth, 1e-300)) +
                   0.01 * std::log(std::max(final_l3_ratio, 1e-300));
    result.elapsed_seconds =
        std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start)
            .count();
    return result;
}

RunResult runCandidateSafely(
    const ns_cascade::PseudospectralSystem& system,
    const Candidate& candidate,
    const Options& options) {
    try {
        return runCandidate(system, candidate, options);
    } catch (const std::exception& error) {
        RunResult failed(candidate, system.gridSize(), system.cutoff());
        failed.failure = error.what();
        return failed;
    }
}

bool betterCandidate(const RunResult& left, const RunResult& right) {
    if (left.completed != right.completed) return left.completed;
    if (left.resolution_ok != right.resolution_ok) return left.resolution_ok;
    if (left.score != right.score) return left.score > right.score;
    return left.candidate.id < right.candidate.id;
}

double relativeDifference(double left, double right) {
    const double scale = std::max(1e-12, std::max(std::abs(left), std::abs(right)));
    return std::abs(left - right) / scale;
}

bool crossResolutionOk(const RunResult& coarse,
                       const RunResult& fine,
                       double& l3_difference,
                       double& vorticity_difference,
                       double& flux_difference) {
    const double coarse_l3_ratio = ratio(coarse.peak_l3, coarse.initial_l3);
    const double fine_l3_ratio = ratio(fine.peak_l3, fine.initial_l3);
    const double coarse_vorticity_ratio =
        ratio(coarse.peak_vorticity, coarse.initial_vorticity);
    const double fine_vorticity_ratio =
        ratio(fine.peak_vorticity, fine.initial_vorticity);
    l3_difference = relativeDifference(coarse_l3_ratio, fine_l3_ratio);
    vorticity_difference = relativeDifference(
        coarse_vorticity_ratio, fine_vorticity_ratio);
    flux_difference = relativeDifference(
        coarse.maximum_positive_forward_flux,
        fine.maximum_positive_forward_flux);
    const bool l3_growth_signal_agrees =
        std::max(coarse_l3_ratio, fine_l3_ratio) <= 1.005 ||
        (coarse_l3_ratio > 1.005 && fine_l3_ratio > 1.005);
    const bool vorticity_growth_signal_agrees =
        std::max(coarse_vorticity_ratio, fine_vorticity_ratio) <= 1.005 ||
        (coarse_vorticity_ratio > 1.005 && fine_vorticity_ratio > 1.005);
    return coarse.completed && fine.completed && coarse.resolution_ok &&
           fine.resolution_ok && l3_difference <= 0.02 &&
           vorticity_difference <= 0.10 && flux_difference <= 0.25 &&
           l3_growth_signal_agrees && vorticity_growth_signal_agrees;
}

std::string csvString(const std::string& value) {
    std::string result = "\"";
    for (std::size_t i = 0; i < value.size(); ++i) {
        if (value[i] == '\"') result += '\"';
        result += value[i];
    }
    result += "\"";
    return result;
}

void writeHeader(std::ostream& output) {
    output
        << "stage,selection_rank,candidate_id,status,grid_size,cutoff,"
        << "core_radius,separation,bend_amplitude,axial_wavenumber,viscosity,"
        << "requested_initial_energy,final_time,steps,min_dt,max_dt,"
        << "max_advective_cfl_upper_bound,max_viscous_stability_number,"
        << "initial_3d_energy_fraction,initial_energy,final_energy,"
        << "initial_l3,final_l3,peak_l3,peak_l3_ratio,"
        << "initial_sampled_vorticity,final_sampled_vorticity,"
        << "peak_sampled_vorticity,peak_vorticity_ratio,"
        << "initial_enstrophy,final_enstrophy,peak_enstrophy,"
        << "peak_enstrophy_ratio,peak_palinstrophy,max_positive_forward_flux,"
        << "peak_cutoff_energy_fraction,max_divergence_defect,"
        << "max_reality_defect,max_energy_balance_residual,resolution_ok,"
        << "ranking_score,cpu_seconds,reference_grid,"
        << "l3_ratio_relative_difference,vorticity_ratio_relative_difference,"
        << "flux_relative_difference,cross_resolution_ok,failure\n";
}

void writeResult(std::ostream& output,
                 const std::string& stage,
                 int selection_rank,
                 const RunResult& result,
                 const Options& options,
                 const RunResult* reference) {
    output << stage << ',' << selection_rank << ',' << result.candidate.id << ','
           << (result.completed ? "completed" : "failed") << ','
           << result.grid_size << ',' << result.cutoff << ','
           << result.candidate.parameters.core_radius << ','
           << result.candidate.parameters.separation << ','
           << result.candidate.parameters.bend_amplitude << ','
           << result.candidate.parameters.axial_wavenumber << ','
           << options.viscosity << ',' << options.initial_energy << ','
           << options.final_time << ',' << result.steps << ','
           << result.minimum_time_step << ',' << result.maximum_time_step << ','
           << result.maximum_cfl_bound << ','
           << result.maximum_viscous_number << ','
           << result.initial_three_dimensional_energy_fraction << ','
           << result.initial_energy << ',' << result.final_energy << ','
           << result.initial_l3 << ',' << result.final_l3 << ','
           << result.peak_l3 << ','
           << ratio(result.peak_l3, result.initial_l3) << ','
           << result.initial_vorticity << ',' << result.final_vorticity << ','
           << result.peak_vorticity << ','
           << ratio(result.peak_vorticity, result.initial_vorticity) << ','
           << result.initial_enstrophy << ',' << result.final_enstrophy << ','
           << result.peak_enstrophy << ','
           << ratio(result.peak_enstrophy, result.initial_enstrophy) << ','
           << result.peak_palinstrophy << ','
           << result.maximum_positive_forward_flux << ','
           << result.peak_cutoff_fraction << ','
           << result.maximum_divergence_defect << ','
           << result.maximum_reality_defect << ','
           << result.maximum_energy_balance_residual << ','
           << (result.resolution_ok ? "true" : "false") << ','
           << result.score << ',' << result.elapsed_seconds << ',';
    if (reference == NULL) {
        output << ",,,,,";
    } else {
        double l3_difference = 0.0;
        double vorticity_difference = 0.0;
        double flux_difference = 0.0;
        const bool consistent = crossResolutionOk(
            *reference,
            result,
            l3_difference,
            vorticity_difference,
            flux_difference);
        output << reference->grid_size << ',' << l3_difference << ','
               << vorticity_difference << ',' << flux_difference << ','
               << (consistent ? "true" : "false") << ',';
    }
    output << csvString(result.failure) << '\n';
}

void printProgress(const std::string& stage,
                   int index,
                   int count,
                   const RunResult& result) {
    std::cout << stage << ' ' << index << '/' << count << ": candidate "
              << result.candidate.id << " core="
              << result.candidate.parameters.core_radius << " separation="
              << result.candidate.parameters.separation << " bend="
              << result.candidate.parameters.bend_amplitude << " axial="
              << result.candidate.parameters.axial_wavenumber;
    if (!result.completed) {
        std::cout << " FAILED: " << result.failure << '\n';
    } else {
        std::cout << " L3 ratio=" << ratio(result.peak_l3, result.initial_l3)
                  << " vorticity ratio="
                  << ratio(result.peak_vorticity, result.initial_vorticity)
                  << " cutoff fraction=" << result.peak_cutoff_fraction
                  << " resolved=" << (result.resolution_ok ? "yes" : "no")
                  << '\n';
    }
    std::cout.flush();
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parseOptions(argc, argv);
        const std::vector<Candidate> candidates = buildCandidates(options);
        const ns_cascade::PseudospectralSystem coarse_system(
            options.coarse_grid, options.viscosity);
        const ns_cascade::PseudospectralSystem fine_system(
            options.fine_grid, options.viscosity);

        std::vector<RunResult> coarse_results;
        coarse_results.reserve(candidates.size());
        for (std::size_t i = 0; i < candidates.size(); ++i) {
            const RunResult result =
                runCandidateSafely(coarse_system, candidates[i], options);
            printProgress("coarse", static_cast<int>(i + 1),
                          static_cast<int>(candidates.size()), result);
            coarse_results.push_back(result);
        }
        std::sort(coarse_results.begin(), coarse_results.end(), betterCandidate);

        const int finalist_count = std::min(
            options.top_count, static_cast<int>(coarse_results.size()));
        std::vector<RunResult> fine_results;
        std::vector<RunResult> references;
        fine_results.reserve(static_cast<std::size_t>(finalist_count));
        references.reserve(static_cast<std::size_t>(finalist_count));
        for (int i = 0; i < finalist_count; ++i) {
            const RunResult fine = runCandidateSafely(
                fine_system, coarse_results[static_cast<std::size_t>(i)].candidate,
                options);
            printProgress("fine", i + 1, finalist_count, fine);
            fine_results.push_back(fine);
            references.push_back(coarse_results[static_cast<std::size_t>(i)]);
        }

        std::ofstream csv(options.output.c_str());
        if (!csv) throw std::runtime_error("Cannot open output file: " + options.output);
        csv << std::setprecision(17);
        writeHeader(csv);
        for (std::size_t i = 0; i < coarse_results.size(); ++i) {
            writeResult(csv,
                        "coarse",
                        static_cast<int>(i + 1),
                        coarse_results[i],
                        options,
                        NULL);
        }
        for (std::size_t i = 0; i < fine_results.size(); ++i) {
            writeResult(csv,
                        "fine",
                        static_cast<int>(i + 1),
                        fine_results[i],
                        options,
                        &references[i]);
        }
        csv.close();
        if (!csv) throw std::runtime_error("Failed while writing " + options.output);

        int coarse_resolved = 0;
        for (std::size_t i = 0; i < coarse_results.size(); ++i) {
            if (coarse_results[i].resolution_ok) ++coarse_resolved;
        }
        int cross_resolved = 0;
        for (std::size_t i = 0; i < fine_results.size(); ++i) {
            double l3_difference = 0.0;
            double vorticity_difference = 0.0;
            double flux_difference = 0.0;
            if (crossResolutionOk(references[i],
                                  fine_results[i],
                                  l3_difference,
                                  vorticity_difference,
                                  flux_difference)) {
                ++cross_resolved;
            }
        }

        std::cout << std::setprecision(8)
                  << "Search complete: " << candidates.size()
                  << " coarse candidates, " << coarse_resolved
                  << " passed cutoff/constraint gates; " << cross_resolved << '/'
                  << finalist_count
                  << " finalists passed the preliminary cross-resolution gate.\n";
        if (!fine_results.empty()) {
            std::size_t best_index = 0;
            for (std::size_t i = 1; i < fine_results.size(); ++i) {
                if (betterCandidate(fine_results[i], fine_results[best_index])) {
                    best_index = i;
                }
            }
            const RunResult& best = fine_results[best_index];
            std::cout << "Best fine-grid candidate " << best.candidate.id
                      << ": peak L3 ratio="
                      << ratio(best.peak_l3, best.initial_l3)
                      << ", peak sampled-vorticity ratio="
                      << ratio(best.peak_vorticity, best.initial_vorticity)
                      << ", peak cutoff fraction=" << best.peak_cutoff_fraction
                      << ".\n";
            if (ratio(best.peak_l3, best.initial_l3) <= 1.005) {
                std::cout
                    << "No material critical-L3 growth was found in this short search.\n";
            } else {
                std::cout
                    << "Critical-L3 growth is only a candidate signal; it is not evidence of singularity without much stronger convergence and analysis.\n";
            }
        }
        std::cout << "Results written to " << options.output << "\n"
                  << "The ranking score is a triage heuristic, not a theorem or statistical significance measure.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
