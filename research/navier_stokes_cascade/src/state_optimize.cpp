#include "ns_cascade/candidate_score.hpp"
#include "ns_cascade/spectral_profile.hpp"
#include "ns_cascade/state_optimizer.hpp"

#include <algorithm>
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

const double kPi = 3.1415926535897932384626433832795;

enum class SeedFamily {
    WavePackets,
    VortexTubes,
    OrthogonalBundle
};

const char* seedFamilyName(SeedFamily family) {
    switch (family) {
        case SeedFamily::WavePackets: return "wave-packets";
        case SeedFamily::VortexTubes: return "vortex-tubes";
        case SeedFamily::OrthogonalBundle: return "orthogonal-bundle";
    }
    return "unknown";
}

SeedFamily parseSeedFamily(const std::string& text) {
    if (text == "wave-packets") return SeedFamily::WavePackets;
    if (text == "vortex-tubes") return SeedFamily::VortexTubes;
    if (text == "orthogonal-bundle") return SeedFamily::OrthogonalBundle;
    throw std::invalid_argument("Unknown state-optimizer seed family: " + text);
}

struct Options {
    int grid_size = 16;
    int fine_grid_size = 32;
    int cutoff = 0;
    int fine_cutoff = 0;
    int seed_bandwidth = 2;
    double viscosity = 0.02;
    double fixed_time_step = 0.0005;
    double fine_maximum_time_step = 0.005;
    double final_time = 0.04;
    double target_cfl = 0.35;
    double diffusion_safety = 2.0;
    int diagnostic_every = 5;
    double initial_energy = 10.0;
    SeedFamily seed_family = SeedFamily::WavePackets;
    int iterations = 2;
    double gradient_check_angle = 0.002;
    double gradient_tolerance = 0.01;
    double initial_line_angle = 0.12;
    double minimum_line_angle = 0.002;
    double armijo_fraction = 1e-4;
    ns_cascade::StateObjectiveWeights objective_weights;
    int profile_bin_count = 32;
    double profile_maximum_coordinate = 4.0;
    double profile_log_scale_window = 0.01;
    double profile_drift_threshold = 1.0;
    double minimum_characteristic_growth = 1.10;
    double minimum_critical_growth = 1.005;
    std::string state_input;
    std::string output = "navier_stokes_state_optimization.csv";
    std::string state_output = "navier_stokes_optimized_state.csv";
};

struct Evaluation {
    bool completed = false;
    bool valid = false;
    bool fixed_step_safe = true;
    std::string failure;
    int steps = 0;
    ns_cascade::StateObjectiveValue objective;
    ns_cascade::OptimizationState final_state;
    std::vector<ns_cascade::OptimizationState> trajectory;
    std::vector<double> time_steps;
    double initial_energy = 0.0;
    double final_energy = 0.0;
    double initial_h_half = 0.0;
    double final_h_half = 0.0;
    double peak_h_half = 0.0;
    double initial_l3 = 0.0;
    double final_l3 = 0.0;
    double peak_l3 = 0.0;
    double initial_vorticity = 0.0;
    double final_vorticity = 0.0;
    double peak_vorticity = 0.0;
    double initial_characteristic_wavenumber = 0.0;
    double final_characteristic_wavenumber = 0.0;
    int profile_scale_windows = 0;
    double first_profile_drift = std::numeric_limits<double>::quiet_NaN();
    double latest_profile_drift = std::numeric_limits<double>::quiet_NaN();
    double minimum_profile_drift = std::numeric_limits<double>::infinity();
    double latest_profile_window_time =
        std::numeric_limits<double>::quiet_NaN();
    bool profile_stationarity_improving = false;
    double peak_cutoff_fraction = 0.0;
    double maximum_divergence_defect = 0.0;
    double maximum_reality_defect = 0.0;
    double maximum_energy_balance_residual = 0.0;
    double maximum_positive_forward_flux = 0.0;
    double maximum_cfl_bound = 0.0;
    double maximum_viscous_number = 0.0;
    ns_cascade::SearchScore score;
};

struct PairAssessment {
    bool preliminary_cross_resolution_ok = false;
    bool profile_drift_comparison_valid = false;
    bool profile_windows_recent = false;
    double l3_ratio_relative_difference = 0.0;
    double h_half_ratio_relative_difference = 0.0;
    double vorticity_ratio_relative_difference = 0.0;
    double flux_relative_difference = 0.0;
    double characteristic_scale_relative_difference = 0.0;
    double profile_drift_relative_difference = 1.0;
    ns_cascade::PairedSearchScore score;
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

void printUsage(const char* program) {
    std::cout
        << "Usage: " << program << " [options]\n\n"
        << "Use the verified reverse discrete adjoint to optimize a low-band, "
           "real, divergence-free initial Fourier state on a fixed-energy "
           "sphere. The accepted state is rerun on a finer grid and remains "
           "subject to cutoff and profile gates.\n\n"
        << "  --grid N                    Coarse optimization grid (default: 16)\n"
        << "  --fine-grid N               Validation grid (default: 32)\n"
        << "  --cutoff K                  Coarse retained cutoff; zero is automatic\n"
        << "  --fine-cutoff K             Fine retained cutoff; zero is automatic\n"
        << "  --seed-bandwidth K          Initial-variable Fourier box (default: 2)\n"
        << "  --viscosity NU              Viscosity (default: 0.02)\n"
        << "  --dt DT                     Fixed differentiated step (default: 0.0005)\n"
        << "  --fine-max-dt DT            Fine adaptive maximum (default: 0.005)\n"
        << "  --final-time T              Objective horizon (default: 0.04)\n"
        << "  --cfl C                     Conservative CFL target (default: 0.35)\n"
        << "  --diffusion-safety S        Viscous RK4 bound (default: 2)\n"
        << "  --diagnostic-every N        Full diagnostic cadence (default: 5)\n"
        << "  --energy E                  Fixed initial energy (default: 10)\n"
        << "  --initial-family NAME       wave-packets, vortex-tubes, or "
           "orthogonal-bundle\n"
        << "  --iterations N              Maximum accepted steps (default: 2)\n"
        << "  --gradient-check-angle A    Geodesic FD check angle (default: 0.002)\n"
        << "  --gradient-tolerance R      Relative adjoint/FD gate (default: 0.01)\n"
        << "  --line-angle A              Initial ascent angle (default: 0.12)\n"
        << "  --minimum-line-angle A      Backtracking stop (default: 0.002)\n"
        << "  --scale-weight W            k_rms-growth reward (default: 0.15)\n"
        << "  --cutoff-weight W           Smooth cutoff cost (default: 0.04)\n"
        << "  --cutoff-threshold F        Hard cutoff gate (default: 0.01)\n"
        << "  --profile-scale-window X    log(k_rms) window (default: 0.01)\n"
        << "  --state-input PATH          Resume or replay a saved coefficient CSV\n"
        << "  --output PATH               Optimization trace CSV\n"
        << "  --state-output PATH         Optimized coefficient CSV\n"
        << "  --help                      Show this message\n";
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
        } else if (flag == "--fine-grid") {
            options.fine_grid_size =
                parseNumber<int>(requireValue(i, argc, argv), flag);
        } else if (flag == "--cutoff") {
            options.cutoff = parseNumber<int>(requireValue(i, argc, argv), flag);
        } else if (flag == "--fine-cutoff") {
            options.fine_cutoff =
                parseNumber<int>(requireValue(i, argc, argv), flag);
        } else if (flag == "--seed-bandwidth") {
            options.seed_bandwidth =
                parseNumber<int>(requireValue(i, argc, argv), flag);
        } else if (flag == "--viscosity") {
            options.viscosity =
                parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--dt") {
            options.fixed_time_step =
                parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--fine-max-dt") {
            options.fine_maximum_time_step =
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
        } else if (flag == "--energy") {
            options.initial_energy =
                parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--initial-family") {
            options.seed_family = parseSeedFamily(requireValue(i, argc, argv));
        } else if (flag == "--iterations") {
            options.iterations = parseNumber<int>(requireValue(i, argc, argv), flag);
        } else if (flag == "--gradient-check-angle") {
            options.gradient_check_angle =
                parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--gradient-tolerance") {
            options.gradient_tolerance =
                parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--line-angle") {
            options.initial_line_angle =
                parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--minimum-line-angle") {
            options.minimum_line_angle =
                parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--scale-weight") {
            options.objective_weights.characteristic_scale_weight =
                parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--cutoff-weight") {
            options.objective_weights.cutoff_penalty_weight =
                parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--cutoff-threshold") {
            options.objective_weights.cutoff_fraction_threshold =
                parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--profile-scale-window") {
            options.profile_log_scale_window =
                parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--state-input") {
            options.state_input = requireValue(i, argc, argv);
        } else if (flag == "--output") {
            options.output = requireValue(i, argc, argv);
        } else if (flag == "--state-output") {
            options.state_output = requireValue(i, argc, argv);
        } else {
            throw std::invalid_argument("Unknown option: " + flag);
        }
    }

    if (options.grid_size < 8 || options.fine_grid_size <= options.grid_size ||
        options.cutoff < 0 || options.fine_cutoff < 0 ||
        options.seed_bandwidth < 1 || options.iterations < 0 ||
        options.diagnostic_every < 1 || options.output.empty() ||
        options.state_output.empty() || options.output == options.state_output ||
        (!options.state_input.empty() &&
         (options.state_input == options.output ||
          options.state_input == options.state_output))) {
        throw std::invalid_argument("State-optimizer integer or path is invalid");
    }
    const double positive[] = {
        options.viscosity,
        options.fixed_time_step,
        options.fine_maximum_time_step,
        options.final_time,
        options.target_cfl,
        options.diffusion_safety,
        options.initial_energy,
        options.gradient_check_angle,
        options.gradient_tolerance,
        options.initial_line_angle,
        options.minimum_line_angle,
        options.armijo_fraction,
        options.profile_maximum_coordinate,
        options.profile_log_scale_window,
        options.profile_drift_threshold,
        options.minimum_characteristic_growth,
        options.minimum_critical_growth};
    for (std::size_t i = 0; i < sizeof(positive) / sizeof(positive[0]); ++i) {
        if (!std::isfinite(positive[i]) || positive[i] <= 0.0) {
            throw std::invalid_argument(
                "Positive state-optimizer option is invalid");
        }
    }
    if (options.minimum_line_angle > options.initial_line_angle ||
        options.gradient_check_angle >= 0.25 ||
        options.initial_line_angle >= 1.0 ||
        options.profile_log_scale_window >= 1.0 ||
        options.minimum_characteristic_growth <= 1.0 ||
        options.minimum_critical_growth <= 1.0) {
        throw std::invalid_argument("State-optimizer scale or gate is invalid");
    }
    ns_cascade::validateStateObjectiveWeights(options.objective_weights);
    return options;
}

bool finiteDiagnostics(
    const ns_cascade::PseudospectralSystem::Diagnostics& values) {
    return std::isfinite(values.energy) &&
           std::isfinite(values.enstrophy) &&
           std::isfinite(values.palinstrophy) &&
           std::isfinite(values.critical_h_half) &&
           std::isfinite(values.critical_l3_sample) &&
           std::isfinite(values.sampled_vorticity_max) &&
           std::isfinite(values.high_shell_energy_fraction) &&
           std::isfinite(values.divergence_defect) &&
           std::isfinite(values.reality_defect) &&
           std::isfinite(values.energy_balance_residual);
}

double positiveForwardFlux(
    const std::vector<ns_cascade::GalerkinSystem::ShellDiagnostics>& shells) {
    double maximum = 0.0;
    for (std::size_t i = 0; i < shells.size(); ++i) {
        maximum = std::max(maximum, shells[i].forward_flux);
    }
    return maximum;
}

double ratio(double value, double reference) {
    return reference == 0.0 ? 0.0 : value / reference;
}

double relativeDifference(double left, double right) {
    return std::abs(left - right) /
           std::max(1e-15, std::max(std::abs(left), std::abs(right)));
}

void recordDiagnostics(
    const ns_cascade::PseudospectralSystem& system,
    const ns_cascade::OptimizationState& state,
    const ns_cascade::PseudospectralSystem::Diagnostics& values,
    Evaluation& result) {
    if (!finiteDiagnostics(values)) {
        throw std::runtime_error("State-optimizer diagnostic became non-finite");
    }
    result.final_energy = values.energy;
    result.final_h_half = values.critical_h_half;
    result.final_l3 = values.critical_l3_sample;
    result.final_vorticity = values.sampled_vorticity_max;
    result.peak_h_half = std::max(result.peak_h_half, values.critical_h_half);
    result.peak_l3 = std::max(result.peak_l3, values.critical_l3_sample);
    result.peak_vorticity =
        std::max(result.peak_vorticity, values.sampled_vorticity_max);
    result.peak_cutoff_fraction = std::max(
        result.peak_cutoff_fraction, values.high_shell_energy_fraction);
    result.maximum_divergence_defect = std::max(
        result.maximum_divergence_defect, values.divergence_defect);
    result.maximum_reality_defect = std::max(
        result.maximum_reality_defect, values.reality_defect);
    result.maximum_energy_balance_residual = std::max(
        result.maximum_energy_balance_residual,
        values.energy_balance_residual);
    result.maximum_positive_forward_flux = std::max(
        result.maximum_positive_forward_flux,
        positiveForwardFlux(system.shellDiagnostics(state)));
}

ns_cascade::SearchScoreEvidence scoreEvidence(
    const Evaluation& result,
    const Options& options) {
    ns_cascade::SearchScoreEvidence evidence;
    evidence.peak_l3_ratio = ratio(result.peak_l3, result.initial_l3);
    evidence.peak_h_half_ratio =
        ratio(result.peak_h_half, result.initial_h_half);
    evidence.final_l3_ratio = ratio(result.final_l3, result.initial_l3);
    evidence.peak_vorticity_ratio =
        ratio(result.peak_vorticity, result.initial_vorticity);
    evidence.peak_cutoff_fraction = result.peak_cutoff_fraction;
    evidence.cutoff_fraction_threshold =
        options.objective_weights.cutoff_fraction_threshold;
    evidence.profile_scale_windows = result.profile_scale_windows;
    evidence.first_profile_drift = result.profile_scale_windows > 0
                                       ? result.first_profile_drift
                                       : 0.0;
    evidence.latest_profile_drift = result.profile_scale_windows > 0
                                        ? result.latest_profile_drift
                                        : 0.0;
    evidence.minimum_profile_drift = result.profile_scale_windows > 0
                                         ? result.minimum_profile_drift
                                         : 0.0;
    evidence.profile_drift_threshold = options.profile_drift_threshold;
    return evidence;
}

Evaluation evaluateTrajectory(
    const ns_cascade::PseudospectralSystem& system,
    const ns_cascade::OptimizationState& initial_state,
    const Options& options,
    bool fixed_time_step,
    bool keep_trajectory) {
    Evaluation result;
    ns_cascade::OptimizationState state = initial_state;
    const ns_cascade::PseudospectralSystem::Diagnostics initial =
        system.diagnostics(state);
    if (!finiteDiagnostics(initial)) {
        throw std::runtime_error("Initial state-optimizer diagnostic is invalid");
    }
    const ns_cascade::SpectrumProfile initial_profile =
        ns_cascade::rescaledSpectrumProfile(
            system,
            state,
            options.profile_bin_count,
            options.profile_maximum_coordinate);
    ns_cascade::SpectrumProfile profile_anchor = initial_profile;
    result.initial_energy = initial.energy;
    result.final_energy = initial.energy;
    result.initial_h_half = initial.critical_h_half;
    result.final_h_half = initial.critical_h_half;
    result.peak_h_half = initial.critical_h_half;
    result.initial_l3 = initial.critical_l3_sample;
    result.final_l3 = initial.critical_l3_sample;
    result.peak_l3 = initial.critical_l3_sample;
    result.initial_vorticity = initial.sampled_vorticity_max;
    result.final_vorticity = initial.sampled_vorticity_max;
    result.peak_vorticity = initial.sampled_vorticity_max;
    result.initial_characteristic_wavenumber =
        initial_profile.characteristic_wavenumber;
    result.final_characteristic_wavenumber =
        initial_profile.characteristic_wavenumber;
    recordDiagnostics(system, state, initial, result);

    const double time_tolerance =
        16.0 * std::numeric_limits<double>::epsilon() *
        std::max(1.0, options.final_time);
    double time = 0.0;
    while (time < options.final_time) {
        if (result.steps >= 1000000) {
            throw std::runtime_error(
                "State optimizer exceeded one million trajectory steps");
        }
        const double remaining = options.final_time - time;
        const double proposed = std::min(
            fixed_time_step ? options.fixed_time_step
                            : options.fine_maximum_time_step,
            remaining);
        const ns_cascade::AdaptiveStepInfo safe =
            system.chooseAdaptiveTimeStep(
                state,
                proposed,
                options.target_cfl,
                options.diffusion_safety);
        double dt = safe.time_step;
        if (fixed_time_step) {
            const double tolerance =
                64.0 * std::numeric_limits<double>::epsilon() * proposed;
            if (safe.time_step + tolerance < proposed) {
                result.fixed_step_safe = false;
                result.failure = "fixed timestep violates conservative bound";
                result.final_state = state;
                return result;
            }
            dt = proposed;
        }
        if (!std::isfinite(dt) || dt <= 0.0 || time + dt == time) {
            throw std::runtime_error("State optimizer selected an invalid step");
        }
        if (keep_trajectory) {
            result.trajectory.push_back(state);
            result.time_steps.push_back(dt);
        }
        system.stepRungeKutta4(state, dt);
        time += dt;
        if (options.final_time - time <= time_tolerance) {
            time = options.final_time;
        }
        ++result.steps;
        result.maximum_cfl_bound = std::max(
            result.maximum_cfl_bound, safe.advective_cfl_upper_bound);
        result.maximum_viscous_number = std::max(
            result.maximum_viscous_number, safe.viscous_stability_number);

        const double energy = system.energy(state);
        const double enstrophy = system.enstrophy(state);
        const double cutoff_fraction =
            system.cutoffShellEnergyFraction(state);
        if (!std::isfinite(energy) || energy <= 0.0 ||
            !std::isfinite(enstrophy) || enstrophy <= 0.0 ||
            !std::isfinite(cutoff_fraction) || cutoff_fraction < 0.0 ||
            cutoff_fraction > 1.0 + 1e-12) {
            throw std::runtime_error("Per-step state spectral gate is invalid");
        }
        result.peak_cutoff_fraction = std::max(
            result.peak_cutoff_fraction, cutoff_fraction);
        const double characteristic_wavenumber =
            std::sqrt(enstrophy / energy);
        const double forward_log_scale = std::log(
            characteristic_wavenumber /
            profile_anchor.characteristic_wavenumber);
        if (forward_log_scale >= options.profile_log_scale_window) {
            const ns_cascade::SpectrumProfile current_profile =
                ns_cascade::rescaledSpectrumProfile(
                    system,
                    state,
                    options.profile_bin_count,
                    options.profile_maximum_coordinate);
            const ns_cascade::SpectrumProfileChange change =
                ns_cascade::spectrumProfileChange(
                    current_profile, profile_anchor);
            if (!change.scale_normalized_drift_valid) {
                throw std::runtime_error(
                    "Completed state profile window has invalid drift");
            }
            if (result.profile_scale_windows == 0) {
                result.first_profile_drift =
                    change.l1_per_log_scale_change;
            }
            ++result.profile_scale_windows;
            result.latest_profile_drift =
                change.l1_per_log_scale_change;
            result.minimum_profile_drift = std::min(
                result.minimum_profile_drift,
                change.l1_per_log_scale_change);
            result.latest_profile_window_time = time;
            result.final_characteristic_wavenumber =
                current_profile.characteristic_wavenumber;
            profile_anchor = current_profile;
        }

        if (result.steps % options.diagnostic_every == 0 ||
            time >= options.final_time) {
            recordDiagnostics(
                system, state, system.diagnostics(state), result);
        }
    }

    result.completed = true;
    result.final_state = state;
    result.objective = ns_cascade::evaluateStateObjective(
        system, initial_state, state, options.objective_weights);
    const ns_cascade::SpectrumProfile final_profile =
        ns_cascade::rescaledSpectrumProfile(
            system,
            state,
            options.profile_bin_count,
            options.profile_maximum_coordinate);
    result.final_characteristic_wavenumber =
        final_profile.characteristic_wavenumber;
    result.valid =
        result.fixed_step_safe &&
        result.peak_cutoff_fraction <=
            options.objective_weights.cutoff_fraction_threshold &&
        result.maximum_divergence_defect <= 1e-9 &&
        result.maximum_reality_defect <= 1e-9 &&
        result.maximum_energy_balance_residual <= 1e-9;
    if (!result.valid) {
        if (result.peak_cutoff_fraction >
            options.objective_weights.cutoff_fraction_threshold) {
            result.failure = "cutoff-shell energy exceeds the hard gate";
        } else if (result.maximum_divergence_defect > 1e-9) {
            result.failure = "divergence defect exceeds the invariant gate";
        } else if (result.maximum_reality_defect > 1e-9) {
            result.failure = "Fourier-reality defect exceeds the invariant gate";
        } else if (result.maximum_energy_balance_residual > 1e-9) {
            result.failure = "energy-balance residual exceeds the invariant gate";
        } else {
            result.failure = "trajectory failed an unspecified hard gate";
        }
    }
    result.score = ns_cascade::scoreSingleResolution(
        scoreEvidence(result, options));
    result.profile_stationarity_improving =
        result.score.profile_stationarity_improving;
    return result;
}

const char* stateCsvHeader() {
    return "family,source_grid,simulation_cutoff,seed_bandwidth,"
           "target_energy,kx,ky,kz,ux_real,ux_imag,uy_real,uy_imag,"
           "uz_real,uz_imag";
}

std::vector<std::string> splitCsvRow(const std::string& line) {
    std::vector<std::string> fields;
    std::istringstream stream(line);
    std::string field;
    while (std::getline(stream, field, ',')) fields.push_back(field);
    if (!line.empty() && line[line.size() - 1U] == ',') fields.push_back("");
    return fields;
}

ns_cascade::OptimizationState readOptimizationState(
    const std::string& path,
    const ns_cascade::PseudospectralSystem& system,
    const Options& options) {
    std::ifstream input(path.c_str());
    if (!input) throw std::runtime_error("Could not open state-input CSV");
    std::string line;
    if (!std::getline(input, line) || line != stateCsvHeader()) {
        throw std::runtime_error("State-input CSV header is invalid");
    }

    ns_cascade::OptimizationState state = system.zeroState();
    std::vector<bool> seen(system.gridPointCount(), false);
    std::size_t rows = 0U;
    int source_grid = 0;
    int source_cutoff = 0;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> fields = splitCsvRow(line);
        if (fields.size() != 14U) {
            throw std::runtime_error("State-input CSV row has the wrong width");
        }
        if (fields[0] != seedFamilyName(options.seed_family)) {
            throw std::runtime_error(
                "State-input family does not match --initial-family");
        }
        const int row_grid = parseNumber<int>(fields[1], "state source grid");
        const int row_cutoff =
            parseNumber<int>(fields[2], "state source cutoff");
        const int row_bandwidth =
            parseNumber<int>(fields[3], "state seed bandwidth");
        const double row_energy =
            parseNumber<double>(fields[4], "state target energy");
        if (rows == 0U) {
            source_grid = row_grid;
            source_cutoff = row_cutoff;
        }
        if (row_grid != source_grid || row_cutoff != source_cutoff ||
            row_bandwidth != options.seed_bandwidth || source_grid < 8 ||
            source_cutoff < options.seed_bandwidth ||
            !std::isfinite(row_energy) ||
            relativeDifference(row_energy, options.initial_energy) > 1e-13) {
            throw std::runtime_error("State-input CSV metadata is inconsistent");
        }
        const ns_cascade::WaveVector wave(
            parseNumber<int>(fields[5], "state kx"),
            parseNumber<int>(fields[6], "state ky"),
            parseNumber<int>(fields[7], "state kz"));
        if (wave.normSquared() == 0 ||
            ns_cascade::stateMaximumComponent(wave) >
                options.seed_bandwidth) {
            throw std::runtime_error("State-input CSV contains a forbidden mode");
        }
        const std::size_t index = system.indexOf(wave);
        if (seen[index]) {
            throw std::runtime_error("State-input CSV repeats a Fourier mode");
        }
        const double ux_real = parseNumber<double>(fields[8], "state ux real");
        const double ux_imag = parseNumber<double>(fields[9], "state ux imag");
        const double uy_real = parseNumber<double>(fields[10], "state uy real");
        const double uy_imag = parseNumber<double>(fields[11], "state uy imag");
        const double uz_real = parseNumber<double>(fields[12], "state uz real");
        const double uz_imag = parseNumber<double>(fields[13], "state uz imag");
        const double components[] = {
            ux_real, ux_imag, uy_real, uy_imag, uz_real, uz_imag};
        for (std::size_t component = 0;
             component < sizeof(components) / sizeof(components[0]);
             ++component) {
            if (!std::isfinite(components[component])) {
                throw std::runtime_error(
                    "State-input CSV contains a non-finite coefficient");
            }
        }
        state[index] = ns_cascade::ComplexVector(
            ns_cascade::Complex(ux_real, ux_imag),
            ns_cascade::Complex(uy_real, uy_imag),
            ns_cascade::Complex(uz_real, uz_imag));
        seen[index] = true;
        ++rows;
    }
    if (!input.eof()) {
        throw std::runtime_error("Failed while reading state-input CSV");
    }
    const std::size_t expected_rows =
        ns_cascade::stateOptimizationDegreesOfFreedom(
            options.seed_bandwidth) /
        2U;
    if (rows != expected_rows ||
        relativeDifference(system.energy(state), options.initial_energy) >
            2e-13 ||
        system.divergenceDefect(state) > 2e-12 ||
        system.realityDefect(state) > 2e-12) {
        throw std::runtime_error(
            "State-input coefficients fail completeness or invariant checks");
    }
    return state;
}

ns_cascade::OptimizationState initialState(const Options& options,
                                           const ns_cascade::PseudospectralSystem& system) {
    if (!options.state_input.empty()) {
        return readOptimizationState(options.state_input, system, options);
    }
    ns_cascade::OptimizationState raw;
    if (options.seed_family == SeedFamily::WavePackets) {
        raw = system.interactingWavePacketState(
            ns_cascade::WavePacketParameters(
                1.0877734104170571,
                1,
                0.7577474832313891,
                0.299291384334349),
            options.initial_energy);
    } else if (options.seed_family == SeedFamily::VortexTubes) {
        raw = system.vortexTubePairState(
            ns_cascade::VortexTubeParameters(0.55, 1.2, 0.3, 2),
            options.initial_energy);
    } else {
        raw = system.vortexBundleState(
            ns_cascade::VortexBundleParameters(
                ns_cascade::VortexTubeParameters(0.70, 1.2, 0.25, 2),
                0.75,
                kPi / 3.0),
            options.initial_energy);
    }
    return ns_cascade::makeBandLimitedOptimizationState(
        system,
        raw,
        options.seed_bandwidth,
        options.initial_energy);
}

ns_cascade::OptimizationState fullInitialGradient(
    const ns_cascade::PseudospectralSystem& system,
    const ns_cascade::OptimizationState& initial,
    const Evaluation& evaluation,
    const Options& options) {
    if (!evaluation.completed ||
        evaluation.trajectory.size() != evaluation.time_steps.size()) {
        throw std::invalid_argument(
            "Adjoint gradient requires a complete stored trajectory");
    }
    ns_cascade::OptimizationState gradient =
        ns_cascade::terminalStateObjectiveGradient(
            system,
            evaluation.final_state,
            options.objective_weights);
    for (std::size_t step = evaluation.trajectory.size(); step-- > 0;) {
        gradient = system.adjointRungeKutta4Step(
            evaluation.trajectory[step],
            gradient,
            evaluation.time_steps[step]);
    }
    return ns_cascade::addOptimizationStates(
        gradient,
        ns_cascade::initialStateObjectiveGradient(
            system, initial, options.objective_weights),
        1.0);
}

PairAssessment assessPair(const Evaluation& coarse,
                          const Evaluation& fine,
                          const Options& options) {
    PairAssessment assessment;
    if (!coarse.completed || !fine.completed) return assessment;
    const double coarse_l3 = ratio(coarse.peak_l3, coarse.initial_l3);
    const double fine_l3 = ratio(fine.peak_l3, fine.initial_l3);
    const double coarse_h = ratio(coarse.peak_h_half, coarse.initial_h_half);
    const double fine_h = ratio(fine.peak_h_half, fine.initial_h_half);
    const double coarse_vorticity =
        ratio(coarse.peak_vorticity, coarse.initial_vorticity);
    const double fine_vorticity =
        ratio(fine.peak_vorticity, fine.initial_vorticity);
    assessment.l3_ratio_relative_difference =
        relativeDifference(coarse_l3, fine_l3);
    assessment.h_half_ratio_relative_difference =
        relativeDifference(coarse_h, fine_h);
    assessment.vorticity_ratio_relative_difference =
        relativeDifference(coarse_vorticity, fine_vorticity);
    assessment.flux_relative_difference = relativeDifference(
        coarse.maximum_positive_forward_flux,
        fine.maximum_positive_forward_flux);
    const bool l3_signal = std::max(coarse_l3, fine_l3) <= 1.005 ||
                           (coarse_l3 > 1.005 && fine_l3 > 1.005);
    const bool h_signal = std::max(coarse_h, fine_h) <= 1.005 ||
                          (coarse_h > 1.005 && fine_h > 1.005);
    const bool vorticity_signal =
        std::max(coarse_vorticity, fine_vorticity) <= 1.005 ||
        (coarse_vorticity > 1.005 && fine_vorticity > 1.005);
    assessment.preliminary_cross_resolution_ok =
        coarse.valid && fine.valid &&
        assessment.l3_ratio_relative_difference <= 0.02 &&
        assessment.h_half_ratio_relative_difference <= 0.02 &&
        assessment.vorticity_ratio_relative_difference <= 0.10 &&
        assessment.flux_relative_difference <= 0.25 &&
        l3_signal && h_signal && vorticity_signal;

    const double coarse_scale = ratio(
        coarse.final_characteristic_wavenumber,
        coarse.initial_characteristic_wavenumber);
    const double fine_scale = ratio(
        fine.final_characteristic_wavenumber,
        fine.initial_characteristic_wavenumber);
    assessment.characteristic_scale_relative_difference =
        relativeDifference(coarse_scale, fine_scale);
    assessment.profile_drift_comparison_valid =
        coarse.profile_scale_windows > 0 && fine.profile_scale_windows > 0;
    assessment.profile_windows_recent =
        assessment.profile_drift_comparison_valid &&
        coarse.latest_profile_window_time >= 0.75 * options.final_time &&
        fine.latest_profile_window_time >= 0.75 * options.final_time;
    if (assessment.profile_drift_comparison_valid) {
        assessment.profile_drift_relative_difference = relativeDifference(
            coarse.latest_profile_drift, fine.latest_profile_drift);
    }

    ns_cascade::PairedScoreEvidence evidence;
    evidence.coarse = scoreEvidence(coarse, options);
    evidence.fine = scoreEvidence(fine, options);
    evidence.preliminary_cross_resolution_ok =
        assessment.preliminary_cross_resolution_ok;
    evidence.l3_ratio_relative_difference =
        assessment.l3_ratio_relative_difference;
    evidence.h_half_ratio_relative_difference =
        assessment.h_half_ratio_relative_difference;
    evidence.characteristic_scale_relative_difference =
        assessment.characteristic_scale_relative_difference;
    evidence.profile_drift_relative_difference =
        assessment.profile_drift_relative_difference;
    evidence.profile_drift_comparison_valid =
        assessment.profile_drift_comparison_valid;
    evidence.coarse_profile_window_recent =
        coarse.latest_profile_window_time >= 0.75 * options.final_time;
    evidence.fine_profile_window_recent =
        fine.latest_profile_window_time >= 0.75 * options.final_time;
    evidence.coarse_characteristic_growth = coarse_scale;
    evidence.fine_characteristic_growth = fine_scale;
    evidence.minimum_characteristic_growth =
        options.minimum_characteristic_growth;
    evidence.minimum_critical_growth = options.minimum_critical_growth;
    assessment.score = ns_cascade::scoreResolutionPair(evidence);
    return assessment;
}

void writeTraceRow(std::ostream& output,
                   int iteration,
                   const char* stage,
                   const char* resolution,
                   const Evaluation& value,
                   const Options& options,
                   double gradient_norm,
                   double gradient_check,
                   double gradient_slope,
                   double line_angle,
                   int trajectory_count,
                   const PairAssessment* pair = nullptr) {
    output << iteration << ',' << stage << ',' << resolution << ','
           << seedFamilyName(options.seed_family) << ','
           << options.seed_bandwidth << ','
           << ns_cascade::stateOptimizationDegreesOfFreedom(
                  options.seed_bandwidth)
           << ',' << value.objective.total << ','
           << (value.valid ? "true" : "false") << ',' << value.steps << ','
           << value.objective.critical_log_growth << ','
           << value.objective.characteristic_log_growth << ','
           << value.objective.cutoff_penalty << ','
           << ratio(value.peak_h_half, value.initial_h_half) << ','
           << ratio(value.final_h_half, value.initial_h_half) << ','
           << ratio(value.peak_l3, value.initial_l3) << ','
           << ratio(value.final_l3, value.initial_l3) << ','
           << ratio(value.peak_vorticity, value.initial_vorticity) << ','
           << ratio(value.final_characteristic_wavenumber,
                    value.initial_characteristic_wavenumber)
           << ',' << value.profile_scale_windows << ','
           << value.first_profile_drift << ','
           << value.latest_profile_drift << ','
           << value.minimum_profile_drift << ','
           << value.latest_profile_window_time << ','
           << (value.profile_stationarity_improving ? "true" : "false")
           << ',' << value.score.total << ',' << value.peak_cutoff_fraction
           << ',' << value.maximum_divergence_defect << ','
           << value.maximum_reality_defect << ','
           << value.maximum_energy_balance_residual << ','
           << value.maximum_cfl_bound << ',' << value.maximum_viscous_number
           << ',' << gradient_norm << ',' << gradient_check << ','
           << gradient_slope << ',' << line_angle << ',' << trajectory_count
           << ',' << (pair == nullptr ?
                          std::numeric_limits<double>::quiet_NaN() :
                          pair->score.total)
           << ',' << (pair != nullptr && pair->preliminary_cross_resolution_ok
                          ? "true" : "false")
           << ',' << (pair != nullptr && pair->score.refinement_eligible
                          ? "true" : "false")
           << '\n';
}

void writeState(const std::string& path,
                const ns_cascade::PseudospectralSystem& system,
                const ns_cascade::OptimizationState& state,
                const Options& options) {
    std::ofstream output(path.c_str());
    if (!output) throw std::runtime_error("Could not open optimized-state CSV");
    output << std::setprecision(17) << stateCsvHeader() << '\n';
    const std::vector<ns_cascade::WaveVector>& modes = system.gridModes();
    for (std::size_t i = 0; i < state.size(); ++i) {
        const ns_cascade::WaveVector& wave = modes[i];
        if (wave.normSquared() == 0 ||
            ns_cascade::stateMaximumComponent(wave) >
                options.seed_bandwidth) {
            continue;
        }
        output << seedFamilyName(options.seed_family) << ','
               << system.gridSize() << ',' << system.cutoff() << ','
               << options.seed_bandwidth << ',' << options.initial_energy
               << ',' << wave.x << ',' << wave.y << ',' << wave.z << ','
               << std::real(state[i].x) << ',' << std::imag(state[i].x)
               << ',' << std::real(state[i].y) << ','
               << std::imag(state[i].y) << ',' << std::real(state[i].z)
               << ',' << std::imag(state[i].z) << '\n';
    }
    if (!output) {
        throw std::runtime_error("Failed while writing optimized-state CSV");
    }
}

int run(const Options& options) {
    const ns_cascade::PseudospectralSystem coarse_system(
        options.grid_size, options.viscosity, options.cutoff);
    const ns_cascade::PseudospectralSystem fine_system(
        options.fine_grid_size, options.viscosity, options.fine_cutoff);
    if (options.seed_bandwidth > coarse_system.cutoff()) {
        throw std::invalid_argument(
            "Seed bandwidth exceeds the coarse retained cutoff");
    }

    ns_cascade::OptimizationState current_state =
        initialState(options, coarse_system);
    std::ofstream trace(options.output.c_str());
    if (!trace) throw std::runtime_error("Could not open state-optimizer CSV");
    trace << std::setprecision(17)
          << "iteration,stage,resolution,initial_family,seed_bandwidth,"
             "real_degrees_of_freedom,objective,valid,steps,"
             "critical_log_growth,characteristic_log_growth,cutoff_penalty,"
             "peak_h_half_ratio,final_h_half_ratio,peak_l3_ratio,"
             "final_l3_ratio,peak_vorticity_ratio,characteristic_ratio,"
             "profile_scale_windows,first_profile_drift,latest_profile_drift,"
             "minimum_profile_drift,latest_profile_window_time,"
             "profile_stationarity_improving,search_score,"
             "peak_cutoff_fraction,max_divergence_defect,max_reality_defect,"
             "max_energy_balance_residual,max_cfl_bound,max_viscous_number,"
             "constrained_gradient_norm,gradient_relative_error,"
             "gradient_slope,line_angle,trajectory_count,pair_search_score,"
             "preliminary_cross_resolution_ok,refinement_eligible\n";

    int trajectory_count = 0;
    int accepted_steps = 0;
    Evaluation current = evaluateTrajectory(
        coarse_system, current_state, options, true, true);
    ++trajectory_count;
    if (!current.valid) {
        throw std::runtime_error(
            "Initial state fails coarse fixed-step or resolution gates: " +
            current.failure);
    }
    const double nan = std::numeric_limits<double>::quiet_NaN();
    writeTraceRow(trace,
                  0,
                  "initial",
                  "coarse",
                  current,
                  options,
                  nan,
                  nan,
                  nan,
                  0.0,
                  trajectory_count);

    const ns_cascade::OptimizationState initial_fine_state =
        ns_cascade::liftOptimizationState(
            coarse_system, current_state, fine_system);
    const Evaluation initial_fine = evaluateTrajectory(
        fine_system, initial_fine_state, options, false, false);
    ++trajectory_count;
    const PairAssessment initial_pair = assessPair(
        current, initial_fine, options);
    writeTraceRow(trace,
                  0,
                  "fine-check",
                  "fine",
                  initial_fine,
                  options,
                  nan,
                  nan,
                  nan,
                  0.0,
                  trajectory_count,
                  &initial_pair);
    ns_cascade::OptimizationState best_state = current_state;
    Evaluation best_coarse = current;
    best_coarse.trajectory.clear();
    best_coarse.time_steps.clear();
    Evaluation best_fine = initial_fine;
    PairAssessment best_pair = initial_pair;
    int best_iteration = 0;

    for (int iteration = 1; iteration <= options.iterations; ++iteration) {
        if (current.trajectory.empty()) {
            current = evaluateTrajectory(
                coarse_system, current_state, options, true, true);
            ++trajectory_count;
        }
        const ns_cascade::OptimizationState full_gradient =
            fullInitialGradient(
                coarse_system, current_state, current, options);
        const ns_cascade::OptimizationState projected_gradient =
            ns_cascade::projectStateObjectiveGradient(
                coarse_system,
                current_state,
                full_gradient,
                options.seed_bandwidth);
        const double gradient_norm =
            ns_cascade::optimizationStateNorm(projected_gradient);
        const ns_cascade::OptimizationState direction =
            ns_cascade::stateEnergySphereDirection(
                coarse_system,
                current_state,
                full_gradient,
                options.seed_bandwidth);
        const double gradient_slope =
            ns_cascade::stateRealInnerProduct(full_gradient, direction);
        if (!std::isfinite(gradient_slope) || gradient_slope <= 0.0) {
            throw std::runtime_error(
                "Projected adjoint direction is not an ascent direction");
        }

        const ns_cascade::OptimizationState check_plus_state =
            ns_cascade::stateEnergySphereStep(
                coarse_system,
                current_state,
                direction,
                options.gradient_check_angle);
        const ns_cascade::OptimizationState check_minus_state =
            ns_cascade::stateEnergySphereStep(
                coarse_system,
                current_state,
                direction,
                -options.gradient_check_angle);
        const Evaluation check_plus = evaluateTrajectory(
            coarse_system, check_plus_state, options, true, false);
        const Evaluation check_minus = evaluateTrajectory(
            coarse_system, check_minus_state, options, true, false);
        trajectory_count += 2;
        double checked_slope = nan;
        double gradient_relative_error =
            std::numeric_limits<double>::infinity();
        if (check_plus.valid && check_minus.valid) {
            checked_slope =
                (check_plus.objective.total - check_minus.objective.total) /
                (2.0 * options.gradient_check_angle);
            gradient_relative_error = relativeDifference(
                checked_slope, gradient_slope);
        }
        writeTraceRow(trace,
                      iteration,
                      "gradient-check",
                      "coarse",
                      current,
                      options,
                      gradient_norm,
                      gradient_relative_error,
                      gradient_slope,
                      0.0,
                      trajectory_count);
        if (!check_plus.valid || !check_minus.valid ||
            gradient_relative_error > options.gradient_tolerance) {
            std::cout << "iteration " << iteration
                      << " gradient rejected: relative error="
                      << gradient_relative_error << '\n';
            break;
        }

        double line_angle = options.initial_line_angle;
        bool accepted = false;
        Evaluation accepted_value;
        ns_cascade::OptimizationState accepted_state;
        while (line_angle >= options.minimum_line_angle) {
            const ns_cascade::OptimizationState trial_state =
                ns_cascade::stateEnergySphereStep(
                    coarse_system,
                    current_state,
                    direction,
                    line_angle);
            const Evaluation trial = evaluateTrajectory(
                coarse_system, trial_state, options, true, false);
            ++trajectory_count;
            if (trial.valid &&
                trial.objective.total >
                    current.objective.total +
                        options.armijo_fraction * line_angle *
                            gradient_slope &&
                trial.score.total > current.score.total) {
                accepted = true;
                accepted_state = trial_state;
                accepted_value = trial;
                break;
            }
            line_angle *= 0.5;
        }
        if (!accepted) {
            writeTraceRow(trace,
                          iteration,
                          "line-search-rejected",
                          "coarse",
                          current,
                          options,
                          gradient_norm,
                          gradient_relative_error,
                          gradient_slope,
                          0.0,
                          trajectory_count);
            std::cout << "iteration " << iteration
                      << " stopped: no valid jointly improving step\n";
            break;
        }
        current_state = accepted_state;
        current = accepted_value;
        ++accepted_steps;
        writeTraceRow(trace,
                      iteration,
                      "accepted",
                      "coarse",
                      current,
                      options,
                      gradient_norm,
                      gradient_relative_error,
                      gradient_slope,
                      line_angle,
                      trajectory_count);
        const ns_cascade::OptimizationState accepted_fine_state =
            ns_cascade::liftOptimizationState(
                coarse_system, current_state, fine_system);
        const Evaluation accepted_fine = evaluateTrajectory(
            fine_system, accepted_fine_state, options, false, false);
        ++trajectory_count;
        const PairAssessment accepted_pair = assessPair(
            current, accepted_fine, options);
        writeTraceRow(trace,
                      iteration,
                      "fine-check",
                      "fine",
                      accepted_fine,
                      options,
                      nan,
                      nan,
                      nan,
                      0.0,
                      trajectory_count,
                      &accepted_pair);
        const bool accepted_pair_preferred =
            (accepted_pair.preliminary_cross_resolution_ok &&
             !best_pair.preliminary_cross_resolution_ok) ||
            (accepted_pair.preliminary_cross_resolution_ok ==
                 best_pair.preliminary_cross_resolution_ok &&
             accepted_pair.score.total > best_pair.score.total);
        if (accepted_pair_preferred) {
            best_state = current_state;
            best_coarse = current;
            best_fine = accepted_fine;
            best_pair = accepted_pair;
            best_iteration = iteration;
        }
        std::cout << std::setprecision(10)
                  << "iteration " << iteration
                  << " objective=" << current.objective.total
                  << " H1/2="
                  << ratio(current.peak_h_half, current.initial_h_half)
                  << " k_rms="
                  << ratio(current.final_characteristic_wavenumber,
                           current.initial_characteristic_wavenumber)
                  << " profile-drift=" << current.latest_profile_drift
                  << " cutoff=" << current.peak_cutoff_fraction
                  << " pair-score=" << accepted_pair.score.total
                  << " selected="
                  << (best_iteration == iteration ? "yes" : "no")
                  << " gradient-error=" << gradient_relative_error << '\n';
    }

    writeTraceRow(trace,
                  best_iteration,
                  "selected",
                  "coarse",
                  best_coarse,
                  options,
                  nan,
                  nan,
                  nan,
                  0.0,
                  trajectory_count,
                  &best_pair);
    writeTraceRow(trace,
                  best_iteration,
                  "fine-validation",
                  "fine",
                  best_fine,
                  options,
                  nan,
                  nan,
                  nan,
                  0.0,
                  trajectory_count,
                  &best_pair);
    if (!trace) {
        throw std::runtime_error("Failed while writing state-optimizer CSV");
    }
    writeState(options.state_output,
               coarse_system,
               best_state,
               options);

    std::cout << std::setprecision(12)
              << "Adjoint Fourier-state optimization complete\n"
              << "  family/bandwidth/real-dofs: "
              << seedFamilyName(options.seed_family) << '/'
              << options.seed_bandwidth << '/'
              << ns_cascade::stateOptimizationDegreesOfFreedom(
                     options.seed_bandwidth)
              << '\n'
              << "  accepted steps/selected iteration/trajectories: "
              << accepted_steps << '/' << best_iteration << '/'
              << trajectory_count << '\n'
              << "  coarse objective/H1/2/k_rms/drift/cutoff: "
              << best_coarse.objective.total << '/'
              << ratio(best_coarse.peak_h_half, best_coarse.initial_h_half)
              << '/'
              << ratio(best_coarse.final_characteristic_wavenumber,
                       best_coarse.initial_characteristic_wavenumber)
              << '/' << best_coarse.latest_profile_drift << '/'
              << best_coarse.peak_cutoff_fraction << '\n'
              << "  fine objective/H1/2/k_rms/drift/cutoff: "
              << best_fine.objective.total << '/'
              << ratio(best_fine.peak_h_half, best_fine.initial_h_half) << '/'
              << ratio(best_fine.final_characteristic_wavenumber,
                       best_fine.initial_characteristic_wavenumber)
              << '/' << best_fine.latest_profile_drift << '/'
              << best_fine.peak_cutoff_fraction << '\n'
              << "  cross-resolution/refinement-eligible: "
              << (best_pair.preliminary_cross_resolution_ok ? "yes" : "no")
              << '/'
              << (best_pair.score.refinement_eligible ? "yes" : "no")
              << '\n'
              << "  outputs: " << options.output << " and "
              << options.state_output << '\n'
              << "This constrained floating-point candidate is not a PDE "
                 "proof and cannot bypass the resolution/profile gates.\n";
    return 0;
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
