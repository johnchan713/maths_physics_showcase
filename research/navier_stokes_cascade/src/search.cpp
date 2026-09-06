#include "ns_cascade/candidate_score.hpp"
#include "ns_cascade/pseudospectral.hpp"
#include "ns_cascade/spectral_profile.hpp"

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

enum class VortexFamily {
    Pair,
    OrthogonalBundle,
    WavePackets
};

const char* vortexFamilyName(VortexFamily family) {
    if (family == VortexFamily::Pair) return "pair";
    if (family == VortexFamily::OrthogonalBundle) return "orthogonal-bundle";
    return "wave-packets";
}

VortexFamily parseVortexFamily(const std::string& value) {
    if (value == "pair") return VortexFamily::Pair;
    if (value == "orthogonal-bundle" || value == "bundle") {
        return VortexFamily::OrthogonalBundle;
    }
    if (value == "wave-packets" || value == "packets") {
        return VortexFamily::WavePackets;
    }
    throw std::invalid_argument(
        "Unknown vortex family: " + value +
        " (expected pair, orthogonal-bundle, or wave-packets)");
}

std::vector<VortexFamily> parseVortexFamilies(const std::string& text) {
    std::vector<VortexFamily> families;
    std::istringstream stream(text);
    std::string item;
    while (std::getline(stream, item, ',')) {
        if (item.empty()) {
            throw std::invalid_argument("Empty entry in --families");
        }
        families.push_back(parseVortexFamily(item));
    }
    if (families.empty()) {
        throw std::invalid_argument("--families cannot be empty");
    }
    return families;
}

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
    std::vector<double> initial_energies = {1.0};
    double cutoff_fraction_threshold = 0.01;
    int profile_bin_count = 32;
    double profile_maximum_coordinate = 4.0;
    double profile_log_scale_window = 0.025;
    double profile_drift_threshold = 1.0;
    double minimum_characteristic_growth = 1.10;
    double minimum_critical_growth = 1.005;
    std::vector<VortexFamily> families = {VortexFamily::Pair};
    std::vector<double> core_radii = {0.55, 0.70};
    std::vector<double> separations = {1.2, 1.8};
    std::vector<double> bend_amplitudes = {0.0, 0.30};
    std::vector<int> axial_wavenumbers = {1, 2};
    std::vector<double> orthogonal_pair_weights = {0.5, 1.0};
    std::vector<double> phase_offsets = {
        0.0, 1.0471975511965977461542144610932};
    std::vector<double> packet_widths = {0.8, 1.1};
    std::vector<int> carrier_wavenumbers = {1, 2};
    std::vector<double> packet_secondary_weights = {0.75, 1.0};
    std::vector<double> packet_phase_offsets = {
        0.0, 1.0471975511965977461542144610932};
    std::string output = "navier_stokes_candidate_search.csv";
};

struct Candidate {
    int id;
    VortexFamily family;
    ns_cascade::VortexTubeParameters tube_parameters;
    double orthogonal_pair_weight;
    double phase_offset;
    ns_cascade::WavePacketParameters packet_parameters;
    double initial_energy;
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
    double initial_h_half;
    double final_h_half;
    double peak_h_half;
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
    double maximum_nonlinear_enstrophy_production;
    double maximum_enstrophy_production_to_dissipation;
    double maximum_net_enstrophy_rate;
    double initial_characteristic_wavenumber;
    double final_characteristic_wavenumber;
    int profile_scale_windows;
    double first_profile_drift;
    double latest_profile_drift;
    double minimum_profile_drift;
    double latest_profile_window_time;
    bool profile_stationarity_improving;
    double profile_drift_cost;
    double profile_trend_reward;
    double profile_rebound_cost;
    double cutoff_cost;
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
          initial_h_half(std::numeric_limits<double>::quiet_NaN()),
          final_h_half(std::numeric_limits<double>::quiet_NaN()),
          peak_h_half(std::numeric_limits<double>::quiet_NaN()),
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
          maximum_nonlinear_enstrophy_production(
              std::numeric_limits<double>::quiet_NaN()),
          maximum_enstrophy_production_to_dissipation(
              std::numeric_limits<double>::quiet_NaN()),
          maximum_net_enstrophy_rate(
              std::numeric_limits<double>::quiet_NaN()),
          initial_characteristic_wavenumber(
              std::numeric_limits<double>::quiet_NaN()),
          final_characteristic_wavenumber(
              std::numeric_limits<double>::quiet_NaN()),
          profile_scale_windows(0),
          first_profile_drift(std::numeric_limits<double>::quiet_NaN()),
          latest_profile_drift(std::numeric_limits<double>::quiet_NaN()),
          minimum_profile_drift(std::numeric_limits<double>::quiet_NaN()),
          latest_profile_window_time(
              std::numeric_limits<double>::quiet_NaN()),
          profile_stationarity_improving(false),
          profile_drift_cost(std::numeric_limits<double>::quiet_NaN()),
          profile_trend_reward(std::numeric_limits<double>::quiet_NaN()),
          profile_rebound_cost(std::numeric_limits<double>::quiet_NaN()),
          cutoff_cost(std::numeric_limits<double>::quiet_NaN()),
          score(-std::numeric_limits<double>::infinity()),
          elapsed_seconds(std::numeric_limits<double>::quiet_NaN()) {}
};

struct PairAssessment {
    bool preliminary_cross_resolution_ok = false;
    double l3_ratio_relative_difference =
        std::numeric_limits<double>::quiet_NaN();
    double h_half_ratio_relative_difference =
        std::numeric_limits<double>::quiet_NaN();
    double vorticity_ratio_relative_difference =
        std::numeric_limits<double>::quiet_NaN();
    double flux_relative_difference =
        std::numeric_limits<double>::quiet_NaN();
    double characteristic_scale_relative_difference =
        std::numeric_limits<double>::quiet_NaN();
    double profile_drift_relative_difference =
        std::numeric_limits<double>::quiet_NaN();
    bool profile_drift_comparison_valid = false;
    bool profile_windows_recent = false;
    ns_cascade::PairedSearchScore score;
};

struct FinalistResult {
    RunResult coarse;
    RunResult fine;
    PairAssessment assessment;

    FinalistResult(const RunResult& coarse_value,
                   const RunResult& fine_value,
                   const PairAssessment& assessment_value)
        : coarse(coarse_value),
          fine(fine_value),
          assessment(assessment_value) {}
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
        << "Coarse-to-fine search over smooth periodic cascade candidates.\n"
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
        << "  --energy E              One normalized initial energy (default: 1)\n"
        << "  --energies LIST         Initial energies to search\n"
        << "  --families LIST         pair,orthogonal-bundle,wave-packets\n"
        << "  --cores LIST            Core radii (default: 0.55,0.70)\n"
        << "  --separations LIST      Pair separations (default: 1.2,1.8)\n"
        << "  --bends LIST            Helical bend amplitudes (default: 0,0.30)\n"
        << "  --axial-modes LIST      Positive axial modes (default: 1,2)\n"
        << "  --orthogonal-weights LIST  Relative x/y-pair weights for bundles\n"
        << "  --phase-offsets LIST    Bundle helical phase offsets in radians\n"
        << "  --packet-widths LIST    Periodic Gaussian envelope widths\n"
        << "  --carrier-modes LIST    Wave-packet carrier wavenumbers\n"
        << "  --packet-weights LIST   Relative weights of packets two and three\n"
        << "  --packet-phases LIST    Wave-packet triad phase offsets\n"
        << "  --cutoff-threshold F    Max cutoff-shell energy fraction (default: 0.01)\n"
        << "  --profile-bins N        Rescaled-spectrum bins before overflow (default: 32)\n"
        << "  --profile-max-xi X      Last finite |k|/k_rms coordinate (default: 4)\n"
        << "  --profile-scale-window W  Forward log(k_rms) per drift window (default: 0.025)\n"
        << "  --profile-drift-threshold D  Max refinement drift (default: 1)\n"
        << "  --minimum-scale-growth G  Min k_rms ratio for refinement (default: 1.10)\n"
        << "  --minimum-critical-growth G  Min critical-norm ratio (default: 1.005)\n"
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
            options.initial_energies.assign(
                1, parseNumber<double>(requireValue(i, argc, argv), flag));
        } else if (flag == "--energies") {
            options.initial_energies =
                parseList<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--families") {
            options.families =
                parseVortexFamilies(requireValue(i, argc, argv));
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
        } else if (flag == "--orthogonal-weights") {
            options.orthogonal_pair_weights =
                parseList<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--phase-offsets") {
            options.phase_offsets =
                parseList<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--packet-widths") {
            options.packet_widths =
                parseList<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--carrier-modes") {
            options.carrier_wavenumbers =
                parseList<int>(requireValue(i, argc, argv), flag);
        } else if (flag == "--packet-weights") {
            options.packet_secondary_weights =
                parseList<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--packet-phases") {
            options.packet_phase_offsets =
                parseList<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--cutoff-threshold") {
            options.cutoff_fraction_threshold =
                parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--profile-bins") {
            options.profile_bin_count =
                parseNumber<int>(requireValue(i, argc, argv), flag);
        } else if (flag == "--profile-max-xi") {
            options.profile_maximum_coordinate =
                parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--profile-scale-window") {
            options.profile_log_scale_window =
                parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--profile-drift-threshold") {
            options.profile_drift_threshold =
                parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--minimum-scale-growth") {
            options.minimum_characteristic_growth =
                parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--minimum-critical-growth") {
            options.minimum_critical_growth =
                parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--output") {
            options.output = requireValue(i, argc, argv);
        } else {
            throw std::invalid_argument("Unknown option: " + flag);
        }
    }

    removeDuplicates(options.families);
    removeDuplicates(options.core_radii);
    removeDuplicates(options.separations);
    removeDuplicates(options.bend_amplitudes);
    removeDuplicates(options.axial_wavenumbers);
    removeDuplicates(options.orthogonal_pair_weights);
    removeDuplicates(options.phase_offsets);
    removeDuplicates(options.packet_widths);
    removeDuplicates(options.carrier_wavenumbers);
    removeDuplicates(options.packet_secondary_weights);
    removeDuplicates(options.packet_phase_offsets);
    removeDuplicates(options.initial_energies);
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
    for (std::size_t i = 0; i < options.initial_energies.size(); ++i) {
        if (!finiteAndPositive(options.initial_energies[i])) {
            throw std::invalid_argument(
                "Every requested initial energy must be finite and positive");
        }
    }
    if (!finiteAndPositive(options.cutoff_fraction_threshold) ||
        options.cutoff_fraction_threshold >= 1.0) {
        throw std::invalid_argument(
            "--cutoff-threshold must be finite and in (0,1)");
    }
    if (options.profile_bin_count < 4 || options.profile_bin_count > 4096) {
        throw std::invalid_argument("--profile-bins must be between 4 and 4096");
    }
    if (!finiteAndPositive(options.profile_maximum_coordinate)) {
        throw std::invalid_argument(
            "--profile-max-xi must be finite and positive");
    }
    if (!finiteAndPositive(options.profile_log_scale_window) ||
        options.profile_log_scale_window >= 1.0) {
        throw std::invalid_argument(
            "--profile-scale-window must be finite and in (0,1)");
    }
    if (!finiteAndPositive(options.profile_drift_threshold)) {
        throw std::invalid_argument(
            "--profile-drift-threshold must be finite and positive");
    }
    if (!std::isfinite(options.minimum_characteristic_growth) ||
        options.minimum_characteristic_growth <= 1.0) {
        throw std::invalid_argument(
            "--minimum-scale-growth must be finite and greater than one");
    }
    if (!std::isfinite(options.minimum_critical_growth) ||
        options.minimum_critical_growth <= 1.0) {
        throw std::invalid_argument(
            "--minimum-critical-growth must be finite and greater than one");
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
    for (std::size_t i = 0; i < options.orthogonal_pair_weights.size(); ++i) {
        if (!finiteAndPositive(options.orthogonal_pair_weights[i])) {
            throw std::invalid_argument(
                "Every orthogonal-pair weight must be finite and positive");
        }
    }
    for (std::size_t i = 0; i < options.phase_offsets.size(); ++i) {
        if (!std::isfinite(options.phase_offsets[i])) {
            throw std::invalid_argument(
                "Every bundle phase offset must be finite");
        }
    }
    const bool uses_wave_packets =
        std::find(options.families.begin(),
                  options.families.end(),
                  VortexFamily::WavePackets) != options.families.end();
    for (std::size_t i = 0; uses_wave_packets &&
         i < options.packet_widths.size(); ++i) {
        if (!std::isfinite(options.packet_widths[i]) ||
            options.packet_widths[i] <= 0.0 || options.packet_widths[i] > pi) {
            throw std::invalid_argument(
                "Every packet width must be finite and in (0,pi]");
        }
    }
    for (std::size_t i = 0; uses_wave_packets &&
         i < options.carrier_wavenumbers.size(); ++i) {
        if (options.carrier_wavenumbers[i] < 1 ||
            options.carrier_wavenumbers[i] > coarse_safe_cutoff / 2) {
            throw std::invalid_argument(
                "Every carrier mode must fit below half the coarse cutoff");
        }
    }
    for (std::size_t i = 0; uses_wave_packets &&
         i < options.packet_secondary_weights.size(); ++i) {
        if (!finiteAndPositive(options.packet_secondary_weights[i])) {
            throw std::invalid_argument(
                "Every packet weight must be finite and positive");
        }
    }
    for (std::size_t i = 0; uses_wave_packets &&
         i < options.packet_phase_offsets.size(); ++i) {
        if (!std::isfinite(options.packet_phase_offsets[i])) {
            throw std::invalid_argument(
                "Every packet phase must be finite");
        }
    }
    return options;
}

std::vector<Candidate> buildCandidates(const Options& options) {
    std::vector<Candidate> candidates;
    int id = 1;
    for (std::size_t energy = 0;
         energy < options.initial_energies.size();
         ++energy) {
        for (std::size_t family = 0;
             family < options.families.size();
             ++family) {
            if (options.families[family] == VortexFamily::WavePackets) {
                for (std::size_t width = 0;
                     width < options.packet_widths.size(); ++width) {
                    for (std::size_t carrier = 0;
                         carrier < options.carrier_wavenumbers.size(); ++carrier) {
                        for (std::size_t weight = 0;
                             weight < options.packet_secondary_weights.size();
                             ++weight) {
                            for (std::size_t phase = 0;
                                 phase < options.packet_phase_offsets.size();
                                 ++phase) {
                                Candidate candidate;
                                candidate.id = id++;
                                candidate.family = VortexFamily::WavePackets;
                                candidate.tube_parameters =
                                    ns_cascade::VortexTubeParameters();
                                candidate.orthogonal_pair_weight = 0.0;
                                candidate.phase_offset = 0.0;
                                candidate.packet_parameters =
                                    ns_cascade::WavePacketParameters(
                                        options.packet_widths[width],
                                        options.carrier_wavenumbers[carrier],
                                        options.packet_secondary_weights[weight],
                                        options.packet_phase_offsets[phase]);
                                candidate.initial_energy =
                                    options.initial_energies[energy];
                                candidates.push_back(candidate);
                            }
                        }
                    }
                }
                continue;
            }
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
                            // Axial mode and phase have no effect on straight tubes.
                            if (options.bend_amplitudes[bend] == 0.0 && mode != 0) {
                                continue;
                            }
                            const std::size_t weight_count =
                                options.families[family] == VortexFamily::Pair
                                    ? 1U
                                    : options.orthogonal_pair_weights.size();
                            const std::size_t phase_count =
                                options.families[family] == VortexFamily::Pair ||
                                        options.bend_amplitudes[bend] == 0.0
                                    ? 1U
                                    : options.phase_offsets.size();
                            for (std::size_t weight = 0;
                                 weight < weight_count;
                                 ++weight) {
                                for (std::size_t phase = 0;
                                     phase < phase_count;
                                     ++phase) {
                                    Candidate candidate;
                                    candidate.id = id++;
                                    candidate.family = options.families[family];
                                    candidate.initial_energy =
                                        options.initial_energies[energy];
                                    candidate.tube_parameters =
                                        ns_cascade::VortexTubeParameters(
                                            options.core_radii[core],
                                            options.separations[separation],
                                            options.bend_amplitudes[bend],
                                            options.axial_wavenumbers[mode]);
                                    candidate.orthogonal_pair_weight =
                                        candidate.family == VortexFamily::Pair
                                            ? 0.0
                                            : options.orthogonal_pair_weights[weight];
                                    candidate.phase_offset =
                                        candidate.family == VortexFamily::Pair
                                            ? 0.0
                                            : options.phase_offsets[phase];
                                    candidate.packet_parameters =
                                        ns_cascade::WavePacketParameters();
                                    candidates.push_back(candidate);
                                }
                            }
                        }
                    }
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
           std::isfinite(values.critical_h_half) &&
           std::isfinite(values.critical_l3_sample) &&
           std::isfinite(values.sampled_vorticity_max) &&
           std::isfinite(values.vorticity_sup_upper_bound) &&
           std::isfinite(values.spectral_centroid) &&
           std::isfinite(values.high_shell_energy_fraction) &&
           std::isfinite(values.divergence_defect) &&
           std::isfinite(values.reality_defect) &&
           std::isfinite(values.nonlinear_enstrophy_production) &&
           std::isfinite(values.viscous_enstrophy_destruction) &&
           std::isfinite(values.net_enstrophy_rate) &&
           std::isfinite(values.enstrophy_production_to_dissipation) &&
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
    result.final_h_half = values.critical_h_half;
    result.final_vorticity = values.sampled_vorticity_max;
    result.final_enstrophy = values.enstrophy;
    if (values.energy <= 0.0 || values.enstrophy <= 0.0) {
        throw std::runtime_error(
            "Candidate lost positive energy or enstrophy");
    }
    result.final_characteristic_wavenumber =
        std::sqrt(values.enstrophy / values.energy);
    result.peak_l3 = std::max(result.peak_l3, values.critical_l3_sample);
    result.peak_h_half =
        std::max(result.peak_h_half, values.critical_h_half);
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
    result.maximum_nonlinear_enstrophy_production = std::max(
        result.maximum_nonlinear_enstrophy_production,
        values.nonlinear_enstrophy_production);
    result.maximum_enstrophy_production_to_dissipation = std::max(
        result.maximum_enstrophy_production_to_dissipation,
        values.enstrophy_production_to_dissipation);
    result.maximum_net_enstrophy_rate = std::max(
        result.maximum_net_enstrophy_rate,
        values.net_enstrophy_rate);
    result.maximum_positive_forward_flux = std::max(
        result.maximum_positive_forward_flux,
        positiveForwardFlux(system.shellDiagnostics(state)));
}

double ratio(double peak, double initial) {
    return initial == 0.0 ? 0.0 : peak / initial;
}

ns_cascade::SearchScoreEvidence scoreEvidence(const RunResult& result,
                                               const Options& options) {
    ns_cascade::SearchScoreEvidence evidence;
    evidence.peak_l3_ratio = ratio(result.peak_l3, result.initial_l3);
    evidence.peak_h_half_ratio =
        ratio(result.peak_h_half, result.initial_h_half);
    evidence.final_l3_ratio = ratio(result.final_l3, result.initial_l3);
    evidence.peak_vorticity_ratio =
        ratio(result.peak_vorticity, result.initial_vorticity);
    evidence.peak_cutoff_fraction = result.peak_cutoff_fraction;
    evidence.cutoff_fraction_threshold = options.cutoff_fraction_threshold;
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

RunResult runCandidate(const ns_cascade::PseudospectralSystem& system,
                       const Candidate& candidate,
                       const Options& options) {
    RunResult result(candidate, system.gridSize(), system.cutoff());
    const std::chrono::steady_clock::time_point start =
        std::chrono::steady_clock::now();
    ns_cascade::PseudospectralSystem::State state;
    if (candidate.family == VortexFamily::Pair) {
        state = system.vortexTubePairState(
            candidate.tube_parameters, candidate.initial_energy);
    } else if (candidate.family == VortexFamily::OrthogonalBundle) {
        state = system.vortexBundleState(
            ns_cascade::VortexBundleParameters(
                candidate.tube_parameters,
                candidate.orthogonal_pair_weight,
                candidate.phase_offset),
            candidate.initial_energy);
    } else {
        state = system.interactingWavePacketState(
            candidate.packet_parameters, candidate.initial_energy);
    }
    const ns_cascade::PseudospectralSystem::Diagnostics initial =
        system.diagnostics(state);
    if (!finiteDiagnostics(initial)) {
        throw std::runtime_error("Initial diagnostics are non-finite");
    }
    const ns_cascade::SpectrumProfile initial_profile =
        ns_cascade::rescaledSpectrumProfile(
            system,
            state,
            options.profile_bin_count,
            options.profile_maximum_coordinate);
    ns_cascade::SpectrumProfile profile_anchor = initial_profile;

    result.initial_three_dimensional_energy_fraction =
        threeDimensionalEnergyFraction(system, state);
    result.initial_energy = initial.energy;
    result.final_energy = initial.energy;
    result.initial_l3 = initial.critical_l3_sample;
    result.final_l3 = initial.critical_l3_sample;
    result.peak_l3 = initial.critical_l3_sample;
    result.initial_h_half = initial.critical_h_half;
    result.final_h_half = initial.critical_h_half;
    result.peak_h_half = initial.critical_h_half;
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
    result.maximum_nonlinear_enstrophy_production =
        -std::numeric_limits<double>::infinity();
    result.maximum_enstrophy_production_to_dissipation =
        -std::numeric_limits<double>::infinity();
    result.maximum_net_enstrophy_rate =
        -std::numeric_limits<double>::infinity();
    result.initial_characteristic_wavenumber =
        initial_profile.characteristic_wavenumber;
    result.final_characteristic_wavenumber =
        initial_profile.characteristic_wavenumber;
    result.minimum_profile_drift =
        std::numeric_limits<double>::infinity();
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

        const double scale_check_energy = system.energy(state);
        const double scale_check_enstrophy = system.enstrophy(state);
        const double scale_check_cutoff_fraction =
            system.cutoffShellEnergyFraction(state);
        if (!finiteAndPositive(scale_check_energy) ||
            !finiteAndPositive(scale_check_enstrophy) ||
            !std::isfinite(scale_check_cutoff_fraction) ||
            scale_check_cutoff_fraction < 0.0 ||
            scale_check_cutoff_fraction > 1.0 + 1e-12) {
            throw std::runtime_error(
                "Per-step spectral gate became invalid");
        }
        result.peak_cutoff_fraction = std::max(
            result.peak_cutoff_fraction, scale_check_cutoff_fraction);
        const double scale_check_wavenumber =
            std::sqrt(scale_check_enstrophy / scale_check_energy);
        const int completed_profile_windows =
            ns_cascade::completedForwardProfileScaleWindows(
                scale_check_wavenumber,
                initial_profile.characteristic_wavenumber,
                options.profile_log_scale_window);
        if (completed_profile_windows > result.profile_scale_windows) {
            const ns_cascade::SpectrumProfile current_profile =
                ns_cascade::rescaledSpectrumProfile(
                    system,
                    state,
                    options.profile_bin_count,
                    options.profile_maximum_coordinate);
            result.final_characteristic_wavenumber =
                current_profile.characteristic_wavenumber;
            const ns_cascade::SpectrumProfileChange change =
                ns_cascade::spectrumProfileChange(
                    current_profile, profile_anchor);
            if (!change.scale_normalized_drift_valid) {
                throw std::runtime_error(
                    "A completed profile scale window has invalid drift");
            }
            if (result.profile_scale_windows == 0) {
                result.first_profile_drift =
                    change.l1_per_log_scale_change;
            }
            result.profile_scale_windows = completed_profile_windows;
            result.latest_profile_drift =
                change.l1_per_log_scale_change;
            result.minimum_profile_drift = std::min(
                result.minimum_profile_drift,
                change.l1_per_log_scale_change);
            result.latest_profile_window_time = time;
            profile_anchor = current_profile;
        }

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
    const ns_cascade::SearchScore score =
        ns_cascade::scoreSingleResolution(scoreEvidence(result, options));
    result.profile_stationarity_improving =
        score.profile_stationarity_improving;
    result.profile_drift_cost = score.profile_drift_cost;
    result.profile_trend_reward = score.profile_trend_reward;
    result.profile_rebound_cost = score.profile_rebound_cost;
    result.cutoff_cost = score.cutoff_cost;
    result.score = score.total;
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
                       double& h_half_difference,
                       double& vorticity_difference,
                       double& flux_difference) {
    const double coarse_l3_ratio = ratio(coarse.peak_l3, coarse.initial_l3);
    const double fine_l3_ratio = ratio(fine.peak_l3, fine.initial_l3);
    const double coarse_h_half_ratio =
        ratio(coarse.peak_h_half, coarse.initial_h_half);
    const double fine_h_half_ratio =
        ratio(fine.peak_h_half, fine.initial_h_half);
    const double coarse_vorticity_ratio =
        ratio(coarse.peak_vorticity, coarse.initial_vorticity);
    const double fine_vorticity_ratio =
        ratio(fine.peak_vorticity, fine.initial_vorticity);
    l3_difference = relativeDifference(coarse_l3_ratio, fine_l3_ratio);
    h_half_difference = relativeDifference(
        coarse_h_half_ratio, fine_h_half_ratio);
    vorticity_difference = relativeDifference(
        coarse_vorticity_ratio, fine_vorticity_ratio);
    flux_difference = relativeDifference(
        coarse.maximum_positive_forward_flux,
        fine.maximum_positive_forward_flux);
    const bool l3_growth_signal_agrees =
        std::max(coarse_l3_ratio, fine_l3_ratio) <= 1.005 ||
        (coarse_l3_ratio > 1.005 && fine_l3_ratio > 1.005);
    const bool h_half_growth_signal_agrees =
        std::max(coarse_h_half_ratio, fine_h_half_ratio) <= 1.005 ||
        (coarse_h_half_ratio > 1.005 && fine_h_half_ratio > 1.005);
    const bool vorticity_growth_signal_agrees =
        std::max(coarse_vorticity_ratio, fine_vorticity_ratio) <= 1.005 ||
        (coarse_vorticity_ratio > 1.005 && fine_vorticity_ratio > 1.005);
    return coarse.completed && fine.completed && coarse.resolution_ok &&
           fine.resolution_ok && l3_difference <= 0.02 &&
           h_half_difference <= 0.02 &&
           vorticity_difference <= 0.10 && flux_difference <= 0.25 &&
           l3_growth_signal_agrees && h_half_growth_signal_agrees &&
           vorticity_growth_signal_agrees;
}

PairAssessment assessResolutionPair(const RunResult& coarse,
                                    const RunResult& fine,
                                    const Options& options) {
    PairAssessment assessment;
    if (!coarse.completed || !fine.completed) return assessment;

    assessment.preliminary_cross_resolution_ok = crossResolutionOk(
        coarse,
        fine,
        assessment.l3_ratio_relative_difference,
        assessment.h_half_ratio_relative_difference,
        assessment.vorticity_ratio_relative_difference,
        assessment.flux_relative_difference);
    const double coarse_characteristic_growth = ratio(
        coarse.final_characteristic_wavenumber,
        coarse.initial_characteristic_wavenumber);
    const double fine_characteristic_growth = ratio(
        fine.final_characteristic_wavenumber,
        fine.initial_characteristic_wavenumber);
    assessment.characteristic_scale_relative_difference = relativeDifference(
        coarse_characteristic_growth, fine_characteristic_growth);
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
        assessment.profile_drift_comparison_valid
            ? assessment.profile_drift_relative_difference
            : 1.0;
    evidence.profile_drift_comparison_valid =
        assessment.profile_drift_comparison_valid;
    evidence.coarse_profile_window_recent =
        coarse.latest_profile_window_time >= 0.75 * options.final_time;
    evidence.fine_profile_window_recent =
        fine.latest_profile_window_time >= 0.75 * options.final_time;
    evidence.coarse_characteristic_growth = coarse_characteristic_growth;
    evidence.fine_characteristic_growth = fine_characteristic_growth;
    evidence.minimum_characteristic_growth =
        options.minimum_characteristic_growth;
    evidence.minimum_critical_growth = options.minimum_critical_growth;
    assessment.score = ns_cascade::scoreResolutionPair(evidence);
    return assessment;
}

bool betterFinalist(const FinalistResult& left,
                    const FinalistResult& right) {
    if (left.fine.completed != right.fine.completed) {
        return left.fine.completed;
    }
    if (left.assessment.score.refinement_eligible !=
        right.assessment.score.refinement_eligible) {
        return left.assessment.score.refinement_eligible;
    }
    if (left.assessment.preliminary_cross_resolution_ok !=
        right.assessment.preliminary_cross_resolution_ok) {
        return left.assessment.preliminary_cross_resolution_ok;
    }
    if (left.assessment.score.total != right.assessment.score.total) {
        return left.assessment.score.total > right.assessment.score.total;
    }
    return left.fine.candidate.id < right.fine.candidate.id;
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
        << "initial_family,orthogonal_pair_weight,phase_offset,"
        << "packet_width,carrier_wavenumber,packet_secondary_weight,packet_phase,"
        << "core_radius,separation,bend_amplitude,axial_wavenumber,viscosity,"
        << "requested_initial_energy,final_time,steps,min_dt,max_dt,"
        << "max_advective_cfl_upper_bound,max_viscous_stability_number,"
        << "initial_3d_energy_fraction,initial_energy,final_energy,"
        << "initial_l3,final_l3,peak_l3,peak_l3_ratio,"
        << "initial_h_half,final_h_half,peak_h_half,peak_h_half_ratio,"
        << "initial_sampled_vorticity,final_sampled_vorticity,"
        << "peak_sampled_vorticity,peak_vorticity_ratio,"
        << "initial_enstrophy,final_enstrophy,peak_enstrophy,"
        << "peak_enstrophy_ratio,peak_palinstrophy,max_positive_forward_flux,"
        << "peak_cutoff_energy_fraction,max_divergence_defect,"
        << "max_reality_defect,max_energy_balance_residual,resolution_ok,"
        << "max_nonlinear_enstrophy_production,"
        << "max_enstrophy_production_to_dissipation,max_net_enstrophy_rate,"
        << "initial_characteristic_wavenumber,final_characteristic_wavenumber,"
        << "characteristic_wavenumber_ratio,profile_scale_windows,"
        << "first_profile_drift,latest_profile_drift,minimum_profile_drift,"
        << "latest_profile_window_time,profile_stationarity_improving,"
        << "profile_drift_cost,profile_trend_reward,profile_rebound_cost,"
        << "cutoff_cost,"
        << "ranking_score,cpu_seconds,reference_grid,"
        << "l3_ratio_relative_difference,h_half_ratio_relative_difference,"
        << "vorticity_ratio_relative_difference,"
        << "flux_relative_difference,characteristic_scale_relative_difference,"
        << "profile_drift_relative_difference,"
        << "profile_drift_comparison_valid,profile_windows_recent,"
        << "multi_resolution_cutoff_cost,"
        << "paired_ranking_score,refinement_eligible,"
        << "cross_resolution_ok,failure\n";
}

void writeResult(std::ostream& output,
                 const std::string& stage,
                 int selection_rank,
                 const RunResult& result,
                 const Options& options,
                 const RunResult* reference,
                 const PairAssessment* assessment) {
    output << stage << ',' << selection_rank << ',' << result.candidate.id << ','
           << (result.completed ? "completed" : "failed") << ','
           << result.grid_size << ',' << result.cutoff << ','
           << vortexFamilyName(result.candidate.family) << ','
           << result.candidate.orthogonal_pair_weight << ','
           << result.candidate.phase_offset << ','
           << result.candidate.packet_parameters.envelope_width << ','
           << result.candidate.packet_parameters.carrier_wavenumber << ','
           << result.candidate.packet_parameters.secondary_weight << ','
           << result.candidate.packet_parameters.phase_offset << ','
           << result.candidate.tube_parameters.core_radius << ','
           << result.candidate.tube_parameters.separation << ','
           << result.candidate.tube_parameters.bend_amplitude << ','
           << result.candidate.tube_parameters.axial_wavenumber << ','
           << options.viscosity << ',' << result.candidate.initial_energy << ','
           << options.final_time << ',' << result.steps << ','
           << result.minimum_time_step << ',' << result.maximum_time_step << ','
           << result.maximum_cfl_bound << ','
           << result.maximum_viscous_number << ','
           << result.initial_three_dimensional_energy_fraction << ','
           << result.initial_energy << ',' << result.final_energy << ','
           << result.initial_l3 << ',' << result.final_l3 << ','
           << result.peak_l3 << ','
           << ratio(result.peak_l3, result.initial_l3) << ','
           << result.initial_h_half << ',' << result.final_h_half << ','
           << result.peak_h_half << ','
           << ratio(result.peak_h_half, result.initial_h_half) << ','
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
           << result.maximum_nonlinear_enstrophy_production << ','
           << result.maximum_enstrophy_production_to_dissipation << ','
           << result.maximum_net_enstrophy_rate << ','
           << result.initial_characteristic_wavenumber << ','
           << result.final_characteristic_wavenumber << ','
           << ratio(result.final_characteristic_wavenumber,
                    result.initial_characteristic_wavenumber)
           << ',' << result.profile_scale_windows << ','
           << result.first_profile_drift << ','
           << result.latest_profile_drift << ','
           << result.minimum_profile_drift << ','
           << result.latest_profile_window_time << ','
           << (result.profile_stationarity_improving ? "true" : "false")
           << ',' << result.profile_drift_cost << ','
           << result.profile_trend_reward << ','
           << result.profile_rebound_cost << ',' << result.cutoff_cost << ','
           << result.score << ',' << result.elapsed_seconds << ',';
    if (reference == NULL) {
        output << ",,,,,,,,,,,,,";
    } else {
        if (assessment == NULL) {
            throw std::logic_error(
                "A fine-grid CSV row requires a pair assessment");
        }
        output << reference->grid_size << ','
               << assessment->l3_ratio_relative_difference << ','
               << assessment->h_half_ratio_relative_difference << ','
               << assessment->vorticity_ratio_relative_difference << ','
               << assessment->flux_relative_difference << ','
               << assessment->characteristic_scale_relative_difference << ','
               << assessment->profile_drift_relative_difference << ','
               << (assessment->profile_drift_comparison_valid ? "true"
                                                               : "false")
               << ',' << (assessment->profile_windows_recent ? "true" : "false")
               << ',' << assessment->score.multi_resolution_cutoff_cost << ','
               << assessment->score.total << ','
               << (assessment->score.refinement_eligible ? "true" : "false")
               << ','
               << (assessment->preliminary_cross_resolution_ok ? "true"
                                                                : "false")
               << ',';
    }
    output << csvString(result.failure) << '\n';
}

void printProgress(const std::string& stage,
                   int index,
                   int count,
                   const RunResult& result) {
    std::cout << stage << ' ' << index << '/' << count << ": candidate "
              << result.candidate.id << " family="
              << vortexFamilyName(result.candidate.family);
    if (result.candidate.family == VortexFamily::WavePackets) {
        std::cout << " width="
                  << result.candidate.packet_parameters.envelope_width
                  << " carrier="
                  << result.candidate.packet_parameters.carrier_wavenumber
                  << " packet-weight="
                  << result.candidate.packet_parameters.secondary_weight
                  << " packet-phase="
                  << result.candidate.packet_parameters.phase_offset;
    } else {
        std::cout << " core="
                  << result.candidate.tube_parameters.core_radius
                  << " separation="
                  << result.candidate.tube_parameters.separation << " bend="
                  << result.candidate.tube_parameters.bend_amplitude << " axial="
                  << result.candidate.tube_parameters.axial_wavenumber;
    }
    if (result.candidate.family == VortexFamily::OrthogonalBundle) {
        std::cout << " orthogonal-weight="
                  << result.candidate.orthogonal_pair_weight
                  << " phase-offset=" << result.candidate.phase_offset;
    }
    std::cout << " energy="
              << result.candidate.initial_energy;
    if (!result.completed) {
        std::cout << " FAILED: " << result.failure << '\n';
    } else {
        std::cout << " L3 ratio=" << ratio(result.peak_l3, result.initial_l3)
                  << " H1/2 ratio="
                  << ratio(result.peak_h_half, result.initial_h_half)
                  << " vorticity ratio="
                  << ratio(result.peak_vorticity, result.initial_vorticity)
                  << " k_rms ratio="
                  << ratio(result.final_characteristic_wavenumber,
                           result.initial_characteristic_wavenumber)
                  << " profile windows=" << result.profile_scale_windows;
        if (result.profile_scale_windows > 0) {
            std::cout << " latest drift=" << result.latest_profile_drift;
        }
        std::cout
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
        std::vector<FinalistResult> finalists;
        finalists.reserve(static_cast<std::size_t>(finalist_count));
        for (int i = 0; i < finalist_count; ++i) {
            const RunResult& coarse =
                coarse_results[static_cast<std::size_t>(i)];
            const RunResult fine = runCandidateSafely(
                fine_system, coarse.candidate, options);
            printProgress("fine", i + 1, finalist_count, fine);
            const PairAssessment assessment =
                assessResolutionPair(coarse, fine, options);
            std::cout << "resolution pair for candidate " << fine.candidate.id
                      << " score=" << assessment.score.total
                      << " cross-resolution="
                      << (assessment.preliminary_cross_resolution_ok ? "yes"
                                                                     : "no")
                      << " refinement-eligible="
                      << (assessment.score.refinement_eligible ? "yes" : "no")
                      << '\n';
            finalists.push_back(FinalistResult(coarse, fine, assessment));
        }
        std::sort(finalists.begin(), finalists.end(), betterFinalist);

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
                        NULL,
                        NULL);
        }
        for (std::size_t i = 0; i < finalists.size(); ++i) {
            writeResult(csv,
                        "fine",
                        static_cast<int>(i + 1),
                        finalists[i].fine,
                        options,
                        &finalists[i].coarse,
                        &finalists[i].assessment);
        }
        csv.close();
        if (!csv) throw std::runtime_error("Failed while writing " + options.output);

        int coarse_resolved = 0;
        for (std::size_t i = 0; i < coarse_results.size(); ++i) {
            if (coarse_results[i].resolution_ok) ++coarse_resolved;
        }
        int cross_resolved = 0;
        int refinement_eligible = 0;
        for (std::size_t i = 0; i < finalists.size(); ++i) {
            if (finalists[i].assessment.preliminary_cross_resolution_ok) {
                ++cross_resolved;
            }
            if (finalists[i].assessment.score.refinement_eligible) {
                ++refinement_eligible;
            }
        }

        std::cout << std::setprecision(8)
                  << "Search complete: " << candidates.size()
                  << " coarse candidates, " << coarse_resolved
                  << " passed cutoff/constraint gates; " << cross_resolved << '/'
                  << finalist_count
                  << " finalists passed the preliminary cross-resolution gate; "
                  << refinement_eligible
                  << " passed the stricter profile-refinement gate.\n";
        if (!finalists.empty()) {
            const FinalistResult& best_pair = finalists.front();
            const RunResult& best = best_pair.fine;
            std::cout << "Best cross-resolution candidate " << best.candidate.id
                      << ": peak L3 ratio="
                      << ratio(best.peak_l3, best.initial_l3)
                      << ", peak H1/2 ratio="
                      << ratio(best.peak_h_half, best.initial_h_half)
                      << ", peak sampled-vorticity ratio="
                      << ratio(best.peak_vorticity, best.initial_vorticity)
                      << ", k_rms ratio="
                      << ratio(best.final_characteristic_wavenumber,
                               best.initial_characteristic_wavenumber)
                      << ", profile windows=" << best.profile_scale_windows;
            if (best.profile_scale_windows > 0) {
                std::cout << ", latest profile drift="
                          << best.latest_profile_drift;
            }
            std::cout << ", paired score="
                      << best_pair.assessment.score.total
                      << ", peak cutoff fraction=" << best.peak_cutoff_fraction
                      << ".\n";
            if (refinement_eligible == 0) {
                std::cout
                    << "No finalist may be promoted automatically: none combined converged critical-norm growth with an improving, low-drift rescaled profile on both grids.\n";
            } else {
                std::cout
                    << "Only refinement-eligible rows should seed the next parameter sweep; eligibility remains a numerical heuristic, not evidence of singularity.\n";
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
