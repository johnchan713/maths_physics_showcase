#include "ns_cascade/fftw_reference.hpp"
#include "ns_cascade/pseudospectral.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {

struct Options {
    int grid_size = 16;
    int cutoff = 0;
    double viscosity = 0.02;
    double maximum_time_step = 0.005;
    double final_time = 0.02;
    double target_cfl = 0.35;
    double diffusion_safety = 2.0;
    int diagnostic_every = 10;
    double initial_energy = 10.0;
    ns_cascade::InitialCondition initial_condition =
        ns_cascade::InitialCondition::VortexTubes;
    ns_cascade::VortexTubeParameters tube_parameters =
        ns_cascade::VortexTubeParameters(0.50, 1.30, 0.35, 2);
    bool use_orthogonal_bundle = false;
    double orthogonal_pair_weight = 0.75;
    double phase_offset = 1.0471975511965977461542144610932;
    bool use_wave_packets = false;
    ns_cascade::WavePacketParameters packet_parameters;
    bool fixed_time_step = false;
    double state_tolerance = 1e-9;
    double diagnostic_tolerance = 1e-9;
    std::string output = "navier_stokes_fftw_comparison.csv";
};

struct Comparison {
    double maximum_coefficient_difference;
    double relative_state_l2_difference;
    double maximum_diagnostic_scaled_difference;
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
    ++index;
    return argv[index];
}

ns_cascade::InitialCondition parseInitialCondition(const std::string& value) {
    if (value == "deterministic") {
        return ns_cascade::InitialCondition::Deterministic;
    }
    if (value == "taylor-green") {
        return ns_cascade::InitialCondition::TaylorGreen;
    }
    if (value == "abc") return ns_cascade::InitialCondition::ABC;
    if (value == "vortex-tubes") {
        return ns_cascade::InitialCondition::VortexTubes;
    }
    throw std::invalid_argument(
        "Unknown initial condition: " + value +
        " (expected deterministic, taylor-green, abc, or vortex-tubes)");
}

bool parseVortexFamily(const std::string& value) {
    if (value == "pair") return false;
    if (value == "orthogonal-bundle" || value == "bundle") return true;
    throw std::invalid_argument(
        "Unknown vortex family: " + value +
        " (expected pair or orthogonal-bundle)");
}

void printUsage(const char* program) {
    std::cout
        << "Usage: " << program << " [options]\n\n"
        << "Evolve the same state with the internal FFT and an independent FFTW "
           "RK4 path.\n\n"
        << "Options:\n"
        << "  --grid N                 FFT grid size (default: 16)\n"
        << "  --cutoff K               Retained cube cutoff; 0 chooses safe maximum\n"
        << "  --viscosity NU           Non-negative viscosity (default: 0.02)\n"
        << "  --dt DT                  Maximum step (default: 0.005)\n"
        << "  --final-time T           Common target time (default: 0.02)\n"
        << "  --fixed-dt               Disable adaptive CFL selection\n"
        << "  --cfl C                  Adaptive CFL target (default: 0.35)\n"
        << "  --diffusion-safety S     RK4 viscous bound (default: 2)\n"
        << "  --diagnostic-every N     CSV sampling interval (default: 10)\n"
        << "  --initial-condition C    deterministic, taylor-green, abc, or "
           "vortex-tubes\n"
        << "  --energy E               Initial normalized energy (default: 10)\n"
        << "  --tube-core A            Vortex core radius (default: 0.50)\n"
        << "  --tube-separation D      Vortex separation (default: 1.30)\n"
        << "  --tube-bend B            Helical bend (default: 0.35)\n"
        << "  --tube-axial-mode M      Helical mode (default: 2)\n"
        << "  --vortex-family F        pair or orthogonal-bundle (default: pair)\n"
        << "  --orthogonal-weight W    Relative x/y-pair weight (default: 0.75)\n"
        << "  --phase-offset P         Bundle helical phase in radians (default: pi/3)\n"
        << "  --wave-packets           Use the interacting spectral-packet family\n"
        << "  --packet-width W         Periodic Gaussian envelope width\n"
        << "  --carrier-mode M         Packet carrier wavenumber\n"
        << "  --packet-weight W        Relative weight of packets two and three\n"
        << "  --packet-phase P         Packet triad phase offset in radians\n"
        << "  --state-tolerance X      Relative trajectory gate (default: 1e-9)\n"
        << "  --diagnostic-tolerance X Scaled diagnostic gate (default: 1e-9)\n"
        << "  --output PATH            Comparison CSV path\n"
        << "  --help                   Show this message\n";
}

Options parseOptions(int argc, char** argv) {
    Options options;
    for (int i = 1; i < argc; ++i) {
        const std::string flag = argv[i];
        if (flag == "--help") {
            printUsage(argv[0]);
            std::exit(0);
        } else if (flag == "--grid") {
            options.grid_size =
                parseNumber<int>(requireValue(i, argc, argv), flag);
        } else if (flag == "--cutoff") {
            options.cutoff =
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
        } else if (flag == "--fixed-dt") {
            options.fixed_time_step = true;
        } else if (flag == "--cfl") {
            options.target_cfl =
                parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--diffusion-safety") {
            options.diffusion_safety =
                parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--diagnostic-every") {
            options.diagnostic_every =
                parseNumber<int>(requireValue(i, argc, argv), flag);
        } else if (flag == "--initial-condition") {
            options.initial_condition =
                parseInitialCondition(requireValue(i, argc, argv));
        } else if (flag == "--energy") {
            options.initial_energy =
                parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--tube-core") {
            options.tube_parameters.core_radius =
                parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--tube-separation") {
            options.tube_parameters.separation =
                parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--tube-bend") {
            options.tube_parameters.bend_amplitude =
                parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--tube-axial-mode") {
            options.tube_parameters.axial_wavenumber =
                parseNumber<int>(requireValue(i, argc, argv), flag);
        } else if (flag == "--vortex-family") {
            options.use_orthogonal_bundle =
                parseVortexFamily(requireValue(i, argc, argv));
        } else if (flag == "--orthogonal-weight") {
            options.orthogonal_pair_weight =
                parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--phase-offset") {
            options.phase_offset =
                parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--wave-packets") {
            options.use_wave_packets = true;
        } else if (flag == "--packet-width") {
            options.packet_parameters.envelope_width =
                parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--carrier-mode") {
            options.packet_parameters.carrier_wavenumber =
                parseNumber<int>(requireValue(i, argc, argv), flag);
        } else if (flag == "--packet-weight") {
            options.packet_parameters.secondary_weight =
                parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--packet-phase") {
            options.packet_parameters.phase_offset =
                parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--state-tolerance") {
            options.state_tolerance =
                parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--diagnostic-tolerance") {
            options.diagnostic_tolerance =
                parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--output") {
            options.output = requireValue(i, argc, argv);
        } else {
            throw std::invalid_argument("Unknown option: " + flag);
        }
    }

    if (options.cutoff < 0) {
        throw std::invalid_argument("--cutoff cannot be negative");
    }
    if (!std::isfinite(options.viscosity) || options.viscosity < 0.0) {
        throw std::invalid_argument("--viscosity must be finite and non-negative");
    }
    if (!std::isfinite(options.maximum_time_step) ||
        options.maximum_time_step <= 0.0) {
        throw std::invalid_argument("--dt must be finite and positive");
    }
    if (!std::isfinite(options.final_time) || options.final_time <= 0.0) {
        throw std::invalid_argument("--final-time must be finite and positive");
    }
    if (!std::isfinite(options.target_cfl) || options.target_cfl <= 0.0 ||
        !std::isfinite(options.diffusion_safety) ||
        options.diffusion_safety <= 0.0) {
        throw std::invalid_argument(
            "--cfl and --diffusion-safety must be finite and positive");
    }
    if (options.diagnostic_every < 1) {
        throw std::invalid_argument("--diagnostic-every must be at least one");
    }
    if (!std::isfinite(options.initial_energy) ||
        options.initial_energy <= 0.0) {
        throw std::invalid_argument("--energy must be finite and positive");
    }
    if (!std::isfinite(options.orthogonal_pair_weight) ||
        options.orthogonal_pair_weight < 0.0) {
        throw std::invalid_argument(
            "--orthogonal-weight must be finite and non-negative");
    }
    if (!std::isfinite(options.phase_offset)) {
        throw std::invalid_argument("--phase-offset must be finite");
    }
    if (!std::isfinite(options.packet_parameters.envelope_width) ||
        options.packet_parameters.envelope_width <= 0.0 ||
        !std::isfinite(options.packet_parameters.secondary_weight) ||
        options.packet_parameters.secondary_weight <= 0.0 ||
        !std::isfinite(options.packet_parameters.phase_offset) ||
        options.packet_parameters.carrier_wavenumber < 1) {
        throw std::invalid_argument("Wave-packet parameters are invalid");
    }
    if (!std::isfinite(options.state_tolerance) ||
        options.state_tolerance <= 0.0 ||
        !std::isfinite(options.diagnostic_tolerance) ||
        options.diagnostic_tolerance <= 0.0) {
        throw std::invalid_argument("Comparison tolerances must be positive");
    }
    if (options.output.empty()) {
        throw std::invalid_argument("--output cannot be empty");
    }
    return options;
}

double scaledDifference(double left, double right) {
    if (!std::isfinite(left) || !std::isfinite(right)) {
        return std::numeric_limits<double>::infinity();
    }
    const double scale = std::max(
        std::max(std::abs(left), std::abs(right)), 1.0);
    return std::abs(left - right) / scale;
}

bool diagnosticsAreFinite(
    const ns_cascade::PseudospectralSystem::Diagnostics& values,
    const ns_cascade::FftwReferenceDiagnostics& reference) {
    return std::isfinite(values.energy) &&
           std::isfinite(values.enstrophy) &&
           std::isfinite(values.critical_h_half) &&
           std::isfinite(values.critical_l3_sample) &&
           std::isfinite(values.sampled_vorticity_max) &&
           std::isfinite(values.high_shell_energy_fraction) &&
           std::isfinite(values.divergence_defect) &&
           std::isfinite(values.reality_defect) &&
           std::isfinite(reference.energy) &&
           std::isfinite(reference.enstrophy) &&
           std::isfinite(reference.critical_h_half) &&
           std::isfinite(reference.critical_l3_sample) &&
           std::isfinite(reference.sampled_vorticity_max) &&
           std::isfinite(reference.cutoff_shell_energy_fraction) &&
           std::isfinite(reference.divergence_defect) &&
           std::isfinite(reference.reality_defect);
}

Comparison compareStates(
    const ns_cascade::PseudospectralSystem::State& internal_state,
    const ns_cascade::FftwReferenceSystem::State& reference_state,
    const ns_cascade::PseudospectralSystem::Diagnostics& internal,
    const ns_cascade::FftwReferenceDiagnostics& reference) {
    if (internal_state.size() != reference_state.size()) {
        throw std::invalid_argument("Cannot compare states of different sizes");
    }
    double difference_squared = 0.0;
    double reference_squared = 0.0;
    Comparison comparison;
    comparison.maximum_coefficient_difference = 0.0;
    for (std::size_t i = 0; i < internal_state.size(); ++i) {
        const ns_cascade::ComplexVector difference =
            internal_state[i] - reference_state[i];
        difference_squared += ns_cascade::normSquared(difference);
        reference_squared += ns_cascade::normSquared(reference_state[i]);
        comparison.maximum_coefficient_difference = std::max(
            comparison.maximum_coefficient_difference,
            ns_cascade::norm(difference));
    }
    comparison.relative_state_l2_difference = std::sqrt(
        difference_squared / std::max(reference_squared, 1e-300));

    comparison.maximum_diagnostic_scaled_difference = 0.0;
    comparison.maximum_diagnostic_scaled_difference = std::max(
        comparison.maximum_diagnostic_scaled_difference,
        scaledDifference(internal.energy, reference.energy));
    comparison.maximum_diagnostic_scaled_difference = std::max(
        comparison.maximum_diagnostic_scaled_difference,
        scaledDifference(internal.enstrophy, reference.enstrophy));
    comparison.maximum_diagnostic_scaled_difference = std::max(
        comparison.maximum_diagnostic_scaled_difference,
        scaledDifference(internal.critical_h_half,
                         reference.critical_h_half));
    comparison.maximum_diagnostic_scaled_difference = std::max(
        comparison.maximum_diagnostic_scaled_difference,
        scaledDifference(internal.critical_l3_sample,
                         reference.critical_l3_sample));
    comparison.maximum_diagnostic_scaled_difference = std::max(
        comparison.maximum_diagnostic_scaled_difference,
        scaledDifference(internal.sampled_vorticity_max,
                         reference.sampled_vorticity_max));
    comparison.maximum_diagnostic_scaled_difference = std::max(
        comparison.maximum_diagnostic_scaled_difference,
        scaledDifference(internal.high_shell_energy_fraction,
                         reference.cutoff_shell_energy_fraction));
    return comparison;
}

void writeHeader(std::ostream& output) {
    output
        << "step,time,time_step,max_coefficient_difference,"
        << "relative_state_l2_difference,max_diagnostic_scaled_difference,"
        << "internal_energy,fftw_energy,internal_enstrophy,fftw_enstrophy,"
        << "internal_critical_h_half,fftw_critical_h_half,"
        << "internal_critical_l3,fftw_critical_l3,"
        << "internal_vorticity_max,fftw_vorticity_max,"
        << "internal_cutoff_fraction,fftw_cutoff_fraction,"
        << "internal_divergence_defect,fftw_divergence_defect,"
        << "internal_reality_defect,fftw_reality_defect,"
        << "internal_cfl_bound,fftw_cfl_bound,viscous_stability_number\n";
}

void writeRow(std::ostream& output,
              std::uint64_t step,
              double time,
              double time_step,
              const Comparison& comparison,
              const ns_cascade::PseudospectralSystem::Diagnostics& internal,
              const ns_cascade::FftwReferenceDiagnostics& reference,
              double internal_cfl,
              double reference_cfl,
              double viscous_number) {
    output << step << ',' << time << ',' << time_step << ','
           << comparison.maximum_coefficient_difference << ','
           << comparison.relative_state_l2_difference << ','
           << comparison.maximum_diagnostic_scaled_difference << ','
           << internal.energy << ',' << reference.energy << ','
           << internal.enstrophy << ',' << reference.enstrophy << ','
           << internal.critical_h_half << ',' << reference.critical_h_half << ','
           << internal.critical_l3_sample << ','
           << reference.critical_l3_sample << ','
           << internal.sampled_vorticity_max << ','
           << reference.sampled_vorticity_max << ','
           << internal.high_shell_energy_fraction << ','
           << reference.cutoff_shell_energy_fraction << ','
           << internal.divergence_defect << ','
           << reference.divergence_defect << ','
           << internal.reality_defect << ',' << reference.reality_defect << ','
           << internal_cfl << ',' << reference_cfl << ',' << viscous_number
           << '\n';
}

ns_cascade::PseudospectralSystem::State makeInitialState(
    const ns_cascade::PseudospectralSystem& system,
    const Options& options) {
    if (options.use_wave_packets) {
        return system.interactingWavePacketState(
            options.packet_parameters, options.initial_energy);
    }
    if (options.initial_condition == ns_cascade::InitialCondition::VortexTubes) {
        if (options.use_orthogonal_bundle) {
            return system.vortexBundleState(
                ns_cascade::VortexBundleParameters(
                    options.tube_parameters,
                    options.orthogonal_pair_weight,
                    options.phase_offset),
                options.initial_energy);
        }
        return system.vortexTubePairState(
            options.tube_parameters, options.initial_energy);
    }
    return system.initialState(options.initial_condition, options.initial_energy);
}

int run(const Options& options) {
    const ns_cascade::PseudospectralSystem internal_system(
        options.grid_size, options.viscosity, options.cutoff);
    const ns_cascade::FftwReferenceSystem reference_system(
        options.grid_size, options.viscosity, options.cutoff);
    ns_cascade::PseudospectralSystem::State internal_state =
        makeInitialState(internal_system, options);
    ns_cascade::FftwReferenceSystem::State reference_state = internal_state;

    std::ofstream csv(options.output.c_str());
    if (!csv) {
        throw std::runtime_error("Could not open comparison CSV: " +
                                 options.output);
    }
    csv << std::setprecision(17);
    writeHeader(csv);

    double time = 0.0;
    std::uint64_t step = 0U;
    double peak_state_difference = 0.0;
    double peak_diagnostic_difference = 0.0;
    double peak_divergence_defect = 0.0;
    double peak_reality_defect = 0.0;

    const ns_cascade::PseudospectralSystem::Diagnostics initial_internal =
        internal_system.diagnostics(internal_state);
    const ns_cascade::FftwReferenceDiagnostics initial_reference =
        reference_system.diagnostics(reference_state);
    if (!diagnosticsAreFinite(initial_internal, initial_reference)) {
        throw std::runtime_error("Initial comparison diagnostics are non-finite");
    }
    const Comparison initial_comparison = compareStates(
        internal_state, reference_state, initial_internal, initial_reference);
    peak_diagnostic_difference =
        initial_comparison.maximum_diagnostic_scaled_difference;
    peak_divergence_defect = std::max(
        initial_internal.divergence_defect,
        initial_reference.divergence_defect);
    peak_reality_defect = std::max(
        initial_internal.reality_defect,
        initial_reference.reality_defect);
    writeRow(csv,
             step,
             time,
             0.0,
             initial_comparison,
             initial_internal,
             initial_reference,
             0.0,
             0.0,
             0.0);

    while (time < options.final_time) {
        const double remaining = options.final_time - time;
        const double proposed_maximum =
            std::min(options.maximum_time_step, remaining);
        double time_step = proposed_maximum;
        const ns_cascade::AdaptiveStepInfo internal_step =
            internal_system.chooseAdaptiveTimeStep(
                internal_state,
                proposed_maximum,
                options.target_cfl,
                options.diffusion_safety);
        const ns_cascade::FftwReferenceStepInfo reference_step =
            reference_system.chooseAdaptiveTimeStep(
                reference_state,
                proposed_maximum,
                options.target_cfl,
                options.diffusion_safety);
        if (!options.fixed_time_step) {
            time_step = std::min(internal_step.time_step,
                                 reference_step.time_step);
        }
        const double internal_velocity_bound =
            internal_step.velocity_supremum_bound;
        const double reference_velocity_bound =
            reference_step.velocity_supremum_bound;
        if (time + time_step == time) {
            throw std::runtime_error(
                "Time step is too small to advance floating-point time");
        }

        internal_system.stepRungeKutta4(internal_state, time_step);
        reference_system.stepRungeKutta4(reference_state, time_step);
        time += time_step;
        if (time > options.final_time ||
            options.final_time - time <=
                8.0 * std::numeric_limits<double>::epsilon() *
                    std::max(1.0, options.final_time)) {
            time = options.final_time;
        }
        ++step;

        double difference_squared = 0.0;
        double reference_squared = 0.0;
        for (std::size_t i = 0; i < internal_state.size(); ++i) {
            difference_squared += ns_cascade::normSquared(
                internal_state[i] - reference_state[i]);
            reference_squared += ns_cascade::normSquared(reference_state[i]);
        }
        const double current_state_difference = std::sqrt(
            difference_squared / std::max(reference_squared, 1e-300));
        if (!std::isfinite(current_state_difference)) {
            throw std::runtime_error("A trajectory comparison became non-finite");
        }
        peak_state_difference = std::max(
            peak_state_difference, current_state_difference);

        const bool diagnostic_step =
            step % static_cast<std::uint64_t>(options.diagnostic_every) == 0U ||
            time >= options.final_time;
        if (!diagnostic_step) continue;

        const ns_cascade::PseudospectralSystem::Diagnostics internal =
            internal_system.diagnostics(internal_state);
        const ns_cascade::FftwReferenceDiagnostics reference =
            reference_system.diagnostics(reference_state);
        if (!diagnosticsAreFinite(internal, reference)) {
            throw std::runtime_error("Comparison diagnostics became non-finite");
        }
        const Comparison comparison = compareStates(
            internal_state, reference_state, internal, reference);
        peak_diagnostic_difference = std::max(
            peak_diagnostic_difference,
            comparison.maximum_diagnostic_scaled_difference);
        peak_divergence_defect = std::max(
            peak_divergence_defect,
            std::max(internal.divergence_defect,
                     reference.divergence_defect));
        peak_reality_defect = std::max(
            peak_reality_defect,
            std::max(internal.reality_defect, reference.reality_defect));

        const double largest_wave_number =
            std::sqrt(3.0) * static_cast<double>(internal_system.cutoff());
        const double internal_cfl =
            time_step * internal_velocity_bound * largest_wave_number;
        const double reference_cfl =
            time_step * reference_velocity_bound * largest_wave_number;
        const double viscous_number =
            time_step * options.viscosity * 3.0 * internal_system.cutoff() *
            internal_system.cutoff();
        writeRow(csv,
                 step,
                 time,
                 time_step,
                 comparison,
                 internal,
                 reference,
                 internal_cfl,
                 reference_cfl,
                 viscous_number);
    }

    if (!csv) throw std::runtime_error("Failed while writing comparison CSV");
    const bool passed =
        std::isfinite(peak_state_difference) &&
        std::isfinite(peak_diagnostic_difference) &&
        peak_state_difference <= options.state_tolerance &&
        peak_diagnostic_difference <= options.diagnostic_tolerance &&
        peak_divergence_defect <= options.diagnostic_tolerance &&
        peak_reality_defect <= options.diagnostic_tolerance;

    std::cout << std::setprecision(12)
              << "Independent FFTW trajectory comparison\n"
              << "  grid/cutoff: " << options.grid_size << '/'
              << internal_system.cutoff() << '\n'
              << "  initial family: "
              << (options.use_wave_packets
                      ? "wave-packets"
                      : (options.use_orthogonal_bundle ? "orthogonal-bundle"
                                                       : "pair"))
              << '\n'
              << "  final time/steps: " << time << '/' << step << '\n'
              << "  peak relative state difference: "
              << peak_state_difference << '\n'
              << "  peak diagnostic scaled difference: "
              << peak_diagnostic_difference << '\n'
              << "  peak divergence/reality defects: "
              << peak_divergence_defect << '/' << peak_reality_defect << '\n'
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
