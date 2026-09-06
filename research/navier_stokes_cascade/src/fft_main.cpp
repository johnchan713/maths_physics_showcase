#include "ns_cascade/checkpoint.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstdint>
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
    int grid_size = 16;
    int cutoff = 0;
    double viscosity = 0.05;
    double time_step = 0.0005;
    double final_time = 0.0;
    bool final_time_set = false;
    int steps = 200;
    int diagnostic_every = 20;
    double initial_energy = 1.0;
    bool adaptive = true;
    double target_cfl = 0.4;
    double diffusion_safety = 2.0;
    ns_cascade::VortexTubeParameters tube_parameters;
    bool custom_tube_parameters = false;
    ns_cascade::InitialCondition initial_condition =
        ns_cascade::InitialCondition::TaylorGreen;
    std::string output = "navier_stokes_fft.csv";
    std::string shell_output = "navier_stokes_fft_shells.csv";
    std::string profile_output = "navier_stokes_fft_profile.csv";
    int profile_bin_count = 32;
    double profile_maximum_coordinate = 4.0;
    std::string restart_path;
    std::string checkpoint_output;
    int checkpoint_every = 100;
    bool scientific_configuration_set = false;
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

void printUsage(const char* program) {
    std::cout
        << "Usage: " << program << " [options]\n\n"
        << "Strictly dealiased, radix-2 FFT Navier-Stokes experiment.\n\n"
        << "Options:\n"
        << "  --grid N                Power-of-two grid size (default: 16)\n"
        << "  --cutoff K              Retained cube cutoff; 0 uses floor((N-1)/3)\n"
        << "  --viscosity NU          Non-negative viscosity (default: 0.05)\n"
        << "  --dt DT                 Maximum adaptive step (default: 0.0005)\n"
        << "  --steps N               Sets final time to N*dt unless --final-time is used\n"
        << "  --final-time T          Physical end time (default: steps*dt)\n"
        << "  --adaptive-dt           Enable adaptive control (default)\n"
        << "  --fixed-dt              Use fixed dt, shortening only the final step\n"
        << "  --cfl C                 Conservative advective CFL target (default: 0.4)\n"
        << "  --diffusion-safety S    Bound on nu*|k|max^2*dt (default: 2)\n"
        << "  --diagnostic-every N    CSV sampling interval (default: 20)\n"
        << "  --energy E              Initial normalized energy (default: 1)\n"
        << "  --initial-condition C   deterministic, taylor-green, abc, or vortex-tubes\n"
        << "                            (default: taylor-green)\n"
        << "  --tube-core R           Vortex-tube core radius (default: 0.55)\n"
        << "  --tube-separation D     Tube-centre separation (default: 1.6)\n"
        << "  --tube-bend A           Helical bend amplitude (default: 0.25)\n"
        << "  --tube-axial-mode M     Positive axial wavenumber (default: 1)\n"
        << "  --output PATH           Time-series CSV path\n"
        << "  --shell-output PATH     Shell-transfer CSV path\n"
        << "  --profile-output PATH   Rescaled-spectrum CSV path\n"
        << "  --profile-bins N        Finite bins before overflow (default: 32)\n"
        << "  --profile-max-xi X      Last finite |k|/k_rms coordinate (default: 4)\n"
        << "  --checkpoint-output P  Atomically replace checkpoint P while running\n"
        << "  --checkpoint-every N   Save every N accepted steps (default: 100)\n"
        << "  --restart PATH          Resume the exact state stored in PATH\n"
        << "  --help                  Show this message\n";
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
            options.scientific_configuration_set = true;
        } else if (flag == "--cutoff") {
            options.cutoff = parseNumber<int>(requireValue(i, argc, argv), flag);
            options.scientific_configuration_set = true;
        } else if (flag == "--viscosity") {
            options.viscosity = parseNumber<double>(requireValue(i, argc, argv), flag);
            options.scientific_configuration_set = true;
        } else if (flag == "--dt") {
            options.time_step = parseNumber<double>(requireValue(i, argc, argv), flag);
            options.scientific_configuration_set = true;
        } else if (flag == "--final-time") {
            options.final_time =
                parseNumber<double>(requireValue(i, argc, argv), flag);
            options.final_time_set = true;
        } else if (flag == "--steps") {
            options.steps = parseNumber<int>(requireValue(i, argc, argv), flag);
            options.scientific_configuration_set = true;
        } else if (flag == "--adaptive-dt") {
            options.adaptive = true;
            options.scientific_configuration_set = true;
        } else if (flag == "--fixed-dt") {
            options.adaptive = false;
            options.scientific_configuration_set = true;
        } else if (flag == "--cfl") {
            options.target_cfl =
                parseNumber<double>(requireValue(i, argc, argv), flag);
            options.scientific_configuration_set = true;
        } else if (flag == "--diffusion-safety") {
            options.diffusion_safety =
                parseNumber<double>(requireValue(i, argc, argv), flag);
            options.scientific_configuration_set = true;
        } else if (flag == "--diagnostic-every") {
            options.diagnostic_every =
                parseNumber<int>(requireValue(i, argc, argv), flag);
        } else if (flag == "--energy") {
            options.initial_energy =
                parseNumber<double>(requireValue(i, argc, argv), flag);
            options.scientific_configuration_set = true;
        } else if (flag == "--initial-condition") {
            options.initial_condition =
                parseInitialCondition(requireValue(i, argc, argv));
            options.scientific_configuration_set = true;
        } else if (flag == "--tube-core") {
            options.tube_parameters.core_radius =
                parseNumber<double>(requireValue(i, argc, argv), flag);
            options.custom_tube_parameters = true;
            options.scientific_configuration_set = true;
        } else if (flag == "--tube-separation") {
            options.tube_parameters.separation =
                parseNumber<double>(requireValue(i, argc, argv), flag);
            options.custom_tube_parameters = true;
            options.scientific_configuration_set = true;
        } else if (flag == "--tube-bend") {
            options.tube_parameters.bend_amplitude =
                parseNumber<double>(requireValue(i, argc, argv), flag);
            options.custom_tube_parameters = true;
            options.scientific_configuration_set = true;
        } else if (flag == "--tube-axial-mode") {
            options.tube_parameters.axial_wavenumber =
                parseNumber<int>(requireValue(i, argc, argv), flag);
            options.custom_tube_parameters = true;
            options.scientific_configuration_set = true;
        } else if (flag == "--output") {
            options.output = requireValue(i, argc, argv);
        } else if (flag == "--shell-output") {
            options.shell_output = requireValue(i, argc, argv);
        } else if (flag == "--profile-output") {
            options.profile_output = requireValue(i, argc, argv);
        } else if (flag == "--profile-bins") {
            options.profile_bin_count =
                parseNumber<int>(requireValue(i, argc, argv), flag);
            options.scientific_configuration_set = true;
        } else if (flag == "--profile-max-xi") {
            options.profile_maximum_coordinate =
                parseNumber<double>(requireValue(i, argc, argv), flag);
            options.scientific_configuration_set = true;
        } else if (flag == "--checkpoint-output") {
            options.checkpoint_output = requireValue(i, argc, argv);
        } else if (flag == "--checkpoint-every") {
            options.checkpoint_every =
                parseNumber<int>(requireValue(i, argc, argv), flag);
        } else if (flag == "--restart") {
            options.restart_path = requireValue(i, argc, argv);
        } else {
            throw std::invalid_argument("Unknown option: " + flag);
        }
    }

    if (options.cutoff < 0) throw std::invalid_argument("--cutoff cannot be negative");
    if (!std::isfinite(options.viscosity) || options.viscosity < 0.0) {
        throw std::invalid_argument("--viscosity must be non-negative");
    }
    if (!std::isfinite(options.time_step) || options.time_step <= 0.0) {
        throw std::invalid_argument("--dt must be finite and positive");
    }
    if (!std::isfinite(options.final_time) ||
        (options.final_time_set && options.final_time <= 0.0)) {
        throw std::invalid_argument("--final-time must be finite and positive");
    }
    if (options.steps < 1) throw std::invalid_argument("--steps must be >= 1");
    if (options.diagnostic_every < 1) {
        throw std::invalid_argument("--diagnostic-every must be >= 1");
    }
    if (!std::isfinite(options.initial_energy) || options.initial_energy <= 0.0) {
        throw std::invalid_argument("--energy must be finite and positive");
    }
    if (!std::isfinite(options.target_cfl) || options.target_cfl <= 0.0) {
        throw std::invalid_argument("--cfl must be finite and positive");
    }
    if (!std::isfinite(options.diffusion_safety) ||
        options.diffusion_safety <= 0.0) {
        throw std::invalid_argument(
            "--diffusion-safety must be finite and positive");
    }
    if (options.custom_tube_parameters &&
        options.initial_condition != ns_cascade::InitialCondition::VortexTubes) {
        throw std::invalid_argument(
            "Tube geometry options require --initial-condition vortex-tubes");
    }
    if (options.profile_bin_count < 4 || options.profile_bin_count > 4096) {
        throw std::invalid_argument("--profile-bins must be between 4 and 4096");
    }
    if (!std::isfinite(options.profile_maximum_coordinate) ||
        options.profile_maximum_coordinate <= 0.0) {
        throw std::invalid_argument("--profile-max-xi must be finite and positive");
    }
    if (options.checkpoint_every < 1) {
        throw std::invalid_argument("--checkpoint-every must be >= 1");
    }
    if (options.output.empty() || options.shell_output.empty() ||
        options.profile_output.empty()) {
        throw std::invalid_argument("Output paths cannot be empty");
    }
    if (options.output == options.shell_output ||
        options.output == options.profile_output ||
        options.shell_output == options.profile_output) {
        throw std::invalid_argument("CSV output paths must be different");
    }
    if (!options.checkpoint_output.empty() &&
        (options.checkpoint_output == options.output ||
         options.checkpoint_output == options.shell_output ||
         options.checkpoint_output == options.profile_output)) {
        throw std::invalid_argument(
            "Checkpoint and CSV output paths must be different");
    }
    if (!options.checkpoint_output.empty()) {
        const std::string temporary_checkpoint_path =
            options.checkpoint_output + ".tmp";
        if (temporary_checkpoint_path == options.output ||
            temporary_checkpoint_path == options.shell_output ||
            temporary_checkpoint_path == options.profile_output ||
            temporary_checkpoint_path == options.restart_path) {
            throw std::invalid_argument(
                "Temporary checkpoint path collides with an input or output path");
        }
    }
    if (!options.restart_path.empty()) {
        if (options.restart_path == options.output ||
            options.restart_path == options.shell_output ||
            options.restart_path == options.profile_output) {
            throw std::invalid_argument(
                "Restart input and CSV output paths must be different");
        }
        if (options.scientific_configuration_set) {
            throw std::invalid_argument(
                "A restart restores all scientific parameters; do not combine "
                "--restart with grid, solver, initial-data, or profile options");
        }
        if (!options.final_time_set) {
            throw std::invalid_argument(
                "--restart requires an absolute --final-time greater than the checkpoint time");
        }
    }
    return options;
}

void writeTimeHeader(std::ostream& output) {
    output
        << "step,time,dt_used,adaptive,advective_cfl_upper_bound,"
        << "viscous_stability_number,velocity_supremum_bound,grid_size,cutoff,"
        << "initial_condition,tube_core_radius,tube_separation,"
        << "tube_bend_amplitude,tube_axial_wavenumber,"
        << "energy,enstrophy,palinstrophy,critical_h_half,"
        << "critical_l3_sample,sampled_vorticity_max,"
        << "vorticity_sup_upper_bound,spectral_centroid,"
        << "high_shell_energy_fraction,divergence_defect,reality_defect,"
        << "nonlinear_enstrophy_production,viscous_enstrophy_destruction,"
        << "net_enstrophy_rate,enstrophy_production_to_dissipation,"
        << "energy_balance_residual,bkm_sampled_integral,"
        << "characteristic_wavenumber,characteristic_wavenumber_ratio,"
        << "rescaled_profile_l1_from_initial,"
        << "rescaled_profile_l1_from_previous,"
        << "rescaled_profile_overlap_from_previous,"
        << "rescaled_profile_absolute_log_scale_change,"
        << "rescaled_profile_l1_per_log_scale_change,"
        << "rescaled_profile_log_scale_change_valid,"
        << "rescaled_profile_mean,rescaled_profile_standard_deviation,"
        << "analyticity_radius_estimate,"
        << "analyticity_radius_times_characteristic_wavenumber,"
        << "analyticity_fit_r_squared,"
        << "analyticity_fit_shells,analyticity_fit_valid\n";
}

void writeTimeRow(std::ostream& output,
                  std::uint64_t step,
                  double time,
                  const Options& options,
                  const ns_cascade::PseudospectralSystem& system,
                  const ns_cascade::PseudospectralSystem::Diagnostics& values,
                  const ns_cascade::AdaptiveStepInfo& step_information,
                  double bkm_integral,
                  const ns_cascade::SpectrumProfile& profile,
                  const ns_cascade::SpectrumProfile& initial_profile,
                  const ns_cascade::SpectrumProfile& previous_profile) {
    const double l1_from_initial =
        ns_cascade::spectrumProfileL1Distance(profile, initial_profile);
    const ns_cascade::SpectrumProfileChange profile_change =
        ns_cascade::spectrumProfileChange(profile, previous_profile);
    output << step << ',' << time << ',' << step_information.time_step << ','
           << (options.adaptive ? "true" : "false") << ','
           << step_information.advective_cfl_upper_bound << ','
           << step_information.viscous_stability_number << ','
           << step_information.velocity_supremum_bound << ','
           << system.gridSize() << ',' << system.cutoff() << ','
           << ns_cascade::initialConditionName(options.initial_condition) << ','
           << options.tube_parameters.core_radius << ','
           << options.tube_parameters.separation << ','
           << options.tube_parameters.bend_amplitude << ','
           << options.tube_parameters.axial_wavenumber << ','
           << values.energy << ',' << values.enstrophy
           << ',' << values.palinstrophy << ',' << values.critical_h_half << ','
           << values.critical_l3_sample << ','
           << values.sampled_vorticity_max << ','
           << values.vorticity_sup_upper_bound << ',' << values.spectral_centroid
           << ',' << values.high_shell_energy_fraction << ','
           << values.divergence_defect << ',' << values.reality_defect << ','
           << values.nonlinear_enstrophy_production << ','
           << values.viscous_enstrophy_destruction << ','
           << values.net_enstrophy_rate << ','
           << values.enstrophy_production_to_dissipation << ','
           << values.energy_balance_residual << ',' << bkm_integral << ','
           << profile.characteristic_wavenumber << ','
           << profile.characteristic_wavenumber /
                  initial_profile.characteristic_wavenumber
           << ',' << l1_from_initial << ',' << profile_change.l1_distance << ','
           << profile_change.overlap << ','
           << profile_change.absolute_log_scale_change << ','
           << profile_change.l1_per_log_scale_change << ','
           << (profile_change.scale_normalized_drift_valid ? "true" : "false")
           << ',' << profile.rescaled_mean << ','
           << profile.rescaled_standard_deviation << ','
           << profile.analyticity_radius_estimate << ','
           << profile.analyticity_radius_estimate *
                  profile.characteristic_wavenumber
           << ','
           << profile.analyticity_fit_r_squared << ','
           << profile.analyticity_fit_shells << ','
           << (profile.analyticity_fit_valid ? "true" : "false") << '\n';
}

void writeShellHeader(std::ostream& output) {
    output
        << "initial_condition,tube_core_radius,tube_separation,"
        << "tube_bend_amplitude,tube_axial_wavenumber,"
        << "step,time,grid_size,cutoff,shell,lower_radius,"
        << "upper_radius,shell_energy,nonlinear_transfer,viscous_dissipation,"
        << "forward_flux\n";
}

double writeShellRows(
    std::ostream& output,
    const Options& options,
    std::uint64_t step,
    double time,
    const ns_cascade::PseudospectralSystem& system,
    const std::vector<ns_cascade::PseudospectralSystem::ShellDiagnostics>& shells) {
    double largest_forward_flux = 0.0;
    for (std::size_t i = 0; i < shells.size(); ++i) {
        const ns_cascade::PseudospectralSystem::ShellDiagnostics& values = shells[i];
        output << ns_cascade::initialConditionName(options.initial_condition) << ','
               << options.tube_parameters.core_radius << ','
               << options.tube_parameters.separation << ','
               << options.tube_parameters.bend_amplitude << ','
               << options.tube_parameters.axial_wavenumber << ',' << step
               << ',' << time << ',' << system.gridSize() << ',' << system.cutoff()
               << ',' << values.shell << ',' << values.lower_radius << ','
               << values.upper_radius << ',' << values.energy << ','
               << values.nonlinear_transfer << ',' << values.viscous_dissipation
               << ',' << values.forward_flux << '\n';
        largest_forward_flux =
            std::max(largest_forward_flux, values.forward_flux);
    }
    return largest_forward_flux;
}

void writeProfileHeader(std::ostream& output) {
    output
        << "step,time,grid_size,cutoff,bin,rescaled_bin_center,"
        << "overflow_node,energy_weight,characteristic_wavenumber,"
        << "analyticity_radius_estimate,"
        << "analyticity_radius_times_characteristic_wavenumber,"
        << "analyticity_fit_r_squared,"
        << "analyticity_fit_shells,analyticity_fit_valid\n";
}

void writeProfileRows(std::ostream& output,
                      std::uint64_t step,
                      double time,
                      const ns_cascade::PseudospectralSystem& system,
                      const ns_cascade::SpectrumProfile& profile) {
    const std::size_t finite_bin_count = profile.energy_fractions.size() - 1;
    const double width =
        profile.maximum_rescaled_coordinate / finite_bin_count;
    for (std::size_t bin = 0; bin < profile.energy_fractions.size(); ++bin) {
        const bool overflow = bin == finite_bin_count;
        output << step << ',' << time << ',' << system.gridSize() << ','
               << system.cutoff() << ',' << bin << ','
               << width * (static_cast<double>(bin) + 0.5) << ','
               << (overflow ? "true" : "false") << ','
               << profile.energy_fractions[bin] << ','
               << profile.characteristic_wavenumber << ','
               << profile.analyticity_radius_estimate << ','
               << profile.analyticity_radius_estimate *
                      profile.characteristic_wavenumber
               << ','
               << profile.analyticity_fit_r_squared << ','
               << profile.analyticity_fit_shells << ','
               << (profile.analyticity_fit_valid ? "true" : "false") << '\n';
    }
}

void restoreCheckpointConfiguration(
    Options& options,
    const ns_cascade::CheckpointConfiguration& configuration) {
    options.grid_size = configuration.grid_size;
    options.cutoff = configuration.cutoff;
    options.viscosity = configuration.viscosity;
    options.initial_condition = configuration.initial_condition;
    options.tube_parameters = configuration.tube_parameters;
    options.initial_energy = configuration.initial_energy;
    options.adaptive = configuration.adaptive;
    options.time_step = configuration.maximum_time_step;
    options.target_cfl = configuration.target_cfl;
    options.diffusion_safety = configuration.diffusion_safety;
    options.profile_bin_count = configuration.profile_bin_count;
    options.profile_maximum_coordinate =
        configuration.profile_maximum_coordinate;
}

ns_cascade::SimulationCheckpoint makeCheckpoint(
    const Options& options,
    const ns_cascade::PseudospectralSystem& system,
    const ns_cascade::PseudospectralSystem::State& state,
    const ns_cascade::CheckpointProgress& progress,
    const ns_cascade::SpectrumProfile& initial_profile,
    const ns_cascade::SpectrumProfile& previous_profile) {
    ns_cascade::SimulationCheckpoint checkpoint;
    checkpoint.configuration.grid_size = system.gridSize();
    checkpoint.configuration.cutoff = system.cutoff();
    checkpoint.configuration.viscosity = system.viscosity();
    checkpoint.configuration.initial_condition = options.initial_condition;
    checkpoint.configuration.tube_parameters = options.tube_parameters;
    checkpoint.configuration.initial_energy = options.initial_energy;
    checkpoint.configuration.adaptive = options.adaptive;
    checkpoint.configuration.maximum_time_step = options.time_step;
    checkpoint.configuration.target_cfl = options.target_cfl;
    checkpoint.configuration.diffusion_safety = options.diffusion_safety;
    checkpoint.configuration.profile_bin_count = options.profile_bin_count;
    checkpoint.configuration.profile_maximum_coordinate =
        options.profile_maximum_coordinate;
    checkpoint.progress = progress;
    checkpoint.initial_profile = initial_profile;
    checkpoint.previous_profile = previous_profile;
    checkpoint.state = state;
    return checkpoint;
}

ns_cascade::AdaptiveStepInfo fixedStepInformation(
    const ns_cascade::PseudospectralSystem& system,
    const ns_cascade::PseudospectralSystem::State& state,
    double time_step) {
    ns_cascade::AdaptiveStepInfo information;
    information.time_step = time_step;
    information.velocity_supremum_bound =
        system.velocitySupremumUpperBound(state);
    const double largest_wave_number =
        std::sqrt(3.0) * static_cast<double>(system.cutoff());
    const double maximum_wave_squared =
        3.0 * static_cast<double>(system.cutoff()) * system.cutoff();
    information.advective_cfl_upper_bound =
        time_step * information.velocity_supremum_bound * largest_wave_number;
    information.viscous_stability_number =
        time_step * system.viscosity() * maximum_wave_squared;
    return information;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        Options options = parseOptions(argc, argv);
        const bool restarted = !options.restart_path.empty();
        ns_cascade::SimulationCheckpoint restored_checkpoint;
        if (restarted) {
            restored_checkpoint =
                ns_cascade::loadSimulationCheckpoint(options.restart_path);
            restoreCheckpointConfiguration(
                options, restored_checkpoint.configuration);
        }
        const ns_cascade::PseudospectralSystem system(
            options.grid_size, options.viscosity, options.cutoff);
        ns_cascade::PseudospectralSystem::State state;
        ns_cascade::CheckpointProgress progress;
        ns_cascade::SpectrumProfile initial_profile;
        ns_cascade::SpectrumProfile previous_profile;
        if (restarted) {
            state = restored_checkpoint.state;
            progress = restored_checkpoint.progress;
            initial_profile = restored_checkpoint.initial_profile;
            previous_profile = restored_checkpoint.previous_profile;
        } else if (options.initial_condition ==
                   ns_cascade::InitialCondition::VortexTubes) {
            state = system.vortexTubePairState(options.tube_parameters,
                                               options.initial_energy);
        } else {
            state = system.initialState(options.initial_condition,
                                        options.initial_energy);
        }
        const double target_time = options.final_time_set
                                       ? options.final_time
                                       : options.steps * options.time_step;
        if (!std::isfinite(target_time) || target_time <= 0.0) {
            throw std::invalid_argument("Requested final time is not finite and positive");
        }
        const double time_tolerance =
            16.0 * std::numeric_limits<double>::epsilon() *
            std::max(1.0, target_time);
        if (restarted && target_time - progress.time <= time_tolerance) {
            throw std::invalid_argument(
                "Restart --final-time must be greater than the checkpoint time");
        }

        std::ofstream csv(options.output.c_str());
        std::ofstream shell_csv(options.shell_output.c_str());
        std::ofstream profile_csv(options.profile_output.c_str());
        if (!csv) throw std::runtime_error("Cannot open output file: " + options.output);
        if (!shell_csv) {
            throw std::runtime_error("Cannot open shell output file: " +
                                     options.shell_output);
        }
        if (!profile_csv) {
            throw std::runtime_error("Cannot open profile output file: " +
                                     options.profile_output);
        }
        csv << std::setprecision(17);
        shell_csv << std::setprecision(17);
        profile_csv << std::setprecision(17);
        writeTimeHeader(csv);
        writeShellHeader(shell_csv);
        writeProfileHeader(profile_csv);

        ns_cascade::PseudospectralSystem::Diagnostics diagnostics =
            system.diagnostics(state);
        if (diagnostics.divergence_defect > 1e-9 ||
            diagnostics.reality_defect > 1e-9) {
            throw std::runtime_error(
                "Initial or restarted state violates Fourier constraints");
        }
        ns_cascade::SpectrumProfile current_profile =
            ns_cascade::rescaledSpectrumProfile(
                system,
                state,
                options.profile_bin_count,
                options.profile_maximum_coordinate);
        if (!restarted) {
            progress.step = 0;
            progress.time = 0.0;
            progress.bkm_sampled_integral = 0.0;
            progress.previous_vorticity_max = diagnostics.sampled_vorticity_max;
            progress.previous_diagnostic_time = 0.0;
            progress.initial_critical_l3 = diagnostics.critical_l3_sample;
            progress.initial_critical_h_half = diagnostics.critical_h_half;
            progress.initial_sampled_vorticity =
                diagnostics.sampled_vorticity_max;
            progress.peak_high_shell_fraction =
                diagnostics.high_shell_energy_fraction;
            progress.peak_critical_l3 = diagnostics.critical_l3_sample;
            progress.peak_critical_h_half = diagnostics.critical_h_half;
            progress.peak_sampled_vorticity =
                diagnostics.sampled_vorticity_max;
            progress.maximum_production_to_dissipation =
                diagnostics.enstrophy_production_to_dissipation;
            progress.peak_forward_flux = 0.0;
            progress.minimum_time_step =
                std::numeric_limits<double>::infinity();
            progress.maximum_time_step = 0.0;
            progress.maximum_cfl_bound = 0.0;
            progress.maximum_viscous_number = 0.0;
            initial_profile = current_profile;
            previous_profile = current_profile;
        } else {
            progress.peak_high_shell_fraction = std::max(
                progress.peak_high_shell_fraction,
                diagnostics.high_shell_energy_fraction);
            progress.peak_critical_l3 = std::max(
                progress.peak_critical_l3, diagnostics.critical_l3_sample);
            progress.peak_critical_h_half = std::max(
                progress.peak_critical_h_half, diagnostics.critical_h_half);
            progress.peak_sampled_vorticity = std::max(
                progress.peak_sampled_vorticity,
                diagnostics.sampled_vorticity_max);
            progress.maximum_production_to_dissipation = std::max(
                progress.maximum_production_to_dissipation,
                diagnostics.enstrophy_production_to_dissipation);
        }
        const ns_cascade::AdaptiveStepInfo initial_step_information =
            fixedStepInformation(system, state, 0.0);
        writeTimeRow(csv,
                     progress.step,
                     progress.time,
                     options,
                     system,
                     diagnostics,
                     initial_step_information,
                     progress.bkm_sampled_integral,
                     current_profile,
                     initial_profile,
                     previous_profile);
        const double starting_forward_flux = writeShellRows(
            shell_csv,
            options,
            progress.step,
            progress.time,
            system,
            system.shellDiagnostics(state));
        progress.peak_forward_flux = std::max(
            progress.peak_forward_flux, starting_forward_flux);
        writeProfileRows(profile_csv,
                         progress.step,
                         progress.time,
                         system,
                         current_profile);
        ns_cascade::SpectrumProfileChange final_profile_change =
            ns_cascade::spectrumProfileChange(
                current_profile, previous_profile);

        const std::chrono::steady_clock::time_point start =
            std::chrono::steady_clock::now();
        while (progress.time < target_time) {
            if (progress.step >= 10000000U) {
                throw std::runtime_error("Adaptive integration exceeded 10000000 steps");
            }
            const double remaining_time = target_time - progress.time;
            const double permitted_time_step =
                std::min(options.time_step, remaining_time);
            ns_cascade::AdaptiveStepInfo step_information;
            if (options.adaptive) {
                step_information = system.chooseAdaptiveTimeStep(
                    state,
                    permitted_time_step,
                    options.target_cfl,
                    options.diffusion_safety);
            } else {
                step_information =
                    fixedStepInformation(system, state, permitted_time_step);
            }
            if (progress.time + step_information.time_step == progress.time) {
                throw std::runtime_error(
                    "Time step is too small to advance floating-point time");
            }
            system.stepRungeKutta4(state, step_information.time_step);
            progress.time += step_information.time_step;
            if (target_time - progress.time <= time_tolerance) {
                progress.time = target_time;
            }
            ++progress.step;

            progress.minimum_time_step = std::min(
                progress.minimum_time_step, step_information.time_step);
            progress.maximum_time_step = std::max(
                progress.maximum_time_step, step_information.time_step);
            progress.maximum_cfl_bound = std::max(
                progress.maximum_cfl_bound,
                step_information.advective_cfl_upper_bound);
            progress.maximum_viscous_number = std::max(
                progress.maximum_viscous_number,
                step_information.viscous_stability_number);
            const bool final_step = progress.time >= target_time;
            const bool scheduled_checkpoint_step =
                !options.checkpoint_output.empty() &&
                progress.step %
                        static_cast<std::uint64_t>(options.checkpoint_every) ==
                    0U;
            const bool diagnostic_step =
                progress.step %
                        static_cast<std::uint64_t>(options.diagnostic_every) ==
                    0U ||
                scheduled_checkpoint_step ||
                final_step;
            if (diagnostic_step) {
                diagnostics = system.diagnostics(state);
                current_profile = ns_cascade::rescaledSpectrumProfile(
                    system,
                    state,
                    options.profile_bin_count,
                    options.profile_maximum_coordinate);
                const double interval =
                    progress.time - progress.previous_diagnostic_time;
                progress.bkm_sampled_integral +=
                    0.5 * interval *
                    (progress.previous_vorticity_max +
                     diagnostics.sampled_vorticity_max);
                progress.previous_vorticity_max =
                    diagnostics.sampled_vorticity_max;
                progress.previous_diagnostic_time = progress.time;
                progress.peak_high_shell_fraction = std::max(
                    progress.peak_high_shell_fraction,
                    diagnostics.high_shell_energy_fraction);
                progress.peak_critical_l3 = std::max(
                    progress.peak_critical_l3,
                    diagnostics.critical_l3_sample);
                progress.peak_critical_h_half = std::max(
                    progress.peak_critical_h_half,
                    diagnostics.critical_h_half);
                progress.peak_sampled_vorticity = std::max(
                    progress.peak_sampled_vorticity,
                    diagnostics.sampled_vorticity_max);
                progress.maximum_production_to_dissipation = std::max(
                    progress.maximum_production_to_dissipation,
                    diagnostics.enstrophy_production_to_dissipation);
                final_profile_change = ns_cascade::spectrumProfileChange(
                    current_profile, previous_profile);
                writeTimeRow(csv,
                             progress.step,
                             progress.time,
                             options,
                             system,
                             diagnostics,
                             step_information,
                             progress.bkm_sampled_integral,
                             current_profile,
                             initial_profile,
                             previous_profile);
                progress.peak_forward_flux = std::max(
                    progress.peak_forward_flux,
                    writeShellRows(shell_csv,
                                   options,
                                   progress.step,
                                   progress.time,
                                   system,
                                   system.shellDiagnostics(state)));
                writeProfileRows(profile_csv,
                                 progress.step,
                                 progress.time,
                                 system,
                                 current_profile);
                previous_profile = current_profile;
            }
            const bool checkpoint_step =
                !options.checkpoint_output.empty() &&
                (scheduled_checkpoint_step || final_step);
            if (checkpoint_step) {
                ns_cascade::saveSimulationCheckpoint(
                    options.checkpoint_output,
                    makeCheckpoint(options,
                                   system,
                                   state,
                                   progress,
                                   initial_profile,
                                   previous_profile));
            }
        }
        const std::chrono::duration<double> elapsed =
            std::chrono::steady_clock::now() - start;

        std::cout << std::setprecision(8)
                  << "Initial condition: "
                  << ns_cascade::initialConditionName(options.initial_condition) << '\n'
                  << "Grid: " << system.gridSize() << "^3, retained cutoff K="
                  << system.cutoff() << " (safe maximum "
                  << system.safeDealiasCutoff() << ")\n"
                  << "Retained non-zero modes: " << system.modeCount() << '\n'
                  << "Time stepping: "
                  << (options.adaptive ? "adaptive" : "fixed") << ", "
                  << progress.step
                  << " accepted steps to t=" << target_time << '\n'
                  << "Accepted dt range: [" << progress.minimum_time_step
                  << ", " << progress.maximum_time_step << "]\n"
                  << "Maximum conservative CFL bound: "
                  << progress.maximum_cfl_bound
                  << " (target " << options.target_cfl << ")\n"
                  << "Maximum viscous stability number: "
                  << progress.maximum_viscous_number << '\n'
                  << "Final normalized energy: " << diagnostics.energy << '\n'
                  << "Peak sampled L3 norm: " << progress.peak_critical_l3 << '\n'
                  << "Peak/initial sampled L3 ratio: "
                  << progress.peak_critical_l3 /
                         progress.initial_critical_l3
                  << '\n'
                  << "Peak/initial H1/2 ratio: "
                  << progress.peak_critical_h_half /
                         progress.initial_critical_h_half
                  << '\n'
                  << "Peak/initial sampled vorticity ratio: "
                  << progress.peak_sampled_vorticity /
                         progress.initial_sampled_vorticity
                  << '\n'
                  << "Maximum enstrophy production/dissipation ratio: "
                  << progress.maximum_production_to_dissipation << '\n'
                  << "Sampled BKM integral: "
                  << progress.bkm_sampled_integral << '\n'
                  << "Peak cutoff-shell energy fraction: "
                  << progress.peak_high_shell_fraction << '\n'
                  << "Peak forward shell flux: "
                  << progress.peak_forward_flux << '\n'
                  << "Characteristic-wavenumber ratio: "
                  << current_profile.characteristic_wavenumber /
                         initial_profile.characteristic_wavenumber
                  << '\n'
                  << "Final rescaled-profile L1 change from initial: "
                  << ns_cascade::spectrumProfileL1Distance(
                         current_profile, initial_profile)
                  << '\n'
                  << "Final rescaled-profile L1 change from prior sample: "
                  << final_profile_change.l1_distance << '\n'
                  << "Final absolute log spectral-scale change: "
                  << final_profile_change.absolute_log_scale_change << '\n'
                  << "Final rescaled-profile L1/log-scale drift: "
                  << final_profile_change.l1_per_log_scale_change
                  << " (valid="
                  << (final_profile_change.scale_normalized_drift_valid
                          ? "true"
                          : "false")
                  << ")\n"
                  << "Final analyticity-radius estimate: "
                  << current_profile.analyticity_radius_estimate << " (R^2="
                  << current_profile.analyticity_fit_r_squared << ", shells="
                  << current_profile.analyticity_fit_shells << ")\n"
                  << "Final rescaled tail-slope product delta*k_rms: "
                  << current_profile.analyticity_radius_estimate *
                         current_profile.characteristic_wavenumber
                  << '\n'
                  << "Evolution time: " << elapsed.count() << " seconds\n"
                  << "Diagnostics written to " << options.output << '\n'
                  << "Shell transfers written to " << options.shell_output << '\n'
                  << "Rescaled profiles written to " << options.profile_output
                  << '\n';
        if (restarted) {
            std::cout << "Restarted from " << options.restart_path << "\n";
        }
        if (!options.checkpoint_output.empty()) {
            std::cout << "Final checkpoint written to "
                      << options.checkpoint_output << '\n';
        }
        if (options.initial_condition == ns_cascade::InitialCondition::VortexTubes) {
            std::cout << "Vortex tubes: core="
                      << options.tube_parameters.core_radius << ", separation="
                      << options.tube_parameters.separation << ", bend="
                      << options.tube_parameters.bend_amplitude
                      << ", axial mode="
                      << options.tube_parameters.axial_wavenumber << '\n';
        }
        if (progress.peak_high_shell_fraction > 0.01) {
            std::cout
                << "Resolution warning: more than 1% of energy reached the cutoff shell.\n";
        }
        const double characteristic_growth =
            current_profile.characteristic_wavenumber /
            initial_profile.characteristic_wavenumber;
        if (characteristic_growth >= 1.1 &&
            final_profile_change.scale_normalized_drift_valid &&
            final_profile_change.l1_per_log_scale_change <= 1.0 &&
            progress.peak_high_shell_fraction <= 0.01) {
            std::cout
                << "Profile candidate heuristic: spectral scale moved at least 10% with low rescaled-shape drift; refinement is still required.\n";
        } else {
            std::cout
                << "No stable moving rescaled-profile signal passed the heuristic gate (L1/|d log k_rms| <= 1).\n";
        }
        std::cout
            << "This dealiased finite grid still cannot prove regularity or blow-up.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
