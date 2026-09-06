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
        } else if (flag == "--cutoff") {
            options.cutoff = parseNumber<int>(requireValue(i, argc, argv), flag);
        } else if (flag == "--viscosity") {
            options.viscosity = parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--dt") {
            options.time_step = parseNumber<double>(requireValue(i, argc, argv), flag);
        } else if (flag == "--final-time") {
            options.final_time =
                parseNumber<double>(requireValue(i, argc, argv), flag);
            options.final_time_set = true;
        } else if (flag == "--steps") {
            options.steps = parseNumber<int>(requireValue(i, argc, argv), flag);
        } else if (flag == "--adaptive-dt") {
            options.adaptive = true;
        } else if (flag == "--fixed-dt") {
            options.adaptive = false;
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
        } else if (flag == "--initial-condition") {
            options.initial_condition =
                parseInitialCondition(requireValue(i, argc, argv));
        } else if (flag == "--tube-core") {
            options.tube_parameters.core_radius =
                parseNumber<double>(requireValue(i, argc, argv), flag);
            options.custom_tube_parameters = true;
        } else if (flag == "--tube-separation") {
            options.tube_parameters.separation =
                parseNumber<double>(requireValue(i, argc, argv), flag);
            options.custom_tube_parameters = true;
        } else if (flag == "--tube-bend") {
            options.tube_parameters.bend_amplitude =
                parseNumber<double>(requireValue(i, argc, argv), flag);
            options.custom_tube_parameters = true;
        } else if (flag == "--tube-axial-mode") {
            options.tube_parameters.axial_wavenumber =
                parseNumber<int>(requireValue(i, argc, argv), flag);
            options.custom_tube_parameters = true;
        } else if (flag == "--output") {
            options.output = requireValue(i, argc, argv);
        } else if (flag == "--shell-output") {
            options.shell_output = requireValue(i, argc, argv);
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
    if (options.output.empty() || options.shell_output.empty()) {
        throw std::invalid_argument("Output paths cannot be empty");
    }
    if (options.output == options.shell_output) {
        throw std::invalid_argument("--output and --shell-output must be different");
    }
    return options;
}

void writeTimeHeader(std::ostream& output) {
    output
        << "step,time,dt_used,adaptive,advective_cfl_upper_bound,"
        << "viscous_stability_number,velocity_supremum_bound,grid_size,cutoff,"
        << "initial_condition,tube_core_radius,tube_separation,"
        << "tube_bend_amplitude,tube_axial_wavenumber,"
        << "energy,enstrophy,palinstrophy,"
        << "critical_l3_sample,sampled_vorticity_max,"
        << "vorticity_sup_upper_bound,spectral_centroid,"
        << "high_shell_energy_fraction,divergence_defect,reality_defect,"
        << "energy_balance_residual,bkm_sampled_integral\n";
}

void writeTimeRow(std::ostream& output,
                  int step,
                  double time,
                  const Options& options,
                  const ns_cascade::PseudospectralSystem& system,
                  const ns_cascade::PseudospectralSystem::Diagnostics& values,
                  const ns_cascade::AdaptiveStepInfo& step_information,
                  double bkm_integral) {
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
           << ',' << values.palinstrophy << ',' << values.critical_l3_sample << ','
           << values.sampled_vorticity_max << ','
           << values.vorticity_sup_upper_bound << ',' << values.spectral_centroid
           << ',' << values.high_shell_energy_fraction << ','
           << values.divergence_defect << ',' << values.reality_defect << ','
           << values.energy_balance_residual << ',' << bkm_integral << '\n';
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
    int step,
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
        const Options options = parseOptions(argc, argv);
        const ns_cascade::PseudospectralSystem system(
            options.grid_size, options.viscosity, options.cutoff);
        ns_cascade::PseudospectralSystem::State state;
        if (options.initial_condition == ns_cascade::InitialCondition::VortexTubes) {
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

        std::ofstream csv(options.output.c_str());
        std::ofstream shell_csv(options.shell_output.c_str());
        if (!csv) throw std::runtime_error("Cannot open output file: " + options.output);
        if (!shell_csv) {
            throw std::runtime_error("Cannot open shell output file: " +
                                     options.shell_output);
        }
        csv << std::setprecision(17);
        shell_csv << std::setprecision(17);
        writeTimeHeader(csv);
        writeShellHeader(shell_csv);

        ns_cascade::PseudospectralSystem::Diagnostics diagnostics =
            system.diagnostics(state);
        double bkm_sampled_integral = 0.0;
        double previous_vorticity_max = diagnostics.sampled_vorticity_max;
        double peak_high_shell_fraction = diagnostics.high_shell_energy_fraction;
        double peak_critical_l3 = diagnostics.critical_l3_sample;
        double peak_sampled_vorticity = diagnostics.sampled_vorticity_max;
        const double initial_critical_l3 = diagnostics.critical_l3_sample;
        const double initial_sampled_vorticity = diagnostics.sampled_vorticity_max;
        double previous_diagnostic_time = 0.0;
        const ns_cascade::AdaptiveStepInfo initial_step_information =
            fixedStepInformation(system, state, 0.0);
        writeTimeRow(csv,
                     0,
                     0.0,
                     options,
                     system,
                     diagnostics,
                     initial_step_information,
                     bkm_sampled_integral);
        double peak_forward_flux = writeShellRows(
            shell_csv,
            options,
            0,
            0.0,
            system,
            system.shellDiagnostics(state));

        const std::chrono::steady_clock::time_point start =
            std::chrono::steady_clock::now();
        int step = 0;
        double time = 0.0;
        double minimum_time_step = std::numeric_limits<double>::infinity();
        double maximum_time_step = 0.0;
        double maximum_cfl_bound = 0.0;
        double maximum_viscous_number = 0.0;
        const double time_tolerance =
            16.0 * std::numeric_limits<double>::epsilon() *
            std::max(1.0, target_time);
        while (time < target_time) {
            if (step >= 10000000) {
                throw std::runtime_error("Adaptive integration exceeded 10000000 steps");
            }
            const double remaining_time = target_time - time;
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
            if (time + step_information.time_step == time) {
                throw std::runtime_error(
                    "Time step is too small to advance floating-point time");
            }
            system.stepRungeKutta4(state, step_information.time_step);
            time += step_information.time_step;
            if (target_time - time <= time_tolerance) time = target_time;
            ++step;

            minimum_time_step =
                std::min(minimum_time_step, step_information.time_step);
            maximum_time_step =
                std::max(maximum_time_step, step_information.time_step);
            maximum_cfl_bound = std::max(
                maximum_cfl_bound,
                step_information.advective_cfl_upper_bound);
            maximum_viscous_number = std::max(
                maximum_viscous_number,
                step_information.viscous_stability_number);
            const bool final_step = time >= target_time;
            if (step % options.diagnostic_every != 0 && !final_step) continue;

            diagnostics = system.diagnostics(state);
            const double interval = time - previous_diagnostic_time;
            bkm_sampled_integral +=
                0.5 * interval *
                (previous_vorticity_max + diagnostics.sampled_vorticity_max);
            previous_vorticity_max = diagnostics.sampled_vorticity_max;
            previous_diagnostic_time = time;
            peak_high_shell_fraction =
                std::max(peak_high_shell_fraction,
                         diagnostics.high_shell_energy_fraction);
            peak_critical_l3 =
                std::max(peak_critical_l3, diagnostics.critical_l3_sample);
            peak_sampled_vorticity = std::max(
                peak_sampled_vorticity, diagnostics.sampled_vorticity_max);
            writeTimeRow(csv,
                         step,
                         time,
                         options,
                         system,
                         diagnostics,
                         step_information,
                         bkm_sampled_integral);
            peak_forward_flux = std::max(
                peak_forward_flux,
                writeShellRows(shell_csv,
                               options,
                               step,
                               time,
                               system,
                               system.shellDiagnostics(state)));
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
                  << (options.adaptive ? "adaptive" : "fixed") << ", " << step
                  << " accepted steps to t=" << target_time << '\n'
                  << "Accepted dt range: [" << minimum_time_step << ", "
                  << maximum_time_step << "]\n"
                  << "Maximum conservative CFL bound: " << maximum_cfl_bound
                  << " (target " << options.target_cfl << ")\n"
                  << "Maximum viscous stability number: "
                  << maximum_viscous_number << '\n'
                  << "Final normalized energy: " << diagnostics.energy << '\n'
                  << "Peak sampled L3 norm: " << peak_critical_l3 << '\n'
                  << "Peak/initial sampled L3 ratio: "
                  << peak_critical_l3 / initial_critical_l3 << '\n'
                  << "Peak/initial sampled vorticity ratio: "
                  << peak_sampled_vorticity / initial_sampled_vorticity << '\n'
                  << "Sampled BKM integral: " << bkm_sampled_integral << '\n'
                  << "Peak cutoff-shell energy fraction: "
                  << peak_high_shell_fraction << '\n'
                  << "Peak forward shell flux: " << peak_forward_flux << '\n'
                  << "Evolution time: " << elapsed.count() << " seconds\n"
                  << "Diagnostics written to " << options.output << '\n'
                  << "Shell transfers written to " << options.shell_output << '\n';
        if (options.initial_condition == ns_cascade::InitialCondition::VortexTubes) {
            std::cout << "Vortex tubes: core="
                      << options.tube_parameters.core_radius << ", separation="
                      << options.tube_parameters.separation << ", bend="
                      << options.tube_parameters.bend_amplitude
                      << ", axial mode="
                      << options.tube_parameters.axial_wavenumber << '\n';
        }
        if (peak_high_shell_fraction > 0.01) {
            std::cout
                << "Resolution warning: more than 1% of energy reached the cutoff shell.\n";
        }
        std::cout
            << "This dealiased finite grid still cannot prove regularity or blow-up.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
